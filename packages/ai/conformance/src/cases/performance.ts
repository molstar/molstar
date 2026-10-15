/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import type { Message } from '@molstar/ai-agent';
import {
  type BeginInput,
  buildMove,
  type CommitBatch,
  createRunner,
  defaultDrive,
  type EntryId,
  type HeadId,
  type HeadRecord,
  type Inspectable,
  MAIN,
  type Move,
  type RetentionSnapshot,
  type Runner,
  type RunnerOptions,
  type RunnerPort,
  SCHEMA_VERSION,
  type ScopeId,
  type Segment,
  type ToolFrame,
} from '@molstar/ai-harness';
import { difference, type StorageCounts, totalReads, totalWrites } from '../meter.js';
import { type ConformanceFixture, createSession, type Session, type SessionSetup } from '../session.js';
import { simulatedIds } from '../testing/ids.js';
import { assistant, declareTools, user } from '../testing/messages.js';
import type { ModelStep } from '../testing/model.js';
import type { ToolRun, ToolSpec } from '../testing/tools.js';
import { expectWith, hashes, stepEvents, TURN, with_ } from './helpers.js';

/** Slack per comparison, as in anyfold's suite. */
export const SLACK = { records: 16, bytes: 16 * 1024 };

/** "A few" segments per moved pointer in one commit (P2 settlement-shape). */
const FEW = 4;

// ---------------------------------------------------------------------------------------------
// Comparing counts
// ---------------------------------------------------------------------------------------------

const METRICS: readonly (readonly [string, (counts: StorageCounts) => number, number])[] = [
  ['transactions', (c) => c.transactions, SLACK.records],
  ['records read', totalReads, SLACK.records],
  ['records written', totalWrites, SLACK.records],
  ['bytes read', (c) => c.bytesRead, SLACK.bytes],
  ['bytes written', (c) => c.bytesWritten, SLACK.bytes],
];

/** Every metered quantity of `larger` is at most `factor` times `smaller`'s, plus slack. */
function expectGrowth(label: string, smaller: StorageCounts, larger: StorageCounts, factor: number): void {
  for (const [name, value, slack] of METRICS)
    expectWith(value(larger), `${label}: ${name}`).toBeLessThanOrEqual(factor * value(smaller) + slack);
}

/**
 * A series of runs at increasing sizes: each is at most (size ratio, or `cap` if lower) times the one
 * before, plus slack. Doubling sizes give the "at most twice" of P1; 4/40/80 and 1/8/64 are bounded
 * by their ratio, which is the linear bound.
 */
function expectSeries(label: string, sizes: readonly number[], costs: readonly StorageCounts[], cap?: number): void {
  for (let i = 1; i < sizes.length; i += 1) {
    const [smaller, larger] = [sizes[i - 1] as number, sizes[i] as number];
    const factor = Math.min(larger / smaller, cap ?? Number.POSITIVE_INFINITY);
    expectGrowth(`${label} ${smaller}->${larger}`, costs[i - 1] as StorageCounts, costs[i] as StorageCounts, factor);
  }
}

/** `other` costs the same as `base`, within slack. */
const expectSame = (label: string, base: StorageCounts, other: StorageCounts) => expectGrowth(label, base, other, 1);

/** Storage reads of the kinds that hold the journal (heads, segments, messages) and blobs. */
const journalReads = (c: StorageCounts) => c.reads.heads + c.reads.segments + c.reads.messages;

/** The cost of each step between consecutive marks. */
const between = (marks: readonly StorageCounts[]): StorageCounts[] =>
  marks.slice(1).map((mark, index) => difference(mark, marks[index] as StorageCounts));

/** Writes of every kind but the Execution, events and blobs: what a nested run must never cause. */
const durableWrites = (c: StorageCounts) =>
  c.writes.execution + c.writes.events + c.writes.messages + c.writes.pointers + c.writes.heads + c.writes.segments;

// ---------------------------------------------------------------------------------------------
// Workload builders
// ---------------------------------------------------------------------------------------------

interface Holder {
  s?: Session;
  runner?: Runner;
}
function need<T>(value: T | undefined): T {
  if (value === undefined) throw new Error('performance case: used before it was set');
  return value;
}
const retentionOf = (runner: Runner): RetentionSnapshot => (runner as unknown as Inspectable).retention();

const replies = (count: number): ModelStep[] => Array.from({ length: count }, (_, i) => ({ text: `reply ${i}` }));
const callsOf = (tool: string, count: number): NonNullable<ModelStep['calls']> =>
  Array.from({ length: count }, (_, i) => ({ tool, arguments: { i } }));
const work: ToolRun = (context) => ({ text: `worked ${JSON.stringify(context.arguments)}` });

/** A session whose first `history` prompts are already answered; `tail` is the model's script after them. */
async function seeded(
  fixture: ConformanceFixture,
  history: number,
  tail: ModelStep[],
  tools: Readonly<Record<string, ToolSpec | ToolRun>> = {},
  setup: SessionSetup = {},
): Promise<Session> {
  const s = await createSession(fixture, { ...setup, model: [...replies(history), ...tail], tools });
  for (let i = 0; i < history; i += 1) await s.send(`past ${i}`);
  return s;
}

/** Submits the prompts as items (the first with the tools' declaration) and returns the launch key. */
async function queue(s: Session, prompts: readonly (string | readonly Message[])[]): Promise<string> {
  let launch: string | undefined;
  for (const [i, prompt] of prompts.entries()) {
    const sent = await s.send(prompt, { key: `item-${i + 1}`, run: false });
    launch ??= sent.launch;
  }
  return need(launch);
}

/**
 * A next-turn drive that begins one waiting item per turn and releases when none waits; `before` runs
 * at each decision. Its boundaries steer nothing (the default steers every waiting item into the turn).
 */
const oneAtATime = (before: () => void) =>
  with_({
    boundary: () => 'proceed',
    'next-turn': (_decision, context) => {
      before();
      const item = context.inbox[0];
      return item === undefined ? 'release' : { begin: { keys: [item.key] } };
    },
  });

interface RoundsRun {
  s: Session;
  runner: Runner;
  /** The storage counts at the start of each round's tool. */
  marks: StorageCounts[];
  /** The runner's retention at the start of each round's tool, when `inspect`. */
  retained: RetentionSnapshot[];
  cost: StorageCounts;
}

/** One turn of `rounds` serial model->tool rounds on one warm runner, then a final answer. */
async function serialRounds(
  fixture: ConformanceFixture,
  rounds: number,
  options: { inspect?: boolean; onTool?: () => void } = {},
): Promise<RoundsRun> {
  const ref: Holder = {};
  const marks: StorageCounts[] = [];
  const retained: RetentionSnapshot[] = [];
  const s = await createSession(fixture, {
    tools: {
      work: (context) => {
        options.onTool?.();
        marks.push(need(ref.s).counts());
        if (options.inspect === true) retained.push(retentionOf(need(ref.runner)));
        return work(context);
      },
    },
    model: [
      ...Array.from({ length: rounds }, (_, round) => ({ calls: [{ tool: 'work', arguments: { round } }] })),
      { text: 'done' },
    ],
  });
  ref.s = s;
  const runner = s.runner();
  ref.runner = runner;
  const before = s.counts();
  await s.send('go', { runner });
  return { s, runner, marks, retained, cost: difference(s.counts(), before) };
}

/** One tools step of `width` calls. */
async function toolWidth(fixture: ConformanceFixture, width: number) {
  const s = await createSession(fixture, {
    tools: { work },
    model: [{ calls: callsOf('work', width) }, { text: 'done' }],
  });
  const before = s.counts();
  await s.send('go');
  return { s, cost: difference(s.counts(), before) };
}

/** One answer streamed as `frames` text deltas. */
async function streamed(fixture: ConformanceFixture, frames: number) {
  const deltas = Array.from({ length: frames }, () => 'x');
  const s = await createSession(fixture, { model: [{ text: deltas.join(''), deltas }] });
  const runner = s.runner();
  const before = s.counts();
  await s.send('go', { runner });
  return { s, runner, cost: difference(s.counts(), before) };
}

/** A tool streaming `frames` progress frames, then its result. */
async function toolFrames(fixture: ConformanceFixture, frames: number) {
  const s = await createSession(fixture, {
    tools: {
      stream: async function* (): AsyncGenerator<ToolFrame> {
        for (let i = 0; i < frames; i += 1) yield { kind: 'progress', value: { i } };
        yield { kind: 'result', status: 'succeeded', parts: [{ kind: 'text', text: 'streamed' }] };
      },
    },
    model: [{ calls: callsOf('stream', 1) }, { text: 'done' }],
  });
  const before = s.counts();
  await s.send('go');
  return { s, cost: difference(s.counts(), before) };
}

/** `turns` prompts answered in a row by one warm runner, one turn at a time; marks at each turn's start and the end. */
async function turnsInARow(fixture: ConformanceFixture, turns: number) {
  const ref: Holder = {};
  const marks: StorageCounts[] = [];
  const s = await createSession(fixture, {
    model: replies(turns),
    drive: oneAtATime(() => marks.push(need(ref.s).counts())),
  });
  ref.s = s;
  const launch = await queue(
    s,
    Array.from({ length: turns }, (_, i) => `prompt ${i}`),
  );
  const before = s.counts();
  await s.run(launch);
  return { s, marks, cost: difference(s.counts(), before) };
}

/** A prompt after `history` answered ones: the cost of opening the turn (until its first tool starts) and the rest. */
async function openTurn(fixture: ConformanceFixture, history: number) {
  const ref: Holder = {};
  const marks: StorageCounts[] = [];
  const s = await seeded(fixture, history, [{ calls: callsOf('work', 1) }, { text: 'final' }], {
    work: (context) => {
      marks.push(need(ref.s).counts());
      return work(context);
    },
  });
  ref.s = s;
  const { launch } = await s.send('measure', { run: false });
  const before = s.counts();
  await s.run(need(launch));
  const after = s.counts();
  const started = marks[0] as StorageCounts;
  return { open: difference(started, before), warm: difference(after, started), total: difference(after, before) };
}

interface Recovered {
  /** From the new runner's claim to the re-run call starting. */
  open: StorageCounts;
  /** From there to the end of the turn. */
  warm: StorageCounts;
  modelsBefore: number;
  modelsAtResume: number;
  toolsBefore: number;
  toolsAtResume: number;
}

/**
 * A runner that dies in the middle of a tool call (aborted locally, nothing committed), its lease
 * expiring, and a cold runner claiming the unowned turn after `history` answered prompts.
 */
async function recover(fixture: ConformanceFixture, history: number): Promise<Recovered> {
  const ref: Holder & { resumed?: { counts: StorageCounts; models: number; tools: number } } = {};
  let crashed = false;
  const s = await seeded(
    fixture,
    history,
    [{ calls: callsOf('work', 1) }, { text: 'final' }],
    {
      work: () => {
        const session = need(ref.s);
        if (!crashed) {
          crashed = true;
          need(ref.runner).abort('crash');
          return { text: 'lost with the runner' };
        }
        ref.resumed = {
          counts: session.counts(),
          models: session.model.invocations.length,
          tools: session.tools.calls.length,
        };
        return { text: 'found again' };
      },
    },
    { host: { leaseMs: 30_000 } },
  );
  ref.s = s;
  const { launch } = await s.send('now', { run: false });
  const first = s.runner();
  ref.runner = first;
  expect((await s.run(need(launch), first)).reason).toBe('aborted');
  await s.clock.advance(31_000);
  const { launch: relaunch } = await s.client.check();
  const before = s.counts();
  const [modelsBefore, toolsBefore] = [s.model.invocations.length, s.tools.calls.length];
  await s.run(need(relaunch));
  const after = s.counts();
  const resumed = need(ref.resumed);
  return {
    open: difference(resumed.counts, before),
    warm: difference(after, resumed.counts),
    modelsBefore,
    modelsAtResume: resumed.models,
    toolsBefore,
    toolsAtResume: resumed.tools,
  };
}

/** A tool `investigate` running a nested run of `rounds` rounds of an `inner` tool. */
async function nestedRounds(fixture: ConformanceFixture, rounds: number) {
  const ref: Holder & { span?: StorageCounts } = {};
  const s = await createSession(fixture, {
    declare: ['investigate'],
    tools: {
      inner: work,
      investigate: async (context) => {
        const session = need(ref.s);
        const before = session.counts();
        const result = await context.run({
          input: [user('look closely')],
          tools: { enabled: ['inner'], preauthorized: ['inner'], definitions: [session.tools.definition('inner')] },
          model: TURN.model,
          limits: {},
        });
        ref.span = difference(session.counts(), before);
        return { text: result.status };
      },
    },
    model: [
      { calls: callsOf('investigate', 1) },
      ...Array.from({ length: rounds }, () => ({ calls: callsOf('inner', 1) })),
      { text: 'inner done' },
      { text: 'final' },
    ],
  });
  ref.s = s;
  const before = s.counts();
  await s.send('go');
  return { s, cost: difference(s.counts(), before), span: need(ref.span) };
}

// ---------------------------------------------------------------------------------------------
// Protocol-level turns (one commit at a time) and the shared head scheme
// ---------------------------------------------------------------------------------------------

/** A runner's view of its heads: one `buildMove` over a local cache of segments, per pointer. */
function headBuilder() {
  const segments = new Map<string, Segment>();
  const heads = new Map<ScopeId, HeadRecord>();
  return {
    async append(pointer: ScopeId, entries: readonly EntryId[], start?: HeadId): Promise<Move> {
      const from = heads.get(pointer);
      const move = await buildMove({
        pointer,
        ...(from !== undefined ? { from } : {}),
        edit: { append: entries },
        read: (ids) => ids.flatMap((id) => segments.get(id) ?? []),
      });
      for (const segment of move.segments) segments.set(segment.id, segment);
      heads.set(pointer, move.head);
      const base = from?.id ?? start;
      return { ...move, ...(base !== undefined ? { from: base } : {}) };
    },
  };
}

/** A claimed lease and a begun turn (`turn-x`, pending `step-2`); `commit` settles one model step at a time. */
async function protocolTurn(s: Session) {
  const hello = user('hello');
  const submitted = await s.client.submit({ key: 'item-1', messages: [hello] });
  const claim = { launch: need(submitted.launch), runnerId: 'r1', schemaVersion: SCHEMA_VERSION };
  const claimed = await s.port.claim(claim);
  const builder = headBuilder();
  const move = await builder.append(MAIN, await hashes([hello]), claimed.execution.main);
  const begin: BeginInput = {
    key: 'begin-1',
    turnId: 'turn-x',
    step: 'step-1',
    items: ['item-1'],
    configuration: TURN,
    messages: [hello],
    move,
    calls: [],
    next: { id: 'step-2', kind: 'model' },
  };
  const begun = await s.port.begin(claimed.lease, begin);
  let lease = begun.lease;
  let step = 2;
  return {
    submitted,
    claim,
    claimed,
    begin,
    begun,
    /** Settles the pending model step with these appends (per pointer); returns the commit's own cost. */
    async commit(appends: Readonly<Record<ScopeId, readonly Message[]>>, end = false) {
      const moves: Move[] = [];
      const messages: Message[] = [];
      for (const [scope, list] of Object.entries(appends)) {
        moves.push(await builder.append(scope, await hashes(list)));
        messages.push(...list);
      }
      const id = step;
      step += 1;
      const batch: CommitBatch = {
        key: `commit-${id}`,
        turnId: 'turn-x',
        step: `step-${id}`,
        ending: 'settled',
        calls: [],
        messages,
        moves,
        consumed: [],
        ...(end
          ? { end: { status: 'completed' as const } }
          : { next: { id: `step-${id + 1}`, kind: 'model' as const } }),
      };
      const before = s.counts();
      const committed = await s.port.commit(lease, batch);
      const cost = difference(s.counts(), before);
      lease = committed.lease;
      return { cost, moves, batch, committed };
    },
  };
}

// ---------------------------------------------------------------------------------------------
// Observing a runner's own messages to the host
// ---------------------------------------------------------------------------------------------

interface PortCall {
  readonly op: string;
  readonly at: number;
  readonly args: readonly unknown[];
}

/** A runner over a recording proxy of the session's port: what it sent, in order, and when (virtual time). */
function observedRunner(s: Session, overrides: Partial<Omit<RunnerOptions, 'port'>> = {}) {
  const log: PortCall[] = [];
  const port = new Proxy(s.port, {
    get(target, prop) {
      const value: unknown = Reflect.get(target, prop, target);
      if (typeof value !== 'function') return value;
      return (...args: unknown[]) => {
        log.push({ op: String(prop), at: s.clock.now(), args });
        return value.apply(target, args);
      };
    },
  }) as RunnerPort;
  const ids = s.tools.definitions.map((definition) => definition.tool.id);
  const runner = createRunner({
    port,
    runnerId: 'observed',
    models: s.model.bindings,
    tools: s.tools.bindings,
    drive: defaultDrive(),
    configure: () => ({ ...TURN, tools: { enabled: ids, preauthorized: ids } }),
    clock: s.clock,
    ids: simulatedIds('observed-'),
    ...overrides,
  });
  return { runner, log };
}

/** Counts `crypto.subtle.digest` calls (every hash the harness takes goes through it) while `run` works. */
async function withDigestSpy<T>(run: (digests: { calls(): number; bytes(): number }) => Promise<T>): Promise<T> {
  const spy = vi.spyOn(globalThis.crypto.subtle, 'digest');
  try {
    return await run({
      calls: () => spy.mock.calls.length,
      bytes: () => spy.mock.calls.reduce((sum, call) => sum + (call[1] as { byteLength: number }).byteLength, 0),
    });
  } finally {
    spy.mockRestore();
  }
}

const ROUNDS = [2, 4, 8, 16, 32];

/** P1–P11 (Performance): bounds on metered counts, never on time. */
export function performanceCases(fixture: ConformanceFixture): void {
  describe('Performance', () => {
    // -------------------------------------------------------------------------------------------
    // P1. Doubling a workload at most doubles host work.
    // -------------------------------------------------------------------------------------------

    it('P1 history: a turn after 4, 40 and 80 turns of history costs at most the size ratio of the smaller', async () => {
      const sizes = [4, 40, 80];
      const costs: StorageCounts[] = [];
      for (const size of sizes) costs.push((await openTurn(fixture, size)).total);
      expectSeries('history', sizes, costs);
    });

    it('P1 serial-rounds: 2 to 32 serial rounds, each doubling costs at most twice the smaller', async () => {
      const costs: StorageCounts[] = [];
      for (const rounds of ROUNDS) costs.push((await serialRounds(fixture, rounds)).cost);
      expectSeries('rounds', ROUNDS, costs);
    });

    it('P1 tool-width: 1, 8 and 64 calls in one step cost at most the width ratio of the smaller', async () => {
      const sizes = [1, 8, 64];
      const costs: StorageCounts[] = [];
      for (const width of sizes) costs.push((await toolWidth(fixture, width)).cost);
      expectSeries('width', sizes, costs);
    });

    it('P1 frames: 1 and 100 stream frames cost at most twice the smaller', async () => {
      const costs = [(await streamed(fixture, 1)).cost, (await streamed(fixture, 100)).cost];
      expectSeries('frames', [1, 100], costs, 2);
    });

    it('P1 sequential-turns: 3, 6 and 12 turns in a row, each doubling costs at most twice the smaller', async () => {
      const sizes = [3, 6, 12];
      const costs: StorageCounts[] = [];
      for (const turns of sizes) costs.push((await turnsInARow(fixture, turns)).cost);
      expectSeries('turns', sizes, costs);
    });

    it('P1 pages: reading a page of 5 and of 20 retained turns costs at most the ratio of the smaller', async () => {
      const { s } = await turnsInARow(fixture, 20);
      const sizes = [5, 20];
      const costs: StorageCounts[] = [];
      for (const limit of sizes) {
        const before = s.counts();
        expect(await s.client.turns({ limit })).toHaveLength(limit);
        costs.push(difference(s.counts(), before));
      }
      expectSeries('pages', sizes, costs);
    });

    it('P1 nested-rounds: a nested run of 2 to 32 rounds, each doubling costs at most twice the smaller', async () => {
      const costs: StorageCounts[] = [];
      for (const rounds of ROUNDS) costs.push((await nestedRounds(fixture, rounds)).cost);
      expectSeries('nested rounds', ROUNDS, costs);
    });

    // -------------------------------------------------------------------------------------------
    // P2. A commit costs a constant plus what it adds.
    // -------------------------------------------------------------------------------------------

    it('P2 settlement-shape: a commit costs the Execution, each moved pointer and what it adds, and reads no old record', async () => {
      const s = await createSession(fixture);
      const t = await protocolTurn(s);
      // The first commit creates the side pointer; the measured one moves it and `main`.
      await t.commit({ [MAIN]: [assistant('a1')], side: [assistant('s1')] });
      const { cost: d } = await t.commit({ [MAIN]: [user('u2'), assistant('a3')], side: [assistant('s2')] });
      expect(d.transactions).toBe(1);
      expect([d.reads.execution, d.writes.execution]).toEqual([1, 1]);
      expect([d.reads.pointers, d.writes.pointers]).toEqual([1, 1]);
      expect(d.writes.events).toBe(1);
      expect(d.writes.heads).toBe(2);
      expect(d.writes.segments).toBeGreaterThanOrEqual(2);
      expect(d.writes.segments).toBeLessThanOrEqual(2 * FEW);
      expect(d.writes.messages).toBe(3);
      expect(d.writes.blobs).toBe(0);
      expect(d.reads.heads + d.reads.segments + d.reads.messages).toBe(0);
    });

    it('P2 append-near-constant: 2,000 single-entry appends write near-constant bytes each, in about log n segments', async () => {
      const n = 2_000;
      const s = await createSession(fixture);
      const t = await protocolTurn(s);
      let [bytes, first, copied, widest, reads] = [0, 0, 0, 0, 0];
      for (let i = 0; i < n; i += 1) {
        const { cost, moves } = await t.commit({ [MAIN]: [user(`append ${i}`)] });
        bytes += cost.bytesWritten;
        if (i === 0) first = cost.bytesWritten;
        reads += journalReads(cost);
        for (const move of moves) {
          widest = Math.max(widest, move.head.segments.length);
          copied += move.segments.reduce((sum, segment) => sum + segment.entries.length, 0);
        }
      }
      const log = Math.ceil(Math.log2(n));
      expect(bytes / n).toBeLessThanOrEqual(first + SLACK.bytes);
      expect(widest).toBeLessThanOrEqual(log + 2);
      expect(copied).toBeLessThanOrEqual(n * (log + 2));
      expect(reads).toBe(0);
    });

    it('P2 dedup-write: appending a message already stored writes no message record', async () => {
      const s = await createSession(fixture);
      const t = await protocolTurn(s);
      const { cost: d } = await t.commit({ [MAIN]: [user('hello')] });
      expect(d.writes.messages).toBe(0);
      expect(d.reads.messages).toBe(1);
      // The same through a runner: a prompt and an answer repeated word for word are stored once.
      const r = await createSession(fixture, { model: [{ text: 'same' }, { text: 'same' }] });
      await r.send('same');
      const before = r.counts();
      await r.send('same');
      expect(difference(r.counts(), before).writes.messages).toBe(0);
    });

    // -------------------------------------------------------------------------------------------
    // P3. A warm runner re-reads nothing it committed.
    // -------------------------------------------------------------------------------------------

    it('P3 serial-rounds: each step of 2 to 32 rounds reads what the first did, and no head or message it wrote', async () => {
      for (const rounds of ROUNDS) {
        const steps = between((await serialRounds(fixture, rounds)).marks);
        const first = steps[0] as StorageCounts;
        for (const step of steps) {
          expectWith(step.reads.heads + step.reads.messages, `${rounds} rounds`).toBe(0);
          expectWith(totalReads(step), `${rounds} rounds`).toBeLessThanOrEqual(totalReads(first) + SLACK.records);
          expectWith(step.bytesRead, `${rounds} rounds`).toBeLessThanOrEqual(first.bytesRead + SLACK.bytes);
        }
      }
    });

    it('P3 sequential-turns: 12 turns on one warm runner cost the same each and read no head or message', async () => {
      const { marks } = await turnsInARow(fixture, 12);
      const turns = between(marks);
      expect(turns).toHaveLength(12);
      const [first, ...rest] = turns as [StorageCounts, ...StorageCounts[]];
      for (const turn of rest) {
        expect(turn.reads.heads + turn.reads.messages).toBe(0);
        expectSame('turn', first, turn);
      }
    });

    // -------------------------------------------------------------------------------------------
    // P4. Cold work is linear in the selected head.
    // -------------------------------------------------------------------------------------------

    it('P4 selected-history: opening a turn (cold begin, and a claim that resumes) is linear in the history, warm steps are not', async () => {
      const sizes = [4, 40, 80];
      for (const mode of ['begin', 'resume'] as const) {
        const runs: { open: StorageCounts; warm: StorageCounts }[] = [];
        for (const size of sizes)
          runs.push(mode === 'begin' ? await openTurn(fixture, size) : await recover(fixture, size));
        const base = runs[0] as (typeof runs)[0];
        for (let i = 1; i < sizes.length; i += 1) {
          const [small, large] = [runs[i - 1] as (typeof runs)[0], runs[i] as (typeof runs)[0]];
          const factor = (sizes[i] as number) / (sizes[i - 1] as number);
          expectWith(totalReads(large.open), `${mode} open reads`).toBeLessThanOrEqual(
            factor * totalReads(small.open) + SLACK.records,
          );
          expectWith(large.open.bytesRead, `${mode} open bytes`).toBeLessThanOrEqual(
            factor * small.open.bytesRead + SLACK.bytes,
          );
          expectSame(`${mode} warm steps at ${sizes[i]}`, base.warm, large.warm);
        }
      }
    });

    // -------------------------------------------------------------------------------------------
    // P5. Work does not depend on retained history it does not select.
    // -------------------------------------------------------------------------------------------

    it("P5 unselected-history: with 20 past turns, a reset main and other pointers, a step costs what a fresh harness's does", async () => {
      const note: ToolRun = (context) => {
        const n = (context.arguments as { n: number }).n;
        context.journal.append(`thread-${n}`, [assistant(`note ${n}`)]);
        return { text: `noted ${n}` };
      };
      const past: ModelStep[] = Array.from({ length: 20 }, (_, n) => n).flatMap((n) => [
        { calls: [{ tool: 'note', arguments: { n } }] },
        { text: `ok ${n}` },
      ]);
      const old = await createSession(fixture, { tools: { note }, model: [...past, { text: 'fresh' }] });
      const empty = (await old.client.execution()).main;
      for (let i = 0; i < 20; i += 1) await old.send(`turn ${i}`);
      await old.client.reset({ key: 'reset-1', pointers: { main: empty } });
      expect((await old.client.pointers(['thread-0', 'thread-19'])).every((p) => p !== undefined)).toBe(true);
      const fresh = await createSession(fixture, { tools: { note }, model: [{ text: 'fresh' }], declare: [] });
      const [oldBefore, freshBefore] = [old.counts(), fresh.counts()];
      await old.send('measure');
      await fresh.send('measure');
      const [a, b] = [difference(old.counts(), oldBefore), difference(fresh.counts(), freshBefore)];
      expect(a.transactions).toBe(b.transactions);
      expect(a.reads).toEqual(b.reads);
      expect(a.writes).toEqual(b.writes);
      expect(Math.abs(a.bytesRead - b.bytesRead)).toBeLessThanOrEqual(SLACK.bytes);
      expect(Math.abs(a.bytesWritten - b.bytesWritten)).toBeLessThanOrEqual(SLACK.bytes);
    });

    it('P5 page-reads: turns, events and steps read only the records they return, through their indexes', async () => {
      const { s } = await turnsInARow(fixture, 20);
      for (const limit of [5, 20]) {
        const before = s.counts();
        expect(await s.client.turns({ limit })).toHaveLength(limit);
        const d = difference(s.counts(), before);
        expectWith(d.reads.events, `turns ${limit}`).toBe(limit);
        // Besides the records returned: the Execution (the active turn) and the tombstone check.
        expectWith(totalReads(d) - d.reads.execution - d.reads.store, `turns ${limit}`).toBe(limit);
      }
      const all = await s.events();
      let before = s.counts();
      const page = await s.client.events({ after: (all[all.length - 11] as { seq: number }).seq, limit: 500 });
      let d = difference(s.counts(), before);
      expect(page).toHaveLength(10);
      expect(d.reads.events).toBe(10);
      expect(totalReads(d)).toBe(10);
      const ids = [...new Set(all.flatMap((event) => (event.kind === 'step' ? [event.step] : [])))].filter(
        (_, i) => i % 5 === 0,
      );
      before = s.counts();
      const found = await s.client.steps(ids);
      d = difference(s.counts(), before);
      expect(found).toHaveLength(ids.length);
      expect(d.reads.events).toBe(ids.length);
      expect(totalReads(d) - d.reads.execution - d.reads.store).toBe(ids.length);
    });

    // -------------------------------------------------------------------------------------------
    // P6. Streaming stays outside transactions.
    // -------------------------------------------------------------------------------------------

    it('P6 frames: 100 frames and 1 frame give the same number of transactions', async () => {
      const measure = async (deltas: string[]) => {
        const s = await createSession(fixture, { model: [{ text: deltas.join(''), deltas }] });
        const before = s.counts();
        await s.send('go');
        return difference(s.counts(), before);
      };
      const one = await measure(['x'.repeat(100)]);
      const many = await measure(Array.from({ length: 100 }, () => 'x'));
      expect(many.transactions).toBe(one.transactions);
    });

    it('P6 frames: the same commits and writes, and only the observation count differs', async () => {
      const [one, many] = [await streamed(fixture, 1), await streamed(fixture, 100)];
      expect(many.cost.transactions).toBe(one.cost.transactions);
      expect(many.cost.writes).toEqual(one.cost.writes);
      expect((await many.s.events()).length).toBe((await one.s.events()).length);
      expect(many.runner.observations.length).toBeGreaterThan(one.runner.observations.length);
    });

    it('P6 frames: 100 tool progress frames and 1 give the same transactions and writes', async () => {
      const [one, many] = [await toolFrames(fixture, 1), await toolFrames(fixture, 100)];
      expect(many.cost.transactions).toBe(one.cost.transactions);
      expect(many.cost.writes).toEqual(one.cost.writes);
      expect((await many.s.events()).length).toBe((await one.s.events()).length);
    });

    // -------------------------------------------------------------------------------------------
    // P7. A wider tools step is still one commit.
    // -------------------------------------------------------------------------------------------

    it('P7 width: widths 1, 8 and 64 each settle in one commit; records written grow linearly, transactions do not', async () => {
      const widths = [1, 8, 64];
      const runs = [];
      for (const width of widths) runs.push(await toolWidth(fixture, width));
      for (const [i, width] of widths.entries()) {
        const run = runs[i] as (typeof runs)[0];
        const tools = stepEvents(await run.s.events()).filter((event) => event.stepKind === 'tools');
        expectWith(tools, `width ${width}`).toHaveLength(1);
        expectWith(tools[0]?.calls, `width ${width}`).toHaveLength(width);
        expectWith(tools[0]?.moves, `width ${width}`).toHaveLength(1);
        expectWith(run.cost.transactions, `width ${width}`).toBe((runs[0] as (typeof runs)[0]).cost.transactions);
      }
      const [w1, w8, w64] = runs.map((run) => totalWrites(run.cost)) as [number, number, number];
      expect(w8).toBeLessThanOrEqual(8 * w1 + SLACK.records);
      expect(w64).toBeLessThanOrEqual(8 * w8 + SLACK.records);
      // Linear, not quadratic: the cost of the 56 extra calls is at most 8 times that of the 7 before them.
      expect(w64 - w8).toBeLessThanOrEqual(8 * (w8 - w1) + SLACK.records);
    });

    // -------------------------------------------------------------------------------------------
    // P8. Nested runs add no transactions.
    // -------------------------------------------------------------------------------------------

    it("P8 nested-rounds: a nested run of 2 to 32 rounds causes no transaction of its own; the parent's commit grows linearly", async () => {
      const runs = [];
      for (const rounds of ROUNDS) runs.push(await nestedRounds(fixture, rounds));
      for (const [i, rounds] of ROUNDS.entries()) {
        const run = runs[i] as (typeof runs)[0];
        expectWith(run.span.transactions, `${rounds} rounds`).toBe(0);
        expectWith(durableWrites(run.span), `${rounds} rounds`).toBe(0);
        expectWith(run.cost.transactions, `${rounds} rounds`).toBe((runs[0] as (typeof runs)[0]).cost.transactions);
        const tools = stepEvents(await run.s.events()).filter((event) => event.stepKind === 'tools');
        expect(tools).toHaveLength(1);
        const nested = tools[0]?.calls.filter((call) => call.parentCallId !== undefined);
        expectWith(nested, `${rounds} rounds`).toHaveLength(2 * rounds + 1);
      }
      expectSeries(
        'nested',
        ROUNDS,
        runs.map((run) => run.cost),
      );
    });

    // -------------------------------------------------------------------------------------------
    // P9. Replay, renewal and recovery are cheap.
    // -------------------------------------------------------------------------------------------

    it('P9 replay: repeating any operation with its key writes nothing and calls no model or tool', async () => {
      // Protocol operations: submit, claim, begin, commit, release and create; then stop and reset.
      const p = await createSession(fixture);
      const t = await protocolTurn(p);
      const done = await t.commit({ [MAIN]: [assistant('hi')] }, true);
      const release = { seenSeq: done.committed.seq };
      const released = await p.port.release(done.committed.lease, release);
      await p.client.submit({ key: 'item-2', messages: [user('later')] });
      const stopped = await p.client.stop({ key: 'stop-1', item: 'item-2' });
      const reset = await p.client.reset({ key: 'reset-1' });
      const [modelCalls, toolCalls] = [p.model.invocations.length, p.tools.calls.length];
      const seq = (await p.client.execution()).seq;
      const before = p.counts();
      expect(await p.client.submit({ key: 'item-1', messages: [user('hello')] })).toEqual(t.submitted);
      // A replayed claim returns the original lease; its Execution is a fresh read (progress.md, phase 3).
      expect((await p.port.claim(t.claim)).lease).toEqual(t.claimed.lease);
      expect(await p.port.begin(t.claimed.lease, t.begin)).toEqual(t.begun);
      expect(await p.port.commit(t.begun.lease, done.batch)).toEqual(done.committed);
      expect(await p.port.release(done.committed.lease, release)).toEqual(released);
      expect(await p.client.stop({ key: 'stop-1', item: 'item-2' })).toEqual(stopped);
      expect(await p.client.reset({ key: 'reset-1' })).toEqual(reset);
      await p.host.create(p.harnessId, { key: `create-${p.harnessId}` });
      expect(totalWrites(difference(p.counts(), before))).toBe(0);
      expect((await p.client.execution()).seq).toBe(seq);
      expect([p.model.invocations.length, p.tools.calls.length]).toEqual([modelCalls, toolCalls]);

      // A prompt answered by a model and a tool, then its submit repeated.
      const s = await createSession(fixture, {
        tools: { work },
        declare: [],
        model: [{ calls: callsOf('work', 1) }, { text: 'done' }],
      });
      // The same input both times (the session adds declarations to its first prompt only).
      const prompt = [declareTools([s.tools.definition('work')]), user('go')];
      const first = await s.send(prompt, { key: 'k1' });
      const [models, tools, answered] = [s.model.invocations.length, s.tools.calls.length, s.counts()];
      const again = await s.send(prompt, { key: 'k1', run: false });
      expect(again.launch).toBe(first.launch);
      expect(totalWrites(difference(s.counts(), answered))).toBe(0);
      expect([s.model.invocations.length, s.tools.calls.length]).toEqual([models, tools]);
    });

    it('P9 renew: 100 renewals write only the runner state: no event, no other record, no seq change', async () => {
      const s = await createSession(fixture);
      const { launch } = await s.client.submit({ key: 'item-1', messages: [user('hi')] });
      let { lease } = await s.port.claim({ launch: need(launch), runnerId: 'r1', schemaVersion: SCHEMA_VERSION });
      const seq = (await s.client.execution()).seq;
      const before = s.counts();
      for (let i = 0; i < 100; i += 1)
        lease = await s.port.renew(lease, { progress: { phase: 'idle', since: s.clock.now() } });
      const d = difference(s.counts(), before);
      expect(d.writes.execution).toBeLessThanOrEqual(100);
      expect(totalWrites(d) - d.writes.execution).toBe(0);
      expect((await s.client.execution()).seq).toBe(seq);
    });

    it('P9 busy-no-heartbeat: a runner committing a round every second for a minute sends no renew', async () => {
      const ref: Holder = {};
      const s = await createSession(fixture, {
        host: { leaseMs: 30_000 },
        tools: {
          tick: async () => {
            await need(ref.s).clock.sleep(1_000);
            return { text: 'tick' };
          },
        },
        model: [...Array.from({ length: 60 }, () => ({ calls: callsOf('tick', 1) })), { text: 'done' }],
      });
      ref.s = s;
      const { launch } = await s.send('go', { run: false });
      const { runner, log } = observedRunner(s, { limits: { signalEveryMs: 10_000 } });
      const started = s.clock.now();
      await s.run(need(launch), runner);
      expect(s.clock.now() - started).toBeGreaterThanOrEqual(60_000);
      expect(log.filter((call) => call.op === 'renew')).toHaveLength(0);
    });

    it('P9 idle-heartbeat: a 60 s retry wait with signalEveryMs 10 s sends about six renewals with the same progress', async () => {
      const retry = with_({ failure: () => ({ retry: { afterMs: 60_000 } }) });
      const s = await createSession(fixture, {
        host: { leaseMs: 30_000 },
        model: [{ fail: { transient: true, message: 'try later' } }, { text: 'done' }],
      });
      const { launch } = await s.send('go', { run: false });
      const { runner, log } = observedRunner(s, { drive: retry, limits: { signalEveryMs: 10_000 } });
      await s.run(need(launch), runner);
      const renewals = log.flatMap((call, index) => (call.op === 'renew' ? [{ call, previous: log[index - 1] }] : []));
      expect(renewals.length).toBeGreaterThanOrEqual(5);
      expect(renewals.length).toBeLessThanOrEqual(7);
      const progress = (call: PortCall) => (call.args[1] as { progress?: unknown } | undefined)?.progress;
      expect(progress(renewals[0]?.call as PortCall)).toBeDefined();
      for (const { call, previous } of renewals) {
        expect(progress(call)).toEqual(progress(renewals[0]?.call as PortCall));
        expect(call.at - (previous as PortCall).at).toBeGreaterThanOrEqual(10_000);
      }
    });

    it('P9 progress-coalesced: 1,000 tool progress frames in 5 s cause no extra message; the next message carries the latest', async () => {
      const ref: Holder & { log?: PortCall[]; from?: number; to?: number; lastFrameAt?: number } = {};
      const s = await createSession(fixture, {
        host: { leaseMs: 30_000 },
        tools: {
          stream: async function* (): AsyncGenerator<ToolFrame> {
            const session = need(ref.s);
            ref.from = need(ref.log).length;
            for (let i = 0; i < 1_000; i += 1) {
              ref.lastFrameAt = session.clock.now();
              yield { kind: 'progress', value: { i } };
              if (i % 10 === 9) await session.clock.sleep(50);
            }
            ref.to = need(ref.log).length;
            yield { kind: 'result', status: 'succeeded', parts: [{ kind: 'text', text: 'streamed' }] };
          },
        },
        model: [{ calls: callsOf('stream', 1) }, { text: 'done' }],
      });
      ref.s = s;
      const { launch } = await s.send('go', { run: false });
      const { runner, log } = observedRunner(s, { limits: { signalEveryMs: 10_000 } });
      ref.log = log;
      await s.run(need(launch), runner);
      expect(log.slice(need(ref.from), need(ref.to))).toHaveLength(0);
      const next = log.slice(need(ref.to)).find((call) => call.op === 'commit');
      const batch = next?.args[1] as CommitBatch | undefined;
      expect(batch?.progress).toBeDefined();
      expect(batch?.progress?.since).toBeGreaterThanOrEqual(need(ref.lastFrameAt));
    });

    it('P9 recovery: a runner claiming an unowned turn calls no model or tool before re-running the pending step, and reads little', async () => {
      const r = await recover(fixture, 3);
      expect(r.modelsAtResume).toBe(r.modelsBefore);
      expect(r.toolsAtResume).toBe(r.toolsBefore + 1);
      expect(r.open.reads.execution).toBeLessThanOrEqual(2);
      expect(r.open.reads.heads).toBeLessThanOrEqual(1);
      expect(r.open.reads.segments).toBeLessThanOrEqual(8);
      expect(r.open.reads.messages).toBeLessThanOrEqual(12);
      expect(r.open.reads.blobs).toBeLessThanOrEqual(2);
      expect(r.open.reads.events).toBeLessThanOrEqual(1);
      expect(r.open.reads.pointers).toBe(0);
      expect(r.open.writes.messages + r.open.writes.heads + r.open.writes.segments + r.open.writes.pointers).toBe(0);
    });

    // -------------------------------------------------------------------------------------------
    // P10. A runner retains only its selected context and the open step.
    // -------------------------------------------------------------------------------------------

    it('P10 cache-no-invalidation: a warm runner reads no message, head, segment or blob, and the Execution only when seq moved', async () => {
      const ref: Holder & { log?: PortCall[]; from?: number; mark?: StorageCounts; outside?: boolean } = {};
      const s = await createSession(fixture, {
        tools: {
          work: async (context) => {
            const session = need(ref.s);
            ref.from ??= need(ref.log).length;
            ref.mark ??= session.counts();
            // One change from outside the runner, mid-turn.
            if (ref.outside !== true && (context.arguments as { i: number }).i === 4) {
              ref.outside = true;
              await session.client.submit({ key: 'outside', messages: [user('from outside')] });
            }
            return { text: 'worked' };
          },
        },
        model: [
          ...Array.from({ length: 8 }, (_, i) => ({ calls: [{ tool: 'work', arguments: { i } }] })),
          { text: 'done' },
          { text: 'answered' },
        ],
      });
      ref.s = s;
      const { launch } = await s.send('go', { run: false });
      const { runner, log } = observedRunner(s);
      ref.log = log;
      await s.run(need(launch), runner);
      const warm = log.slice(need(ref.from));
      const calls = (...ops: string[]) => warm.filter((call) => ops.includes(call.op)).length;
      expect(calls('heads', 'entries', 'blobs', 'pointers')).toBe(0);
      expect(calls('execution')).toBeLessThanOrEqual(1);
      const served = difference(s.counts(), need(ref.mark));
      expect(journalReads(served) + served.reads.blobs).toBe(0);
    });

    it('P10 evict-on-edit: a compaction edit removing 500 entries leaves the runner holding none of them, with no reads', async () => {
      const ref: Holder = {};
      const marks: StorageCounts[] = [];
      const held: RetentionSnapshot[] = [];
      const s = await createSession(fixture, { model: [{ text: 'one' }, { text: 'two' }] });
      ref.s = s;
      // The second turn begins with an edit keeping a summary and the last 100 entries of the 601 before it.
      let turns = 0;
      const runner = s.runner({
        drive: with_({
          boundary: () => 'proceed',
          'next-turn': (_decision, context) => {
            marks.push(s.counts());
            held.push(retentionOf(need(ref.runner)));
            const item = context.inbox[0];
            if (item === undefined) return 'release';
            turns += 1;
            if (turns < 2) return { begin: { keys: [item.key] } };
            const summary = user('summary of the earlier entries');
            const compact = { entries: [summary, ...context.head.entries.slice(-100)] };
            return { begin: { keys: [item.key], compact } };
          },
        }),
      });
      ref.runner = runner;
      const launch = await queue(s, [Array.from({ length: 600 }, (_, i) => user(`entry ${i}`)), 'small']);
      await s.run(launch, runner);
      const [, beforeEdit, afterEdit] = held as [RetentionSnapshot, RetentionSnapshot, RetentionSnapshot];
      expect(beforeEdit.entries).toBeGreaterThanOrEqual(600);
      // After: the summary, 100 kept, the second prompt and its answer.
      expect(afterEdit.entries).toBeLessThanOrEqual(103 + SLACK.records);
      expect(afterEdit.segments).toBeLessThanOrEqual(Math.ceil(Math.log2(103)) + SLACK.records);
      const turn = difference(marks[2] as StorageCounts, marks[1] as StorageCounts);
      expect(journalReads(turn) + turn.reads.blobs).toBe(0);
    });

    it('P10 spill-bounded: evicted items never exceed cacheBytes; a reopened pointer is served from it, or read cold with 0', async () => {
      const reopen = async (cacheBytes: number) => {
        const ref: Holder = {};
        const marks: StorageCounts[] = [];
        const spills: number[] = [];
        const s = await createSession(fixture, {
          declare: ['investigate'],
          tools: {
            investigate: async (context) => {
              const result = await context.run({
                input: [user('look closely')],
                tools: { enabled: [], preauthorized: [], definitions: [] },
                model: TURN.model,
                pointer: 'inv',
                limits: {},
              });
              return { text: result.status };
            },
          },
          model: [
            { calls: callsOf('investigate', 1) },
            { text: 'found 1' },
            { text: 'done 1' },
            { calls: callsOf('investigate', 1) },
            { text: 'found 2' },
            { text: 'done 2' },
          ],
          drive: oneAtATime(() => {
            marks.push(need(ref.s).counts());
            spills.push(retentionOf(need(ref.runner)).spillBytes);
          }),
        });
        ref.s = s;
        const runner = s.runner({ limits: { cacheBytes } });
        ref.runner = runner;
        await s.run(await queue(s, ['first', 'second']), runner);
        const second = difference(marks[2] as StorageCounts, marks[1] as StorageCounts);
        return { spills, reads: journalReads(second) };
      };
      const [roomy, none, tiny] = [await reopen(1 << 20), await reopen(0), await reopen(64)];
      for (const spill of roomy.spills) expect(spill).toBeLessThanOrEqual(1 << 20);
      for (const spill of none.spills) expect(spill).toBe(0);
      for (const spill of tiny.spills) expect(spill).toBeLessThanOrEqual(64);
      expect(roomy.reads).toBe(0);
      expect(none.reads).toBeGreaterThan(0);
    });

    it('P10 retention: over 2 to 32 rounds the runner holds the selected head plus a constant, and a finished nested run holds only its records', async () => {
      for (const rounds of ROUNDS) {
        const { retained } = await serialRounds(fixture, rounds, { inspect: true });
        expect(retained).toHaveLength(rounds);
        for (const [k, held] of retained.entries()) {
          // The head holds the declaration, the prompt and two messages per finished round; the open step adds the call.
          const selected = 3 + 2 * k;
          expectWith(held.entries, `${rounds} rounds, round ${k}`).toBeLessThanOrEqual(selected + SLACK.records);
          expectWith(held.segments, `${rounds} rounds, round ${k}`).toBeLessThanOrEqual(
            Math.ceil(Math.log2(selected + 1)) + SLACK.records,
          );
          if (k > 0)
            expectWith(held.bytes - (retained[k - 1] as RetentionSnapshot).bytes, `round ${k}`).toBeLessThanOrEqual(
              SLACK.bytes,
            );
        }
      }
      // A nested run of 32 rounds finishes in round 1; a later round's tool finds none of its entries retained.
      const ref: Holder & { held: RetentionSnapshot[] } = { held: [] };
      const s = await createSession(fixture, {
        declare: ['investigate', 'probe'],
        tools: {
          inner: work,
          probe: () => {
            ref.held.push(retentionOf(need(ref.runner)));
            return { text: 'probed' };
          },
          investigate: async (context) => {
            ref.held.push(retentionOf(need(ref.runner)));
            const result = await context.run({
              input: [user('look closely')],
              tools: {
                enabled: ['inner'],
                preauthorized: ['inner'],
                definitions: [need(ref.s).tools.definition('inner')],
              },
              model: TURN.model,
              limits: {},
            });
            return { text: result.status };
          },
        },
        model: [
          { calls: callsOf('investigate', 1) },
          ...Array.from({ length: 32 }, () => ({ calls: callsOf('inner', 1) })),
          { text: 'inner done' },
          { calls: callsOf('probe', 1) },
          { text: 'final' },
        ],
      });
      ref.s = s;
      ref.runner = s.runner();
      await s.send('go', { runner: ref.runner });
      const [start, later] = ref.held as [RetentionSnapshot, RetentionSnapshot];
      // The head then holds: the declaration, the prompt, the first call and its result, and the open step adds one call.
      expect(later.entries).toBeLessThanOrEqual(start.entries + 2 + SLACK.records);
    });

    // -------------------------------------------------------------------------------------------
    // P11. Hashing is linear and done once.
    // -------------------------------------------------------------------------------------------

    it('P11 hash-once: no round hashes more than the first, however long the head', async () => {
      await withDigestSpy(async (digests) => {
        const calls: number[] = [];
        await serialRounds(fixture, 16, { onTool: () => calls.push(digests.calls()) });
        const rounds = calls.slice(1).map((count, i) => count - (calls[i] as number));
        expect(digests.calls()).toBeGreaterThan(0);
        for (const count of rounds) expect(count).toBeLessThanOrEqual((rounds[0] as number) + 2);
      });
    });

    it('P11 hash-once: a turn hashes the same after 2 turns of history and after 40', async () => {
      const cost = async (history: number) => {
        const s = await seeded(fixture, history, [{ text: 'measure' }]);
        const { launch } = await s.send('measure', { run: false });
        return await withDigestSpy(async (digests) => {
          await s.run(need(launch));
          return digests.calls();
        });
      };
      const [short, long] = [await cost(2), await cost(40)];
      expect(short).toBeGreaterThan(0);
      expect(long).toBeLessThanOrEqual(short + 2);
    });

    it('P11 hash-once: hashing a new message is linear in its size', async () => {
      const bytes = async (size: number) => {
        const s = await createSession(fixture, { model: [{ text: 'ok' }] });
        const { launch } = await s.send(`${'x'.repeat(size)}`, { run: false });
        return await withDigestSpy(async (digests) => {
          await s.run(need(launch));
          return digests.bytes();
        });
      };
      const sizes = [64 * 1024, 128 * 1024, 256 * 1024];
      const hashed: number[] = [];
      for (const size of sizes) hashed.push(await bytes(size));
      expect(hashed[0]).toBeGreaterThan(0);
      for (let i = 1; i < hashed.length; i += 1)
        expect(hashed[i] as number).toBeLessThanOrEqual(2 * (hashed[i - 1] as number) + SLACK.bytes);
    });
  });
}
