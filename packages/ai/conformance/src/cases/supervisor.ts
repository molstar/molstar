/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  createRunner,
  defaultDrive,
  type Execution,
  type LaunchKey,
  type Progress,
  type Runner,
  type RunnerPort,
  type RunnerState,
  SCHEMA_VERSION,
  type ToolFrame,
} from '@molstar/ai-harness';
import * as fc from 'fast-check';
import { describe, expect, it } from '@jest/globals';
import { type ConformanceFixture, createSession, type Session } from '../session.js';
import { settle } from '../testing/clock.js';
import { simulatedIds } from '../testing/ids.js';
import { user } from '../testing/messages.js';
import { simulatedSupervisor } from '../testing/supervisor.js';
import { controls, delegate, expectWith, TURN } from './helpers.js';

/** Known lease and launch windows, so the cases can expire them (harness.md §4). */
const WINDOWS = { leaseMs: 30_000, launchMs: 10_000 };

/** A runner over `port` with the session's scripted model and tools and no tools enabled in its turns. */
const runnerOver = (s: Session, port: RunnerPort, runnerId: string, handOffAfterMs?: number): Runner =>
  createRunner({
    port,
    runnerId,
    models: s.model.bindings,
    tools: s.tools.bindings,
    drive: defaultDrive(),
    configure: () => TURN,
    clock: s.clock,
    ids: simulatedIds(`${runnerId}-`),
    ...(handOffAfterMs !== undefined ? { limits: { handOffAfterMs } } : {}),
  });

/** What one runner did through its port, in the order of a shared counter (S1). */
interface Hold {
  claimDone?: number;
  lastOkStart: number;
  leaseUntil: number;
  released: boolean;
  /** Set when this runner's claim succeeded while another runner still held a live lease. */
  violation?: string;
}

/** Records, per runner, when its successful operations started and what its lease was. */
function traced(port: RunnerPort, hold: Hold, holds: readonly Hold[], s: Session, next: () => number): RunnerPort {
  const ok = (start: number) => {
    hold.lastOkStart = start;
  };
  return delegate(port, {
    async claim(input) {
      const before = holds
        .filter((h) => h !== hold && h.claimDone !== undefined)
        .map((h) => ({ h, until: h.leaseUntil }));
      const result = await port.claim(input);
      hold.claimDone = next();
      hold.leaseUntil = result.lease.until;
      for (const { h, until } of before)
        if (!h.released && until > s.clock.now()) hold.violation = `claimed while a lease ran until ${until}`;
      return result;
    },
    async putBlobs(lease, blobs) {
      const start = next();
      await port.putBlobs(lease, blobs);
      ok(start);
    },
    async begin(lease, input) {
      const start = next();
      const result = await port.begin(lease, input);
      ok(start);
      hold.leaseUntil = result.lease.until;
      return result;
    },
    async commit(lease, batch) {
      const start = next();
      const result = await port.commit(lease, batch);
      ok(start);
      hold.leaseUntil = result.lease.until;
      return result;
    },
    async renew(lease, input) {
      const start = next();
      const result = await port.renew(lease, input);
      ok(start);
      hold.leaseUntil = result.until;
      return result;
    },
    async release(lease, input) {
      const start = next();
      const result = await port.release(lease, input);
      ok(start);
      hold.released = true;
      return result;
    },
  });
}

type Action =
  | { kind: 'submit' }
  | { kind: 'withdraw'; pick: number }
  | { kind: 'start'; mode: 'normal' | 'fail' | 'handoff' }
  | { kind: 'crash'; pick: number }
  | { kind: 'advance'; ms: number }
  | { kind: 'check' };

const actions: fc.Arbitrary<Action> = fc.oneof(
  { weight: 4, arbitrary: fc.constant<Action>({ kind: 'submit' }) },
  { weight: 1, arbitrary: fc.record({ kind: fc.constant('withdraw' as const), pick: fc.nat(4) }) },
  {
    weight: 5,
    arbitrary: fc.record({
      kind: fc.constant('start' as const),
      mode: fc.constantFrom('normal' as const, 'normal' as const, 'fail' as const, 'handoff' as const),
    }),
  },
  { weight: 2, arbitrary: fc.record({ kind: fc.constant('crash' as const), pick: fc.nat(4) }) },
  {
    weight: 4,
    arbitrary: fc.record({
      kind: fc.constant('advance' as const),
      ms: fc.constantFrom(1_000, 6_000, 12_000, 31_000, 61_000),
    }),
  },
  { weight: 2, arbitrary: fc.constant<Action>({ kind: 'check' }) },
);

/** Work exists (harness.md §4): a turn active and neither parked nor blocked, or no turn and an item waiting. */
const hasWork = (e: Execution): boolean =>
  e.turn !== undefined ? e.turn.park === undefined && e.turn.blocked === undefined : e.inbox.length > 0;

const live = (state: RunnerState): Extract<RunnerState, { phase: 'live' }> => {
  expect(state.phase).toBe('live');
  return state as Extract<RunnerState, { phase: 'live' }>;
};

/** One random sequence of S1's actions against one session; throws on the first broken invariant. */
async function liveness(fixture: ConformanceFixture, sequence: readonly Action[]): Promise<void> {
  const s = await createSession(fixture, {
    host: WINDOWS,
    model: Array.from({ length: 300 }, () => ({ text: 'ok', waitMs: 5_000 })),
  });
  const holds: Hold[] = [];
  const running: { runner: Runner; done: boolean }[] = [];
  const pending: LaunchKey[] = [];
  const known = new Set<LaunchKey>();
  let tick = 0;
  let items = 0;
  let mode: 'normal' | 'fail' | 'handoff' = 'normal';
  const next = () => ++tick;
  const learn = (launch: LaunchKey | undefined) => {
    if (launch === undefined || known.has(launch)) return;
    known.add(launch);
    pending.push(launch);
  };
  const supervisor = simulatedSupervisor({
    failLaunch: () => mode === 'fail',
    start: (launch) => {
      const hold: Hold = { lastOkStart: 0, leaseUntil: 0, released: false };
      holds.push(hold);
      const runner = runnerOver(
        s,
        traced(s.port, hold, holds, s, next),
        `runner-${holds.length}`,
        mode === 'handoff' ? 0 : undefined,
      );
      const entry = { runner, done: false };
      running.push(entry);
      return runner.run(launch).finally(() => {
        entry.done = true;
      });
    },
  });
  const startNext = (as: typeof mode) => {
    const key = pending.shift();
    if (key === undefined) return;
    mode = as;
    supervisor.ensure(s.harnessId, key);
    mode = 'normal';
  };
  /** Reads the Execution; a `launching` state's key is a launch the supervisor has not started yet. */
  const observe = async () => {
    const execution = await s.client.execution();
    if (execution.runner.phase === 'launching') learn(execution.runner.key);
    return execution;
  };

  const invariants = async () => {
    const execution = await observe();
    if (hasWork(execution)) {
      // While work exists a runner is on it or launching (harness.md §4)...
      expectWith(execution.runner.phase, 'work exists but no runner state').not.toBe('none');
      const until =
        execution.runner.phase === 'launching'
          ? execution.runner.until
          : execution.runner.phase === 'live'
            ? execution.runner.leaseUntil
            : Number.POSITIVE_INFINITY;
      // ...or, once the launch or the lease has expired, the next check launches (a failed launch or a crash).
      if (until < s.clock.now()) {
        const { launch } = await s.client.check();
        expectWith(launch, 'work exists, the runner state expired, and check launched nothing').toBeDefined();
        learn(launch);
      }
    }
    // Never two live runners: a claim never overlaps a live lease, and a runner makes no successful
    // call that starts after a later claimant's claim completed.
    expect(holds.flatMap((h) => (h.violation === undefined ? [] : [h.violation]))).toEqual([]);
    for (const a of holds)
      for (const b of holds)
        if (a !== b && a.claimDone !== undefined && b.claimDone !== undefined && a.claimDone < b.claimDone)
          expectWith(a.lastOkStart, 'a superseded runner still acted').toBeLessThan(b.claimDone);
  };

  for (const action of sequence) {
    switch (action.kind) {
      case 'submit': {
        const { launch } = await s.client.submit({ key: `item-${++items}`, messages: [user(`message ${items}`)] });
        learn(launch);
        break;
      }
      case 'withdraw': {
        const { inbox } = await s.client.execution();
        await s.client.withdraw({ key: inbox[action.pick % Math.max(1, inbox.length)]?.key ?? 'item-none' });
        break;
      }
      case 'start':
        startNext(action.mode);
        break;
      case 'crash': {
        // A crash: the runner stops with no commit and no release; its lease just runs out.
        const alive = running.filter((r) => !r.done);
        alive[action.pick % Math.max(1, alive.length)]?.runner.abort('crash');
        break;
      }
      case 'advance':
        await s.clock.advance(action.ms);
        break;
      case 'check':
        learn((await s.client.check()).launch);
        break;
    }
    await settle(16);
    await invariants();
  }

  // Convergence: with no more failures, every launch that was lost and every crash is recovered by a
  // relaunch, and the work is done.
  for (let round = 0; round < 20; round += 1) {
    await s.clock.advance(61_000);
    await settle(16);
    if (!hasWork(await observe())) break;
    learn((await s.client.check()).launch);
    while (pending.length > 0) startNext('normal');
    await settle(16);
    await invariants();
  }
  const end = await observe();
  expectWith(hasWork(end), 'work remained after every launch and lease had expired').toBe(false);
  expect(end.inbox).toEqual([]);
  expect(end.turn).toBeUndefined();
  await invariants();
}

/** S1–S3 (Supervisor and leases). */
export function supervisorCases(fixture: ConformanceFixture): void {
  describe('Supervisor and leases', () => {
    // S1 (harness.md §4)
    it('S1 liveness: while work exists a runner is on it or launching, never two at once, and a failed launch or crash is relaunched', async () => {
      await fc.assert(
        fc.asyncProperty(fc.array(actions, { minLength: 1, maxLength: 20 }), (sequence) => liveness(fixture, sequence)),
        { numRuns: 50 },
      );
    }, 120_000);

    // S2 (harness.md §4, Claim and Release)
    it('S2 stale-claim: a duplicate launch claims after another, gets stale, and exits', async () => {
      const s = await createSession(fixture);
      const { launch } = await s.client.submit({ key: 'item-1', messages: [user('hello')] });
      await s.port.claim({ launch: launch as string, runnerId: 'r1', schemaVersion: SCHEMA_VERSION });
      await expect(
        s.port.claim({ launch: launch as string, runnerId: 'r2', schemaVersion: SCHEMA_VERSION }),
      ).rejects.toMatchObject({ code: 'stale' });
      expect((await s.run(launch as string)).reason).toBe('stale-launch');
      expect((await s.client.execution()).runner).toMatchObject({ phase: 'live', runnerId: 'r1' });
    });

    it("S2 release-race: an item submitted after the runner's last read makes its release fail inbox-changed", async () => {
      const s = await createSession(fixture);
      const { launch } = await s.client.submit({ key: 'item-1', messages: [user('hello')] });
      await s.client.withdraw({ key: 'item-1' });
      const { lease, execution } = await s.port.claim({
        launch: launch as string,
        runnerId: 'r1',
        schemaVersion: SCHEMA_VERSION,
      });
      await s.client.submit({ key: 'item-2', messages: [user('again')] });
      await expect(s.port.release(lease, { seenSeq: execution.seq })).rejects.toMatchObject({
        code: 'inbox-changed',
      });
    });

    // S3 (harness.md §4, Faults are not losses; §4.3, stall policy)
    it("S3 watch-fault: a watch error and one failed read do not end the runner's hold (defect 3)", async () => {
      const s = await createSession(fixture, { model: [{ text: 'ok' }] });
      let claimed = false;
      let readFailed = false;
      let watchFailed = false;
      const once =
        <A extends unknown[], T>(read: (...args: A) => Promise<T>) =>
        async (...args: A): Promise<T> => {
          if (claimed && !readFailed) {
            readFailed = true;
            throw new Error('read failed');
          }
          return await read(...args);
        };
      const port = delegate(s.port, {
        async claim(input) {
          const result = await s.port.claim(input);
          claimed = true;
          return result;
        },
        execution: once(() => s.port.execution()),
        heads: once((ids) => s.port.heads(ids)),
        pointers: once((scopes) => s.port.pointers(scopes)),
        entries: once((ids) => s.port.entries(ids)),
        events: once((input) => s.port.events(input)),
        steps: once((ids) => s.port.steps(ids)),
        turns: once((input) => s.port.turns(input)),
        blobs: once((refs) => s.port.blobs(refs)),
        watch: () => ({
          subscribe(onNext, onError, onComplete) {
            if (watchFailed) return s.port.watch().subscribe(onNext, onError, onComplete);
            watchFailed = true;
            queueMicrotask(() => onError?.(new Error('watch failed')));
            return () => undefined;
          },
        }),
      });
      const { launch } = await s.send('go', { run: false });
      const exit = await s.run(launch as string, runnerOver(s, port, 'runner-1'));
      expect(exit.reason).toBe('released');
      const events = await s.events();
      expect(controls(events, 'claimed')).toHaveLength(1);
      expect(controls(events, 'launched')).toHaveLength(1);
      expect(controls(events, 'stopped')).toHaveLength(0);
      expect((await s.turns())[0]?.status).toBe('completed');
      expect(s.model.invocations).toHaveLength(1);
    });

    it('S3 stalled: a hung tool keeps its lease live while progress.since stands still', async () => {
      let clock: Session['clock'] | undefined;
      let started = false;
      const s = await createSession(fixture, {
        host: WINDOWS,
        tools: {
          hang: async (context) => {
            started = true;
            await clock?.sleep(600_000, context.signal);
            return { text: 'late' };
          },
        },
        model: [{ calls: [{ tool: 'hang', arguments: {} }] }, { text: 'final' }],
      });
      clock = s.clock;
      const { launch } = await s.send('go', { run: false });
      const finished = s.runner().run(launch as string);
      // Waits for the tool itself, not a fixed number of rounds: on a slow backend (fake IndexedDB) the runner may
      // still be answering permission after 40.
      for (let i = 0; i < 500 && !started; i += 1) await settle();
      await s.clock.advance(25_000);
      const first = live((await s.client.execution()).runner);
      await s.clock.advance(60_000);
      const second = live((await s.client.execution()).runner);
      expect(first.progress?.phase).toBe('tools');
      expect(second.progress?.since).toBe((first.progress as Progress).since);
      expect(second.leaseUntil).toBeGreaterThan(first.leaseUntil);
      expect(second.leaseUntil).toBeGreaterThan(s.clock.now());
      await s.clock.advance(600_000);
      await s.clock.runUntil(finished);
    });

    it('S3 stalled: a tool emitting progress frames moves progress.since', async () => {
      let clock: Session['clock'] | undefined;
      let started = false;
      async function* ticking(signal: AbortSignal): AsyncGenerator<ToolFrame> {
        started = true;
        for (let i = 0; i < 24; i += 1) {
          await clock?.sleep(5_000, signal);
          yield { kind: 'progress', value: i };
        }
        yield { kind: 'result', status: 'succeeded', parts: [{ kind: 'text', text: 'done' }] };
      }
      const s = await createSession(fixture, {
        host: WINDOWS,
        tools: { work: (context) => ticking(context.signal) },
        model: [{ calls: [{ tool: 'work', arguments: {} }] }, { text: 'final' }],
      });
      clock = s.clock;
      const { launch } = await s.send('go', { run: false });
      const finished = s.runner().run(launch as string);
      for (let i = 0; i < 500 && !started; i += 1) await settle();
      await s.clock.advance(25_000);
      const first = live((await s.client.execution()).runner);
      await s.clock.advance(30_000);
      const second = live((await s.client.execution()).runner);
      expect(second.progress?.since).toBeGreaterThan((first.progress as Progress).since);
      await s.clock.runUntil(finished);
    });

    it.todo(
      "S3 stalled: onStalled is called once after stalledAfterMs('tools') and never for a progressing tool; the policy hooks belong to the supervisor, which the simulated supervisor does not implement and the suite may not import (createInProcessSupervisor is in ai-components)",
    );
  });
}
