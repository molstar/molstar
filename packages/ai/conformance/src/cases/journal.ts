/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashMessage, type Message } from '@molstar/ai-agent';
import {
  buildMove,
  createHarnessHost,
  type EntryId,
  type HeadRecord,
  type Lease,
  MAIN,
  type Move,
  SCHEMA_VERSION,
  type Segment,
  type StorageTx,
} from '@molstar/ai-harness';
import { describe, expect, it } from '@jest/globals';
import { difference, totalWrites } from '../meter.js';
import { type ConformanceFixture, createSession, type Session } from '../session.js';
import { settle } from '../testing/clock.js';
import { assistant, textOf, toolDefinition, user } from '../testing/messages.js';
import { appendMove, stepEvents, TURN, with_ } from './helpers.js';

/** Polls (letting the runner and the simulated clock's promises settle) until the condition holds. */
async function until(condition: () => boolean, rounds = 200): Promise<void> {
  for (let i = 0; i < rounds; i += 1) {
    if (condition()) return;
    await settle();
  }
  throw new Error('until: the condition never held');
}

/** "About log n" segments per head (harness.md §14), with room for the constant. */
const segmentBound = (entries: number): number => 2 * Math.ceil(Math.log2(Math.max(entries, 2))) + 2;

/**
 * A turn begun by hand, then settled step by step the way a runner commits: each commit builds its
 * move with the shared head scheme (`buildMove`) from a local segment cache and moves `main` once.
 */
interface Chain {
  readonly entries: EntryId[];
  head(): HeadRecord;
  /** One settled step that appends `messages` to `main`. */
  append(messages: readonly Message[]): Promise<void>;
  /** One settled step that replaces `main`'s head with `kept` followed by `added`. */
  edit(kept: readonly EntryId[], added: readonly Message[]): Promise<void>;
  /** One settled step with arbitrary moves (non-`main` pointers, for example). */
  commit(messages: readonly Message[], moves: readonly Move[]): Promise<void>;
  /** A move that creates `pointer` holding `messages`. */
  create(pointer: string, messages: readonly Message[]): Promise<Move>;
  /** A move that appends `messages` to `pointer`, from the head `created` made. */
  extend(pointer: string, created: Move, messages: readonly Message[]): Promise<Move>;
}

async function openChain(s: Session): Promise<Chain> {
  const { launch } = await s.client.submit({ key: 'item-1', messages: [user('hello')] });
  const claimed = await s.port.claim({ launch: launch as string, runnerId: 'r1', schemaVersion: SCHEMA_VERSION });
  const segments = new Map<string, Segment>();
  const remember = (move: Move): Move => {
    for (const segment of move.segments) segments.set(segment.id, segment);
    return move;
  };
  const read = (ids: readonly string[]): Segment[] => ids.map((id) => segments.get(id) as Segment);
  const prompt = user('hello');
  const promptHash = await hashMessage(prompt);
  const first = remember(await buildMove({ pointer: MAIN, edit: { entries: [promptHash] }, read }));
  const begun = await s.port.begin(claimed.lease, {
    key: 'begin-1',
    turnId: 'turn-x',
    step: 'step-1',
    items: ['item-1'],
    configuration: TURN,
    messages: [prompt],
    move: { ...first, from: claimed.execution.main },
    calls: [],
    next: { id: 'step-2', kind: 'model' },
  });
  let lease: Lease = begun.lease;
  let head = first.head;
  let entries: EntryId[] = [promptHash];
  let step = 2;

  const commit = async (messages: readonly Message[], moves: readonly Move[]): Promise<void> => {
    const committed = await s.port.commit(lease, {
      key: `commit-${step}`,
      turnId: 'turn-x',
      step: `step-${step}`,
      ending: 'settled',
      calls: [],
      messages,
      moves,
      consumed: [],
      next: { id: `step-${step + 1}`, kind: 'model' },
    });
    lease = committed.lease;
    step += 1;
  };
  const moveMain = async (messages: readonly Message[], move: Move, next: EntryId[]): Promise<void> => {
    await commit(messages, [{ ...remember(move), from: head.id }]);
    head = move.head;
    entries = next;
  };

  return {
    get entries() {
      return entries;
    },
    head: () => head,
    async append(messages) {
      const added = await Promise.all(messages.map(hashMessage));
      const move = await buildMove({ pointer: MAIN, from: head, edit: { append: added }, read });
      await moveMain(messages, move, [...entries, ...added]);
    },
    async edit(kept, added) {
      const next = [...kept, ...(await Promise.all(added.map(hashMessage)))];
      await moveMain(added, await buildMove({ pointer: MAIN, from: head, edit: { entries: next }, read }), next);
    },
    commit,
    async create(pointer, messages) {
      const hashes = await Promise.all(messages.map(hashMessage));
      return remember(await buildMove({ pointer, edit: { entries: hashes }, read }));
    },
    async extend(pointer, created, messages) {
      const hashes = await Promise.all(messages.map(hashMessage));
      const move = await buildMove({ pointer, from: created.head, edit: { append: hashes }, read });
      return { ...remember(move), from: created.head.id };
    },
  };
}

/** Every way to change one part of a value: a leaf, a dropped or added element, a reordering. */
function variants(value: unknown): unknown[] {
  if (typeof value === 'string') return [`${value}!`];
  if (typeof value === 'number') return [value + 1];
  if (typeof value === 'boolean') return [!value];
  if (value === null) return ['not null'];
  if (Array.isArray(value)) {
    const reversed = [...value].reverse();
    return [
      ...value.flatMap((item, i) => variants(item).map((changed) => value.with(i, changed))),
      value.slice(0, -1),
      [...value, 'extra'],
      ...(JSON.stringify(reversed) === JSON.stringify(value) ? [] : [reversed]),
    ];
  }
  const record = value as Record<string, unknown>;
  return [
    ...Object.entries(record).flatMap(([key, item]) =>
      variants(item).map((changed) => ({ ...record, [key]: changed })),
    ),
    ...Object.keys(record).map((key) => Object.fromEntries(Object.entries(record).filter(([k]) => k !== key))),
    { ...record, extra: 1 },
  ];
}

/** The same value with every object's keys in reverse order. */
const reorder = (value: unknown): unknown =>
  Array.isArray(value)
    ? value.map(reorder)
    : value !== null && typeof value === 'object'
      ? Object.fromEntries(
          Object.entries(value)
            .reverse()
            .map(([key, item]) => [key, reorder(item)]),
        )
      : value;

/** J1–J4 (Journal). */
export function journalCases(fixture: ConformanceFixture): void {
  describe('Journal', () => {
    it('J1 rewrite-refused: a backend write that changes a stored entry or head fails integrity', async () => {
      const storage = await fixture.storage();
      await createHarnessHost(storage).create('harness-j1', { key: 'create-j1' });
      const write = (run: (tx: StorageTx) => void) => storage.transact('harness-j1', async (tx) => run(tx));
      const original = user('original');
      const id = await hashMessage(original);
      const head: HeadRecord = { id: 'head-j1', segments: ['segment-j1'], length: 1 };
      const segment: Segment = { id: 'segment-j1', entries: [id] };

      await write((tx) => {
        tx.messages.put(id, original);
        tx.heads.put(head.id, head);
        tx.segments.put(segment.id, segment);
      });
      await expect(write((tx) => tx.messages.put(id, user('changed')))).rejects.toMatchObject({ code: 'integrity' });
      await expect(write((tx) => tx.heads.put(head.id, { ...head, length: 2 }))).rejects.toMatchObject({
        code: 'integrity',
      });
      await expect(write((tx) => tx.segments.put(segment.id, { ...segment, entries: [id, id] }))).rejects.toMatchObject(
        { code: 'integrity' },
      );

      const stored = await storage.transact('harness-j1', async (tx) => ({
        message: (await tx.messages.get([id]))[0],
        head: (await tx.heads.get([head.id]))[0],
        segment: (await tx.segments.get([segment.id]))[0],
      }));
      expect(stored).toEqual({ message: original, head, segment });
    });

    it('J1 reset-keeps-everything: after a reset, every earlier entry, head and event still reads', async () => {
      const s = await createSession(fixture, { model: [{ text: 'one' }, { text: 'two' }] });
      await s.send('first');
      const first = (await s.client.execution()).main;
      await s.send('second');
      const second = (await s.client.execution()).main;
      expect(second).not.toBe(first);
      const before = { events: await s.events(), main: await s.main() };

      await s.client.reset({ key: 'reset-1', pointers: { main: first } });

      expect((await s.client.execution()).main).toBe(first);
      const [oldHead, newHead] = await s.client.heads([first, second]);
      expect(oldHead?.entries.length).toBeGreaterThan(0);
      expect(newHead?.entries.slice(0, oldHead?.entries.length)).toEqual(oldHead?.entries);
      expect(await s.client.entries(newHead?.entries ?? [])).toEqual(before.main);
      const after = await s.events();
      expect(after.slice(0, before.events.length)).toEqual(before.events);
      expect(after.length).toBeGreaterThan(before.events.length);
    });

    it('J2 stale-move: a move from an older head fails head-changed and writes nothing', async () => {
      const s = await createSession(fixture);
      const { launch } = await s.client.submit({ key: 'item-1', messages: [user('hello')] });
      const { lease, execution } = await s.port.claim({
        launch: launch as string,
        runnerId: 'r1',
        schemaVersion: SCHEMA_VERSION,
      });
      const begin = await appendMove(s.port, MAIN, execution.main, [user('hello')]);
      const begun = await s.port.begin(lease, {
        key: 'begin-1',
        turnId: 'turn-x',
        step: 'step-1',
        items: ['item-1'],
        configuration: TURN,
        messages: [user('hello')],
        move: begin,
        calls: [],
        next: { id: 'step-2', kind: 'model' },
      });
      const stale = await appendMove(s.port, MAIN, execution.main, [assistant('hi')]);
      const before = s.counts();
      await expect(
        s.port.commit(begun.lease, {
          key: 'commit-1',
          turnId: 'turn-x',
          step: 'step-2',
          ending: 'settled',
          calls: [],
          messages: [assistant('hi')],
          moves: [stale],
          consumed: [],
          end: { status: 'completed' },
        }),
      ).rejects.toMatchObject({ code: 'head-changed' });
      expect(totalWrites(difference(s.counts(), before))).toBe(0);
    });

    it('J2 two-moves-one-pointer: a batch with two moves of one pointer fails invalid', async () => {
      const s = await createSession(fixture);
      const chain = await openChain(s);
      const one = [assistant('one')];
      const two = [assistant('two')];
      const first = await appendMove(s.port, MAIN, chain.head().id, one);
      const second = await appendMove(s.port, MAIN, chain.head().id, two);
      await expect(chain.commit([...one, ...two], [first, second])).rejects.toMatchObject({ code: 'invalid' });
    });

    it('J2 many-pointers: with 1,000 pointers, a commit moving one reads and writes only that pointer and the Execution', async () => {
      const s = await createSession(fixture);
      const chain = await openChain(s);
      const created: Move[] = [];
      const messages: Message[] = [];
      for (let i = 0; i < 1_000; i += 1) {
        messages.push(user(`pointer ${i}`));
        created.push(await chain.create(`pointer-${i}`, [messages[i] as Message]));
      }
      await chain.commit(messages, created);

      const more = assistant('one more entry');
      const move = await chain.extend('pointer-0', created[0] as Move, [more]);
      const before = s.counts();
      await chain.commit([more], [move]);
      const used = difference(s.counts(), before);

      expect(used.reads.pointers).toBe(1);
      expect(used.writes.pointers).toBe(1);
      expect(used.reads.execution).toBe(1);
      expect(used.writes.execution).toBe(1);
      expect(used.reads.heads).toBe(0);
      expect(used.reads.segments).toBe(0);
    });

    it('J3 stop-during-tools: Stop while tools run leaves main as it was before the round', async () => {
      let clock: { sleep(ms: number, signal?: AbortSignal): Promise<void> } | undefined;
      const s = await createSession(fixture, {
        tools: {
          slow: async (context) => {
            await clock?.sleep(60_000, context.signal);
            return { text: 'done' };
          },
        },
        model: [{ calls: [{ tool: 'slow', arguments: {} }] }, { text: 'final' }],
      });
      clock = s.clock;
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.tools.callsOf('slow').length === 1);
      const before = await s.main();
      expect(before.some((message) => message.role === 'assistant')).toBe(false);

      await s.client.stop({ key: 'stop-1', turnId: (await s.client.execution()).turn?.id as string });
      await s.clock.runUntil(running);

      expect(await s.main()).toEqual(before);
      expect((await s.turns())[0]?.status).toBe('stopped');
    });

    it('J3 crash-after-model-step: a crash after a model step with calls settled leaves no assistant message in main', async () => {
      let clock: { sleep(ms: number, signal?: AbortSignal): Promise<void> } | undefined;
      const s = await createSession(fixture, {
        tools: {
          slow: async (context) => {
            await clock?.sleep(60_000, context.signal);
            return { text: 'done' };
          },
        },
        model: [{ calls: [{ tool: 'slow', arguments: {} }] }, { text: 'final' }],
      });
      clock = s.clock;
      const { launch } = await s.send('go', { run: false });
      // The runner is stuck in its tool for good and the simulated time never moves: its crash.
      void s.runner().run(launch as string);
      await until(() => s.tools.callsOf('slow').length === 1);

      const modelStep = stepEvents(await s.events()).find((event) => event.stepKind === 'model');
      expect(modelStep).toMatchObject({ ending: 'settled', moves: [], next: { kind: 'tools' } });
      expect((await s.main()).some((message) => message.role === 'assistant')).toBe(false);
    });

    it('J3 interrupt-model-step: an interrupt answered abandon gives main only the steered items', async () => {
      const s: Session = await createSession(fixture, {
        model: [
          { text: 'slow', waitMs: 60_000 },
          { text: 'after', waitMs: 60_000 },
        ],
        drive: with_({
          interrupt: () => 'abandon',
          boundary: (_decision, context) =>
            context.inbox.length > 0 ? { steer: context.inbox.map((item) => item.key) } : 'proceed',
        }),
      });
      const { launch } = await s.send('first', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);
      const turnId = (await s.client.execution()).turn?.id as string;

      await s.send('interrupt me', { key: 'interrupt-1', interrupt: { turnId }, run: false });
      await until(() => s.model.invocations.length === 2);

      const main = await s.main();
      expect(main.map((message) => message.role)).toEqual(['user', 'user']);
      expect(main.map(textOf)).toEqual(['first', 'interrupt me']);
      await s.clock.runUntil(running);
    });

    it('J4 long-history: 2,000 appends, then a read stays within the bound', async () => {
      const s = await createSession(fixture);
      const chain = await openChain(s);
      for (let i = 0; i < 2_000; i += 1) await chain.append([user(`message ${i}`)]);

      const before = s.counts();
      const [head] = await s.client.heads([chain.head().id]);
      const used = difference(s.counts(), before);

      expect(head?.entries).toEqual(chain.entries);
      expect(used.reads.heads).toBeLessThanOrEqual(1);
      expect(used.reads.segments).toBeLessThanOrEqual(segmentBound(chain.entries.length));
      expect(used.reads.messages).toBe(0);
    });

    it('J4 edit-then-append: an edit, then 500 appends, still reads within the bound', async () => {
      const s = await createSession(fixture);
      const chain = await openChain(s);
      for (let i = 0; i < 600; i += 1) await chain.append([user(`message ${i}`)]);
      await chain.edit(
        chain.entries.filter((_, i) => i % 2 === 0),
        [assistant('a summary of the other half')],
      );
      for (let i = 0; i < 500; i += 1) await chain.append([user(`later ${i}`)]);

      const before = s.counts();
      const [head] = await s.client.heads([chain.head().id]);
      const used = difference(s.counts(), before);

      expect(head?.entries).toEqual(chain.entries);
      expect(used.reads.heads).toBeLessThanOrEqual(1);
      expect(used.reads.segments).toBeLessThanOrEqual(segmentBound(chain.entries.length));
      expect(used.reads.messages).toBe(0);
    });

    it('J4 hash-stable: hashMessage ignores key order and changes with any part', async () => {
      const messages: Message[] = [
        {
          role: 'system',
          parts: [
            { kind: 'text', text: 'be brief' },
            {
              kind: 'tool-definition',
              definition: toolDefinition('lookup', { q: { type: 'string' } }),
              deferred: false,
            },
          ],
        },
        {
          role: 'assistant',
          parts: [
            { kind: 'text', text: 'looking' },
            {
              kind: 'tool-call',
              requestId: 'req-1',
              tool: { id: 'lookup' },
              arguments: { q: 'x', options: { a: 1, b: [true, null, 's'] } },
              providerOptions: { signature: 'abc' },
            },
          ],
        },
        {
          role: 'user',
          parts: [
            {
              kind: 'tool-result',
              requestId: 'req-1',
              tool: { id: 'lookup' },
              status: 'succeeded',
              parts: [{ kind: 'text', text: 'found' }],
            },
            { kind: 'reference', reference: { uri: 'mem://a', mediaType: 'text/plain', byteLength: 3 } },
          ],
        },
      ];

      for (const message of messages) {
        const hash = await hashMessage(message);
        expect(await hashMessage(structuredClone(message))).toBe(hash);
        expect(await hashMessage(reorder(message) as Message)).toBe(hash);

        const changed = variants(message);
        expect(changed.length).toBeGreaterThan(10);
        for (const other of changed) expect(await hashMessage(other as Message)).not.toBe(hash);
      }
    });
  });
}
