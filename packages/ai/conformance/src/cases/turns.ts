/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type BeginInput,
  type CommitBatch,
  type HarnessEvent,
  type HarnessStorage,
  type LaunchKey,
  type Lease,
  MAIN,
  type RecordStore,
  SCHEMA_VERSION,
  type StorageTx,
} from '@molstar/ai-harness';
import { describe, expect, it } from '@jest/globals';
import { difference, totalWrites } from '../meter.js';
import { type ConformanceFixture, createSession, type Session } from '../session.js';
import { settle } from '../testing/clock.js';
import { assistant, textOf, user } from '../testing/messages.js';
import { appendMove, controls, lastStep, stepEvents, TURN, with_ } from './helpers.js';

/** Polls (letting the runner and the simulated clock's promises settle) until the condition holds. */
async function until(condition: () => boolean, rounds = 200): Promise<void> {
  for (let i = 0; i < rounds; i += 1) {
    if (condition()) return;
    await settle();
  }
  throw new Error('until: the condition never held');
}

/** Moves the simulated clock just past a lease. */
const expire = (s: Session, lease: Lease): Promise<void> =>
  s.clock.advance(Math.max(0, lease.until - s.clock.now()) + 1);

const statuses = async (s: Session): Promise<string[]> => (await s.turns()).map((turn) => turn.status);

/** A turn begun by hand through the runner port: its lease, ids and the head `main` has after `begin`. */
interface Started {
  readonly lease: Lease;
  readonly turnId: string;
  /** The pending step. */
  readonly step: string;
  readonly head: string;
  readonly begin: BeginInput;
}

type Wrap = <T>(call: () => Promise<T>) => Promise<T>;

/**
 * Submits `item-<tag>` (unless a launch is given), claims as `r-<tag>` and begins `hand-turn-<tag>`
 * with that item. `wrap` lets a case run each keyed call twice.
 */
async function startTurn(s: Session, tag: string, launch?: LaunchKey, wrap: Wrap = (call) => call()): Promise<Started> {
  const prompt = user(`hello ${tag}`);
  const given = launch ?? (await wrap(() => s.client.submit({ key: `item-${tag}`, messages: [prompt] }))).launch;
  if (given === undefined) throw new Error('startTurn: nothing launched');
  const { lease, execution } = await wrap(() =>
    s.port.claim({ launch: given, runnerId: `r-${tag}`, schemaVersion: SCHEMA_VERSION }),
  );
  const move = await appendMove(s.port, MAIN, execution.main, [prompt]);
  const begin: BeginInput = {
    key: `begin-${tag}`,
    turnId: `hand-turn-${tag}`,
    step: `hand-step-${tag}-1`,
    items: [`item-${tag}`],
    configuration: TURN,
    messages: [prompt],
    move,
    calls: [],
    next: { id: `hand-step-${tag}-2`, kind: 'model' },
  };
  const begun = await wrap(() => s.port.begin(lease, begin));
  return { lease: begun.lease, turnId: begin.turnId, step: `hand-step-${tag}-2`, head: move.head.id, begin };
}

/** A commit of the pending step that appends `text` to `main`: ends the turn, or goes on with a model step. */
async function settleBatch(
  s: Session,
  started: Started,
  key: string,
  options: { readonly final?: boolean; readonly text?: string } = {},
): Promise<CommitBatch> {
  const message = assistant(options.text ?? 'done');
  const move = await appendMove(s.port, MAIN, started.head, [message]);
  return {
    key,
    turnId: started.turnId,
    step: started.step,
    ending: 'settled',
    calls: [],
    messages: [message],
    moves: [move],
    consumed: [],
    ...(options.final === false
      ? { next: { id: `${started.step}-next`, kind: 'model' as const } }
      : { end: { status: 'completed' as const } }),
  };
}

const abandonBatch = (started: Started, key: string): CommitBatch => ({
  key,
  turnId: started.turnId,
  step: started.step,
  ending: 'abandoned',
  calls: [],
  messages: [],
  moves: [],
  consumed: [],
  end: { status: 'stopped' },
});

const parkBatch = (started: Started, key: string): CommitBatch => ({
  key,
  turnId: started.turnId,
  step: started.step,
  ending: 'parked',
  calls: [],
  messages: [],
  moves: [],
  consumed: [],
  park: { path: [] },
});

/** Backend faults, armed by a case: one fires once and disarms. */
interface Faults {
  /** After the next transaction commits, the caller gets an error anyway. */
  afterTransaction: boolean;
  /** The nth write of the next transaction throws (1-based; 0: never). */
  failOnWrite: number;
  /** The next `drop` throws (a crash after the tombstone). */
  failDrop: boolean;
  fired: number;
}
const noFaults = (): Faults => ({ afterTransaction: false, failOnWrite: 0, failDrop: false, fired: 0 });

function guarded(tx: StorageTx, tick: () => void): StorageTx {
  const store = <T>(inner: RecordStore<T>): RecordStore<T> => ({
    get: (ids) => inner.get(ids),
    put: (id, value) => {
      tick();
      inner.put(id, value);
    },
  });
  return {
    execution: store(tx.execution),
    messages: store(tx.messages),
    pointers: store(tx.pointers),
    heads: store(tx.heads),
    segments: store(tx.segments),
    events: {
      append: (event) => {
        tick();
        tx.events.append(event);
      },
      after: (seq, limit) => tx.events.after(seq, limit),
      byStep: (ids) => tx.events.byStep(ids),
      byTurn: (turnId, after, limit) => tx.events.byTurn(turnId, after, limit),
      endedTurns: (before, limit) => tx.events.endedTurns(before, limit),
      byKey: (key) => tx.events.byKey(key),
      lastCompaction: (pointer) => tx.events.lastCompaction(pointer),
    },
    store: {
      tombstone: (id) => tx.store.tombstone(id),
      setTombstone: (id) => {
        tick();
        tx.store.setTombstone(id);
      },
    },
  };
}

function withFaults(inner: HarnessStorage, faults: Faults): HarnessStorage {
  return {
    putBlobs: (id, blobs) => inner.putBlobs(id, blobs),
    getBlobs: (id, refs) => inner.getBlobs(id, refs),
    watch: (id) => inner.watch(id),
    async drop(id) {
      if (faults.failDrop) {
        faults.failDrop = false;
        faults.fired += 1;
        throw new Error('fault: crash after the tombstone');
      }
      return await inner.drop(id);
    },
    async transact(id, run) {
      let writes = 0;
      const tick = () => {
        writes += 1;
        if (faults.failOnWrite > 0 && writes === faults.failOnWrite) {
          faults.failOnWrite = 0;
          faults.fired += 1;
          throw new Error('fault: in the middle of a transaction');
        }
      };
      const result = await inner.transact(id, (tx) => run(guarded(tx, tick)));
      if (faults.afterTransaction) {
        faults.afterTransaction = false;
        faults.fired += 1;
        throw new Error('fault: after the write');
      }
      return result;
    },
  };
}

const faulty = (fixture: ConformanceFixture, faults: Faults): ConformanceFixture => ({
  name: fixture.name,
  storage: async () => withFaults(await fixture.storage(), faults),
});

/** T1–T6 (Turns and operations). */
export function turnCases(fixture: ConformanceFixture): void {
  describe('Turns and operations', () => {
    // T1. Operations run only under their preconditions (harness.md §12).

    it('T1 begin-while-active: begin with a turn active fails invalid', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      await s.client.submit({ key: 'item-2', messages: [user('hello 2')] });
      const move = await appendMove(s.port, MAIN, started.head, [user('hello 2')]);
      await expect(
        s.port.begin(started.lease, {
          ...started.begin,
          key: 'begin-2',
          turnId: 'hand-turn-2',
          step: 'hand-step-2-1',
          items: ['item-2'],
          messages: [user('hello 2')],
          move,
          next: { id: 'hand-step-2-2', kind: 'model' },
        }),
      ).rejects.toMatchObject({ code: 'invalid' });
    });

    it('T1 schema-differs: a claim of another schema version fails stale-schema and blocks the turn', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      await expire(s, started.lease);
      const { launch } = await s.client.check();
      expect(launch).toBeDefined();

      await expect(
        s.port.claim({ launch: launch as string, runnerId: 'r-other', schemaVersion: SCHEMA_VERSION + 1 }),
      ).rejects.toMatchObject({ code: 'stale-schema' });

      const execution = await s.client.execution();
      expect(execution.turn?.blocked).toBe('stale-schema');
      expect(execution.runner).toEqual({ phase: 'none' });
      expect(controls(await s.events(), 'blocked')).toHaveLength(1);
      // A blocked turn is not work: neither a check nor a new item launches anything.
      expect((await s.client.check()).launch).toBeUndefined();
      expect((await s.client.submit({ key: 'item-2', messages: [user('later')] })).launch).toBeUndefined();
      // A stop ends it; with an item waiting, the stop launches.
      const stopped = await s.client.stop({ key: 'stop-1', turnId: started.turnId });
      expect((await s.client.execution()).turn).toBeUndefined();
      expect(stopped.launch).toBeDefined();
    });

    it('T1 expired-lease: a commit after the lease expired fails stale', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      const batch = await settleBatch(s, started, 'commit-1');
      await expire(s, started.lease);
      await expect(s.port.commit(started.lease, batch)).rejects.toMatchObject({ code: 'stale' });
    });

    it('T1 wrong-step: a commit naming another step fails stale', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      const batch = await settleBatch(s, started, 'commit-1');
      await expect(s.port.commit(started.lease, { ...batch, step: 'hand-step-other' })).rejects.toMatchObject({
        code: 'stale',
      });
    });

    it('T1 zombie: a runner that lost its lease cannot commit after another settled', async () => {
      const s = await createSession(fixture);
      const a = await startTurn(s, '1');
      await expire(s, a.lease);
      const { launch } = await s.client.check();
      const b = await s.port.claim({ launch: launch as string, runnerId: 'r-b', schemaVersion: SCHEMA_VERSION });
      expect(b.execution.turn?.pending.id).toBe(a.step);
      await s.port.commit(b.lease, await settleBatch(s, a, 'commit-b', { text: 'from B' }));
      const fromA = await settleBatch(s, a, 'commit-a', { text: 'from A' });

      const before = s.counts();
      await expect(s.port.commit(a.lease, fromA)).rejects.toMatchObject({ code: 'stale' });
      expect(totalWrites(difference(s.counts(), before))).toBe(0);
      expect(await statuses(s)).toEqual(['completed']);
    });

    // T2. Every operation is atomic and idempotent by key.

    it('T2 same-key-same-result: each keyed operation repeated with its key returns the first result and writes no event', async () => {
      const s = await createSession(fixture);
      const again: Wrap = async (call) => {
        const first = await call();
        const count = (await s.events()).length;
        expect(await call()).toEqual(first);
        expect((await s.events()).length).toBe(count);
        return first;
      };

      // submit, claim, begin
      const started = await startTurn(s, '1', undefined, again);
      // commit
      const batch = await settleBatch(s, started, 'commit-1', { final: false });
      const committed = await again(() => s.port.commit(started.lease, batch));
      // stop, release
      const stopped = await again(() => s.client.stop({ key: 'stop-1', turnId: started.turnId }));
      await again(() => s.port.release(committed.lease, { seenSeq: stopped.seq }));
      // reset
      await again(() => s.client.reset({ key: 'reset-1' }));
      // resume (of a parked turn)
      const second = await startTurn(s, '2');
      await s.port.commit(second.lease, parkBatch(second, 'park-2'));
      await again(() => s.client.resume({ key: 'resume-1', turnId: second.turnId }));
    });

    it('T2 same-key-other-input: a key reused with different input fails invalid', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      await expect(s.client.submit({ key: 'item-1', messages: [user('something else')] })).rejects.toMatchObject({
        code: 'invalid',
      });
      await expect(
        s.port.begin(started.lease, { ...started.begin, messages: [user('something else')] }),
      ).rejects.toMatchObject({ code: 'invalid' });

      const committed = await s.port.commit(started.lease, await settleBatch(s, started, 'commit-1', { final: false }));
      const other = await settleBatch(s, started, 'commit-1', { final: false, text: 'other' });
      await expect(s.port.commit(committed.lease, other)).rejects.toMatchObject({ code: 'invalid' });

      await s.client.stop({ key: 'stop-1', turnId: started.turnId });
      await expect(s.client.stop({ key: 'stop-1', turnId: started.turnId, reason: 'other' })).rejects.toMatchObject({
        code: 'invalid',
      });
      await s.client.reset({ key: 'reset-1' });
      await expect(s.client.reset({ key: 'reset-1', reason: 'other' })).rejects.toMatchObject({ code: 'invalid' });
    });

    it('T2 unknown-outcome: a backend fault after the write, then a retry with the same key, applies once', async () => {
      const faults = noFaults();
      const s = await createSession(faulty(fixture, faults));
      const item = { key: 'item-1', messages: [user('hello 1')] };

      faults.afterTransaction = true;
      await expect(s.client.submit(item)).rejects.toThrow();
      expect(faults.fired).toBe(1);
      const retried = await s.client.submit(item);
      expect(retried.launch).toBeDefined();
      expect((await s.client.execution()).inbox.map((queued) => queued.key)).toEqual(['item-1']);
      expect(controls(await s.events(), 'submitted')).toHaveLength(1);

      const started = await startTurn(s, '1', retried.launch);
      const batch = await settleBatch(s, started, 'commit-1', { final: false });
      faults.afterTransaction = true;
      await expect(s.port.commit(started.lease, batch)).rejects.toThrow();
      expect(faults.fired).toBe(2);
      await s.port.commit(started.lease, batch);
      expect(stepEvents(await s.events()).filter((event) => event.commitKey === 'commit-1')).toHaveLength(1);
      expect((await s.client.execution()).turn?.pending.id).toBe(`${started.step}-next`);
    });

    it('T2 submit-replay-after-consumed: a submit retried after its item began a turn is found by its key and adds no item', async () => {
      const s = await createSession(fixture, { model: [{ text: 'hello back' }] });
      const item = { key: 'prompt-1', messages: [user('go')] };
      const first = await s.client.submit(item);
      await s.run(first.launch as string);
      expect((await s.client.execution()).inbox).toEqual([]);
      const events = (await s.events()).length;

      const before = s.counts();
      const replay = await s.client.submit(item);
      const used = difference(s.counts(), before);

      expect(replay).toEqual(first);
      expect((await s.client.execution()).inbox).toEqual([]);
      expect((await s.events()).length).toBe(events);
      expect(totalWrites(used)).toBe(0);
      // Found through the key index: no journal read, and a handful of event records at most.
      expect(used.reads.messages + used.reads.heads + used.reads.segments).toBe(0);
      expect(used.reads.events).toBeLessThanOrEqual(3);
    });

    it('T2 late-fault: a fault in the middle of a commit leaves nothing written', async () => {
      const faults = noFaults();
      const s = await createSession(faulty(fixture, faults));
      const started = await startTurn(s, '1');
      const batch = await settleBatch(s, started, 'commit-1', { final: false });
      const before = { execution: await s.client.execution(), events: await s.events() };

      faults.failOnWrite = 3;
      await expect(s.port.commit(started.lease, batch)).rejects.toThrow();
      expect(faults.fired).toBe(1);

      expect(await s.client.execution()).toEqual(before.execution);
      expect(await s.events()).toEqual(before.events);
      // Nothing of it stayed behind: the same commit applies now, once.
      await s.port.commit(started.lease, batch);
      expect(stepEvents(await s.events()).filter((event) => event.commitKey === 'commit-1')).toHaveLength(1);
      expect((await s.client.execution()).seq).toBe(before.execution.seq + 1);
    });

    it('T2 remove-tombstone: a removed id cannot be created again, and a crash after the tombstone leaves every operation failing invalid', async () => {
      const faults = noFaults();
      const removed = await createSession(faulty(fixture, faults));
      await removed.host.remove(removed.harnessId);
      await expect(removed.host.create(removed.harnessId, { key: 'create-again' })).rejects.toMatchObject({
        code: 'invalid',
      });
      await expect(removed.client.submit({ key: 'item-1', messages: [user('hi')] })).rejects.toMatchObject({
        code: 'invalid',
      });

      const crashed = await createSession(faulty(fixture, faults));
      faults.failDrop = true;
      await expect(crashed.host.remove(crashed.harnessId)).rejects.toThrow();
      expect(faults.fired).toBe(1);
      await expect(crashed.host.create(crashed.harnessId, { key: 'create-again' })).rejects.toMatchObject({
        code: 'invalid',
      });
      await expect(crashed.client.execution()).rejects.toMatchObject({ code: 'invalid' });
      await expect(crashed.client.submit({ key: 'item-1', messages: [user('hi')] })).rejects.toMatchObject({
        code: 'invalid',
      });
      await expect(crashed.client.check()).rejects.toMatchObject({ code: 'invalid' });
    });

    // T3. Each item is consumed exactly once, or withdrawn.

    it('T3 steer-across-crash: a steered item, a crash before the commit, and the re-run consumes it once', async () => {
      let asked = 0;
      const s: Session = await createSession(fixture, {
        model: [{ text: 'one' }, { text: 'two' }, { text: 'three' }],
        host: { leaseMs: 10_000, launchMs: 10_000 },
        drive: with_({
          'next-turn': (_decision, context) => ({ begin: { keys: context.inbox.slice(0, 1).map((item) => item.key) } }),
          boundary: (_decision, context) => {
            asked += 1;
            // The first runner never answers its first boundary: it crashes here.
            if (asked === 1) return new Promise<never>(() => undefined);
            return context.inbox.length > 0 ? { steer: context.inbox.map((item) => item.key) } : 'proceed';
          },
        }),
      });
      const { launch } = await s.send('first', { key: 'first', run: false });
      await s.client.submit({ key: 'second', messages: [user('second')] });
      const crashing = s.runner();
      void crashing.run(launch as string);
      await until(() => asked === 1);
      crashing.abort('crash');

      await s.clock.advance(11_000);
      const relaunch = await s.client.check();
      expect(relaunch.launch).toBeDefined();
      await s.run(relaunch.launch as string);

      const main = await s.main();
      expect(main.filter((message) => textOf(message) === 'second')).toHaveLength(1);
      const consumed = stepEvents(await s.events()).flatMap((event) => event.consumed);
      expect(consumed.filter((key) => key === 'second')).toHaveLength(1);
      expect((await s.client.execution()).inbox).toEqual([]);
      expect(await statuses(s)).toEqual(['completed']);
    });

    it('T3 withdrawn-then-steered: steering an item withdrawn meanwhile fails inbox-changed and the drive is asked again without it', async () => {
      const asks: string[][] = [];
      const s: Session = await createSession(fixture, {
        model: [{ text: 'one' }],
        drive: with_({
          'next-turn': (_decision, context) => ({ begin: { keys: context.inbox.slice(0, 1).map((item) => item.key) } }),
          boundary: async (_decision, context) => {
            asks.push(context.inbox.map((item) => item.key));
            if (asks.length > 1) return 'proceed';
            await s.client.withdraw({ key: 'second' });
            return { steer: ['second'] };
          },
        }),
      });
      const { launch } = await s.send('first', { key: 'first', run: false });
      await s.client.submit({ key: 'second', messages: [user('second')] });
      await s.run(launch as string);

      expect(asks).toEqual([['second'], []]);
      expect((await s.main()).map(textOf)).not.toContain('second');
      expect(stepEvents(await s.events()).flatMap((event) => event.consumed)).not.toContain('second');
      expect((await s.client.execution()).inbox).toEqual([]);
      expect(await statuses(s)).toEqual(['completed']);
    });

    it('T3 begin-order: begin with three items appends their messages in the given order', async () => {
      const s = await createSession(fixture, {
        model: [{ text: 'done' }],
        drive: with_({ 'next-turn': () => ({ begin: { keys: ['c', 'a', 'b'] } }) }),
      });
      const { launch } = await s.send('a', { key: 'a', run: false });
      await s.client.submit({ key: 'b', messages: [user('b')] });
      await s.client.submit({ key: 'c', messages: [user('c')] });
      await s.run(launch as string);

      expect((await s.main()).slice(0, 3).map(textOf)).toEqual(['c', 'a', 'b']);
    });

    // T4. Stop ends the turn and touches nothing else (harness.md §5).

    it('T4 stop-no-runner: Stop with no live lease ends the turn stopped at once', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      await expire(s, started.lease);

      const stopped = await s.client.stop({ key: 'stop-1', turnId: started.turnId });

      expect((await s.client.execution()).turn).toBeUndefined();
      expect(stopped.launch).toBeUndefined();
      expect(await statuses(s)).toEqual(['stopped']);
      const events = await s.events();
      expect(controls(events, 'stopped')).toHaveLength(1);
      expect(controls(events, 'stop-requested')).toHaveLength(0);
    });

    it("T4 stop-live-runner: Stop with a live runner commits abandoned with the ended attempts' usage and end stopped", async () => {
      const s = await createSession(fixture, { model: [{ text: 'slow', waitMs: 60_000 }] });
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);

      await s.client.stop({ key: 'stop-1', turnId: (await s.client.execution()).turn?.id as string });
      await s.clock.runUntil(running);

      const last = lastStep(await s.events());
      expect(last).toMatchObject({ ending: 'abandoned', end: { status: 'stopped' } });
      const attempts = last.calls.flatMap((call) => call.attempts);
      expect(attempts).toHaveLength(1);
      expect(attempts[0]).toMatchObject({ outcome: 'aborted' });
      expect(attempts[0]?.usage).toBeDefined();
      expect(last.summary?.usage.modelCalls).toBeGreaterThanOrEqual(1);
      expect(await statuses(s)).toEqual(['stopped']);
    });

    it('T4 stop-keeps-inbox: two items wait; after Stop both still wait', async () => {
      let begun = 0;
      const s = await createSession(fixture, {
        model: [{ text: 'slow', waitMs: 60_000 }],
        // The first item begins a turn; after it, the runner releases rather than starting the others.
        drive: with_({
          'next-turn': (_decision, context) =>
            begun++ === 0 ? { begin: { keys: [context.inbox[0]?.key as string] } } : 'release',
        }),
      });
      const { launch } = await s.send('first', { run: false });
      const running = s.runner().run(launch as string);
      // The model is invoked once the turn has begun (advancing time alone lets nothing run).
      await until(() => s.model.invocations.length === 1);
      await s.client.submit({ key: 'second', messages: [user('second')] });
      await s.client.submit({ key: 'third', messages: [user('third')] });
      const turn = (await s.client.execution()).turn;
      expect(turn).toBeDefined();
      await s.client.stop({ key: 'stop-1', turnId: turn?.id as string });
      await s.clock.runUntil(running);
      const after = await s.client.execution();
      expect(after.turn).toBeUndefined();
      expect(after.inbox.map((item) => item.key)).toEqual(['second', 'third']);
    });

    it("T4 stop-beats-settlement: Stop requested before a settlement's commit refuses the commit and abandons the step", async () => {
      const s: Session = await createSession(fixture, {
        model: [{ text: 'final' }],
        drive: with_({
          // The last decision before the settling commit: Stop lands just before it.
          boundary: async () => {
            await s.client.stop({ key: 'stop-1', turnId: (await s.client.execution()).turn?.id as string });
            return 'proceed' as const;
          },
        }),
      });
      await s.send('go');

      expect(await statuses(s)).toEqual(['stopped']);
      expect(lastStep(await s.events())).toMatchObject({ ending: 'abandoned', end: { status: 'stopped' } });
      expect((await s.main()).map((message) => message.role)).toEqual(['user']);
    });

    it('T4 stop-then-expiry: Stop requested, the runner dies, and a check after the lease expires ends the turn stopped', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      await s.client.stop({ key: 'stop-1', turnId: started.turnId });
      expect((await s.client.execution()).turn?.stopRequested).toBeDefined();
      await expire(s, started.lease);

      const checked = await s.client.check();

      expect(checked.launch).toBeUndefined();
      expect((await s.client.execution()).turn).toBeUndefined();
      expect(await statuses(s)).toEqual(['stopped']);
      expect(controls(await s.events(), 'stopped')).toHaveLength(1);
    });

    it('T4 stop-parked-with-items: Stop on a parked turn with items waiting and no runner ends it and returns a launch key', async () => {
      const s = await createSession(fixture, {
        tools: { guarded: () => ({ text: 'ran' }) },
        configuration: { tools: { enabled: ['guarded'], preauthorized: [] } },
        model: [{ calls: [{ tool: 'guarded', arguments: {} }] }],
        drive: with_({ permission: () => ({ park: {} }) }),
      });
      const sent = await s.send('go');
      expect(sent.exit?.reason).toBe('released');
      const turn = (await s.client.execution()).turn;
      expect(turn?.park).toBeDefined();
      // Items behind a park launch nothing (harness.md §5).
      expect((await s.client.submit({ key: 'waiting', messages: [user('more')] })).launch).toBeUndefined();

      const stopped = await s.client.stop({ key: 'stop-1', turnId: turn?.id as string });

      expect(stopped.launch).toBeDefined();
      const after = await s.client.execution();
      expect(after.turn).toBeUndefined();
      expect(after.inbox.map((item) => item.key)).toEqual(['waiting']);
      expect(await statuses(s)).toEqual(['stopped']);
      expect(s.tools.calls).toHaveLength(0);
    });

    it('T4 stop-item-waiting: stop({ item }) while the item waits withdraws it and no turn begins', async () => {
      const s = await createSession(fixture);
      await s.send('go', { key: 'it1', run: false });

      const stopped = await s.client.stop({ key: 'stop-1', item: 'it1' });

      expect(stopped.withdrawn).toBe(true);
      const execution = await s.client.execution();
      expect(execution.inbox).toEqual([]);
      expect(execution.turn).toBeUndefined();
      expect(await s.turns()).toEqual([]);
      expect(controls(await s.events(), 'withdrawn')).toHaveLength(1);
    });

    it('T4 stop-item-begun: stop({ item }) just after a runner began a turn with the item stops that turn', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');

      const stopped = await s.client.stop({ key: 'stop-1', item: 'item-1' });

      expect(stopped.withdrawn ?? false).toBe(false);
      expect((await s.client.execution()).turn?.stopRequested).toBeDefined();
      await s.port.commit(started.lease, abandonBatch(started, 'abandon-1'));
      expect(await statuses(s)).toEqual(['stopped']);
      const events = await s.events();
      expect(controls(events, 'stop-requested')).toHaveLength(1);
      expect(controls(events, 'withdrawn')).toHaveLength(0);
    });

    it('T4 stop-while-launching: Stop on an unowned turn while a runner is launching ends it at once, and the runner finds no turn and releases', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      await expire(s, started.lease);
      const { launch } = await s.client.check();
      expect(launch).toBeDefined();

      await s.client.stop({ key: 'stop-1', turnId: started.turnId });
      expect((await s.client.execution()).turn).toBeUndefined();

      const exit = await s.run(launch as string);

      expect(exit.reason).toBe('released');
      expect(await statuses(s)).toEqual(['stopped']);
    });

    it("T4 stop-right-after-claim: Stop just after a claim; the runner's first commit is abandoned with end stopped", async () => {
      const s = await createSession(fixture, { model: [{ text: 'slow', waitMs: 60_000 }] });
      const started = await startTurn(s, '1');
      await expire(s, started.lease);
      const { launch } = await s.client.check();
      // A runner claims the unowned turn and re-runs its pending step.
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);

      await s.client.stop({ key: 'stop-1', turnId: started.turnId });
      await s.clock.runUntil(running);

      const steps = stepEvents(await s.events());
      expect(steps.map((event) => event.ending)).toEqual(['settled', 'abandoned']);
      expect(steps.at(-1)).toMatchObject({ end: { status: 'stopped' } });
      expect(await statuses(s)).toEqual(['stopped']);
    });

    it("T4 stop-before-release: Stop between a park commit and the runner's release makes the release end the turn stopped", async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      const parked = await s.port.commit(started.lease, parkBatch(started, 'park-1'));

      const stopped = await s.client.stop({ key: 'stop-1', turnId: started.turnId });
      expect((await s.client.execution()).turn).toBeDefined();
      await s.port.release(parked.lease, { seenSeq: stopped.seq });

      expect((await s.client.execution()).turn).toBeUndefined();
      expect(await statuses(s)).toEqual(['stopped']);
      const events = await s.events();
      expect(controls(events, 'stopped')).toHaveLength(1);
      expect(controls(events, 'released')).toHaveLength(1);
    });

    it('T4 stop-listed: a turn ended by Stop with no runner and one ended by reset both appear in turns() with their summaries', async () => {
      const s = await createSession(fixture);
      const first = await startTurn(s, '1');
      await expire(s, first.lease);
      await s.client.stop({ key: 'stop-1', turnId: first.turnId });
      const second = await startTurn(s, '2');
      await s.client.reset({ key: 'reset-1' });

      const turns = await s.turns();

      expect(turns).toMatchObject([
        { id: second.turnId, status: 'reset', startedBy: ['item-2'], firstStep: 'hand-step-2-1' },
        { id: first.turnId, status: 'stopped', startedBy: ['item-1'], firstStep: 'hand-step-1-1' },
      ]);
      for (const turn of turns) {
        expect(turn.usage).toBeDefined();
        expect(turn.endedAt).toBeDefined();
      }
      const events = await s.events();
      expect(controls(events, 'reset')[0]?.summary).toEqual(turns[0]);
      expect(controls(events, 'stopped')[0]?.summary).toEqual(turns[1]);
    });

    it('T4 only-abandon-after-request: while stopRequested, a settled commit fails stale and the abandoned one is accepted', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      await s.client.stop({ key: 'stop-1', turnId: started.turnId });
      expect((await s.client.execution()).turn?.stopRequested).toBeDefined();

      await expect(s.port.commit(started.lease, await settleBatch(s, started, 'commit-1'))).rejects.toMatchObject({
        code: 'stale',
      });
      await s.port.commit(started.lease, abandonBatch(started, 'abandon-1'));

      expect(await statuses(s)).toEqual(['stopped']);
    });

    // T5. Reset ends any turn and moves only the heads it names (harness.md §5).

    it("T5 reset-busy-turn: reset while a runner works ends the turn reset and the runner's next commit fails stale", async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      const batch = await settleBatch(s, started, 'commit-1');

      await s.client.reset({ key: 'reset-1' });

      expect(await statuses(s)).toEqual(['reset']);
      expect((await s.client.execution()).turn).toBeUndefined();
      await expect(s.port.commit(started.lease, batch)).rejects.toMatchObject({ code: 'stale' });
    });

    /** A completed turn, then a second one begun by hand, whose `begin` moved `main` past the first. */
    const completedThenBegun = async (s: Session) => {
      await s.send('first');
      const completed = (await s.client.execution()).main;
      const completedMessages = await s.main();
      await startTurn(s, '2');
      const during = (await s.client.execution()).main;
      const duringMessages = await s.main();
      expect(during).not.toBe(completed);
      return { completed, completedMessages, during, duringMessages };
    };

    it('T5 reset-default: reset without heads moves main to the head at the end of the last completed turn', async () => {
      const s = await createSession(fixture, { model: [{ text: 'one' }] });
      const { completed, completedMessages } = await completedThenBegun(s);

      await s.client.reset({ key: 'reset-1' });

      expect((await s.client.execution()).main).toBe(completed);
      expect(await s.main()).toEqual(completedMessages);
      expect(await statuses(s)).toEqual(['reset', 'completed']);
    });

    it("T5 reset-undo: a second reset to the first one's starting head restores it", async () => {
      const s = await createSession(fixture, { model: [{ text: 'one' }] });
      const { during, duringMessages } = await completedThenBegun(s);
      await s.client.reset({ key: 'reset-1' });

      await s.client.reset({ key: 'reset-2', pointers: { main: during } });

      expect((await s.client.execution()).main).toBe(during);
      expect(await s.main()).toEqual(duringMessages);
    });

    // T6. The event log is complete and ordered (harness.md §3, §12, §13).

    it("T6 dense-seq: a scripted session's events have consecutive seqs", async () => {
      const s = await createSession(fixture, {
        tools: { echo: () => ({ text: 'echo' }) },
        model: [{ calls: [{ tool: 'echo', arguments: {} }] }, { text: 'one' }, { text: 'two' }],
      });
      await s.send('first');
      await s.send('second');
      await s.send('third', { key: 'third', run: false });
      await s.client.withdraw({ key: 'third' });

      const events = await s.events();

      expect(events.some((event) => event.kind === 'step')).toBe(true);
      expect(events.some((event) => event.kind === 'control')).toBe(true);
      expect(events.length).toBeGreaterThan(5);
      const first = events[0]?.seq as number;
      expect(events.map((event) => event.seq)).toEqual(events.map((_, i) => first + i));
      expect((await s.client.execution()).seq).toBe(events.at(-1)?.seq);
    });

    it('T6 renew-silent: 100 renewals write no event and leave seq unchanged; a check that changes nothing writes none', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      const before = { events: await s.events(), seq: (await s.client.execution()).seq };

      let lease = started.lease;
      for (let i = 0; i < 100; i += 1) {
        await s.clock.advance(1);
        lease = await s.port.renew(lease);
      }
      expect(await s.events()).toEqual(before.events);
      expect((await s.client.execution()).seq).toBe(before.seq);

      const checked = await s.client.check();
      expect(checked.launch).toBeUndefined();
      expect(await s.events()).toEqual(before.events);
      expect((await s.client.execution()).seq).toBe(before.seq);
    });

    it('T6 two-events: a claim on a stop-requested turn whose lease expired writes stopped, then claimed', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      await s.client.stop({ key: 'stop-1', turnId: started.turnId });
      await expire(s, started.lease);
      const submitted = await s.client.submit({ key: 'item-2', messages: [user('next')] });
      expect(submitted.launch).toBeDefined();
      const before = (await s.events()).length;

      await s.port.claim({ launch: submitted.launch as string, runnerId: 'r-2', schemaVersion: SCHEMA_VERSION });

      const added = (await s.events()).slice(before);
      expect(added.map((event) => (event.kind === 'control' ? event.control : event.kind))).toEqual([
        'stopped',
        'claimed',
      ]);
      expect(added[1]?.seq).toBe((added[0]?.seq as number) + 1);
      expect(await statuses(s)).toEqual(['stopped']);
    });

    it('T6 catch-up: a subscriber joining mid-turn subscribes, then reads events(after), and sees each event once', async () => {
      let clock: { sleep(ms: number, signal?: AbortSignal): Promise<void> } | undefined;
      const s = await createSession(fixture, {
        tools: {
          slow: async (context) => {
            await clock?.sleep(5_000, context.signal);
            return { text: 'ok' };
          },
        },
        model: [{ calls: [{ tool: 'slow', arguments: {} }] }, { text: 'final' }],
      });
      clock = s.clock;
      const { launch } = await s.send('go', { run: false });
      const runner = s.runner();
      const running = runner.run(launch as string);
      await until(() => s.tools.callsOf('slow').length === 1);

      const seen: number[] = [];
      const observed: number[] = [];
      let last = 0;
      const take = (event: HarnessEvent) => {
        if (event.seq <= last) return;
        last = event.seq;
        seen.push(event.seq);
      };
      runner.observe().subscribe((observation) => {
        if (observation.kind !== 'event') return;
        observed.push(observation.event.seq);
        take(observation.event);
      });
      for (const event of await s.client.events({ after: last })) take(event);
      await s.clock.runUntil(running);
      for (const event of await s.client.events({ after: last })) take(event);

      expect(observed.length).toBeGreaterThan(0);
      expect(seen).toEqual((await s.events()).map((event) => event.seq));
    });
  });
}
