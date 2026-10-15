/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Message } from '@molstar/ai-agent';
import {
  type ControlEvent,
  type Execution,
  HarnessError,
  type HarnessEvent,
  type StepEvent,
  type StorageTx,
} from '@molstar/ai-harness';
import { createMemoryStorage } from '../storage/memory.js';

const H = 'h1';
const usage = { modelCalls: 0, toolCalls: 0, uncachedInput: 0, cacheRead: 0, cacheWrite: 0, output: 0, reasoning: 0 };

const execution = (seq: number, main = 'm0'): Execution => ({
  harnessId: H,
  seq,
  schemaVersion: 1,
  main,
  inbox: [],
  runner: { phase: 'none' },
});
const step = (seq: number, o: Partial<StepEvent> = {}): StepEvent => ({
  kind: 'step',
  seq,
  at: seq,
  turnId: 't1',
  step: `s${seq}`,
  stepKind: 'model',
  ending: 'settled',
  calls: [],
  moves: [],
  consumed: [],
  commitKey: `k${seq}`,
  ...o,
});
const control = (seq: number, o: Partial<ControlEvent> = {}): ControlEvent => ({
  kind: 'control',
  seq,
  at: seq,
  control: 'submitted',
  detail: null,
  ...o,
});
const summary = { id: 't1', status: 'completed' as const, startedBy: [], firstStep: 's1', usage };
const message = (text: string): Message => ({ role: 'user', content: [{ type: 'text', text }] }) as unknown as Message;

const seed = async (storage: ReturnType<typeof createMemoryStorage>, events: HarnessEvent[]) =>
  storage.transact(H, async (tx) => {
    for (const e of events) tx.events.append(e);
  });
const seqs = (events: readonly HarnessEvent[]) => events.map((e) => e.seq);

describe('memory storage: transactions', () => {
  it('reads its own writes and commits them', async () => {
    const s = createMemoryStorage();
    await s.transact(H, async (tx) => {
      tx.execution.put(H, execution(1));
      tx.events.append(step(1));
      expect((await tx.execution.get([H]))[0]?.seq).toBe(1);
      expect(seqs(await tx.events.after(0, 10))).toEqual([1]);
    });
    await s.transact(H, async (tx) => {
      expect((await tx.execution.get([H]))[0]?.seq).toBe(1);
      expect(seqs(await tx.events.after(0, 10))).toEqual([1]);
    });
  });

  it('applies nothing when run throws, and propagates the error', async () => {
    const s = createMemoryStorage();
    await expect(
      s.transact(H, async (tx) => {
        tx.execution.put(H, execution(1));
        tx.events.append(step(1));
        tx.messages.put('a', message('a'));
        tx.store.setTombstone('x');
        throw new Error('boom');
      }),
    ).rejects.toThrow('boom');
    await s.transact(H, async (tx) => {
      expect(await tx.execution.get([H])).toEqual([undefined]);
      expect(await tx.events.after(0, 10)).toEqual([]);
      expect(await tx.messages.get(['a'])).toEqual([undefined]);
      expect(await tx.store.tombstone('x')).toBe(false);
    });
  });

  it('runs transactions of one harness one at a time, in call order', async () => {
    const s = createMemoryStorage();
    const log: string[] = [];
    const slow = (name: string, ms: number) =>
      s.transact(H, async (tx) => {
        log.push(`${name}:start`);
        const [e] = await tx.execution.get([H]);
        await new Promise((r) => setTimeout(r, ms));
        tx.execution.put(H, execution((e?.seq ?? 0) + 1));
        log.push(`${name}:end`);
      });
    await Promise.all([slow('a', 20), slow('b', 1), slow('c', 1)]);
    expect(log).toEqual(['a:start', 'a:end', 'b:start', 'b:end', 'c:start', 'c:end']);
    await s.transact(H, async (tx) => expect((await tx.execution.get([H]))[0]?.seq).toBe(3));
  });

  it('keeps running after a failed transaction, and lets other harnesses interleave', async () => {
    const s = createMemoryStorage();
    const fail = s.transact(H, async () => {
      throw new Error('x');
    });
    await expect(fail).rejects.toThrow('x');
    let release: () => void = () => {};
    const gate = new Promise<void>((r) => (release = r));
    const blocked = s.transact('other', async () => gate);
    await s.transact(H, async (tx) => tx.execution.put(H, execution(1)));
    release();
    await blocked;
  });

  it('stores copies, so callers cannot mutate records', async () => {
    const s = createMemoryStorage();
    const e = execution(1);
    await s.transact(H, async (tx) => tx.execution.put(H, e));
    await s.transact(H, async (tx) => {
      const [got] = await tx.execution.get([H]);
      (got as { main: string }).main = 'changed';
    });
    await s.transact(H, async (tx) => expect((await tx.execution.get([H]))[0]?.main).toBe('m0'));
  });
});

describe('memory storage: immutability', () => {
  const integrity = (f: () => void) => {
    try {
      f();
    } catch (e) {
      expect(e).toBeInstanceOf(HarnessError);
      expect((e as HarnessError).code).toBe('integrity');
      return;
    }
    throw new Error('expected integrity');
  };

  it('refuses a changed put of a message, head record or segment, at the put call', async () => {
    const s = createMemoryStorage();
    await s.transact(H, async (tx) => {
      tx.messages.put('a', message('a'));
      tx.heads.put('h', { id: 'h', segments: ['g'], length: 1 });
      tx.segments.put('g', { id: 'g', entries: ['a'] });
    });
    await s.transact(H, async (tx) => {
      tx.messages.put('a', message('a'));
      tx.heads.put('h', { id: 'h', segments: ['g'], length: 1 });
      tx.segments.put('g', { id: 'g', entries: ['a'] });
      integrity(() => tx.messages.put('a', message('b')));
      integrity(() => tx.heads.put('h', { id: 'h', segments: ['g'], length: 2 }));
      integrity(() => tx.segments.put('g', { id: 'g', entries: ['a', 'b'] }));
    });
  });

  it('refuses a change against a record written earlier in the same transaction', async () => {
    const s = createMemoryStorage();
    await s.transact(H, async (tx) => {
      tx.messages.put('a', message('a'));
      integrity(() => tx.messages.put('a', message('b')));
    });
  });

  it('lets the Execution and pointers be overwritten', async () => {
    const s = createMemoryStorage();
    await s.transact(H, async (tx) => {
      tx.execution.put(H, execution(1));
      tx.pointers.put('p', { scope: 'p', head: 'a' });
    });
    await s.transact(H, async (tx) => {
      tx.execution.put(H, execution(2));
      tx.pointers.put('p', { scope: 'p', head: 'b' });
    });
    await s.transact(H, async (tx) => {
      expect((await tx.execution.get([H]))[0]?.seq).toBe(2);
      expect((await tx.pointers.get(['p']))[0]?.head).toBe('b');
    });
  });

  it('refuses an event whose seq exists', async () => {
    const s = createMemoryStorage();
    await seed(s, [step(1)]);
    await s.transact(H, async (tx) => {
      integrity(() => tx.events.append(step(1)));
      tx.events.append(step(2));
      integrity(() => tx.events.append(control(2)));
    });
  });
});

describe('memory storage: event log', () => {
  const events: HarnessEvent[] = [
    step(1, { turnId: 't1', step: 'a', commitKey: 'k1' }),
    control(2, { turnId: 't1', key: 'k1' }),
    step(3, {
      turnId: 't2',
      step: 'a',
      commitKey: 'k3',
      compaction: [{ pointer: 'main', sizeAfter: { inputTokens: 1, reservedOutput: 1, window: 9 } }],
    }),
    step(4, { turnId: 't1', step: 'b', end: { status: 'completed' }, summary }),
    control(5, { turnId: 't2', control: 'stopped', summary: { ...summary, id: 't2' } }),
    step(6, {
      turnId: 't3',
      step: 'c',
      compaction: [
        { pointer: 'side', sizeAfter: { inputTokens: 1, reservedOutput: 1, window: 9 } },
        { pointer: 'main', sizeAfter: { inputTokens: 2, reservedOutput: 1, window: 9 } },
      ],
    }),
    step(7, { turnId: 't3', step: 'd', end: { status: 'failed' } }),
  ];

  const both = async (check: (tx: StorageTx) => Promise<void>) => {
    const s = createMemoryStorage();
    // Same answers from committed records and from a transaction's own writes.
    await s.transact(H, async (tx) => {
      for (const e of events) tx.events.append(e);
      await check(tx);
    });
    await s.transact(H, check);
    // And split across the commit boundary.
    const split = createMemoryStorage();
    await seed(split, events.slice(0, 3));
    await split.transact(H, async (tx) => {
      for (const e of events.slice(3)) tx.events.append(e);
      await check(tx);
    });
  };

  it('after: ascending, past seq, at most limit', () =>
    both(async (tx) => {
      expect(seqs(await tx.events.after(0, 100))).toEqual([1, 2, 3, 4, 5, 6, 7]);
      expect(seqs(await tx.events.after(2, 3))).toEqual([3, 4, 5]);
      expect(seqs(await tx.events.after(7, 3))).toEqual([]);
    }));

  it('byStep: step events of the given ids in seq order', () =>
    both(async (tx) => {
      expect(seqs(await tx.events.byStep(['b', 'a']))).toEqual([1, 3, 4]);
      expect(seqs(await tx.events.byStep(['zzz']))).toEqual([]);
    }));

  it('byTurn: step and control events of the turn after a seq', () =>
    both(async (tx) => {
      expect(seqs(await tx.events.byTurn('t1', undefined, 10))).toEqual([1, 2, 4]);
      expect(seqs(await tx.events.byTurn('t1', 1, 10))).toEqual([2, 4]);
      expect(seqs(await tx.events.byTurn('t1', undefined, 2))).toEqual([1, 2]);
      expect(seqs(await tx.events.byTurn('none', undefined, 2))).toEqual([]);
    }));

  it('endedTurns: turn-ending events, newest first, before a seq', () =>
    both(async (tx) => {
      expect(seqs(await tx.events.endedTurns(undefined, 10))).toEqual([7, 5, 4]);
      expect(seqs(await tx.events.endedTurns(7, 10))).toEqual([5, 4]);
      expect(seqs(await tx.events.endedTurns(5, 10))).toEqual([4]);
      expect(seqs(await tx.events.endedTurns(undefined, 2))).toEqual([7, 5]);
      expect(seqs(await tx.events.endedTurns(4, 10))).toEqual([]);
    }));

  it('byKey: step commit keys and control keys in seq order', () =>
    both(async (tx) => {
      expect(seqs(await tx.events.byKey('k1'))).toEqual([1, 2]);
      expect(seqs(await tx.events.byKey('k3'))).toEqual([3]);
      expect(await tx.events.byKey('nope')).toEqual([]);
    }));

  it('lastCompaction: the newest step compacting the pointer', () =>
    both(async (tx) => {
      expect((await tx.events.lastCompaction('main'))?.seq).toBe(6);
      expect((await tx.events.lastCompaction('side'))?.seq).toBe(6);
      expect(await tx.events.lastCompaction('other')).toBeUndefined();
    }));

  it('lastCompaction finds the older event when the pointer was compacted once', async () => {
    const s = createMemoryStorage();
    await seed(s, events.slice(0, 5));
    await s.transact(H, async (tx) => expect((await tx.events.lastCompaction('main'))?.seq).toBe(3));
  });
});

describe('memory storage: watch', () => {
  it('emits the current seq on subscribe, then each seq change, and stops on teardown', async () => {
    const s = createMemoryStorage();
    const seen: number[] = [];
    const stop = s.watch(H).subscribe((n) => seen.push(n));
    expect(seen).toEqual([0]);
    await s.transact(H, async (tx) => tx.execution.put(H, execution(1)));
    await s.transact(H, async (tx) => tx.execution.put(H, execution(1, 'm1')));
    await s.transact(H, async (tx) => tx.events.append(step(1)));
    await s
      .transact(H, async (tx) => {
        tx.execution.put(H, execution(2));
        throw new Error('rolled back');
      })
      .catch(() => {});
    expect(seen).toEqual([0, 1]);
    await s.transact(H, async (tx) => tx.execution.put(H, execution(2)));
    expect(seen).toEqual([0, 1, 2]);
    stop();
    await s.transact(H, async (tx) => tx.execution.put(H, execution(3)));
    expect(seen).toEqual([0, 1, 2]);
    const late: number[] = [];
    s.watch(H).subscribe((n) => late.push(n));
    expect(late).toEqual([3]);
  });
});

describe('memory storage: drop, tombstones, blobs', () => {
  it("drop removes a harness's records but not its tombstone", async () => {
    const s = createMemoryStorage();
    await s.transact(H, async (tx) => {
      tx.execution.put(H, execution(1));
      tx.events.append(step(1));
      tx.messages.put('a', message('a'));
      tx.pointers.put('p', { scope: 'p', head: 'a' });
      tx.heads.put('h', { id: 'h', segments: [], length: 0 });
      tx.segments.put('g', { id: 'g', entries: [] });
      tx.store.setTombstone(H);
    });
    await s.transact('keep', async (tx) => tx.messages.put('a', message('a')));
    await s.putBlobs(H, [{ ref: 'r', body: 1 }]);
    await s.drop(H);
    await s.transact(H, async (tx) => {
      expect(await tx.execution.get([H])).toEqual([undefined]);
      expect(await tx.events.after(0, 10)).toEqual([]);
      expect(await tx.messages.get(['a'])).toEqual([undefined]);
      expect(await tx.pointers.get(['p'])).toEqual([undefined]);
      expect(await tx.heads.get(['h'])).toEqual([undefined]);
      expect(await tx.segments.get(['g'])).toEqual([undefined]);
      expect(await tx.store.tombstone(H)).toBe(true);
    });
    expect(await s.getBlobs(H, ['r'])).toEqual([undefined]);
    await s.transact('keep', async (tx) => expect((await tx.messages.get(['a']))[0]).toBeDefined());
    // A dropped harness can be reused: events start over.
    await seed(s, [step(1)]);
  });

  it('tombstones are store-wide and transactional', async () => {
    const s = createMemoryStorage();
    await s.transact('a', async (tx) => {
      expect(await tx.store.tombstone('z')).toBe(false);
      tx.store.setTombstone('z');
      expect(await tx.store.tombstone('z')).toBe(true);
    });
    await s.transact('b', async (tx) => expect(await tx.store.tombstone('z')).toBe(true));
  });

  it('blobs: put, get by ref, undefined for a missing ref, per harness, copied', async () => {
    const s = createMemoryStorage();
    const body = { a: [1] };
    await s.putBlobs(H, [
      { ref: 'r1', body },
      { ref: 'r2', body: 'text' },
    ]);
    body.a.push(2);
    const [r1, missing, r2] = await s.getBlobs(H, ['r1', 'r3', 'r2']);
    expect(r1).toEqual({ a: [1] });
    expect(missing).toBeUndefined();
    expect(r2).toBe('text');
    expect(await s.getBlobs('other', ['r1'])).toEqual([undefined]);
  });
});
