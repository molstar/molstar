/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { HistoryManager, type HistoryPort } from '../history.js';

/** A state of named values; checkpoints are copies. */
function store(initial: Record<string, number> = {}) {
  let state = { ...initial };
  const port: HistoryPort<Record<string, number>> = {
    checkpoint: () => ({ ...state }),
    restore: async (c) => {
      state = { ...c };
    },
  };
  return {
    port,
    get: () => state,
    set: (k: string, v: number) => {
      state = { ...state, [k]: v };
    },
  };
}

const person = { kind: 'person' };

describe('HistoryManager', () => {
  it('records a step per run and undoes the latest', async () => {
    const s = store({ a: 0 });
    const h = new HistoryManager(s.port);
    expect(h.canUndo).toBe(false);
    await h.run('Set a to 1', person, async () => s.set('a', 1));
    await h.run('Set a to 2', person, async () => s.set('a', 2));
    expect(h.steps.value.map((x) => x.label)).toEqual(['Set a to 2', 'Set a to 1']);
    const r = await h.undo();
    expect(r.ok).toBe(true);
    expect(s.get().a).toBe(1);
    await h.undo();
    expect(s.get().a).toBe(0);
    expect(h.undoBlocked.value).toBe('Nothing to undo');
  });

  it('restores the checkpoint of a failed run and records nothing', async () => {
    const s = store({ a: 0 });
    const h = new HistoryManager(s.port);
    const out = await h.run(
      'Half done',
      person,
      async () => {
        s.set('a', 5);
        return { ok: false };
      },
      { failed: (v) => !v.ok },
    );
    expect(out.restored).toBe(true);
    expect(s.get().a).toBe(0);
    expect(h.steps.value).toHaveLength(0);
    await expect(
      h.run('Throws', person, async () => {
        s.set('a', 7);
        throw new Error('no');
      }),
    ).rejects.toThrow('no');
    expect(s.get().a).toBe(0);
  });

  it('runs one at a time', async () => {
    const s = store({ a: 0 });
    const h = new HistoryManager(s.port);
    const order: string[] = [];
    let resolveFirst!: () => void;
    const first = h.run('First', person, () => {
      order.push('first start');
      return new Promise<void>((r) => (resolveFirst = r)).then(() => {
        order.push('first end');
      });
    });
    const second = h.run('Second', person, async () => {
      order.push('second');
    });
    expect(h.undoBlocked.value).toMatch(/running/);
    await new Promise((r) => setTimeout(r, 0));
    resolveFirst();
    await Promise.all([first, second]);
    expect(order).toEqual(['first start', 'first end', 'second']);
  });

  it('makes an atomic group one step, and restores all of it on failure', async () => {
    const s = store({ a: 0, b: 0 });
    const h = new HistoryManager(s.port);
    await h.group('Prepare', { kind: 'script' }, async () => {
      await h.run('a', person, async () => s.set('a', 1));
      await h.run('b', person, async () => s.set('b', 1));
    });
    expect(h.steps.value).toHaveLength(1);
    expect(h.steps.value[0].runs).toEqual(['a', 'b']);
    await expect(
      h.group('Broken', { kind: 'script' }, async () => {
        await h.run('a', person, async () => s.set('a', 2));
        throw new Error('stop');
      }),
    ).rejects.toThrow('stop');
    expect(s.get()).toEqual({ a: 1, b: 1 });
    await h.undo();
    expect(s.get()).toEqual({ a: 0, b: 0 });
  });

  it('joins a turn into one entry until another origin runs', async () => {
    const s = store({ a: 0 });
    const h = new HistoryManager(s.port);
    const agent = { kind: 'agent', turn: 't1' };
    const end = h.beginGroup('Make it pretty', agent);
    await h.run('one', agent, async () => s.set('a', 1));
    await h.run('two', agent, async () => s.set('a', 2));
    expect(h.steps.value).toHaveLength(1);
    expect(h.steps.value[0].runs).toEqual(['one', 'two']);
    await h.run('person', person, async () => s.set('a', 3));
    await h.run('three', agent, async () => s.set('a', 4));
    end();
    expect(h.steps.value.map((x) => x.label)).toEqual(['Make it pretty', 'person', 'Make it pretty']);

    const asked = await h.undoTurn('t1');
    expect(asked).toMatchObject({ ok: false, laterSteps: 1 });
    const r = await h.undoTurn('t1', { discardLater: true });
    expect(r.ok).toBe(true);
    expect(s.get().a).toBe(0);
    expect(h.steps.value).toHaveLength(0);
  });

  it('marks earlier steps unsafe after a change outside it', async () => {
    const s = store({ a: 0 });
    const h = new HistoryManager(s.port);
    await h.run('a', person, async () => s.set('a', 1));
    h.invalidate('Changed outside the viewer');
    expect(h.steps.value[0].unsafe).toBe('Changed outside the viewer');
    const r = await h.undo();
    expect(r).toMatchObject({ ok: false });
    expect(s.get().a).toBe(1);
  });

  it('turns checkpoints off over the budget, and says so once', async () => {
    const s = store({ a: 0 });
    const port: HistoryPort<Record<string, number>> = { ...s.port, size: () => 10_000 };
    const h = new HistoryManager(port, { budget: { ms: 1e9, bytes: 1000 }, hardLimit: { ms: 1e9, bytes: 1e9 } });
    const offs: unknown[] = [];
    h.checkpointingOff.subscribe((m) => offs.push(m));
    await h.run('first', person, async () => s.set('a', 1));
    await h.run('second', person, async () => s.set('a', 2));
    expect(offs).toHaveLength(1);
    expect(h.checkpointing.value).toBe(false);
    // The first had its checkpoint taken; the second has none.
    expect(h.steps.value[0].unsafe).toMatch(/large/);
    expect(h.steps.value[1].unsafe).toBeUndefined();
    h.clear();
    expect(h.checkpointing.value).toBe(true);
  });

  it('records a step that cannot be undone', async () => {
    const s = store({ a: 0 });
    const h = new HistoryManager(s.port);
    await h.run('Upload', person, async () => s.set('a', 1), { undoable: false });
    expect(h.steps.value[0].unsafe).toBe('It cannot be undone');
    expect(h.canUndo).toBe(false);
  });
});
