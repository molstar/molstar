/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { HistoryModel, type HistoryPort } from '../history-model.js';

// The history model in plain Node (history §1): steps, groups, turns, undo and redo, checkpoints that switch off.
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
const script = { kind: 'script' };
/** Asserts an undo or redo went. */
const expectOk = (r: { ok: boolean }) => expect(r).toMatchObject({ ok: true });
const tick = () => new Promise((r) => setTimeout(r, 0));

/** A promise and its resolver. */
function gate() {
  let open!: () => void;
  const p = new Promise<void>((r) => (open = r));
  return { p, open };
}

describe('HistoryModel', () => {
  it('records a step per run and undoes the latest', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    expect(h.canUndo.value).toBe(false);
    await h.run('Set a to 1', person, async () => s.set('a', 1));
    await h.run('Set a to 2', person, async () => s.set('a', 2));
    expect(h.steps.value.map((x) => x.label)).toEqual(['Set a to 2', 'Set a to 1']);
    expect(h.canUndo.value).toBe(true);
    const r = await h.undo();
    expect(r.ok).toBe(true);
    expect(s.get().a).toBe(1);
    await h.undo();
    expect(s.get().a).toBe(0);
    expect(h.undoBlocked.value).toBe('Nothing to undo');
  });

  it('redoes what undo took back, in order, until a new step', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    await h.run('Set a to 1', person, async () => s.set('a', 1));
    await h.run('Set a to 2', person, async () => s.set('a', 2));
    expect((await h.redo()).ok).toBe(false);
    await h.undo();
    await h.undo();
    expect(s.get().a).toBe(0);
    expect(h.redoable.value.map((x) => x.label)).toEqual(['Set a to 1', 'Set a to 2']);
    expect((await h.redo()).ok).toBe(true);
    expect(s.get().a).toBe(1);
    expect(h.steps.value.map((x) => x.label)).toEqual(['Set a to 1']);
    // Undo again works on the redone step; Redo goes on.
    await h.redo();
    expect(s.get().a).toBe(2);
    expect(h.steps.value.map((x) => x.label)).toEqual(['Set a to 2', 'Set a to 1']);
    await h.undo();
    // A new step: nothing to redo.
    await h.run('Set a to 5', person, async () => s.set('a', 5));
    expect(h.redoable.value).toEqual([]);
    expect((await h.redo()).ok).toBe(false);
    // A change outside the history empties it too.
    await h.undo();
    expect(h.redoable.value).toHaveLength(1);
    h.invalidate('changed elsewhere');
    expect(h.redoable.value).toEqual([]);
  });

  it('restores the checkpoint of a failed run and records nothing', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    const out = await h.run(
      'Half done',
      person,
      async () => {
        s.set('a', 5);
        return { ok: false };
      },
      { failed: (v) => !v.ok },
    );
    expect(out).toMatchObject({ failed: true, restored: true });
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

  it('runs one at a time, and Undo waits for queued runs', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    const order: string[] = [];
    const g = gate();
    const first = h.run('First', person, async () => {
      order.push('first start');
      await g.p;
      order.push('first end');
    });
    const second = h.run('Second', person, async () => {
      order.push('second');
    });
    expect(h.undoBlocked.value).toMatch(/running/);
    expect(h.undoBlockedCode.value).toBe('running');
    await tick();
    expect(await h.undo()).toMatchObject({ ok: false, code: 'running' });
    g.open();
    await Promise.all([first, second]);
    expect(order).toEqual(['first start', 'first end', 'second']);
  });

  it('runs work that calls back in with its context as part of the step, without waiting', async () => {
    const s = store({ a: 0, b: 0 });
    const h = new HistoryModel(s.port);
    await h.run('Outer', person, async (ctx) => {
      s.set('a', 1);
      await h.run('Inner', person, async () => s.set('b', 1), { within: ctx });
    });
    expect(h.steps.value).toHaveLength(1);
    expect(h.steps.value[0].runs).toEqual(['Outer', 'Inner']);
    await h.undo();
    expect(s.get()).toEqual({ a: 0, b: 0 });
  });

  it('makes an atomic group one step, all or nothing', async () => {
    const s = store({ a: 0, b: 0 });
    const h = new HistoryModel(s.port);
    await h.group('Prepare', script, async (ctx) => {
      await h.run('a', script, async () => s.set('a', 1), { within: ctx });
      await h.run('b', script, async () => s.set('b', 1), { within: ctx });
    });
    expect(h.steps.value).toHaveLength(1);
    expect(h.steps.value[0].runs).toEqual(['a', 'b']);
    await expect(
      h.group('Broken', script, async (ctx) => {
        await h.run('a', script, async () => s.set('a', 2), { within: ctx });
        throw new Error('stop');
      }),
    ).rejects.toThrow('stop');
    expect(s.get()).toEqual({ a: 1, b: 1 });
    await h.undo();
    expect(s.get()).toEqual({ a: 0, b: 0 });
  });

  it('fails a group when a run inside failed, unless it continues', async () => {
    const s = store({ a: 0, b: 0 });
    const h = new HistoryModel(s.port);
    const failing = async (ctx: { id: number }) =>
      h.run(
        'bad',
        script,
        async () => {
          s.set('b', 9);
          return false;
        },
        { within: ctx, failed: (v) => v === false },
      );
    const g = await h.group('Ignores a failure', script, async (ctx) => {
      await h.run('a', script, async () => s.set('a', 1), { within: ctx });
      await failing(ctx);
    });
    expect(g).toMatchObject({ failed: true, restored: true, failures: 1 });
    expect(s.get()).toEqual({ a: 0, b: 0 });

    const c = await h.group(
      'Goes on',
      script,
      async (ctx) => {
        await h.run('a', script, async () => s.set('a', 1), { within: ctx });
        await failing(ctx);
      },
      { onFailure: 'continue' },
    );
    expect(c).toMatchObject({ failed: false, failures: 1, runs: ['a'] });
    // The failed run alone was undone.
    expect(s.get()).toEqual({ a: 1, b: 0 });
  });

  it("keeps another origin's call out of a group in flight: it runs after, as its own step", async () => {
    const s = store({ a: 0, p: 0 });
    const h = new HistoryModel(s.port);
    const g = gate();
    const group = h.group(
      'Slow script',
      script,
      async (ctx) => {
        await h.run('a', script, async () => s.set('a', 1), { within: ctx });
        await g.p;
        return 'failed';
      },
      { failed: (v) => v === 'failed' },
    );
    const click = h.run('Click', person, async () => s.set('p', 1));
    await tick();
    expect(s.get().p).toBe(0);
    g.open();
    await group;
    await click;
    // The group failed and was undone; the click ran after it and stands.
    expect(s.get()).toEqual({ a: 0, p: 1 });
    expect(h.steps.value.map((x) => x.label)).toEqual(['Click']);
  });

  it('joins a turn into one entry until another origin runs; a group joins it too', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    const agent = { kind: 'agent', turn: 't1' };
    const end = h.beginGroup('Make it pretty', agent);
    await h.run('one', agent, async () => s.set('a', 1));
    await h.group('snippet', agent, async (ctx) => {
      await h.run('two', agent, async () => s.set('a', 2), { within: ctx });
    });
    expect(h.steps.value).toHaveLength(1);
    expect(h.steps.value[0].runs).toEqual(['one', 'two']);
    await h.run('person', person, async () => s.set('a', 3));
    await h.run('three', agent, async () => s.set('a', 4));
    end();
    expect(h.steps.value.map((x) => x.label)).toEqual(['Make it pretty', 'person', 'Make it pretty']);

    const asked = await h.undoTurn('t1');
    expect(asked).toMatchObject({
      ok: false,
      code: 'later-steps',
      laterSteps: 1,
      reason: '1 later step would be undone too',
    });
    const r = await h.undoTurn(agent, { discardLater: true });
    expect(r.ok).toBe(true);
    expect(s.get().a).toBe(0);
    expect(h.steps.value).toHaveLength(0);
  });

  it('records nothing for a run that changed nothing, or that cannot be undone', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    await h.run('Set a', person, async () => s.set('a', 1));
    await h.run('Already so', person, async () => 'same', { record: () => false });
    await h.run('Upload', person, async () => undefined, { undoable: false });
    expect(h.steps.value.map((x) => x.label)).toEqual(['Set a']);
    expect((await h.undo()).ok).toBe(true);
    expect(s.get().a).toBe(0);
  });

  it('marks earlier steps unsafe after a change outside it', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    await h.run('a', person, async () => s.set('a', 1));
    h.invalidate('Changed outside the viewer');
    expect(h.steps.value[0].unsafe).toBe('Changed outside the viewer');
    const r = await h.undo();
    expect(r).toMatchObject({ ok: false, code: 'unsafe', reason: 'Changed outside the viewer' });
    expect(h.steps.value[0].unsafeCode).toBe('changed');
    expect(s.get().a).toBe(1);
  });

  it('turns checkpoints off over the budget and says so once; a group keeps its checkpoint', async () => {
    const s = store({ a: 0 });
    const port: HistoryPort<Record<string, number>> = { ...s.port, size: () => 10_000 };
    const h = new HistoryModel(port, { budget: { ms: 1e9, bytes: 1000 }, hardLimit: { ms: 1e9, bytes: 1e9 } });
    const offs: unknown[] = [];
    const stop = h.checkpointsOff.listen!((m) => offs.push(m));
    await h.run('first', person, async () => s.set('a', 1));
    await h.run('second', person, async () => s.set('a', 2));
    stop();
    expect(offs).toHaveLength(1);
    expect(h.checkpointsOff.value).toMatchObject({ bytes: 10_000 });
    expect(h.checkpointing.value).toBe(false);
    expect(h.steps.value[0].unsafe).toMatch(/large/);
    expect(h.steps.value[1].unsafe).toBeUndefined();

    const g = await h
      .group('Group', script, async () => {
        s.set('a', 3);
        throw new Error('x');
      })
      .catch(() => 'thrown');
    expect(g).toBe('thrown');
    expect(s.get().a).toBe(2);
    h.clear();
    expect(h.checkpointing.value).toBe(true);
  });

  it('takes no checkpoint over the hard limit, and says it is too large', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel({ ...s.port, size: () => 1e9 }, { hardLimit: { ms: 1e9, bytes: 1000 } });
    await h.run('big', person, async () => s.set('a', 1));
    expect(h.steps.value[0]).toMatchObject({ unsafe: 'Too large to undo', unsafeCode: 'too-large' });
    expect(h.canUndo.value).toBe(false);
  });

  it('past the hard limit, takes no snapshot for a while, then tries again; a clear or undo turned on lets it at once', async () => {
    const s = store({ a: 0 });
    let taken = 0;
    let time = 0;
    let big = true;
    const port: HistoryPort<Record<string, number>> = {
      ...s.port,
      checkpoint: () => {
        taken++;
        return s.port.checkpoint();
      },
      size: () => (big ? 1e9 : 10),
    };
    const h = new HistoryModel(port, { hardLimit: { ms: 1e9, bytes: 1000 }, now: () => time });
    const end = h.beginGroup('Turn', { kind: 'agent', turn: 't' });
    for (let i = 0; i < 5; i++) await h.run(`call ${i}`, { kind: 'agent', turn: 't' }, async () => s.set('a', i));
    end();
    expect(taken).toBe(1);
    expect(h.steps.value[0].unsafeCode).toBe('too-large');
    time += 60_000;
    await h.group('Later', script, async (ctx) => {
      await h.run('x', script, async () => s.set('a', 9), { within: ctx });
    });
    expect(taken).toBe(2);
    // Undo turned back on (Settings) tries again at once, as a clear does.
    h.setCheckpointing(true);
    await h.group('Again', script, async (ctx) => {
      await h.run('y', script, async () => s.set('a', 8), { within: ctx });
    });
    expect(taken).toBe(3);
    big = false;
    h.clear();
    await h.run('Small', script, async () => s.set('a', 10));
    expect(taken).toBe(4);
    expect(h.steps.value[0].unsafe).toBeUndefined();
  });

  it('counts the size pass in the time budget', async () => {
    const s = store({ a: 0 });
    const slowSize = () => {
      const t0 = performance.now();
      while (performance.now() - t0 < 3);
      return 10;
    };
    const h = new HistoryModel({ ...s.port, size: slowSize }, { budget: { ms: 1, bytes: 1e9 } });
    for (let i = 0; i < 3; i++) await h.run(`r${i}`, person, async () => s.set('a', i));
    expect(h.checkpointing.value).toBe(false);
  });

  it('keeps checkpoints for the latest steps only; older ones stay listed, too old to undo', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port, { capacity: 6, checkpoints: 3 });
    for (let i = 1; i <= 8; i++) await h.run(`Set a to ${i}`, person, async () => s.set('a', i));
    expect(h.steps.value.map((x) => x.label)).toEqual([8, 7, 6, 5, 4, 3].map((i) => `Set a to ${i}`));
    expect(h.steps.value.map((x) => !!x.unsafe)).toEqual([false, false, false, true, true, true]);
    expect(h.steps.value[3].unsafe).toMatch(/Too old to undo: the last 3 steps/);
    for (const a of [7, 6, 5]) {
      expectOk(await h.undo());
      expect(s.get().a).toBe(a);
    }
    expect(await h.undo()).toMatchObject({ ok: false, reason: expect.stringMatching(/Too old/) });
    expect(s.get().a).toBe(5);
    // Redo puts the three back, each undoable again.
    for (let i = 0; i < 3; i++) expectOk(await h.redo());
    expect(s.get().a).toBe(8);
    expect(h.steps.value.map((x) => !!x.unsafe)).toEqual([false, false, false, true, true, true]);
    expectOk(await h.undo());
    expect(s.get().a).toBe(7);
  });

  it('records no step for a group of only work that changes nothing undoable', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    const g = await h.group('Show', script, async (ctx) => {
      await h.run('Apply a snapshot', script, async () => s.set('a', 1), { within: ctx, undoable: false });
    });
    expect(g.runs).toEqual([]);
    expect(h.steps.value).toEqual([]);
    const both = await h.group('Both', script, async (ctx) => {
      await h.run('Apply a snapshot', script, async () => s.set('a', 2), { within: ctx, undoable: false });
      await h.run('Set b', script, async () => s.set('b', 1), { within: ctx });
    });
    expect(both.runs).toEqual(['Set b']);
    expect(h.steps.value.map((x) => x.runs)).toEqual([['Set b']]);
  });

  it('releases the lock when taking a checkpoint throws', async () => {
    let fail = true;
    const s = store({ a: 0 });
    const h = new HistoryModel({
      ...s.port,
      checkpoint: () => {
        if (fail) throw new Error('cannot');
        return { ...s.get() };
      },
    });
    await expect(h.run('x', person, async () => s.set('a', 1))).rejects.toThrow('cannot');
    await expect(h.group('g', script, async () => undefined)).rejects.toThrow('cannot');
    fail = false;
    await h.run('y', person, async () => s.set('a', 2));
    expect(h.steps.value.map((x) => x.label)).toEqual(['y']);
    expect(h.running.value).toBe(0);
  });

  it('keeps at most `capacity` entries, dropping the oldest', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port, { capacity: 2 });
    for (let i = 1; i <= 4; i++) await h.run(`set ${i}`, person, async () => s.set('a', i));
    expect(h.steps.value.map((x) => x.label)).toEqual(['set 4', 'set 3']);
  });

  it('a turn’s runs join its entry in place, and a step once published keeps its runs', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    const agent = { kind: 'agent', turn: 't1' };
    const end = h.beginGroup('Turn', agent);
    await h.run('one', agent, async () => s.set('a', 1));
    const first = h.steps.value[0];
    await h.run('two', agent, async () => s.set('a', 2));
    const second = h.steps.value[0];
    await h.group('snippet', agent, async (ctx) => {
      await h.run('three', agent, async () => s.set('a', 3), { within: ctx });
      await h.run('four', agent, async () => s.set('a', 4), { within: ctx });
    });
    end();
    // Older snapshots never see later runs; the same snapshot reads the same list.
    expect(first.runs).toEqual(['one']);
    expect(second.runs).toEqual(['one', 'two']);
    expect(h.steps.value[0].runs).toEqual(['one', 'two', 'three', 'four']);
    expect(h.steps.value[0].runs).toBe(h.steps.value[0].runs);
    expect({ ...h.steps.value[0] }).toMatchObject({
      id: first.id,
      label: 'Turn',
      runs: ['one', 'two', 'three', 'four'],
    });
    h.invalidate('changed outside');
    expect(h.steps.value[0]).toMatchObject({ unsafe: 'changed outside', runs: ['one', 'two', 'three', 'four'] });
    expect(second.unsafe).toBeUndefined();
  });
});

describe('HistoryModel, to a step', () => {
  it('undoes a step and every later one, and redoes up to a step undone', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    for (let i = 1; i <= 4; i++) await h.run(`Set a to ${i}`, person, async () => s.set('a', i));
    const second = h.steps.value[2];
    expectOk(await h.undoTo(second.id));
    expect(s.get().a).toBe(1);
    expect(h.steps.value.map((x) => x.label)).toEqual(['Set a to 1']);
    expect(h.redoable.value.map((x) => x.label)).toEqual(['Set a to 2', 'Set a to 3', 'Set a to 4']);
    expect(h.canRedo.value).toBe(true);
    expect(h.lastMove.value).toMatchObject({ kind: 'undo', steps: [{ label: 'Set a to 2' }] });
    expectOk(await h.redoTo(h.redoable.value[1].id));
    expect(s.get().a).toBe(3);
    expect(h.redoable.value.map((x) => x.label)).toEqual(['Set a to 4']);
    expect(h.lastMove.value).toMatchObject({ kind: 'redo', steps: [{ label: 'Set a to 3' }] });
    expect(await h.redoTo(12345)).toMatchObject({ ok: false });
    expect(await h.undoTo(12345)).toMatchObject({ ok: false });
  });

  it('stops undoing at a step it cannot undo, and says why', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    await h.run('one', person, async () => s.set('a', 1));
    h.invalidate('Changed elsewhere');
    await h.run('two', person, async () => s.set('a', 2));
    const r = await h.undoTo(h.steps.value[1].id);
    expect(r).toMatchObject({ ok: false, reason: 'Changed elsewhere' });
    expect(s.get().a).toBe(1);
    expect(h.undoBlocked.value).toBe('Changed elsewhere');
  });

  it('stamps steps with the injected clock', async () => {
    const s = store({ a: 0 });
    let time = 1000;
    const h = new HistoryModel(s.port, { now: () => time });
    await h.run('one', person, async () => s.set('a', 1));
    time = 2000;
    await h.run('two', person, async () => s.set('a', 2));
    expect(h.steps.value.map((x) => x.time)).toEqual([2000, 1000]);
  });
});

describe('HistoryModel.moveTo', () => {
  it('moves the position after a step, before every step, and forward again', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    for (let i = 1; i <= 3; i++) await h.run(`Set a to ${i}`, person, async () => s.set('a', i));
    const [three, two, one] = h.steps.value;
    expect(await h.moveTo(three.id)).toMatchObject({ ok: false });
    expectOk(await h.moveTo(one.id));
    expect(s.get().a).toBe(1);
    expectOk(await h.moveTo(undefined));
    expect(s.get().a).toBe(0);
    expect(h.steps.value).toEqual([]);
    expectOk(await h.moveTo(two.id));
    expect(s.get().a).toBe(2);
    expect(h.redoable.value.map((x) => x.id)).toEqual([three.id]);
  });
});

describe('HistoryModel, its words', () => {
  it('refuses with a code and the labels’ text; the defaults name no domain', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port, {
      labels: { nothingToUndo: 'Rien à annuler', laterSteps: (n) => `${n} de plus` },
    });
    expect(h.undoBlocked.value).toBe('Rien à annuler');
    expect(h.undoBlockedCode.value).toBe('nothing-to-undo');
    expect(await h.redo()).toEqual({ ok: false, code: 'nothing-to-redo', reason: 'Nothing to redo' });
    expect(await h.undoTo(42)).toMatchObject({ code: 'not-in-history' });
    await h.run('one', person, async () => s.set('a', 1));
    expect(await h.moveTo(h.steps.value[0].id)).toMatchObject({ code: 'latest' });
    expect(await h.undoTurn('nope')).toMatchObject({ code: 'turn-empty' });
    const words = Object.values(h.labels).map((v) => (typeof v === 'function' ? v(2) : v));
    expect(words.join(' ')).not.toMatch(/scene|viewer|structure/i);
  });
});

describe('HistoryModel, its limits (history §5)', () => {
  it('keeps at most `maxRedo` redo steps, the farthest dropped', async () => {
    const s = store({ a: 0 });
    let taken = 0;
    const h = new HistoryModel({ ...s.port, checkpoint: () => (taken++, s.port.checkpoint()) }, { maxRedo: 2 });
    for (let i = 1; i <= 4; i++) await h.run(`Set a to ${i}`, person, async () => s.set('a', i));
    for (let i = 0; i < 4; i++) expectOk(await h.undo());
    expect(h.redoable.value.map((x) => x.label)).toEqual(['Set a to 1', 'Set a to 2']);
    expectOk(await h.redo());
    expectOk(await h.redo());
    expect(await h.redo()).toMatchObject({ ok: false, code: 'nothing-to-redo' });
    expect(s.get().a).toBe(2);
    // No Redo at all: Undo takes no checkpoint of the state it leaves.
    h.setLimits({ maxRedo: 0 });
    expect(h.redoable.value).toEqual([]);
    const before = taken;
    expectOk(await h.undo());
    expect(taken).toBe(before);
    expect(h.canRedo.value).toBe(false);
  });

  it('lists at most `maxRuns` runs of a step, the latest, and counts them all', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port, { maxRuns: 3 });
    const agent = { kind: 'agent', turn: 't' };
    const end = h.beginGroup('Turn', agent);
    for (let i = 0; i < 8; i++) await h.run(`r${i}`, agent, async () => s.set('a', i));
    end();
    const [step] = h.steps.value;
    expect(step).toMatchObject({ runCount: 8, droppedRuns: 5, runs: ['r5', 'r6', 'r7'] });
    expect(step.lastRuns(2)).toEqual(['r6', 'r7']);
    // The checkpoint is the turn's: Undo still restores before its first run.
    expectOk(await h.undo());
    expect(s.get().a).toBe(0);
    const g = await h.group('Many', script, async (ctx) => {
      for (let i = 0; i < 5; i++) await h.run(`g${i}`, script, async () => s.set('a', i), { within: ctx });
    });
    expect(g.step).toMatchObject({ runCount: 5, droppedRuns: 2, runs: ['g2', 'g3', 'g4'] });
  });

  it('changes its limits at runtime: past a lowered one the oldest go at once', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    expect(h.limits).toEqual({ capacity: 100, checkpoints: 10, maxRedo: 100, maxRuns: 1000 });
    for (let i = 1; i <= 6; i++) await h.run(`Set a to ${i}`, person, async () => s.set('a', i));
    h.setLimits({ capacity: 4, checkpoints: 2 });
    expect(h.steps.value.map((x) => x.label)).toEqual([6, 5, 4, 3].map((i) => `Set a to ${i}`));
    expect(h.steps.value.map((x) => x.unsafeCode)).toEqual([undefined, undefined, 'too-old', 'too-old']);
    expect(h.steps.value[2].unsafe).toBe('Too old to undo: the last 2 steps can be undone');
    // The checkpoints count never passes the capacity.
    h.setLimits({ capacity: 1, checkpoints: 5 });
    expect(h.limits).toMatchObject({ capacity: 1, checkpoints: 1 });
    expect(h.steps.value).toHaveLength(1);
  });
});

describe('HistoryModel, not recording (history §1)', () => {
  it('runs in order but records no step, and makes the steps before it unsafe', async () => {
    const s = store({ a: 0 });
    let taken = 0;
    const h = new HistoryModel({ ...s.port, checkpoint: () => (taken++, s.port.checkpoint()) });
    await h.run('one', person, async () => s.set('a', 1));
    h.setRecording(false);
    expect(h.recording.value).toBe(false);
    const before = taken;
    const out = await h.run('two', person, async () => s.set('a', 2));
    expect(out.step).toBeUndefined();
    expect(taken).toBe(before);
    expect(s.get().a).toBe(2);
    expect(h.steps.value.map((x) => x.label)).toEqual(['one']);
    expect(h.steps.value[0]).toMatchObject({ unsafeCode: 'not-recorded' });
    expect(await h.undo()).toMatchObject({ ok: false, code: 'unsafe' });
    // A group too; work that changed nothing (`record` false, `undoable: false`) leaves the steps as they are.
    h.setRecording(true);
    await h.run('three', person, async () => s.set('a', 3));
    h.setRecording(false);
    await h.run('same', person, async () => 'same', { record: () => false });
    await h.run('upload', person, async () => undefined, { undoable: false });
    expect(h.steps.value[0].unsafe).toBeUndefined();
    await h.group('g', script, async (ctx) => {
      await h.run('x', script, async () => s.set('a', 4), { within: ctx });
    });
    expect(h.steps.value[0].unsafeCode).toBe('not-recorded');
  });

  it('records by origin and label through `recordIf`', async () => {
    const s = store({ a: 0 });
    const asked: [string, string | undefined][] = [];
    const h = new HistoryModel(s.port, {
      recordIf: (origin, label) => {
        asked.push([origin.kind, label]);
        return origin.kind !== 'script';
      },
    });
    await h.run('mine', person, async () => s.set('a', 1));
    await h.run(
      (v: number) => `computed ${v}`,
      person,
      async () => 2,
    );
    await h.run('scripted', script, async () => s.set('a', 3));
    expect(asked).toEqual([
      ['person', 'mine'],
      ['person', undefined],
      ['script', 'scripted'],
    ]);
    expect(h.steps.value.map((x) => x.label)).toEqual(['computed 2', 'mine']);
    expect(h.steps.value.every((x) => x.unsafeCode === 'not-recorded')).toBe(true);
  });
});

describe('HistoryModel, a clear or a change during a run (history §1)', () => {
  it('records a run that a clear overtook as nothing, and restores nothing when it fails', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    await h.run('old', person, async () => s.set('a', 1));
    const g = gate();
    const run = h.run('slow', person, async () => {
      await g.p;
      s.set('a', 2);
    });
    await tick();
    // A new scene while it runs.
    s.set('a', 100);
    h.clear();
    g.open();
    expect((await run).step).toBeUndefined();
    expect(h.steps.value).toEqual([]);
    expect(h.canUndo.value).toBe(false);

    const g2 = gate();
    const failing = h.run(
      'fails',
      person,
      async () => {
        await g2.p;
        return false;
      },
      { failed: (v) => !v },
    );
    await tick();
    s.set('a', 7);
    h.clear();
    g2.open();
    expect(await failing).toMatchObject({ failed: true, restored: false });
    expect(s.get().a).toBe(7);
  });

  it('records a run that an invalidate overtook as unsafe', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel(s.port);
    const g = gate();
    const run = h.run('slow', person, async () => {
      await g.p;
      s.set('a', 2);
    });
    await tick();
    h.invalidate('Changed elsewhere');
    g.open();
    const out = await run;
    expect(out.step).toMatchObject({ unsafe: 'Changed elsewhere', unsafeCode: 'changed' });
    expect(await h.undo()).toMatchObject({ ok: false, code: 'unsafe' });
    expect(s.get().a).toBe(2);
  });

  it('clears what turned checkpoints off; an open turn goes on in a new entry', async () => {
    const s = store({ a: 0 });
    const h = new HistoryModel({ ...s.port, size: () => 10_000 }, { budget: { ms: 1e9, bytes: 1000 } });
    const agent = { kind: 'agent', turn: 't' };
    h.beginGroup('Turn', agent);
    await h.run('one', agent, async () => s.set('a', 1));
    expect(h.checkpointsOff.value).toBeDefined();
    h.clear();
    expect(h.checkpointsOff.value).toBeUndefined();
    expect(h.checkpointing.value).toBe(true);
    await h.run('two', agent, async () => s.set('a', 2));
    expect(h.steps.value).toHaveLength(1);
    expect(h.steps.value[0]).toMatchObject({ label: 'Turn', runs: ['two'] });
  });
});

describe('HistoryModel.undoTurn', () => {
  it('undoes the turn of one origin, not another kind’s with the same turn id', async () => {
    const s = store({ a: 0, b: 0 });
    const h = new HistoryModel(s.port);
    const one = { kind: 'agent', turn: 't1' };
    const other = { kind: 'helper', turn: 't1' };
    const endOne = h.beginGroup('Agent', one);
    await h.run('a', one, async () => s.set('a', 1));
    endOne();
    const endOther = h.beginGroup('Helper', other);
    await h.run('b', other, async () => s.set('b', 1));
    endOther();
    expectOk(await h.undoTurn(other));
    expect(s.get()).toEqual({ a: 1, b: 0 });
    expect(h.steps.value.map((x) => x.label)).toEqual(['Agent']);
  });
});
