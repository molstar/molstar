/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import 'fake-indexeddb/auto';
import { createIndexedDbRows } from '@molstar/ai-components';
import { createFiberHost, createMemoryRows, type FiberTable, fiberId } from '@molstar/ai-fibers';
import { describe, expect, it, vi } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import {
  createRecordIndex,
  createSessionRecord,
  type RecordIndex,
  type RunEntry,
  summarizeValues,
  type TurnEntry,
} from '../record.js';
import { fakeApi } from './runs.js';

// The session's record (journal.md, plan 4.3).

describe('the session record', () => {
  it('records each change in order, whoever made it, and not the reads', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const record = await createSessionRecord({ api: v.api });
    const person = v.api.as({ kind: 'person' });
    expectOk(await person.styles.apply({ style: 'surface' }));
    expectOk(await v.api.structure.get());
    expectOk(await v.api.as({ kind: 'agent', turn: 't1' }).styles.apply({ style: 'cartoon' }));
    await record.flush();
    const { items } = await record.read(10);
    const runs = items.map((i) => i.entry).filter((e): e is RunEntry => e.kind === 'run');
    expect(runs.map((r) => [r.operation, r.origin.kind])).toEqual([
      ['styles.apply', 'agent'],
      ['styles.apply', 'person'],
    ]);
    expect(runs[1].result).toMatchObject({ ok: true });
    record.dispose();
  });

  it('records an undo as an entry naming the step, and a turn begun then settled', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const record = await createSessionRecord({ api: v.api, undone: v.managers.history.undone });
    const turn = record.beginTurn('t1', 'make it a surface');
    expectOk(await v.api.as({ kind: 'agent', turn: 't1' }).styles.apply({ style: 'surface' }));
    await record.flush();
    expect((await record.read(10)).items.find((i) => i.entry.kind === 'turn')?.open).toBe(true);
    await turn.settle({ reply: 'Done.', outcome: 'done' });
    expectOk(await v.api.as({ kind: 'person' }).history.undo());
    await record.flush();
    const items = (await record.read(10)).items;
    expect(items.map((i) => i.entry.kind)).toEqual(['undo', 'run', 'turn']);
    // The undo names the run's step.
    const run = items[1].entry as RunEntry;
    expect(items[0].entry).toMatchObject({ kind: 'undo', step: run.step, turn: 't1' });
    const turnItem = (await record.read(10)).items.find((i) => i.entry.kind === 'turn')!;
    expect(turnItem.open).toBe(false);
    expect(turnItem.entry).toMatchObject({ reply: 'Done.', outcome: 'done' });
    record.dispose();
  });

  it('is kept on purpose: copied in order and listed in the index, and goes on there', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const record = await createSessionRecord({ api: v.api });
    expectOk(await v.api.styles.apply({ style: 'surface' }));
    const rows = createMemoryRows<FiberTable>();
    const index = await createRecordIndex(rows);
    await record.keep(rows, index, '1CBS');
    expect(record.kept).toBe(true);
    expectOk(await v.api.styles.apply({ style: 'cartoon' }));
    await record.flush();
    expect((await record.read(10)).items).toHaveLength(2);
    const listed = await index.list();
    expect(listed.map((r) => [r.id, r.title])).toEqual([[record.id, '1CBS']]);
    await index.forget(record.id);
    expect(await index.list()).toEqual([]);
    record.dispose();
  });

  it('keeps values as small strict JSON', () => {
    expect(summarizeValues({ target: { loci: {} }, n: Number.NaN, s: 'x'.repeat(300) })).toEqual({
      target: '(atoms picked in the view)',
      n: 'NaN',
      s: `${'x'.repeat(199)}…`,
    });
  });

  it('says functions in words and keeps no query of a URL', () => {
    expect(
      summarizeValues({
        f: () => 'secret source',
        url: 'https://data.example/1cbs.cif?X-Amz-Signature=abc&token=t',
        frag: 'https://data.example/a#key=1',
        plain: 'a?b',
      }),
    ).toEqual({
      f: '(function)',
      url: 'https://data.example/1cbs.cif?…',
      frag: 'https://data.example/a#…',
      plain: 'a?b',
    });
  });
});

/** A clock the test moves, and a heartbeat it ticks by hand. */
function controlled() {
  let t = 1_000_000;
  const beats: (() => void)[] = [];
  let stopped = 0;
  return {
    clock: () => t,
    advance: (ms: number) => {
      t += ms;
    },
    every: (fn: () => void) => {
      beats.push(fn);
      return () => {
        stopped++;
        beats.splice(beats.indexOf(fn), 1);
      };
    },
    tick: () => {
      for (const b of beats) b();
    },
    get beating() {
      return beats.length;
    },
    get stopped() {
      return stopped;
    },
  };
}

const minutes = (n: number) => n * 60_000;
const turnOf = async (record: Awaited<ReturnType<typeof createSessionRecord>>) => {
  await record.flush();
  return (await record.read(50)).items.find((i) => i.entry.kind === 'turn')!;
};

describe('the session record at G6', () => {
  it('keeps a turn under way alive while it waits, and stops the heartbeat when it settles', async () => {
    const { api } = fakeApi();
    const c = controlled();
    const record = await createSessionRecord({ api, clock: c.clock, every: c.every });
    const turn = record.beginTurn('t1', 'remove the ligand');
    await record.flush();
    expect(c.beating).toBe(1);
    // Parked for three hours (the fiber aborts after one hour without activity), alive by the heartbeat.
    for (let i = 0; i < 6; i++) {
      c.advance(minutes(30));
      c.tick();
      await record.flush();
    }
    expect((await turnOf(record)).open).toBe(true);
    await turn.settle({ reply: 'Removed.', outcome: 'done' });
    expect(c.beating).toBe(0);
    expect((await turnOf(record)).entry).toMatchObject({ outcome: 'done', reply: 'Removed.' });
    await record.dispose();
  });

  it('reads a turn the fiber aborted as stopped, not as gone', async () => {
    const { api, run } = fakeApi();
    const c = controlled();
    const record = await createSessionRecord({ api, clock: c.clock, every: c.every });
    const error = vi.spyOn(console, 'error').mockImplementation(() => undefined);
    record.beginTurn('t1', 'remove the ligand');
    run('Agent change', { kind: 'agent', turn: 't1' });
    await record.flush();
    // No heartbeat: past the inactivity timeout.
    c.advance(minutes(61));
    const item = await turnOf(record);
    expect(item.open).toBe(false);
    expect(item.entry).toMatchObject({ kind: 'turn', turn: 't1', prompt: 'remove the ligand', outcome: 'stopped' });
    await record.dispose();
    error.mockRestore();
  });

  it('settles a turn still open when it is disposed, flushes, and forgets its listeners', async () => {
    const { api, run } = fakeApi();
    const record = await createSessionRecord({ api });
    const heard = vi.fn();
    record.onChange(heard);
    record.beginTurn('t1', 'something long');
    await record.dispose();
    heard.mockClear();
    run('After the end');
    await record.flush();
    expect(heard).not.toHaveBeenCalled();
    const { items } = await record.read(10);
    expect(items.map((i) => i.entry.kind)).toEqual(['turn']);
    expect(items[0]).toMatchObject({ open: false, entry: { outcome: 'stopped' } });
  });

  it('goes on writing when a listener throws', async () => {
    const { api, run } = fakeApi();
    const record = await createSessionRecord({ api });
    const error = vi.spyOn(console, 'error').mockImplementation(() => undefined);
    record.onChange(() => {
      throw new Error('a broken view');
    });
    const after = vi.fn();
    record.onChange(after);
    run('One');
    run('Two');
    await record.flush();
    expect((await record.read(10)).items).toHaveLength(2);
    expect(after).toHaveBeenCalledTimes(2);
    await record.dispose();
    error.mockRestore();
  });

  it('reads more than the fiber gives in one window, in windows', async () => {
    const { api, run } = fakeApi();
    const record = await createSessionRecord({ api });
    for (let i = 0; i < 600; i++) run();
    const page = await record.read(500);
    expect(page.items).toHaveLength(500);
    expect(page).toMatchObject({ more: true, start: 100, end: 600 });
    const rest = await record.read(500, page.start);
    expect(rest.items).toHaveLength(100);
    expect(rest.more).toBe(false);
    await record.dispose();
  });

  it('reads what came after a position, newest first, in windows', async () => {
    const { api, run } = fakeApi();
    const record = await createSessionRecord({ api });
    for (let i = 0; i < 600; i++) run();
    const after = await record.readAfter(590, 50);
    expect(after.items.map((i) => i.position)).toEqual([599, 598, 597, 596, 595, 594, 593, 592, 591, 590]);
    expect(after).toMatchObject({ more: false, start: 590, end: 600 });
    // Past one window: the oldest `limit`, and later ones left out.
    const long = await record.readAfter(100, 300);
    expect(long.items).toHaveLength(300);
    expect(long).toMatchObject({ more: true, start: 100, end: 400 });
    expect(long.items[0].position).toBe(399);
    // Nothing new.
    expect(await record.readAfter(600, 10)).toMatchObject({ items: [], more: false, start: 600, end: 600 });
    await record.dispose();
  });

  it('keeps in the order of writes: what comes while it copies goes to the kept record; one keep at a time', async () => {
    const { api, run } = fakeApi();
    const c = controlled();
    const record = await createSessionRecord({ api, clock: c.clock, every: c.every });
    for (let i = 0; i < 300; i++) run();
    const open = record.beginTurn('t1', 'under way at the keep');
    const rows = createMemoryRows<FiberTable>();
    const index = await createRecordIndex(rows);
    const touches = vi.spyOn(index, 'touch');
    const keeping = record.keep(rows, index, '1CBS');
    expect(record.keep(rows, index, '1CBS')).toBe(keeping);
    // Written while the copy runs.
    run('During the keep');
    const late = record.beginTurn('t2', 'begun during the keep');
    await keeping;
    await late.settle({ reply: 'Done.', outcome: 'done' });
    await open.settle({ reply: 'Done too.', outcome: 'done' });
    expect(record.kept).toBe(true);
    const all = await record.read(400);
    expect(all.items).toHaveLength(303);
    expect(all.items.filter((i) => i.open)).toEqual([]);
    expect(all.items.map((i) => i.entry.kind === 'run' && i.entry.result.ok && i.entry.result.did)).toContain(
      'During the keep',
    );

    // The index: on keep, at a turn's end, and at most once a minute while appended to.
    touches.mockClear();
    for (let i = 0; i < 20; i++) run();
    await record.flush();
    expect(touches).not.toHaveBeenCalled();
    c.advance(minutes(2));
    run();
    await record.flush();
    expect(touches).toHaveBeenCalledTimes(1);
    await record.dispose();
  });

  it('is read back from IndexedDB after a reload: its entries, a turn under way at the reload as stopped, and the index', async () => {
    const name = `prism-viewer-records-test-${Math.random().toString(36).slice(2)}`;
    const { api, run } = fakeApi();
    const record = await createSessionRecord({ api });
    run('Before the keep');
    const rows = createIndexedDbRows<FiberTable>(name);
    const index: RecordIndex = await createRecordIndex(rows);
    await record.keep(rows, index, '1CBS');
    run('After the keep');
    const turn = record.beginTurn('t1', 'remove the ligand');
    await turn.settle({ reply: 'Removed.', outcome: 'done' });
    record.beginTurn('t2', 'still under way');
    await record.flush();
    rows.close();

    // A new page: the same database, opened again, read with no record of its own (as Home would).
    const reopened = createIndexedDbRows<FiberTable>(name);
    const listed = await (await createRecordIndex(reopened)).list();
    expect(listed.map((r) => [r.id, r.title])).toEqual([[record.id, '1CBS']]);
    const host = createFiberHost({ rows: reopened, clock: { now: () => Date.now() + minutes(61) } });
    const page = await host.window({ fiberId: fiberId(record.id), limit: 50, direction: 'forward' });
    if ('status' in page) throw new Error(page.reason);
    const values = page.entries.map((e) =>
      e.state.kind === 'settled'
        ? e.state.payload.value
        : { ...(e.descriptor.value as TurnEntry), state: e.state.kind },
    );
    expect(values.map((v) => (v as { kind: string }).kind)).toEqual(['run', 'run', 'turn', 'turn']);
    expect(values[1]).toMatchObject({ kind: 'run', result: { did: 'After the keep' } });
    expect(values[2]).toMatchObject({ turn: 't1', outcome: 'done', reply: 'Removed.' });
    // The page that wrote it went away mid-turn: past the inactivity timeout, the fiber has aborted it.
    expect(values[3]).toMatchObject({ turn: 't2', state: 'aborted' });
    reopened.close();
    await record.dispose();
  });
});
