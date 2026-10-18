/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import { createHistoryModel, type HistoryRunRow, type HistoryTurnRow } from '../history-model.js';
import { createSessionRecord } from '../record.js';
import { fakeApi } from './runs.js';

// What History shows of the record (journal.md §3).

async function setup() {
  const v = await headlessViewer({ files: [cbs()] });
  const record = await createSessionRecord({
    api: v.api,
    undone: v.managers.history.undone,
    redone: v.managers.history.redone,
  });
  const model = createHistoryModel(record, { window: 200 });
  return { v, record, model };
}

describe('the History model', () => {
  it('shows a person run as one row, newest first, and marks it undone', async () => {
    const { v, record, model } = await setup();
    expectOk(await v.api.as({ kind: 'person' }).styles.apply({ style: 'spacefill' }));
    expectOk(await v.api.as({ kind: 'person' }).styles.apply({ style: 'cartoon' }));
    await record.flush();
    model.setActive(true);
    await model.idle();
    const rows = model.rows.value as HistoryRunRow[];
    expect(rows.map((r) => [r.kind, r.origin, r.undone, r.failed])).toEqual([
      ['run', 'person', false, false],
      ['run', 'person', false, false],
    ]);
    expect(rows[0].text).toMatch(/cartoon/i);

    expectOk(await v.api.history.undo());
    await record.flush();
    await model.idle();
    const after = model.rows.value as HistoryRunRow[];
    // The undo is no row; it marks the latest run.
    expect(after.map((r) => r.undone)).toEqual([true, false]);

    // Redo puts it back: done again (the redo itself is a run row); undone once more after another Undo.
    const cartoon = () =>
      (model.rows.value as HistoryRunRow[]).find((r) => r.kind === 'run' && r.text.startsWith('Applied the Cartoon'));
    expectOk(await v.api.history.redo());
    await record.flush();
    await model.idle();
    expect(cartoon()?.undone).toBe(false);
    expectOk(await v.api.history.undo());
    await record.flush();
    await model.idle();
    expect(cartoon()?.undone).toBe(true);
    model.dispose();
    record.dispose();
  });

  it('marks a failed run, with the failure message', async () => {
    const { v, record, model } = await setup();
    const r = await v.api.as({ kind: 'person' }).styles.apply({ style: 'no-such-style' } as never);
    expect(r.ok).toBe(false);
    await record.flush();
    model.setActive(true);
    await model.idle();
    const row = model.rows.value[0] as HistoryRunRow;
    expect(row).toMatchObject({ kind: 'run', failed: true, undone: false });
    expect(row.text).toBe(r.ok ? '' : r.message);
    model.dispose();
    record.dispose();
  });

  it('groups an agent turn: prompt, reply, its runs under it, undone with its steps', async () => {
    const { v, record, model } = await setup();
    model.setActive(true);
    const turn = record.beginTurn('t1', 'make it a surface');
    const agent = v.api.as({ kind: 'agent', turn: 't1' });
    expectOk(await agent.styles.apply({ style: 'spacefill' }));
    await record.flush();
    await model.idle();
    let rows = model.rows.value;
    const open = rows.find((r) => r.kind === 'turn') as HistoryTurnRow;
    expect(open).toMatchObject({ prompt: 'make it a surface', open: true, undone: false });
    expect(open.runs).toHaveLength(1);
    // The turn's run is under the turn, not a row of its own.
    expect(rows.filter((r) => r.kind === 'run' && r.origin === 'agent')).toHaveLength(0);

    await turn.settle({ reply: 'Done.\nIt is a surface.', outcome: 'done' });
    await record.flush();
    await model.idle();
    const settled = model.rows.value.find((r) => r.kind === 'turn') as HistoryTurnRow;
    expect(settled).toMatchObject({ open: false, reply: 'Done.\nIt is a surface.', outcome: 'done', undone: false });

    expectOk(await v.api.history.undoTurn({ turn: 't1' }));
    await record.flush();
    await model.idle();
    rows = model.rows.value;
    const undone = rows.find((r) => r.kind === 'turn') as HistoryTurnRow;
    expect(undone.undone).toBe(true);
    expect(undone.runs.every((r) => r.undone)).toBe(true);
    model.dispose();
    record.dispose();
  });

  it('keeps a row object while it does not change', async () => {
    const { v, record, model } = await setup();
    model.setActive(true);
    expectOk(await v.api.as({ kind: 'person' }).styles.apply({ style: 'spacefill' }));
    await record.flush();
    await model.idle();
    const first = model.rows.value[0];
    expectOk(await v.api.as({ kind: 'person' }).styles.apply({ style: 'cartoon' }));
    await record.flush();
    await model.idle();
    expect(model.rows.value[1]).toBe(first);
    model.dispose();
    record.dispose();
  });

  it('shows session events as separator rows', async () => {
    const { v, record, model } = await setup();
    expectOk(await v.api.as({ kind: 'person' }).session.clear());
    await record.flush();
    model.setActive(true);
    await model.idle();
    expect(model.rows.value[0]).toMatchObject({ kind: 'session', text: expect.stringMatching(/Cleared/) });
    model.dispose();
    record.dispose();
  });

  it('keeps a window loaded while it shows: new steps push the oldest out, "Show earlier" brings them back', async () => {
    const { api, run } = fakeApi();
    const record = await createSessionRecord({ api });
    const model = createHistoryModel(record, { window: 3, schedule: queueMicrotask });
    model.setActive(true);
    await model.idle();
    for (let i = 1; i <= 8; i++) {
      run(`Change ${i}`);
      await record.flush();
      await model.idle();
    }
    const texts = () => model.rows.value.map((r) => (r.kind === 'run' ? r.text : r.kind));
    expect(texts()).toEqual(['Change 8', 'Change 7', 'Change 6']);
    expect(model.more.value).toBe(true);
    await model.showEarlier();
    expect(texts()).toEqual(['Change 8', 'Change 7', 'Change 6', 'Change 5', 'Change 4', 'Change 3']);
    model.dispose();
    record.dispose();
  });

  it('reads in windows: the newest first, earlier ones by "Show earlier"', async () => {
    const { v, record, model: _ } = await setup();
    _.dispose();
    const model = createHistoryModel(record, { window: 2 });
    const person = v.api.as({ kind: 'person' });
    for (const style of ['spacefill', 'cartoon', 'spacefill', 'cartoon', 'spacefill'] as const) {
      expectOk(await person.styles.apply({ style }));
    }
    await record.flush();
    model.setActive(true);
    await model.idle();
    expect(model.rows.value).toHaveLength(2);
    expect(model.more.value).toBe(true);
    await model.showEarlier();
    expect(model.rows.value).toHaveLength(4);
    expect(model.more.value).toBe(true);
    await model.showEarlier();
    expect(model.rows.value).toHaveLength(5);
    expect(model.more.value).toBe(false);
    // A new run joins at the top without dropping the earlier windows.
    expectOk(await person.styles.apply({ style: 'cartoon' }));
    await record.flush();
    await model.idle();
    expect(model.rows.value).toHaveLength(6);
    model.dispose();
    record.dispose();
  });

  it('coalesces changes to one read per frame, and reads nothing while inactive', async () => {
    const { v, record, model } = await setup();
    model.dispose();
    const frames: (() => void)[] = [];
    const m = createHistoryModel(record, { schedule: (fn) => frames.push(fn) });
    const read = vi.spyOn(record, 'read');
    const readAfter = vi.spyOn(record, 'readAfter');
    const person = v.api.as({ kind: 'person' });

    // Inactive: writes read nothing, and nothing is subscribed.
    expectOk(await person.styles.apply({ style: 'spacefill' }));
    await record.flush();
    expect(read).not.toHaveBeenCalled();
    expect(frames).toHaveLength(0);

    m.setActive(true);
    await m.idle();
    expect(read).toHaveBeenCalledTimes(1);
    expect(m.rows.value).toHaveLength(1);

    // Three writes before the frame: one read, of what is new only.
    read.mockClear();
    expectOk(await person.styles.apply({ style: 'cartoon' }));
    expectOk(await person.styles.apply({ style: 'spacefill' }));
    expectOk(await person.styles.apply({ style: 'cartoon' }));
    await record.flush();
    expect(frames).toHaveLength(1);
    expect(read).not.toHaveBeenCalled();
    frames.shift()!();
    await m.idle();
    expect(read).not.toHaveBeenCalled();
    expect(readAfter).toHaveBeenCalledTimes(1);
    expect((await readAfter.mock.results[0].value).items).toHaveLength(3);
    expect(m.rows.value).toHaveLength(4);

    // Inactive again: no more reads.
    m.setActive(false);
    read.mockClear();
    readAfter.mockClear();
    expectOk(await person.styles.apply({ style: 'spacefill' }));
    await record.flush();
    expect(frames).toHaveLength(0);
    expect(read).not.toHaveBeenCalled();
    expect(readAfter).not.toHaveBeenCalled();
    m.dispose();
    record.dispose();
  });

  it("shows earlier past the thread's window with a turn under way, and its settling", async () => {
    const { api, run } = fakeApi();
    const record = await createSessionRecord({ api });
    const model = createHistoryModel(record, { window: 200, schedule: queueMicrotask });
    const turn = record.beginTurn('t1', 'a long turn');
    for (let i = 0; i < 300; i++) run();
    await record.flush();
    model.setActive(true);
    await model.idle();
    expect(model.rows.value).toHaveLength(200);
    await model.showEarlier();
    // All of it: 300 runs and the turn.
    expect(model.rows.value).toHaveLength(301);
    expect(model.more.value).toBe(false);
    // A change while the turn, older than the newest window, is open: everything loaded is read again (over 256).
    run();
    await record.flush();
    await model.idle();
    expect(model.rows.value).toHaveLength(302);
    await turn.settle({ reply: 'Done.', outcome: 'done' });
    await model.idle();
    const last = model.rows.value[model.rows.value.length - 1] as HistoryTurnRow;
    expect(last).toMatchObject({ kind: 'turn', open: false, outcome: 'done' });
    expect(new Set(model.rows.value.map((r) => r.key)).size).toBe(302);
    model.dispose();
    await record.dispose();
  });

  it('keeps "Show earlier" in order with a read that starts over', async () => {
    const { api, run } = fakeApi();
    const record = await createSessionRecord({ api });
    const model = createHistoryModel(record, { window: 5, schedule: queueMicrotask });
    for (let i = 0; i < 12; i++) run();
    await record.flush();
    model.setActive(true);
    await model.idle();
    // More new entries than a window (a read that starts over), and "Show earlier", at once.
    for (let i = 0; i < 8; i++) run();
    const earlier = model.showEarlier();
    await record.flush();
    await earlier;
    await model.idle();
    const keys = model.rows.value.map((r) => r.key);
    expect(new Set(keys).size).toBe(keys.length);
    model.dispose();
    await record.dispose();
  });

  it('shows a turn the thread aborted as stopped, its runs under it', async () => {
    const { api, run } = fakeApi();
    let now = 1_000_000;
    const record = await createSessionRecord({ api, clock: () => now, every: () => () => undefined });
    const error = vi.spyOn(console, 'error').mockImplementation(() => undefined);
    record.beginTurn('t1', 'a parked turn');
    run('Its change', { kind: 'agent', turn: 't1' });
    await record.flush();
    now += 61 * 60_000;
    const model = createHistoryModel(record, { schedule: queueMicrotask });
    model.setActive(true);
    await model.idle();
    const [row] = model.rows.value as HistoryTurnRow[];
    expect(row).toMatchObject({ kind: 'turn', open: false, outcome: 'stopped', prompt: 'a parked turn' });
    expect(row.runs.map((r) => r.text)).toEqual(['Its change']);
    model.dispose();
    await record.dispose();
    error.mockRestore();
  });
});
