/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createMemoryRows, type ThreadTable } from '@molstar/ai';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import { createHistoryModel, type HistoryRunRow, type HistoryTurnRow, metaOf } from '../history-model.js';
import {
  changedOf,
  ChangedNodesMax,
  createRecordIndex,
  createSessionRecord,
  IndexFoldAt,
  recordFromLines,
  recordLines,
  type RunEntry,
  targetsOf,
  type TurnEntry,
  type UndoEntry,
} from '../record.js';
import { fakeApi } from './runs.js';

// What the record holds beyond the run and its words (journal.md §2): a run's targets and `changed`, the turn's model
// and streamed reply, who undid a step, one entry per Undo turn, and the kept records' index folded.

async function setup() {
  const v = await headlessViewer({ files: [cbs()] });
  const record = await createSessionRecord({
    api: v.api,
    undone: v.managers.history.undone,
    redone: v.managers.history.redone,
  });
  return { v, record };
}

describe('a run in the record: its targets and what it changed', () => {
  it('says what a run was pointed at, in order, each once and cut', () => {
    expect(targetsOf({ structure: '1CBS', node: '', target: 'A:10-20', language: 'pymol', n: 3 })).toEqual([
      '1CBS',
      'A:10-20',
    ]);
    expect(targetsOf({ structure: 'A', target: 'A' })).toEqual(['A']);
    expect(targetsOf({ target: { loci: {} } })).toEqual(['atoms picked in the view']);
    expect(targetsOf({ target: { molql: {} } })).toEqual(['(a selector)']);
    expect(targetsOf({ target: 'x'.repeat(200) })[0]).toHaveLength(80);
    // No credentials or query of an address.
    expect(targetsOf({ node: 'https://u:t@data.example/a.cif?sig=1' })).toEqual(['https://data.example/a.cif?…']);
    expect(targetsOf({})).toEqual([]);
  });

  it('keeps `changed` small: a few node refs, the rest counted, the flags; nothing when nothing changed', () => {
    const nodes = Array.from({ length: 12 }, (_, i) => `n${i}`);
    expect(changedOf({ nodes, view: true })).toEqual({
      nodes: nodes.slice(0, ChangedNodesMax),
      more: 12 - ChangedNodesMax,
      view: true,
    });
    expect(changedOf({ selection: true, sequence: true })).toEqual({ selection: true, sequence: true });
    expect(changedOf({ nodes: [] })).toBeUndefined();
    expect(changedOf({})).toBeUndefined();
    expect(changedOf(undefined)).toBeUndefined();
  });

  it('writes both into the run entry, and History shows them in a line', async () => {
    const { v, record } = await setup();
    expectOk(await v.api.as({ kind: 'person' }).select.set({ target: 'ligand' }));
    expectOk(await v.api.as({ kind: 'person' }).styles.apply({ style: 'surface' }));
    await record.flush();
    const runs = (await record.read(10)).items.map((i) => i.entry as RunEntry);
    const select = runs.find((r) => r.operation === 'select.set')!;
    expect(select.targets).toEqual(['ligand']);
    expect(select.changed).toMatchObject({ selection: true });
    const style = runs.find((r) => r.operation === 'styles.apply')!;
    expect(style.changed?.nodes?.length).toBeGreaterThan(0);

    const model = createHistoryModel(record);
    model.setActive(true);
    await model.idle();
    const rows = model.rows.value as HistoryRunRow[];
    expect(rows.find((r) => r.operation === 'select.set')?.meta).toMatch(/^ligand · .*selection/);
    expect(rows.find((r) => r.operation === 'styles.apply')?.meta).toMatch(/\d+ nodes?/);
    model.dispose();
    await record.dispose();
  });

  it('writes what a failed run was pointed at, and no changes', async () => {
    const { api } = fakeApi();
    const record = await createSessionRecord({ api });
    (api.events.finished as unknown as { next(r: unknown): void }).next({
      id: 1,
      operation: 'select.set',
      title: 'Select',
      values: { target: 'resn XYZ' },
      origin: { kind: 'person' },
      started: 1,
      result: { ok: false, code: 'not-found', message: 'Nothing matches.' },
      ms: 1,
    });
    await record.flush();
    const run = (await record.read(5)).items[0].entry as RunEntry;
    expect(run.targets).toEqual(['resn XYZ']);
    expect(run.changed).toBeUndefined();
    expect(metaOf(run)).toBe('resn XYZ');
    await record.dispose();
  });

  it('keeps them through a session file', async () => {
    const { v, record } = await setup();
    expectOk(await v.api.as({ kind: 'person' }).select.set({ target: 'ligand' }));
    const t = record.beginTurn('t1', 'hello', undefined, 'anthropic:claude');
    await t.settle({ reply: 'Hi.', outcome: 'done' });
    await record.flush();
    const file = recordFromLines(await recordLines(record));
    const entries = (await file!.reader.read(10)).items.map((i) => i.entry);
    expect(entries.find((e): e is RunEntry => e.kind === 'run')).toMatchObject({
      targets: ['ligand'],
      changed: { selection: true },
    });
    expect(entries.find((e): e is TurnEntry => e.kind === 'turn')).toMatchObject({ model: 'anthropic:claude' });
    await record.dispose();
  });
});

describe('a turn in the record: its model and its reply as it streams', () => {
  /** A clock for the reply's coalescing that the test fires by hand. */
  function manual() {
    const timers = new Set<() => void>();
    return {
      after: (fn: () => void) => {
        timers.add(fn);
        return () => timers.delete(fn);
      },
      fire: () => {
        for (const fn of [...timers]) {
          timers.delete(fn);
          fn();
        }
      },
      get waiting() {
        return timers.size;
      },
    };
  }

  it('has the model from the prompt on, and keeps it when settled', async () => {
    const { api } = fakeApi();
    const record = await createSessionRecord({ api });
    const turn = record.beginTurn('t1', 'hello', undefined, 'anthropic:claude-sonnet');
    await record.flush();
    let item = (await record.read(5)).items[0];
    expect(item.open).toBe(true);
    expect(item.entry).toMatchObject({ model: 'anthropic:claude-sonnet' });
    await turn.settle({ reply: 'Hi.', outcome: 'done' });
    item = (await record.read(5)).items[0];
    expect(item.entry).toMatchObject({ model: 'anthropic:claude-sonnet', reply: 'Hi.' });
    await record.dispose();
  });

  it('shows the reply as it arrives, coalesced, and the final text replaces it', async () => {
    const { api } = fakeApi();
    const m = manual();
    const record = await createSessionRecord({ api, after: m.after });
    const model = createHistoryModel(record, { schedule: (fn) => queueMicrotask(fn) });
    model.setActive(true);
    const turn = record.beginTurn('t1', 'explain');
    await record.flush();
    // Many updates before the timer: one delta with the newest text.
    turn.reply('The');
    turn.reply('The ligand');
    turn.reply('The ligand is HEM');
    expect(m.waiting).toBe(1);
    m.fire();
    await record.flush();
    await model.idle();
    let row = model.rows.value.find((r): r is HistoryTurnRow => r.kind === 'turn')!;
    expect(row).toMatchObject({ open: true, reply: 'The ligand is HEM' });

    // Appended: only the new part is written; the row grows.
    turn.reply('The ligand is HEM, a heme group.');
    m.fire();
    await record.flush();
    await model.idle();
    row = model.rows.value.find((r): r is HistoryTurnRow => r.kind === 'turn')!;
    expect(row.reply).toBe('The ligand is HEM, a heme group.');

    // A text that does not extend what was sent (a retried call) replaces it.
    turn.reply('Retried.');
    m.fire();
    await record.flush();
    await model.idle();
    row = model.rows.value.find((r): r is HistoryTurnRow => r.kind === 'turn')!;
    expect(row.reply).toBe('Retried.');

    // Settled: the final text is the reply, and nothing more is streamed.
    await turn.settle({ reply: 'Final answer.', outcome: 'done' });
    turn.reply('late');
    expect(m.waiting).toBe(0);
    await model.idle();
    row = model.rows.value.find((r): r is HistoryTurnRow => r.kind === 'turn')!;
    expect(row).toMatchObject({ open: false, reply: 'Final answer.' });
    model.dispose();
    await record.dispose();
  });

  it('shows a replaced reply of the same length (a retried call)', async () => {
    const { api } = fakeApi();
    const m = manual();
    const record = await createSessionRecord({ api, after: m.after });
    const model = createHistoryModel(record, { schedule: (fn) => queueMicrotask(fn) });
    model.setActive(true);
    const turn = record.beginTurn('t1', 'explain');
    const shown = async (text: string) => {
      turn.reply(text);
      m.fire();
      await record.flush();
      await model.idle();
      return model.rows.value.find((r): r is HistoryTurnRow => r.kind === 'turn')!.reply;
    };
    expect(await shown('First try.')).toBe('First try.');
    expect(await shown('Other text')).toBe('Other text');
    await turn.settle({ reply: 'Done.', outcome: 'done' });
    model.dispose();
    await record.dispose();
  });

  it('keeps what had streamed when it settles with nothing to say, and writes nothing for a reply too long to preview', async () => {
    const { api } = fakeApi();
    const m = manual();
    const record = await createSessionRecord({ api, after: m.after });
    const a = record.beginTurn('t1', 'a');
    a.reply('Half an answ');
    await a.settle({ reply: '', outcome: 'stopped' });
    const b = record.beginTurn('t2', 'b');
    b.reply('x'.repeat(40_000));
    m.fire();
    await record.flush();
    const items = (await record.read(5)).items;
    expect((items.find((i) => (i.entry as TurnEntry).turn === 't1')!.entry as TurnEntry).reply).toBe('Half an answ');
    // No preview of a very long text: the turn is still open, with no reply until it settles.
    const open = items.find((i) => (i.entry as TurnEntry).turn === 't2')!;
    expect(open.open).toBe(true);
    expect((open.entry as TurnEntry).reply).toBeUndefined();
    await b.settle({ reply: 'Done.', outcome: 'done' });
    await record.dispose();
  });
});

describe('Undo in the record: who, and one entry per Undo turn', () => {
  it('names who undid a step, and who redid it', async () => {
    const { v, record } = await setup();
    expectOk(await v.api.as({ kind: 'person' }).styles.apply({ style: 'surface' }));
    expectOk(await v.api.as({ kind: 'script' }).history.undo());
    expectOk(await v.api.as({ kind: 'person' }).history.redo());
    expectOk(await v.api.as({ kind: 'person' }).history.undo());
    await record.flush();
    const undos = (await record.read(20)).items
      .map((i) => i.entry)
      .filter((e) => e.kind === 'undo' || e.kind === 'redo');
    expect(undos.map((e) => [e.kind, (e as UndoEntry).by?.kind])).toEqual([
      ['undo', 'person'],
      ['redo', 'person'],
      ['undo', 'script'],
    ]);
    await record.dispose();
  });

  it('names the agent and its turn when it undid a turn', async () => {
    const { v, record } = await setup();
    const earlier = v.api.as({ kind: 'agent', turn: 't1' });
    expectOk(await earlier.styles.apply({ style: 'surface' }));
    expectOk(await v.api.as({ kind: 'agent', turn: 't2' }).history.undoTurn({ turn: 't1' }));
    await record.flush();
    const undo = (await record.read(20)).items.map((i) => i.entry).find((e) => e.kind === 'undo') as UndoEntry;
    expect(undo.by).toEqual({ kind: 'agent', turn: 't2' });
    expect(undo.turn).toBe('t1');
    await record.dispose();
  });

  it('writes one entry for an Undo turn, however many steps it took back, and History marks them all', async () => {
    const { v, record } = await setup();
    const agent = v.api.as({ kind: 'agent', turn: 't1' });
    const end = v.managers.history.beginGroup('Assistant: styles', { kind: 'agent', turn: 't1' });
    expectOk(await agent.styles.apply({ style: 'surface' }));
    // The person's step between two of the turn's: the group closes and opens again, so the turn has two steps.
    expectOk(await v.api.as({ kind: 'person' }).styles.apply({ style: 'cartoon' }));
    expectOk(await agent.styles.apply({ style: 'spacefill' }));
    end();
    const steps = v.managers.history.steps.value.filter((s) => s.origin.turn === 't1');
    expect(steps.length).toBeGreaterThan(1);

    const r = await v.api.as({ kind: 'person' }).history.undoTurn({ turn: 't1', discardLater: true });
    expectOk(r);
    await record.flush();
    const items = (await record.read(30)).items;
    const undos = items.map((i) => i.entry).filter((e): e is UndoEntry => e.kind === 'undo');
    expect(undos).toHaveLength(1);
    expect(undos[0].by).toEqual({ kind: 'person' });
    expect(undos[0].steps?.length).toBeGreaterThanOrEqual(steps.length);
    expect(undos[0].steps?.[0]).toBe(undos[0].step);

    const model = createHistoryModel(record);
    model.setActive(true);
    await model.idle();
    const rows = model.rows.value;
    const runRows = rows.filter((x): x is HistoryRunRow | HistoryTurnRow => x.kind === 'run' || x.kind === 'turn');
    expect(runRows.length).toBeGreaterThan(0);
    const agentRuns = rows.flatMap((x) => (x.kind === 'run' && x.origin === 'agent' ? [x] : []));
    expect(agentRuns.every((x) => x.undone && x.undoneBy === 'person')).toBe(true);
    model.dispose();
    await record.dispose();
  });

  it('shows in History who undid a turn', async () => {
    const { v, record } = await setup();
    const model = createHistoryModel(record);
    model.setActive(true);
    const turn = record.beginTurn('t1', 'make it a surface');
    const end = v.managers.history.beginGroup('Assistant', { kind: 'agent', turn: 't1' });
    expectOk(await v.api.as({ kind: 'agent', turn: 't1' }).styles.apply({ style: 'surface' }));
    end();
    await turn.settle({ reply: 'Done.', outcome: 'done' });
    expectOk(await v.api.as({ kind: 'agent', turn: 't2' }).history.undoTurn({ turn: 't1' }));
    await record.flush();
    await model.idle();
    const row = model.rows.value.find((r): r is HistoryTurnRow => r.kind === 'turn')!;
    expect(row).toMatchObject({ undone: true, undoneBy: 'agent' });
    model.dispose();
    await record.dispose();
  });
});

describe('the kept records index: folded', () => {
  const entries = async (rows: ReturnType<typeof createMemoryRows<ThreadTable>>) => {
    const loaded = await rows.load('thread:prism-viewer-records');
    return loaded ? loaded.rows.filter((r) => r.table === 'entry').length : 0;
  };

  it('folds a long index to one entry per record, keeping titles, times and the first time', async () => {
    const rows = createMemoryRows<ThreadTable>();
    const index = await createRecordIndex(rows);
    await index.touch('a', 'First', 100);
    for (let i = 1; i <= IndexFoldAt; i++) await index.touch('b', 'Second', 100 + i);
    await index.touch('a', 'First, renamed', 5000);
    await index.touch('gone', 'Forgotten', 1);
    await index.forget('gone');
    await index.idle();
    const before = await entries(rows);
    expect(before).toBeLessThan(IndexFoldAt);
    expect(before).toBeGreaterThan(0);

    const listed = await index.list();
    await index.idle();
    expect(listed.map((r) => [r.id, r.title, r.first, r.last])).toEqual([
      ['a', 'First, renamed', 100, 5000],
      ['b', 'Second', 100 + 1, 100 + IndexFoldAt],
    ]);
    // Read again by a fresh index (as after a reload): the same, from a thread folded on the way (a few entries).
    const again = await (await createRecordIndex(rows)).list();
    expect(again).toEqual(listed);
    expect(await entries(rows)).toBeLessThan(IndexFoldAt / 2);
  });

  it('keeps growing no more than its records: folds as it is appended to', async () => {
    const rows = createMemoryRows<ThreadTable>();
    const index = await createRecordIndex(rows);
    for (let i = 0; i < IndexFoldAt * 3; i++) {
      await index.touch(`r${i % 2}`, `Record ${i % 2}`, i + 1);
      await index.idle();
    }
    expect(await entries(rows)).toBeLessThanOrEqual(IndexFoldAt + 1);
    const listed = await index.list();
    expect(listed.map((r) => r.id).sort()).toEqual(['r0', 'r1']);
    expect(listed.find((r) => r.id === 'r0')!.first).toBe(1);
  });

  it('never drops a record: many records past the fold point stay, and are not refolded at every read', async () => {
    const rows = createMemoryRows<ThreadTable>();
    const index = await createRecordIndex(rows);
    const n = IndexFoldAt + 10;
    for (let i = 0; i < n; i++) await index.touch(`r${i}`, `Record ${i}`, i + 1);
    await index.idle();
    expect((await index.list()).length).toBe(n);
    await index.idle();
    const size = await entries(rows);
    expect((await index.list()).length).toBe(n);
    await index.idle();
    expect(await entries(rows)).toBe(size);
  });
});
