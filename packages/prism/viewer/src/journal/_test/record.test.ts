/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createMemoryRows, type FiberTable } from '@molstar/ai-fibers';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import { createRecordIndex, createSessionRecord, type RunEntry, summarizeValues } from '../record.js';

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
    const record = await createSessionRecord({ api: v.api });
    const turn = record.beginTurn('t1', 'make it a surface');
    expectOk(await v.api.as({ kind: 'agent', turn: 't1' }).styles.apply({ style: 'surface' }));
    await record.flush();
    expect((await record.read(10)).items.find((i) => i.entry.kind === 'turn')?.open).toBe(true);
    await turn.settle({ reply: 'Done.', outcome: 'done' });
    expectOk(await v.api.as({ kind: 'person' }).history.undo());
    await record.flush();
    const kinds = (await record.read(10)).items.map((i) => i.entry.kind);
    expect(kinds).toEqual(['undo', 'run', 'turn']);
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
});
