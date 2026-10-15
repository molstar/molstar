/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createMemoryRows, type FiberTable } from '@molstar/ai-fibers';
import { describe, expect, it } from 'vitest';
import { createHistoryModel } from '../history-model.js';
import { createKeptRecordsModel } from '../kept-records.js';
import { createRecordIndex, createSessionRecord, type RunEntry } from '../record.js';
import { fakeApi } from './runs.js';

// The kept records after a reload (journal.md §3, §5): listed without the session's own, one opened read-only for
// History, forgotten, and at most `max` of them.

/** A record of an earlier session, kept in `rows`, with `n` changes. */
async function earlier(rows: ReturnType<typeof createMemoryRows<FiberTable>>, title: string, n: number) {
  const { api, run } = fakeApi();
  const record = await createSessionRecord({ api });
  for (let i = 0; i < n; i++) run(`${title} change ${i + 1}`);
  await record.keep(rows, await createRecordIndex(rows), title);
  await record.dispose();
  return record.id;
}

describe('the kept records', () => {
  it("lists the earlier ones, not this session's own, and opens one read-only", async () => {
    const rows = createMemoryRows<FiberTable>();
    const a = await earlier(rows, '1CBS', 2);
    await earlier(rows, '1TQN', 1);
    const { api } = fakeApi();
    const record = await createSessionRecord({ api });
    await record.keep(rows, await createRecordIndex(rows), 'This one');
    const model = createKeptRecordsModel({ record, rows });
    const stop = model.watch();
    await model.refresh();
    expect(model.list.value.map((k) => k.title).sort()).toEqual(['1CBS', '1TQN']);

    expect(await model.open(a)).toBe(true);
    const opened = model.opened.value!;
    expect(opened.record.title).toBe('1CBS');
    expect(opened.reader.kept).toBe(true);
    const runs = (await opened.reader.read(10)).items.map((i) => i.entry).filter((e) => e.kind === 'run');
    expect((runs as RunEntry[]).map((r) => r.result.ok && r.result.did)).toEqual(['1CBS change 2', '1CBS change 1']);
    // History reads it as it reads the session's.
    const history = createHistoryModel(opened.reader, { schedule: (fn) => fn() });
    history.setActive(true);
    await history.idle();
    expect(history.rows.value.map((r) => (r.kind === 'run' ? r.text : r.kind))).toEqual([
      '1CBS change 2',
      '1CBS change 1',
    ]);
    history.dispose();
    model.close();
    expect(model.opened.value).toBeUndefined();
    stop();
    model.dispose();
  });

  it('forgets one (closing it when open), and forgets this session’s kept copy by unkeeping it', async () => {
    const rows = createMemoryRows<FiberTable>();
    const a = await earlier(rows, '1CBS', 1);
    const { api } = fakeApi();
    const record = await createSessionRecord({ api });
    const model = createKeptRecordsModel({ record, rows });
    await model.open(a);
    await model.forget(a);
    expect(model.opened.value).toBeUndefined();
    expect(model.list.value).toEqual([]);
    expect(await (await createRecordIndex(rows)).list()).toEqual([]);

    await record.keep(rows, await createRecordIndex(rows), 'This one');
    await model.forget(record.id);
    expect(record.kept).toBe(false);
    model.dispose();
  });

  it('names the oldest kept record when keeping one more passes the bound', async () => {
    const rows = createMemoryRows<FiberTable>();
    const first = await earlier(rows, 'First', 1);
    await new Promise((r) => setTimeout(r, 5));
    await earlier(rows, 'Second', 1);
    const { api } = fakeApi();
    const record = await createSessionRecord({ api });
    expect(await createKeptRecordsModel({ record, rows, max: 3 }).overBound()).toBeUndefined();
    const model = createKeptRecordsModel({ record, rows, max: 2 });
    expect((await model.overBound())?.id).toBe(first);
    model.dispose();
  });
});
