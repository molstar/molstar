/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createMemoryRows, type ThreadRows, type ThreadTable } from '@molstar/ai';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { createRecordIndex, createSessionRecord, IndexFoldAt, type RecordRows } from '../record.js';
import { fakeApi } from './runs.js';

// The kept records' index fenced per database (journal.md §1): one viewer makes several indexes on one database (Keep,
// Home, Storage), and unkeyed viewers share one; a fold rewrites the thread, and a page can close in the middle of it.

const IndexScope = 'thread:prism-viewer-records';

/**
 * A handle to shared rows, as each `createIndexedDbRows(database)` is: `pause()` holds the next drop of the index thread
 * right after it is done (a fold's gap: the thread gone, not yet written again); `crash()` makes every write after that
 * drop fail (the page closed mid-fold).
 */
function handle(base: ThreadRows, database = 'db') {
  let hold: { reached: () => void; gate: Promise<void> } | undefined;
  let crashing = false;
  let crashed = false;
  const rows: RecordRows = {
    database,
    load: (scope) => base.load(scope),
    generation: (scope) => base.generation(scope),
    apply: async (scope, expected, writes) => {
      if (crashed) throw new Error('the page is gone');
      return base.apply(scope, expected, writes);
    },
    drop: async (scope) => {
      if (crashed) throw new Error('the page is gone');
      await base.drop(scope);
      if (scope !== IndexScope) return;
      if (crashing) crashed = true;
      const h = hold;
      hold = undefined;
      if (h) {
        h.reached();
        await h.gate;
      }
    },
  };
  return {
    rows,
    /** Resolves once the gap is reached; the returned function ends it. */
    pause() {
      let release!: () => void;
      let reached!: () => void;
      const at = new Promise<void>((r) => {
        reached = r;
      });
      hold = {
        reached,
        gate: new Promise<void>((r) => {
          release = r;
        }),
      };
      return { at, release: () => release() };
    },
    crash() {
      crashing = true;
    },
  };
}

/** An index past its fold point: records a, b, c, and enough touches that the next list folds. */
async function longIndex(rows: RecordRows) {
  const index = await createRecordIndex(rows);
  await index.touch('a', 'A', 1);
  await index.touch('b', 'B', 2);
  await index.touch('c', 'C', 3);
  for (let i = 0; i < IndexFoldAt - 3; i++) await index.touch('a', 'A', 10 + i);
  await index.idle();
  return index;
}

const ids = (list: readonly { id: string }[]) => list.map((r) => r.id).sort();

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('the kept records index, fenced per database', () => {
  it('another index on the database lists everything while a fold rewrites the thread', async () => {
    const base = createMemoryRows<ThreadTable>();
    const one = handle(base);
    const index = await longIndex(one.rows);
    const other = await createRecordIndex(handle(base).rows);

    const paused = one.pause();
    await index.touch('a', 'A', 999); // past the fold point: a fold starts
    await paused.at;
    const listing = other.list();
    let done = false;
    void listing.then(() => {
      done = true;
    });
    await new Promise((r) => setTimeout(r, 10));
    expect(done).toBe(false); // waits for the fold
    paused.release();
    expect(ids(await listing)).toEqual(['a', 'b', 'c']);
    await index.idle();
  });

  it('a forget from another index during a fold is neither lost nor undone by it', async () => {
    const base = createMemoryRows<ThreadTable>();
    const one = handle(base);
    const index = await longIndex(one.rows);
    const other = await createRecordIndex(handle(base).rows);

    const paused = one.pause();
    await index.touch('a', 'A', 999);
    await paused.at;
    const forgetting = other.forget('b');
    paused.release();
    await forgetting;
    await index.idle();
    expect(ids(await other.list())).toEqual(['a', 'c']);
    expect(ids(await (await createRecordIndex(handle(base).rows)).list())).toEqual(['a', 'c']);
  });

  it('a page closed mid-fold leaves the index whole: the next use finishes the fold', async () => {
    const base = createMemoryRows<ThreadTable>();
    const one = handle(base);
    const index = await longIndex(one.rows);
    one.crash();
    const err = vi.spyOn(console, 'error').mockImplementation(() => undefined);
    await index.touch('a', 'A', 999).catch(() => undefined);
    await index.idle();
    err.mockRestore();
    // A new page on the same database.
    const again = await createRecordIndex(handle(base).rows);
    const listed = await again.list();
    expect(ids(listed)).toEqual(['a', 'b', 'c']);
    expect(listed.find((r) => r.id === 'a')).toMatchObject({ first: 1 });
    // The fold's copy is gone once finished.
    expect(await base.load('thread:prism-viewer-records.fold')).toBeUndefined();
  });

  it('takes a Web Lock per database where the browser has them', async () => {
    const names: string[] = [];
    let tail: Promise<unknown> = Promise.resolve();
    vi.stubGlobal('navigator', {
      locks: {
        request<T>(name: string, work: () => Promise<T>): Promise<T> {
          names.push(name);
          const run = tail.then(work, work);
          tail = run.catch(() => undefined);
          return run;
        },
      },
    });
    const index = await createRecordIndex(handle(createMemoryRows<ThreadTable>(), 'prism-default-records').rows);
    await index.touch('a', 'A', 1);
    await index.list();
    expect(names.length).toBeGreaterThanOrEqual(3);
    expect(new Set(names)).toEqual(new Set(['prism-viewer-records-index:prism-default-records']));
  });

  it("lets go of the kept rows only after its index's fold is done", async () => {
    const base = createMemoryRows<ThreadTable>();
    const one = handle(base);
    const index = await createRecordIndex(one.rows);
    for (let i = 0; i < IndexFoldAt; i++) await index.touch('old', 'Old', i + 1);
    await index.idle();
    const release = vi.fn();
    const { api } = fakeApi();
    const record = await createSessionRecord({ api });
    const paused = one.pause();
    await record.keep(one.rows, index, 'This one', release); // its touch passes the fold point
    await paused.at;
    const disposing = record.dispose();
    await new Promise((r) => setTimeout(r, 10));
    expect(release).not.toHaveBeenCalled();
    paused.release();
    await disposing;
    expect(release).toHaveBeenCalledTimes(1);
    expect(ids(await index.list())).toEqual(['old', record.id].sort());
  });
});
