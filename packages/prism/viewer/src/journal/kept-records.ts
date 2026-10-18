/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ThreadRows } from '@molstar/ai';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import {
  createRecordIndex,
  forgetKept,
  type KeptRecord,
  openKeptRecord,
  openRecordsRows,
  type RecordIndex,
  type RecordReader,
  type SessionRecord,
} from './record.js';

// The kept records after a reload (journal.md §3, §5), DOM-free: Home lists them (the session's own aside), one opens
// read-only in History beside whatever is loaded, and each can be forgotten. At most `max` are kept: keeping one more
// offers to forget the oldest. The rows are the spec's `persist.rows`, else the browser's IndexedDB, opened only when its
// database exists (listing never makes an empty one) and only while something lists or opens.

/** Kept records at most, by default (journal.md §5: a setting). */
export const KeptRecordsMax = 20;

export interface OpenedRecord {
  readonly record: KeptRecord;
  readonly reader: RecordReader;
  /** From a session file (journal.md §5): read-only too, but not kept, so nothing to forget. */
  readonly fromFile?: true;
}

export interface KeptRecordsModel {
  /** The kept records but this session's, newest first; read while watched. */
  readonly list: ReadonlyState<readonly KeptRecord[]>;
  /** The one opened in History, read-only. */
  readonly opened: ReadonlyState<OpenedRecord | undefined>;
  readonly max: number;
  /** Lists now and after every keep or forget, until the returned function is called. */
  watch(): () => void;
  refresh(): Promise<void>;
  open(id: string): Promise<boolean>;
  /** A record read from a session file, shown in History as a kept one is (journal.md §5). */
  show(record: KeptRecord, reader: RecordReader): void;
  close(): void;
  forget(id: string): Promise<void>;
  /** Keeping this session's record would pass the bound: the oldest kept one, which would be forgotten. */
  overBound(): Promise<KeptRecord | undefined>;
  dispose(): void;
}

export interface KeptRecordsOptions {
  readonly record: SessionRecord;
  /** The spec's `persist.rows`; the browser's IndexedDB when absent. */
  readonly rows?: ThreadRows;
  /** The IndexedDB database to read when there are no rows (`recordsDatabase(persist.key)`). */
  readonly database?: string;
  readonly max?: number;
}

export function createKeptRecordsModel(o: KeptRecordsOptions): KeptRecordsModel {
  const max = Math.max(1, o.max ?? KeptRecordsMax);
  const list = observableState<readonly KeptRecord[]>([]);
  const opened = observableState<OpenedRecord | undefined>(undefined);
  type Opened = { rows: ThreadRows & { close?(): void }; index: RecordIndex };
  let opening: Promise<Opened | undefined> | undefined;
  let watchers = 0;
  let disposed = false;

  // One opening at a time; none (no database yet: nothing was kept) is tried again on the next read.
  const open = () => {
    opening ??= (async () => {
      const rows = o.rows ?? (await openRecordsRows(o.database));
      if (!rows) return undefined;
      return { rows, index: await createRecordIndex(rows) };
    })().then(
      (r) => {
        if (!r) opening = undefined;
        return r;
      },
      (e) => {
        opening = undefined;
        throw e;
      },
    );
    return opening;
  };
  const own = (k: KeptRecord) => o.record.kept && k.id === o.record.id;
  const all = async () => {
    const at = await open();
    return at ? await at.index.list() : [];
  };
  const refresh = async () => {
    const kept = await all().catch((e) => {
      console.error('[prism-viewer] kept records:', e);
      return [] as KeptRecord[];
    });
    if (disposed) return;
    list.set(kept.filter((k) => !own(k)));
  };
  // A keep or an unkeep of this session's record changes what is listed.
  let wasKept = o.record.kept;
  const unlisten = o.record.onChange(() => {
    if (o.record.kept === wasKept) return;
    wasKept = o.record.kept;
    if (watchers > 0) void refresh();
  });

  return {
    list: list.asReadonly,
    opened: opened.asReadonly,
    max,
    watch() {
      watchers++;
      void refresh();
      let stopped = false;
      return () => {
        if (stopped) return;
        stopped = true;
        watchers--;
      };
    },
    refresh,
    async open(id) {
      const at = await open();
      const record = at && (await at.index.list()).find((k) => k.id === id);
      if (!at || !record || disposed) return false;
      opened.set({ record, reader: openKeptRecord(at.rows, id) });
      return true;
    },
    show(record, reader) {
      opened.set({ record, reader, fromFile: true });
    },
    close() {
      opened.set(undefined);
    },
    async forget(id) {
      const at = await open();
      if (!at) return;
      if (own({ id, title: '', first: 0, last: 0 })) await o.record.unkeep();
      else await forgetKept(at.rows, at.index, id);
      if (opened.value?.record.id === id) opened.set(undefined);
      await refresh();
    },
    async overBound() {
      if (o.record.kept) return undefined;
      const kept = await all().catch(() => [] as KeptRecord[]);
      return kept.length >= max ? kept[kept.length - 1] : undefined;
    },
    dispose() {
      disposed = true;
      unlisten();
      void opening?.then((at) => {
        if (at && !o.rows) at.rows.close?.();
      });
    },
  };
}
