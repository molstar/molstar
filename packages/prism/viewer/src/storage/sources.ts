/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ThreadRows } from '@molstar/ai';
import {
  createRecordIndex,
  forgetKept,
  keptSize,
  openRecordsRows,
  type RecordIndex,
  type SessionRecord,
} from '../journal/record.js';
import { isListedStore, type ListedStore } from '../persist/stores.js';
import type { PrismViewerPersist } from '../spec.js';
import { formatAge, type StorageEntry, type StorageSource } from './model.js';

// The viewer's own storage sources (journal.md §7): kept records (IndexedDB, or the spec's `persist.rows`), the
// preferences (localStorage) and this tab's light restore (sessionStorage), each with Forget.

/** Kept records not written to for this long are offered for removal (journal.md §7: 90 days by default). */
export const KeptRecordMaxAgeDays = 90;
const DayMs = 24 * 60 * 60_000;

export interface KeptRecordsSourceOptions {
  /** Where kept records are: the spec's `persist.rows`; else the browser's IndexedDB, opened only when it exists. */
  readonly rows?: ThreadRows;
  /** This session's record: forgetting it while it is kept takes it back into memory (`unkeep`). */
  readonly record?: SessionRecord;
  readonly now?: () => number;
  /** Older than this (by its last write) is offered for removal; 90 days by default. */
  readonly maxAgeDays?: number;
}

export interface KeptRecordsSource extends StorageSource {
  /** Closes the IndexedDB rows it opened. */
  dispose(): void;
}

export function keptRecordsSource(o: KeptRecordsSourceOptions = {}): KeptRecordsSource {
  const now = o.now ?? Date.now;
  const maxAge = (o.maxAgeDays ?? KeptRecordMaxAgeDays) * DayMs;
  type Opened = { rows: ThreadRows & { close?(): void }; index: Promise<RecordIndex> };
  // One opening at a time; none (no database yet) is tried again on the next listing.
  let opening: Promise<Opened | undefined> | undefined;
  let disposed = false;

  const open = () => {
    opening ??= (async () => {
      const rows: Opened['rows'] | undefined = o.rows ?? (await openRecordsRows());
      if (!rows) return undefined;
      if (disposed) {
        if (!o.rows) rows.close?.();
        return undefined;
      }
      return { rows, index: createRecordIndex(rows) };
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
  const current = (id: string) => !!o.record?.kept && o.record.id === id;

  return {
    id: 'kept-records',
    title: 'Kept records',
    where: o.rows ? 'the host’s store' : 'IndexedDB',
    async list() {
      const at = await open();
      if (!at) return [];
      const kept = await (await at.index).list();
      const t = now();
      return Promise.all(
        kept.map(async (k): Promise<StorageEntry> => {
          const bytes = await keptSize(at.rows, k.id).catch(() => undefined);
          const stale = t - k.last > maxAge;
          const detail = [
            current(k.id) ? 'This session' : undefined,
            `last written ${formatAge(t - k.last)}`,
            stale ? `older than ${Math.round(maxAge / DayMs)} days: forget it?` : undefined,
          ]
            .filter(Boolean)
            .join(' · ');
          return {
            id: k.id,
            label: k.title,
            detail: detail.charAt(0).toUpperCase() + detail.slice(1),
            ...(bytes !== undefined ? { bytes } : {}),
            ...(stale ? { stale: true } : {}),
          };
        }),
      );
    },
    async forget(id) {
      // This session's record goes on in memory; any other is removed with its index entry.
      if (current(id) && o.record) return o.record.unkeep();
      const at = await open();
      if (!at) return;
      await forgetKept(at.rows, await at.index, id);
    },
    dispose() {
      disposed = true;
      if (!o.rows) void opening?.then((r) => r?.rows.close?.()).catch(() => undefined);
      opening = undefined;
    },
  };
}

export interface StoreSourceOptions {
  readonly id: string;
  readonly title: string;
  readonly where: string;
  readonly store: ListedStore;
  /** Each key's label ('theme' → 'Theme'); the key itself otherwise. */
  readonly labels?: Readonly<Record<string, string>>;
  /** Said of every entry. */
  readonly detail?: string;
}

/** The values a light store keeps (journal.md §6), each with its size; Forget removes it. */
export function storeSource(o: StoreSourceOptions): StorageSource {
  return {
    id: o.id,
    title: o.title,
    where: o.where,
    async list() {
      return o.store.entries().map(
        (e): StorageEntry => ({
          id: e.key,
          label: o.labels?.[e.key] ?? e.key,
          bytes: e.bytes,
          ...(o.detail ? { detail: o.detail } : {}),
        }),
      );
    },
    async forget(key) {
      o.store.remove(key);
    },
  };
}

/** What the preferences hold, in words. */
const PreferenceLabels: Readonly<Record<string, string>> = { theme: 'Theme', layout: 'Layout' };

/**
 * The viewer's own sources: kept records always (empty when there are none), the preferences and the light restore when
 * the spec gives stores that list what they keep (`webStore`, `memoryStore`). Each source with a `dispose` is disposed
 * by `dispose` here.
 */
export function viewerStorageSources(o: { persist?: PrismViewerPersist; record?: SessionRecord }): {
  readonly sources: readonly StorageSource[];
  dispose(): void;
} {
  const rows = o.persist?.rows as ThreadRows | undefined;
  const records = keptRecordsSource({ ...(rows ? { rows } : {}), ...(o.record ? { record: o.record } : {}) });
  const sources: StorageSource[] = [records];
  const preferences = o.persist?.preferences;
  if (isListedStore(preferences)) {
    sources.push(
      storeSource({
        id: 'preferences',
        title: 'Preferences',
        where: 'localStorage',
        store: preferences,
        labels: PreferenceLabels,
      }),
    );
  }
  const restore = o.persist?.restore;
  if (isListedStore(restore)) {
    sources.push(
      storeSource({
        id: 'light-restore',
        title: 'Light restore',
        where: 'sessionStorage',
        store: restore,
        labels: { restore: 'This tab’s scene, without its data' },
        detail: 'Ends with the tab; written again after the next change',
      }),
    );
  }
  return { sources, dispose: () => records.dispose() };
}
