/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';

// What the viewer keeps in the browser, and Forget (journal.md §7, README decision 32), DOM-free. Each part that keeps
// something registers a source (the viewer: kept records, preferences, the light restore; the assistant: stored keys);
// the Storage view lists them while it shows (`watch`) and not otherwise: a hidden view reads nothing.

/** One thing kept: a kept record, a preference, a key (never its value). */
export interface StorageEntry {
  readonly id: string;
  readonly label: string;
  /** What else to say: its age, its scope, why it is offered for removal. */
  readonly detail?: string;
  /** Its size where it is known (about: characters as the store keeps them). */
  readonly bytes?: number;
  /** Offered for removal (journal.md §7: a kept record older than 90 days). */
  readonly stale?: boolean;
}

export interface StorageSource {
  readonly id: string;
  readonly title: string;
  /** Where it is kept, in words ('IndexedDB', 'localStorage'). */
  readonly where?: string;
  list(): Promise<readonly StorageEntry[]>;
  forget(entryId: string): Promise<void>;
  /** Forgets every entry at once; without it, Forget everything forgets each listed one. */
  forgetAll?(): Promise<void>;
}

/** A source as the view shows it. */
export interface StorageSourceState {
  readonly id: string;
  readonly title: string;
  readonly where?: string;
  readonly entries: readonly StorageEntry[];
  /** The entries' sizes added up, when any is known. */
  readonly bytes?: number;
  /** Listed at least once. */
  readonly listed: boolean;
  /** Why it could not be listed (or forgotten). */
  readonly error?: string;
}

export interface StorageModel {
  /** The registered sources, in the order they came, with what they listed last. */
  readonly sources: ReadonlyState<readonly StorageSourceState[]>;
  /** The origin that holds the data (`location.origin` where there is one). */
  readonly origin?: string;
  /** Adds a source (its id unique); the returned function removes it. Listed at once while watched. */
  register(source: StorageSource): () => void;
  /** Lists every source again. */
  refresh(): Promise<void>;
  /** Forgets one entry, then lists its source again. */
  forget(sourceId: string, entryId: string): Promise<void>;
  /** Forgets every entry of every source, then lists them again. */
  forgetEverything(): Promise<void>;
  /** The view shows: lists now, and a source registered meanwhile at once. The returned function ends it. */
  watch(): () => void;
  dispose(): void;
}

const messageOf = (e: unknown) => (e instanceof Error ? e.message : String(e));

export function createStorageModel(o: { origin?: string } = {}): StorageModel {
  const registered = new Map<string, StorageSource>();
  const states = new Map<string, StorageSourceState>();
  const sources = observableState<readonly StorageSourceState[]>([]);
  // Each source's latest listing: one that finishes after a newer one started is dropped.
  const generation = new Map<string, number>();
  let watchers = 0;
  let disposed = false;

  const publish = () => {
    if (!disposed) sources.set([...registered.keys()].map((id) => states.get(id)!));
  };
  const blank = (s: StorageSource): StorageSourceState => ({
    id: s.id,
    title: s.title,
    ...(s.where ? { where: s.where } : {}),
    entries: [],
    listed: false,
  });

  const list = async (source: StorageSource) => {
    const n = (generation.get(source.id) ?? 0) + 1;
    generation.set(source.id, n);
    let next: StorageSourceState;
    try {
      const entries = await source.list();
      const sized = entries.filter((e) => e.bytes !== undefined);
      next = {
        ...blank(source),
        entries,
        listed: true,
        ...(sized.length ? { bytes: sized.reduce((sum, e) => sum + (e.bytes ?? 0), 0) } : {}),
      };
    } catch (e) {
      next = { ...blank(source), listed: true, error: `Could not be read: ${messageOf(e)}` };
    }
    if (generation.get(source.id) !== n || registered.get(source.id) !== source) return;
    states.set(source.id, next);
    publish();
  };

  const failed = (source: StorageSource, e: unknown) => {
    const state = states.get(source.id);
    if (state) states.set(source.id, { ...state, error: `Could not be forgotten: ${messageOf(e)}` });
    publish();
  };

  const model: StorageModel = {
    sources: sources.asReadonly,
    ...(o.origin ? { origin: o.origin } : {}),
    register(source) {
      if (registered.has(source.id)) throw new Error(`A storage source '${source.id}' is registered already.`);
      registered.set(source.id, source);
      states.set(source.id, blank(source));
      publish();
      if (watchers > 0) void list(source);
      return () => {
        if (registered.get(source.id) !== source) return;
        registered.delete(source.id);
        states.delete(source.id);
        generation.delete(source.id);
        publish();
      };
    },
    async refresh() {
      await Promise.all([...registered.values()].map(list));
    },
    async forget(sourceId, entryId) {
      const source = registered.get(sourceId);
      if (!source) return;
      try {
        await source.forget(entryId);
      } catch (e) {
        failed(source, e);
        return;
      }
      await list(source);
    },
    async forgetEverything() {
      await Promise.all(
        [...registered.values()].map(async (source) => {
          try {
            if (source.forgetAll) await source.forgetAll();
            else for (const e of await source.list()) await source.forget(e.id);
          } catch (e) {
            failed(source, e);
            return;
          }
          await list(source);
        }),
      );
    },
    watch() {
      watchers++;
      if (watchers === 1) void model.refresh();
      let done = false;
      return () => {
        if (done) return;
        done = true;
        watchers--;
      };
    },
    dispose() {
      disposed = true;
      registered.clear();
      states.clear();
      generation.clear();
    },
  };
  return model;
}

/** A size in words: "312 B", "4.2 KB", "1.3 MB". */
export function formatBytes(bytes: number): string {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(bytes < 10 * 1024 ? 1 : 0)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}

/** An age in words: "just now", "5 minutes ago", "3 days ago". */
export function formatAge(ms: number): string {
  const minutes = Math.floor(ms / 60_000);
  if (minutes < 1) return 'just now';
  if (minutes < 60) return `${minutes} minute${minutes === 1 ? '' : 's'} ago`;
  const hours = Math.floor(minutes / 60);
  if (hours < 24) return `${hours} hour${hours === 1 ? '' : 's'} ago`;
  const days = Math.floor(hours / 24);
  return `${days} day${days === 1 ? '' : 's'} ago`;
}
