/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PrismStore } from '../spec.js';

// The light stores (journal.md §6–§7): localStorage for preferences, sessionStorage for this tab's light restore.
// Every value carries the version it was written by and is dropped when it does not match (§7); each store has a size
// budget, and a value over it is not kept (said once). A private window or a blocked store: nothing is kept, nothing
// fails.

export interface WebStoreOptions {
  /** Every key's prefix: the viewer's `persist.key` (several viewers on a page keep apart). */
  readonly prefix: string;
  /** What the values are stamped with; another version's are dropped when read. */
  readonly version: string;
  /** The largest value kept, in characters of its JSON (localStorage ~64 KB, sessionStorage ~256 KB: §6). */
  readonly budget: number;
}

interface Envelope {
  readonly v: string;
  readonly d: unknown;
}

/** A value kept, as the Storage view lists it (journal.md §7): its key under the prefix, and its size. */
export interface StoredValue {
  readonly key: string;
  /** About: the characters of its key and its JSON, two bytes each (as the browser counts them). */
  readonly bytes: number;
}

/** A store that lists what it keeps; another version's values are dropped while listing (journal.md §7). */
export interface ListedStore extends PrismStore {
  entries(): readonly StoredValue[];
}

/** Whether a store (the spec's, `persist.preferences` or `persist.restore`) lists what it keeps. */
export function isListedStore(store: PrismStore | undefined): store is ListedStore {
  return typeof (store as ListedStore | undefined)?.entries === 'function';
}

const isEnvelope = (e: unknown): e is Envelope => !!e && typeof e === 'object' && typeof (e as Envelope).v === 'string';

function storageOf(kind: 'local' | 'session'): Storage | undefined {
  try {
    return kind === 'local' ? globalThis.localStorage : globalThis.sessionStorage;
  } catch {
    return undefined;
  }
}

/** A store over the window's localStorage or sessionStorage (absent in Node: a store that keeps nothing). */
export function webStore(kind: 'local' | 'session', o: WebStoreOptions): ListedStore {
  const key = (k: string) => `${o.prefix}.${k}`;
  let warned = false;
  return {
    get(k) {
      const storage = storageOf(kind);
      try {
        const text = storage?.getItem(key(k));
        if (!text) return undefined;
        const e = JSON.parse(text) as Envelope;
        if (e?.v !== o.version) {
          storage?.removeItem(key(k));
          return undefined;
        }
        return e.d;
      } catch {
        return undefined;
      }
    },
    set(k, value) {
      const storage = storageOf(kind);
      if (!storage) return;
      try {
        const text = JSON.stringify({ v: o.version, d: value } satisfies Envelope);
        if (text.length > o.budget) {
          storage.removeItem(key(k));
          if (!warned) console.info(`[prism-viewer] '${k}' is over the ${kind} store's budget: not kept.`);
          warned = true;
          return;
        }
        storage.setItem(key(k), text);
      } catch {
        // Full, or not allowed: not kept.
      }
    },
    remove(k) {
      try {
        storageOf(kind)?.removeItem(key(k));
      } catch {
        // Not allowed: nothing was kept.
      }
    },
    entries() {
      const storage = storageOf(kind);
      const out: StoredValue[] = [];
      const dropped: string[] = [];
      const start = `${o.prefix}.`;
      try {
        if (!storage) return out;
        for (let i = 0; i < storage.length; i++) {
          const name = storage.key(i);
          if (!name?.startsWith(start)) continue;
          const text = storage.getItem(name) ?? '';
          let e: unknown;
          try {
            e = JSON.parse(text);
          } catch {
            continue; // not a value of ours
          }
          if (!isEnvelope(e)) continue;
          if (e.v !== o.version) dropped.push(name);
          else out.push({ key: name.slice(start.length), bytes: (name.length + text.length) * 2 });
        }
        // Collected (§7): another version's values go as they are listed, not only when read.
        for (const name of dropped) storage.removeItem(name);
      } catch {
        // Not allowed: nothing is listed.
      }
      return out;
    },
  };
}

/** A store in memory (tests, headless hosts), with the same budget and version rule. */
export function memoryStore(o: Omit<WebStoreOptions, 'prefix'>): ListedStore & { readonly size: number } {
  const values = new Map<string, string>();
  return {
    get size() {
      return values.size;
    },
    get(k) {
      const text = values.get(k);
      if (!text) return undefined;
      const e = JSON.parse(text) as Envelope;
      return e.v === o.version ? e.d : undefined;
    },
    set(k, value) {
      const text = JSON.stringify({ v: o.version, d: value } satisfies Envelope);
      if (text.length > o.budget) values.delete(k);
      else values.set(k, text);
    },
    remove(k) {
      values.delete(k);
    },
    entries() {
      const out: StoredValue[] = [];
      for (const [k, text] of [...values]) {
        const e = JSON.parse(text) as Envelope;
        if (e.v !== o.version) values.delete(k);
        else out.push({ key: k, bytes: (k.length + text.length) * 2 });
      }
      return out;
    },
  };
}
