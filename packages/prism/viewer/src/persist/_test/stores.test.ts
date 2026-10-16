/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { webStore } from '../stores.js';

// The browser stores (journal.md §6): localStorage, else sessionStorage, else nothing; a storage that throws never
// breaks the viewer.

/** A Storage over a map; `refuse` makes its writes throw (full, or a private mode). */
class MemoryStorage {
  readonly #items = new Map<string, string>();
  constructor(readonly refuse = false) {}
  get length() {
    return this.#items.size;
  }
  key(i: number) {
    return [...this.#items.keys()][i] ?? null;
  }
  getItem(k: string) {
    return this.#items.get(k) ?? null;
  }
  setItem(k: string, v: string) {
    if (this.refuse) throw new Error('QuotaExceededError');
    this.#items.set(k, v);
  }
  removeItem(k: string) {
    this.#items.delete(k);
  }
}

const options = { prefix: 'test', version: '1', budget: 1024 };

afterEach(() => vi.unstubAllGlobals());

describe('webStore', () => {
  it('keeps in localStorage when it works', () => {
    const local = new MemoryStorage();
    vi.stubGlobal('localStorage', local);
    vi.stubGlobal('sessionStorage', new MemoryStorage());
    const store = webStore('local', options);
    store.set('theme', 'dark');
    expect(store.get('theme')).toBe('dark');
    expect(local.getItem('test.theme')).not.toBeNull();
  });

  it('falls back to sessionStorage where localStorage refuses writes', () => {
    const session = new MemoryStorage();
    vi.stubGlobal('localStorage', new MemoryStorage(true));
    vi.stubGlobal('sessionStorage', session);
    const store = webStore('local', options);
    store.set('theme', 'dark');
    expect(store.get('theme')).toBe('dark');
    expect(session.getItem('test.theme')).not.toBeNull();
  });

  it('keeps nothing, and throws nothing, where reading the storages throws', () => {
    const blocked = () => {
      throw new Error('SecurityError');
    };
    Object.defineProperty(globalThis, 'localStorage', { configurable: true, get: blocked });
    Object.defineProperty(globalThis, 'sessionStorage', { configurable: true, get: blocked });
    try {
      const store = webStore('local', options);
      expect(() => store.set('theme', 'dark')).not.toThrow();
      expect(store.get('theme')).toBeUndefined();
      expect(store.entries()).toEqual([]);
      expect(() => store.remove('theme')).not.toThrow();
    } finally {
      delete (globalThis as { localStorage?: unknown }).localStorage;
      delete (globalThis as { sessionStorage?: unknown }).sessionStorage;
    }
  });
});
