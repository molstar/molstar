/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { IDBFactory } from 'fake-indexeddb';
import { describe, expect, it } from 'vitest';
import { createKeyStore, type KeyStoreOptions } from '../keys.js';

class MemoryStorage {
  private readonly map = new Map<string, string>();
  get length() {
    return this.map.size;
  }
  key(i: number) {
    return [...this.map.keys()][i] ?? null;
  }
  getItem(k: string) {
    return this.map.get(k) ?? null;
  }
  setItem(k: string, v: string) {
    this.map.set(k, v);
  }
  removeItem(k: string) {
    this.map.delete(k);
  }
  /** Everything stored, as one text (what a read of the storage shows). */
  dump() {
    return JSON.stringify([...this.map]);
  }
}

class ThrowingStorage {
  get length(): number {
    throw new Error('blocked');
  }
  key(): string | null {
    throw new Error('blocked');
  }
  getItem(): string | null {
    throw new Error('blocked');
  }
  setItem(): void {
    throw new Error('blocked');
  }
  removeItem(): void {
    throw new Error('blocked');
  }
}

const secret = 'sk-test-0123456789-SECRET';

function world() {
  const session = new MemoryStorage();
  const local = new MemoryStorage();
  const indexedDB = new IDBFactory();
  const base = { session, local, indexedDB, crypto: globalThis.crypto } satisfies KeyStoreOptions;
  return { session, local, indexedDB, base };
}

describe('key store', () => {
  it.each(['reload', 'session', 'browser'] as const)('round-trips in the %s scope', async (scope) => {
    const { base } = world();
    const store = createKeyStore(base);
    expect(await store.get('anthropic')).toBeUndefined();
    await store.set('anthropic', secret, scope);
    expect(await store.get('anthropic')).toBe(secret);
    // a new page (a new store over the same storage) finds what was stored, not what was in memory only
    const reopened = createKeyStore(base);
    expect(await reopened.get('anthropic')).toBe(scope === 'reload' ? undefined : secret);
  });

  it('never stores the key as it is, in any row', async () => {
    const rows: Partial<KeyStoreOptions>[] = [
      {}, // WebCrypto and IndexedDB
      { indexedDB: undefined }, // WebCrypto only
      { indexedDB: undefined, crypto: undefined }, // neither
    ];
    for (const row of rows) {
      const { session, local, base } = world();
      const store = createKeyStore({ ...base, ...row });
      await store.set('a', secret, 'session');
      await store.set('b', secret, 'browser');
      for (const text of [session.dump(), local.dump()]) {
        expect(text).not.toContain(secret);
        expect(text).not.toContain(btoa(secret));
        expect(text.length).toBeGreaterThan(20);
      }
      expect(await store.get('a')).toBe(secret);
      expect(await createKeyStore({ ...base, ...row }).get('b')).toBe(secret);
    }
  });

  it('reports the protection of each row', async () => {
    const { base } = world();
    expect(await createKeyStore(base).protection()).toBe('per-browser');
    expect(await createKeyStore({ ...base, indexedDB: undefined }).protection()).toBe('app-key');
    expect(await createKeyStore({ ...base, indexedDB: undefined, crypto: undefined }).protection()).toBe('scramble');
    // no WebCrypto, whatever IndexedDB offers
    expect(await createKeyStore({ ...base, crypto: undefined }).protection()).toBe('scramble');
  });

  it('keeps one per-browser key, generated once', async () => {
    const { session, base } = world();
    const first = createKeyStore(base);
    await first.set('a', secret, 'session');
    expect(JSON.parse(session.getItem('prism-assistant.key.a')!).kind).toBe('aes-browser');
    // another page of the same browser opens it
    expect(await createKeyStore(base).get('a')).toBe(secret);
  });

  it('opens a value written under the app key once IndexedDB is available', async () => {
    const { local, base } = world();
    await createKeyStore({ ...base, indexedDB: undefined }).set('a', secret, 'browser');
    expect(JSON.parse(local.getItem('prism-assistant.key.a')!).kind).toBe('aes-app');
    expect(await createKeyStore(base).get('a')).toBe(secret);
  });

  it('opens a scrambled value where the other forms are in use', async () => {
    const { local, base } = world();
    // XOR form, read where WebCrypto is available
    await createKeyStore({ ...base, indexedDB: undefined, crypto: undefined }).set('x', secret, 'browser');
    expect(JSON.parse(local.getItem('prism-assistant.key.x')!).kind).toBe('xor');
    expect(await createKeyStore(base).get('x')).toBe(secret);
    // XOR form, read where there is no WebCrypto
    expect(await createKeyStore({ ...base, indexedDB: undefined, crypto: undefined }).get('x')).toBe(secret);
    // app key, read where IndexedDB was withdrawn
    await createKeyStore({ ...base, indexedDB: undefined }).set('y', secret, 'browser');
    expect(await createKeyStore({ ...base, indexedDB: undefined }).get('y')).toBe(secret);
  });

  it('does not open a per-browser value without the per-browser key', async () => {
    const { base } = world();
    await createKeyStore(base).set('a', secret, 'browser');
    expect(await createKeyStore({ ...base, indexedDB: undefined }).get('a')).toBeUndefined();
    // copied to another browser: its own per-browser key does not open it
    const other = world();
    other.local.setItem('prism-assistant.key.a', base.local.getItem('prism-assistant.key.a')!);
    expect(await createKeyStore(other.base).get('a')).toBeUndefined();
  });

  it('keeps a key in one scope', async () => {
    const { session, local, base } = world();
    const store = createKeyStore(base);
    await store.set('a', secret, 'browser');
    await store.set('a', secret, 'session');
    expect(local.length).toBe(0);
    expect(session.length).toBe(1);
    expect(await store.list()).toEqual([{ provider: 'a', scope: 'session' }]);
    await store.set('a', secret, 'reload');
    expect(session.length).toBe(0);
    expect(await store.list()).toEqual([{ provider: 'a', scope: 'reload' }]);
    expect(await createKeyStore(base).get('a')).toBeUndefined();
  });

  it('removes a key from every scope', async () => {
    const { session, local, base } = world();
    const store = createKeyStore(base);
    await store.set('a', secret, 'reload');
    await store.set('b', secret, 'session');
    await store.set('c', secret, 'browser');
    session.setItem('prism-assistant.key.c', session.getItem('prism-assistant.key.b')!); // stale copy elsewhere
    for (const p of ['a', 'b', 'c']) await store.remove(p);
    expect(session.length).toBe(0);
    expect(local.length).toBe(0);
    expect(await store.list()).toEqual([]);
    for (const p of ['a', 'b', 'c']) expect(await store.get(p)).toBeUndefined();
  });

  it('lists the providers with their scopes, including ones stored by an earlier page', async () => {
    const { session, base } = world();
    const earlier = createKeyStore(base);
    await earlier.set('anthropic', secret, 'browser');
    await earlier.set('openai', secret, 'session');
    session.setItem('unrelated', 'x');
    const store = createKeyStore(base);
    await store.set('google', secret, 'reload');
    const list = [...(await store.list())].sort((a, b) => a.provider.localeCompare(b.provider));
    expect(list).toEqual([
      { provider: 'anthropic', scope: 'browser' },
      { provider: 'google', scope: 'reload' },
      { provider: 'openai', scope: 'session' },
    ]);
  });

  it('uses the prefix given', async () => {
    const { local, base } = world();
    await createKeyStore({ ...base, prefix: 'x.' }).set('a', secret, 'browser');
    expect(local.getItem('x.a')).not.toBeNull();
    expect(await createKeyStore(base).get('a')).toBeUndefined();
  });

  it('falls back to memory when storage throws, without losing the key', async () => {
    const { indexedDB } = world();
    const store = createKeyStore({
      session: new ThrowingStorage(),
      local: new ThrowingStorage(),
      indexedDB,
      crypto: globalThis.crypto,
    });
    await store.set('a', secret, 'session');
    await store.set('b', secret, 'browser');
    expect(await store.get('a')).toBe(secret);
    expect(await store.get('b')).toBe(secret);
    expect(await store.list()).toEqual([
      { provider: 'a', scope: 'session' },
      { provider: 'b', scope: 'browser' },
    ]);
    await store.remove('a');
    expect(await store.get('a')).toBeUndefined();
  });

  it('falls back when IndexedDB fails, and keeps working without any storage', async () => {
    const failing = {
      open() {
        throw new Error('denied');
      },
    } as unknown as IDBFactory;
    const { local, base } = world();
    const store = createKeyStore({ ...base, indexedDB: failing });
    expect(await store.protection()).toBe('app-key');
    await store.set('a', secret, 'browser');
    expect(JSON.parse(local.getItem('prism-assistant.key.a')!).kind).toBe('aes-app');
    const bare = createKeyStore({ session: undefined, local: undefined, indexedDB: undefined, crypto: undefined });
    await bare.set('a', secret, 'browser');
    expect(await bare.get('a')).toBe(secret);
  });

  it('ignores stored values it cannot read', async () => {
    const { local, base } = world();
    local.setItem('prism-assistant.key.a', 'not json');
    local.setItem('prism-assistant.key.b', JSON.stringify({ v: 2, kind: 'xor', data: '' }));
    const store = createKeyStore(base);
    expect(await store.get('a')).toBeUndefined();
    expect(await store.get('b')).toBeUndefined();
    expect(await store.list()).toEqual([]);
  });

  it('does not write a key back that was removed while it was being sealed', async () => {
    const { local, base } = world();
    const store = createKeyStore(base);
    // Set (it seals, which awaits), and remove before that is done.
    const setting = store.set('anthropic', secret, 'browser');
    await store.remove('anthropic');
    await setting;
    expect(await store.get('anthropic')).toBeUndefined();
    expect(local.length).toBe(0);
    expect(await createKeyStore(base).list()).toEqual([]);
  });

  it('keeps the later of two sets that race', async () => {
    const { local, session, base } = world();
    const store = createKeyStore(base);
    const first = store.set('anthropic', 'sk-first', 'browser');
    await store.set('anthropic', 'sk-second', 'session');
    await first;
    expect(local.length).toBe(0);
    expect(await createKeyStore(base).get('anthropic')).toBe('sk-second');
    expect(session.length).toBe(1);
  });

  it('does not list a stored key it cannot open', async () => {
    const { local, base } = world();
    await createKeyStore(base).set('anthropic', secret, 'browser');
    // The per-browser key is gone (IndexedDB cleared): the stored value no longer opens.
    const store = createKeyStore({ ...base, indexedDB: new IDBFactory() });
    expect(local.length).toBe(1);
    expect(await store.list()).toEqual([]);
  });
});
