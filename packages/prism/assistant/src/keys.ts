/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The assistant's key store (spec/agent.md §4, "Keys"). Keys are the person's, kept for the scope they chose; stored
// ones are always scrambled, so the storage read as it is shows no key. Nothing here guards against code running in
// the page. A key is never logged.

import type { KeyScope } from './model.js';

export type { KeyScope } from './model.js';

/** How stored keys are kept here: 'per-browser' protects, the others only scramble. */
export type KeyProtection = 'per-browser' | 'app-key' | 'scramble';

type StorageLike = Pick<Storage, 'getItem' | 'setItem' | 'removeItem' | 'key' | 'length'>;

export interface KeyStore {
  /** The key for a provider, opened (kept in memory while in use); undefined when there is none. */
  get(provider: string): Promise<string | undefined>;
  set(provider: string, key: string, scope: KeyScope): Promise<void>;
  remove(provider: string): Promise<void>;
  /** The providers with a key, and the scope each is kept for. */
  list(): Promise<readonly { provider: string; scope: KeyScope }[]>;
  /** How stored keys are kept here (what the Assistant tab says: 'per-browser' protects; the others only scramble). */
  protection(): Promise<KeyProtection>;
}

export interface KeyStoreOptions {
  /** Default: `globalThis.sessionStorage` when present. */
  readonly session?: StorageLike;
  /** Default: `globalThis.localStorage` when present. */
  readonly local?: StorageLike;
  /** Default: `globalThis.indexedDB` when present. */
  readonly indexedDB?: IDBFactory;
  /** Default: `globalThis.crypto` when `crypto.subtle` is present (needs a secure context). */
  readonly crypto?: Crypto;
  /** The storage key prefix, default 'prism-assistant.key.'. */
  readonly prefix?: string;
}

type Kind = 'aes-browser' | 'aes-app' | 'xor';

interface Envelope {
  v: 1;
  kind: Kind;
  iv?: string;
  data: string;
}

const DefaultPrefix = 'prism-assistant.key.';
const DbName = 'prism-assistant-keys';
const StoreName = 'keys';
const RecordName = 'scrambling-key';

// The app key: in the source, so it scrambles against a naive read and nothing more.
const AppKeyText = 'prism-assistant/app-key/v1/0x5A3C';
const AppKeyBytes = new TextEncoder().encode(AppKeyText.padEnd(32, '.').slice(0, 32));

function toBase64(bytes: Uint8Array): string {
  let s = '';
  for (let i = 0; i < bytes.length; i++) s += String.fromCharCode(bytes[i]);
  return btoa(s);
}

function fromBase64(text: string): Uint8Array<ArrayBuffer> {
  const s = atob(text);
  const bytes = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) bytes[i] = s.charCodeAt(i);
  return bytes;
}

function xor(bytes: Uint8Array): Uint8Array<ArrayBuffer> {
  const out = new Uint8Array(bytes.length);
  for (let i = 0; i < bytes.length; i++) out[i] = bytes[i] ^ AppKeyBytes[i % AppKeyBytes.length];
  return out;
}

function defaultStorage(name: 'sessionStorage' | 'localStorage'): StorageLike | undefined {
  try {
    return (globalThis as { [k: string]: unknown })[name] as StorageLike | undefined;
  } catch {
    return undefined; // access itself can throw (blocked site data)
  }
}

function defaultCrypto(): Crypto | undefined {
  try {
    const c = globalThis.crypto;
    return c?.subtle ? c : undefined;
  } catch {
    return undefined;
  }
}

function defaultIndexedDB(): IDBFactory | undefined {
  try {
    return globalThis.indexedDB ?? undefined;
  } catch {
    return undefined;
  }
}

function request<T>(req: IDBRequest<T>): Promise<T> {
  return new Promise((resolve, reject) => {
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error ?? new Error('IndexedDB request failed'));
  });
}

/** Opens the database, reads the scrambling key or generates and keeps one (the first writer wins). */
async function loadBrowserKey(idb: IDBFactory, crypto: Crypto): Promise<CryptoKey> {
  const open = idb.open(DbName, 1);
  open.onupgradeneeded = () => {
    if (!open.result.objectStoreNames.contains(StoreName)) open.result.createObjectStore(StoreName);
  };
  const db = await request(open);
  try {
    const existing = (await request(db.transaction(StoreName, 'readonly').objectStore(StoreName).get(RecordName))) as
      | CryptoKey
      | undefined;
    if (existing) return existing;
    const created = await crypto.subtle.generateKey({ name: 'AES-GCM', length: 256 }, false, ['encrypt', 'decrypt']);
    const tx = db.transaction(StoreName, 'readwrite');
    const store = tx.objectStore(StoreName);
    const raced = (await request(store.get(RecordName))) as CryptoKey | undefined;
    if (raced) return raced;
    await request(store.put(created, RecordName));
    await new Promise<void>((resolve, reject) => {
      tx.oncomplete = () => resolve();
      tx.onerror = () => reject(tx.error ?? new Error('IndexedDB transaction failed'));
      tx.onabort = () => reject(tx.error ?? new Error('IndexedDB transaction aborted'));
    });
    return created;
  } finally {
    db.close();
  }
}

export function createKeyStore(options: KeyStoreOptions = {}): KeyStore {
  const session = 'session' in options ? options.session : defaultStorage('sessionStorage');
  const local = 'local' in options ? options.local : defaultStorage('localStorage');
  const idb = 'indexedDB' in options ? options.indexedDB : defaultIndexedDB();
  const crypto = 'crypto' in options ? options.crypto : defaultCrypto();
  const prefix = options.prefix ?? DefaultPrefix;
  const subtle = crypto?.subtle;

  /** What is in use: the keys of every scope (a stored one once opened or set), so storage that fails loses nothing. */
  const memory = new Map<string, { key: string; scope: KeyScope }>();

  let browserKey: Promise<CryptoKey> | undefined;
  /** The per-browser key; a failure is not kept, so IndexedDB granted later is tried again. */
  async function getBrowserKey(): Promise<CryptoKey | undefined> {
    if (!idb || !crypto) return undefined;
    if (!browserKey) {
      const attempt = loadBrowserKey(idb, crypto);
      browserKey = attempt;
      attempt.catch(() => {
        if (browserKey === attempt) browserKey = undefined;
      });
    }
    try {
      return await browserKey;
    } catch {
      return undefined;
    }
  }

  let appKey: Promise<CryptoKey> | undefined;
  function getAppKey(): Promise<CryptoKey> {
    return (appKey ??= subtle!.importKey('raw', AppKeyBytes, { name: 'AES-GCM' }, false, ['encrypt', 'decrypt']));
  }

  async function aesEncrypt(key: CryptoKey, text: string): Promise<{ iv: string; data: string }> {
    const iv = crypto!.getRandomValues(new Uint8Array(12));
    const data = await subtle!.encrypt({ name: 'AES-GCM', iv }, key, new TextEncoder().encode(text));
    return { iv: toBase64(iv), data: toBase64(new Uint8Array(data)) };
  }

  async function aesDecrypt(key: CryptoKey, envelope: Envelope): Promise<string> {
    if (!envelope.iv) throw new Error('missing iv');
    const data = await subtle!.decrypt(
      { name: 'AES-GCM', iv: fromBase64(envelope.iv) },
      key,
      fromBase64(envelope.data),
    );
    return new TextDecoder().decode(data);
  }

  async function seal(text: string): Promise<Envelope> {
    if (subtle) {
      try {
        const bk = await getBrowserKey();
        if (bk) return { v: 1, kind: 'aes-browser', ...(await aesEncrypt(bk, text)) };
      } catch {
        // fall to the next row
      }
      try {
        return { v: 1, kind: 'aes-app', ...(await aesEncrypt(await getAppKey(), text)) };
      } catch {
        // fall to the next row
      }
    }
    return { v: 1, kind: 'xor', data: toBase64(xor(new TextEncoder().encode(text))) };
  }

  /** Tries the per-browser key, then the app key (and the XOR form), the stored kind first. */
  async function open(envelope: Envelope): Promise<string | undefined> {
    if (envelope.kind === 'xor') {
      try {
        return new TextDecoder().decode(xor(fromBase64(envelope.data)));
      } catch {
        return undefined;
      }
    }
    if (!subtle) return undefined;
    const order = envelope.kind === 'aes-app' ? ['app', 'browser'] : ['browser', 'app'];
    for (const which of order) {
      try {
        const key = which === 'app' ? await getAppKey() : await getBrowserKey();
        if (key) return await aesDecrypt(key, envelope);
      } catch {
        // try the next one
      }
    }
    return undefined;
  }

  function parse(text: string | null): Envelope | undefined {
    if (!text) return undefined;
    try {
      const e = JSON.parse(text);
      if (
        e &&
        e.v === 1 &&
        typeof e.data === 'string' &&
        (e.kind === 'aes-browser' || e.kind === 'aes-app' || e.kind === 'xor')
      )
        return e as Envelope;
    } catch {
      // not ours
    }
    return undefined;
  }

  function storageOf(scope: KeyScope): StorageLike | undefined {
    return scope === 'session' ? session : scope === 'browser' ? local : undefined;
  }

  function tryRemove(storage: StorageLike | undefined, name: string) {
    try {
      storage?.removeItem(name);
    } catch {
      // storage can throw (a private window)
    }
  }

  function tryRead(storage: StorageLike | undefined, name: string): string | null {
    try {
      return storage?.getItem(name) ?? null;
    } catch {
      return null;
    }
  }

  function storedProviders(storage: StorageLike | undefined): string[] {
    const out: string[] = [];
    try {
      if (!storage) return out;
      for (let i = 0; i < storage.length; i++) {
        const name = storage.key(i);
        if (name?.startsWith(prefix)) out.push(name.slice(prefix.length));
      }
    } catch {
      // storage can throw
    }
    return out;
  }

  return {
    async get(provider) {
      const held = memory.get(provider);
      if (held) return held.key;
      const name = prefix + provider;
      for (const scope of ['session', 'browser'] as const) {
        const envelope = parse(tryRead(storageOf(scope), name));
        if (!envelope) continue;
        const key = await open(envelope);
        if (key === undefined) continue;
        memory.set(provider, { key, scope });
        return key;
      }
      return undefined;
    },

    async set(provider, key, scope) {
      const name = prefix + provider;
      memory.set(provider, { key, scope });
      const target = storageOf(scope);
      // A key lives in one scope: the others lose their copy.
      for (const other of ['session', 'browser'] as const) {
        if (other !== scope) tryRemove(storageOf(other), name);
      }
      if (!target) return; // memory only ('reload'), or no such storage here: kept in memory for the session
      try {
        target.setItem(name, JSON.stringify(await seal(key)));
      } catch {
        // storage refused (a private window, full): the key stays in memory for the session
      }
    },

    async remove(provider) {
      memory.delete(provider);
      const name = prefix + provider;
      tryRemove(session, name);
      tryRemove(local, name);
    },

    async list() {
      const out = new Map<string, KeyScope>();
      for (const [provider, { scope }] of memory) out.set(provider, scope);
      for (const scope of ['session', 'browser'] as const) {
        const storage = storageOf(scope);
        for (const provider of storedProviders(storage)) {
          if (out.has(provider) || !parse(tryRead(storage, prefix + provider))) continue;
          out.set(provider, scope);
        }
      }
      return [...out].map(([provider, scope]) => ({ provider, scope }));
    },

    async protection() {
      if (!subtle) return 'scramble';
      return (await getBrowserKey()) ? 'per-browser' : 'app-key';
    },
  };
}
