/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { IDBFactory } from 'fake-indexeddb';
import { createStorageModel } from '@molstar/prism-viewer/storage/model';
import { describe, expect, it } from 'vitest';
import { createKeyStore, type KeyStoreOptions } from '../../keys.js';
import { createAISettings, type SettingsPrefs } from '../settings-model.js';
import { keysStorageSource } from '../storage.js';

// The assistant's keys in the viewer's Storage view (journal.md §7): by provider and scope, never the key; Forget
// removes it from the store and from the setup.

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
}

const secret = 'sk-test-0123456789-SECRET';

function world() {
  const session = new MemoryStorage();
  const local = new MemoryStorage();
  const base = { session, local, indexedDB: new IDBFactory(), crypto: globalThis.crypto } satisfies KeyStoreOptions;
  return { session, local, base };
}

const memoryPrefs = (): SettingsPrefs => {
  const values = new Map<string, string>();
  return { get: (k) => values.get(k), set: (k, v) => void values.set(k, v) };
};

describe("the assistant's keys in Storage", () => {
  it('lists each stored key by provider and scope, and never the key', async () => {
    const { base } = world();
    const keys = createKeyStore(base);
    await keys.set('anthropic', secret, 'browser');
    await keys.set('google', secret, 'session');
    const settings = createAISettings(keys, memoryPrefs());
    await settings.load();
    const entries = await keysStorageSource(keys, settings).list();
    expect(entries.map((e) => [e.id, e.label])).toEqual([
      ['anthropic', 'Anthropic'],
      ['google', 'Google'],
    ]);
    expect(entries[0].detail).toMatch(/this browser/);
    expect(entries[1].detail).toMatch(/this tab/);
    expect(JSON.stringify(entries)).not.toContain(secret);
  });

  it('forgets a key through the setup: gone from the store and from the providers', async () => {
    const { base, local } = world();
    const keys = createKeyStore(base);
    await keys.set('anthropic', secret, 'browser');
    const settings = createAISettings(keys, memoryPrefs());
    await settings.load();
    expect(settings.providers.value.find((p) => p.id === 'anthropic')?.key).toBeTruthy();
    const model = createStorageModel();
    model.register(keysStorageSource(keys, settings));
    await model.refresh();
    expect(model.sources.value[0].entries).toHaveLength(1);
    await model.forget('assistant-keys', 'anthropic');
    expect(model.sources.value[0].entries).toEqual([]);
    expect(settings.providers.value.find((p) => p.id === 'anthropic')?.key).toBeUndefined();
    expect(local.length).toBe(0);
    // A new page finds none.
    expect(await createKeyStore(base).list()).toEqual([]);
  });

  it('forgets every key with Forget everything, without a setup too', async () => {
    const { base, local, session } = world();
    const keys = createKeyStore(base);
    await keys.set('anthropic', secret, 'browser');
    await keys.set('openai', secret, 'session');
    const model = createStorageModel();
    model.register(keysStorageSource(keys));
    await model.forgetEverything();
    expect(model.sources.value[0].entries).toEqual([]);
    expect(local.length + session.length).toBe(0);
  });
});
