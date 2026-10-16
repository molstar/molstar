/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ModelAdapter } from '@molstar/ai';
import { describe, expect, it, vi } from 'vitest';
import type { KeyScope, KeyStore } from '../../keys.js';
import { createAISettings, type SettingsPrefs } from '../settings-model.js';

// An own endpoint added with no key (a local server): usable, with no empty key kept (agent.md §4).

const made: { resolveApiKey: () => string | Promise<string> }[] = [];
vi.mock('@molstar/ai/providers/openai', () => ({
  createOpenAIModelAdapter: (options: { resolveApiKey: () => string | Promise<string> }) => {
    made.push(options);
    return {
      capabilitiesFor: async () => ({}),
      project: async () => ({}),
      async *invoke() {},
    } as unknown as ModelAdapter;
  },
}));

function memoryKeys(): KeyStore & { stored: Map<string, string> } {
  const stored = new Map<string, string>();
  const scopes = new Map<string, KeyScope>();
  return {
    stored,
    get: async (p) => stored.get(p),
    async set(p, key, scope) {
      stored.set(p, key);
      scopes.set(p, scope);
    },
    async remove(p) {
      stored.delete(p);
      scopes.delete(p);
    },
    list: async () => [...scopes].map(([provider, scope]) => ({ provider, scope })),
    protection: async () => 'per-browser',
  };
}

function memoryPrefs(): SettingsPrefs {
  const values = new Map<string, string>();
  return { get: (k) => values.get(k), set: (k, v) => void values.set(k, v) };
}

const add = (settings: ReturnType<typeof createAISettings>, key: string) =>
  settings.addProvider({
    label: 'Ollama',
    baseURL: 'http://localhost:11434/v1',
    api: 'chat',
    key,
    scope: 'browser',
    models: ['llama'],
  });

describe('an own endpoint with no key', () => {
  it('is in effect, keeps no empty key, and runs turns with a placeholder', async () => {
    const keys = memoryKeys();
    const prefs = memoryPrefs();
    const settings = createAISettings(keys, prefs);
    await settings.load();
    add(settings, '');
    await Promise.resolve();
    expect(keys.stored.size).toBe(0);
    expect(settings.providers.value.find((p) => p.id === 'own-ollama')?.key).toEqual({
      from: 'here',
      scope: 'browser',
    });
    expect(settings.setup.value).toEqual({ provider: 'own-ollama', model: 'llama' });

    await settings.adapterFor('own-ollama')!.capabilitiesFor({} as never);
    expect(await made.at(-1)!.resolveApiKey()).toBeTruthy();

    // After a reload it is still there.
    const again = createAISettings(keys, prefs);
    await again.load();
    expect(again.providers.value.find((p) => p.id === 'own-ollama')?.key).toEqual({ from: 'here', scope: 'browser' });
    expect(again.adapterFor('own-ollama')).toBeDefined();
  });

  it('is without a key once the person removes it, and uses a key given later', async () => {
    const keys = memoryKeys();
    const settings = createAISettings(keys, memoryPrefs());
    await settings.load();
    add(settings, ' ');
    settings.removeKey('own-ollama');
    expect(settings.providers.value.find((p) => p.id === 'own-ollama')?.key).toBeUndefined();
    expect(settings.adapterFor('own-ollama')).toBeUndefined();

    add(settings, '');
    settings.setKey('own-ollama', 'sk-real', 'browser');
    await settings.adapterFor('own-ollama')!.capabilitiesFor({} as never);
    await vi.waitFor(async () => expect(await made.at(-1)!.resolveApiKey()).toBe('sk-real'));
  });
});
