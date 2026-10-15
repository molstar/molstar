/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { KeyScope, KeyStore } from '../../keys.js';
import { checkEndpoint, createAISettings, type SettingsPrefs } from '../settings-model.js';

// The assistant's model setup (agent.md §4): the choice kept apart from the one in effect, failures said in words.

function memoryKeys(initial: Record<string, KeyScope> = {}, fail?: 'set' | 'list'): KeyStore {
  const keys = new Map<string, { key: string; scope: KeyScope }>(
    Object.entries(initial).map(([p, scope]) => [p, { key: `key-${p}`, scope }]),
  );
  return {
    get: async (p) => keys.get(p)?.key,
    async set(p, key, scope) {
      if (fail === 'set') throw new Error('storage refused');
      keys.set(p, { key, scope });
    },
    remove: async (p) => void keys.delete(p),
    async list() {
      if (fail === 'list') throw new Error('storage refused');
      return [...keys].map(([provider, { scope }]) => ({ provider, scope }));
    },
    protection: async () => 'per-browser',
  };
}

function memoryPrefs(initial: Record<string, string> = {}): SettingsPrefs & { values: Map<string, string> } {
  const values = new Map(Object.entries(initial));
  return { values, get: (k) => values.get(k), set: (k, v) => void values.set(k, v) };
}

const until = async (check: () => boolean, ms = 2000) => {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error(`timed out after ${ms} ms`);
    await new Promise((r) => setTimeout(r, 5));
  }
};

describe('the assistant setup', () => {
  it('keeps the chosen provider and model across a load', async () => {
    const choice = { provider: 'openai', model: 'gpt-chosen' };
    const prefs = memoryPrefs({ setup: JSON.stringify(choice) });
    const settings = createAISettings(memoryKeys({ anthropic: 'browser', openai: 'browser' }), prefs);
    await settings.load();
    expect(settings.setup.value).toEqual(choice);
  });

  it('stands in for a choice whose provider has no key, without forgetting it', async () => {
    const choice = { provider: 'openai', model: 'gpt-chosen' };
    const prefs = memoryPrefs({ setup: JSON.stringify(choice) });
    const settings = createAISettings(memoryKeys({ anthropic: 'browser' }), prefs);
    await settings.load();
    expect(settings.setup.value?.provider).toBe('anthropic');
    expect(JSON.parse(prefs.values.get('setup')!)).toEqual(choice);
    // Its key comes: the choice is back.
    settings.setKey('openai', 'sk-x', 'browser');
    expect(settings.setup.value).toEqual(choice);
  });

  it('says when a key cannot be kept, and shows what is kept', async () => {
    const said: string[] = [];
    const settings = createAISettings(memoryKeys({}, 'set'), memoryPrefs(), { report: (t) => said.push(t) });
    await settings.load();
    settings.setKey('anthropic', 'sk-x', 'browser');
    await until(() => said.length > 0);
    expect(said[0]).toMatch(/key for Anthropic could not be kept/);
    await until(() => !settings.providers.value.find((p) => p.id === 'anthropic')?.key);
    expect(said.join()).not.toContain('sk-x');
  });

  it('says when the stored keys cannot be read', async () => {
    const said: string[] = [];
    const settings = createAISettings(memoryKeys({}, 'list'), memoryPrefs(), { report: (t) => said.push(t) });
    await settings.load();
    expect(said[0]).toMatch(/could not be read/);
    expect(settings.setup.value).toBeUndefined();
  });

  it('adds an own endpoint only over https or local http, and keeps two alike names apart', async () => {
    const said: string[] = [];
    const settings = createAISettings(memoryKeys(), memoryPrefs(), { report: (t) => said.push(t) });
    await settings.load();
    const add = (label: string, baseURL: string) =>
      settings.addProvider({ label, baseURL, api: 'chat', key: 'k', scope: 'reload', models: ['m'] });
    add('Remote', 'http://example.org/v1');
    expect(said[0]).toMatch(/must be https/);
    add('Login', 'https://user:pass@example.org/v1');
    expect(said[1]).toMatch(/not in the URL/);
    add('My API', 'https://a.example/v1');
    add('my-api', 'http://localhost:11434/v1');
    const own = settings.providers.value.filter((p) => p.baseURL);
    expect(own.map((p) => [p.id, p.baseURL])).toEqual([
      ['own-my-api', 'https://a.example/v1'],
      ['own-my-api-2', 'http://localhost:11434/v1'],
    ]);
    expect(checkEndpoint('http://[::1]:8080/v1')).toBeUndefined();
  });
});
