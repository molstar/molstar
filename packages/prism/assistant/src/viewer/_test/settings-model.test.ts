/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Provider, ModelFrame } from '@molstar/ai';
import { describe, expect, it, vi } from 'vitest';
import type { KeyScope, KeyStore } from '../../keys.js';
import { checkEndpoint, createAISettings, lazyProvider, type SettingsPrefs } from '../settings-model.js';

// The assistant's model setup (_spec/README.md, "The setup"): the choice kept apart from the one in effect, failures
// said in words.

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

  it("takes the host's preferred model while the person has chosen none and its provider has a key", async () => {
    const preferred = { provider: 'google', model: 'gemini-3.8-flash' };
    const settings = createAISettings(memoryKeys({ anthropic: 'browser', google: 'browser' }), memoryPrefs(), {
      preferred,
    });
    await settings.load();
    expect(settings.setup.value).toEqual(preferred);
    // No key for it: the first provider with one stands in.
    const without = createAISettings(memoryKeys({ anthropic: 'browser' }), memoryPrefs(), { preferred });
    await without.load();
    expect(without.setup.value?.provider).toBe('anthropic');
    // The person's choice wins.
    settings.choose({ provider: 'anthropic', model: 'claude-haiku-5-5' });
    expect(settings.setup.value?.provider).toBe('anthropic');
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

  it("takes a host's models with their provider, in effect with no key of the person's (a story's scripted model)", async () => {
    const given = {} as never;
    const settings = createAISettings(memoryKeys(), memoryPrefs(), {
      provided: [{ id: 'scripted', label: 'Scripted', models: ['scripted'], provider: given }],
    });
    await settings.load();
    expect(settings.providers.value[0]).toMatchObject({ id: 'scripted', key: { from: 'server' } });
    expect(settings.setup.value).toEqual({ provider: 'scripted', model: 'scripted' });
    expect(settings.providerFor('scripted')).toBe(given);
    expect(settings.providerFor('anthropic')).toBeUndefined();
  });

  it('gives a provider without importing its SDK until it is used', async () => {
    const settings = createAISettings(memoryKeys({ anthropic: 'browser' }), memoryPrefs());
    await settings.load();
    const found = settings.providerFor('anthropic')!;
    expect(settings.providerFor('anthropic')).toBe(found);
    const capabilities = await found.capabilities({
      target: { provider: 'anthropic', model: 'claude-sonnet-4-5' },
      parameters: {},
    } as never);
    expect(capabilities.contextWindowTokens).toBeGreaterThan(0);
  });

  it('loads a lazy provider once, on its first use, and again after a failed load', async () => {
    const frame = { kind: 'text-delta', text: 'hi' } as ModelFrame;
    const inner: Provider = {
      capabilities: async () => ({ contextWindowTokens: 1000 }) as never,
      project: async () => ({ payload: {} }) as never,
      async *invoke() {
        yield frame;
      },
    };
    const load = vi.fn<() => Promise<Provider>>().mockRejectedValueOnce(new Error('chunk failed'));
    load.mockResolvedValue(inner);
    const lazy = lazyProvider(load);
    expect(load).not.toHaveBeenCalled();
    await expect(lazy.capabilities({} as never)).rejects.toThrow('chunk failed');
    expect((await lazy.capabilities({} as never)).contextWindowTokens).toBe(1000);
    await lazy.project({} as never);
    const frames: ModelFrame[] = [];
    for await (const f of lazy.invoke({} as never)) frames.push(f);
    expect(frames).toEqual([frame]);
    expect(load).toHaveBeenCalledTimes(2);
  });
});
