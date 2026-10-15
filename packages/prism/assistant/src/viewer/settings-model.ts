/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ModelAdapter } from '@molstar/ai-agent';
import {
  createAnthropicModelAdapter,
  createGoogleModelAdapter,
  createOpenAIModelAdapter,
} from '@molstar/ai-model-adapters';
import { listedModels } from '@molstar/ai-model-adapters/model-catalog';
import { observableState } from '@molstar/prism-ui/core';
import type { KeyStore } from '../keys.js';
import type { AIProvider, AISettingsModel, AISetup, KeyScope } from '../model.js';

// The assistant's model setup (agent.md §4): the providers, their keys (the person's, kept as they chose, always
// scrambled: keys.ts), the model chosen, one's own OpenAI-compatible endpoints. It gives the agent its adapters: a lookup,
// so a key added while the conversation is open is found (`ModelSource.adapterFor`).

const Builtin: readonly Omit<AIProvider, 'models' | 'key'>[] = [
  { id: 'anthropic', label: 'Anthropic', family: 'anthropic' },
  { id: 'google', label: 'Google', family: 'google' },
  { id: 'openai', label: 'OpenAI', family: 'openai' },
];

interface OwnProvider {
  readonly id: string;
  readonly label: string;
  readonly baseURL: string;
  readonly api: 'chat' | 'responses';
  readonly models: readonly string[];
}

export interface ViewerAISettings extends AISettingsModel {
  /** The adapter for a provider with a key (the agent's model source). */
  adapterFor(provider: string): ModelAdapter | undefined;
  /** Reads the stored keys (once, at start: which providers have one). */
  load(): Promise<void>;
}

/** Where the person's choices (not keys) are kept: the model chosen, their own endpoints. */
export interface SettingsPrefs {
  get(key: string): string | undefined;
  set(key: string, value: string): void;
}

export function localPrefs(prefix = 'prism-assistant.'): SettingsPrefs {
  return {
    get: (key) => {
      try {
        return globalThis.localStorage?.getItem(prefix + key) ?? undefined;
      } catch {
        return undefined;
      }
    },
    set: (key, value) => {
      try {
        globalThis.localStorage?.setItem(prefix + key, value);
      } catch {
        // A private window: the choice lasts for the page.
      }
    },
  };
}

export function createAISettings(keys: KeyStore, prefs: SettingsPrefs = localPrefs()): ViewerAISettings {
  const scopes = new Map<string, KeyScope>();
  const own = new Map<string, OwnProvider>();
  for (const p of parse<OwnProvider[]>(prefs.get('providers')) ?? []) own.set(p.id, p);
  const models = new Map<string, readonly string[]>();

  const providers = observableState<readonly AIProvider[]>([]);
  const setup = observableState<AISetup | undefined>(parse<AISetup>(prefs.get('setup')));
  const publish = () => {
    providers.set([
      ...Builtin.map(
        (b): AIProvider => ({
          ...b,
          models: models.get(b.id) ?? listedModels(b.id).map((m) => m.id),
          ...(scopes.has(b.id) ? { key: { from: 'here' as const, scope: scopes.get(b.id)! } } : {}),
        }),
      ),
      ...[...own.values()].map(
        (p): AIProvider => ({
          id: p.id,
          label: p.label,
          family: 'openai',
          baseURL: p.baseURL,
          models: p.models,
          ...(scopes.has(p.id) ? { key: { from: 'here' as const, scope: scopes.get(p.id)! } } : {}),
        }),
      ),
    ]);
    // Without a choice, the first provider with a key and its first model.
    if (!setup.value || !scopes.has(setup.value.provider)) {
      const first = providers.value.find((p) => p.key && p.models.length);
      setup.set(first ? { provider: first.id, model: first.models[0] } : undefined);
    }
  };
  publish();

  const adapters = new Map<string, ModelAdapter>();
  const resolveApiKey = (provider: string) => async () => {
    const key = await keys.get(provider);
    if (!key) throw new Error(`No key for ${provider}.`);
    return key;
  };

  return {
    providers: providers.asReadonly,
    setup: setup.asReadonly,
    async load() {
      for (const { provider, scope } of await keys.list()) scopes.set(provider, scope);
      publish();
    },
    adapterFor(provider) {
      if (!scopes.has(provider)) return undefined;
      let adapter = adapters.get(provider);
      if (!adapter) {
        const o = own.get(provider);
        adapter =
          provider === 'anthropic'
            ? createAnthropicModelAdapter({ resolveApiKey: resolveApiKey(provider) })
            : provider === 'google'
              ? createGoogleModelAdapter({ resolveApiKey: resolveApiKey(provider) })
              : createOpenAIModelAdapter({
                  resolveApiKey: resolveApiKey(provider),
                  ...(o ? { baseURL: o.baseURL, api: o.api, provider: o.id, label: o.label } : {}),
                });
        adapters.set(provider, adapter);
      }
      return adapter;
    },
    setKey(provider, key, scope) {
      scopes.set(provider, scope);
      adapters.delete(provider);
      void keys.set(provider, key, scope);
      publish();
    },
    removeKey(provider) {
      scopes.delete(provider);
      adapters.delete(provider);
      void keys.remove(provider);
      publish();
    },
    choose(next) {
      setup.set(next);
      prefs.set('setup', JSON.stringify(next));
    },
    async refreshModels() {
      // The catalog's models for now; asking each provider's `/models` comes with their adapters.
      publish();
    },
    async listModels() {
      return [];
    },
    addProvider(p) {
      const id = `own-${p.label.toLowerCase().replace(/[^a-z0-9]+/g, '-')}`;
      own.set(id, { id, label: p.label, baseURL: p.baseURL, api: p.api, models: p.models });
      prefs.set('providers', JSON.stringify([...own.values()]));
      this.setKey(id, p.key, p.scope);
    },
  };
}

function parse<T>(text: string | undefined): T | undefined {
  if (!text) return undefined;
  try {
    return JSON.parse(text) as T;
  } catch {
    return undefined;
  }
}
