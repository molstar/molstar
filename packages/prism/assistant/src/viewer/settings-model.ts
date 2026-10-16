/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ModelAdapter } from '@molstar/ai-agent';
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
  /** Added with no key (a local server takes none): the scope chosen, since no key is kept to say it. */
  readonly keyless?: KeyScope;
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

export interface AISettingsOptions {
  /** Says what went wrong in words (a key that could not be kept, an endpoint refused); the viewer's notifications. */
  readonly report?: (text: string) => void;
  /** Models the host brings with their own adapters (a story's scripted model): in effect without a key of the person's. */
  readonly provided?: readonly ProvidedModels[];
}

/** A provider the host brings, its adapter made already (no key asked for: its key is the host's). */
export interface ProvidedModels {
  readonly id: string;
  readonly label: string;
  /** Its models, the default first. */
  readonly models: readonly string[];
  readonly adapter: ModelAdapter;
}

/** An own endpoint's URL is https, or http on this machine (model.ts, `AIProvider.baseURL`); never with a login in it. */
export function checkEndpoint(baseURL: string): string | undefined {
  let url: URL;
  try {
    url = new URL(baseURL);
  } catch {
    return `“${baseURL}” is not a URL.`;
  }
  if (url.username || url.password) return 'Give the key in its own field, not in the URL.';
  const local = url.hostname === 'localhost' || url.hostname === '127.0.0.1' || url.hostname === '[::1]';
  if (url.protocol === 'https:' || (url.protocol === 'http:' && local)) return undefined;
  return 'The endpoint must be https (or http on localhost).';
}

export function createAISettings(
  keys: KeyStore,
  prefs: SettingsPrefs = localPrefs(),
  options: AISettingsOptions = {},
): ViewerAISettings {
  const report = options.report ?? ((text: string) => console.warn(`[prism-assistant] ${text}`));
  const scopes = new Map<string, KeyScope>();
  const own = new Map<string, OwnProvider>();
  for (const p of parse<OwnProvider[]>(prefs.get('providers')) ?? []) own.set(p.id, p);
  const provided = new Map((options.provided ?? []).map((p) => [p.id, p] as const));

  const providers = observableState<readonly AIProvider[]>([]);
  // The person's choice (kept in prefs) apart from the one in effect: a choice whose provider has no key (yet) is not
  // lost, the first provider with a key stands in for it, and the choice comes back once its key is there.
  let chosen = parse<AISetup>(prefs.get('setup'));
  let loaded = false;
  const setup = observableState<AISetup | undefined>(undefined);
  const publish = () => {
    providers.set([
      // The host's first: its key is its own (from the server, as far as the person is concerned).
      ...[...provided.values()].map(
        (p): AIProvider => ({ id: p.id, label: p.label, family: 'openai', models: p.models, key: { from: 'server' } }),
      ),
      ...Builtin.map(
        (b): AIProvider => ({
          ...b,
          // The catalog's models; asking each provider's `/models` (refreshModels) is not built yet.
          models: listedModels(b.id).map((m) => m.id),
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
    // Which keys there are is known only once loaded: until then, nothing is in effect.
    if (!loaded) return;
    if (chosen && (scopes.has(chosen.provider) || provided.has(chosen.provider))) setup.set(chosen);
    else {
      const first = providers.value.find((p) => p.key && p.models.length);
      setup.set(first ? { provider: first.id, model: first.models[0] } : undefined);
    }
  };
  publish();

  /** The scopes as the store has them (after a failure, so what shows is what is kept). */
  const resync = async () => {
    try {
      const listed = await keys.list();
      scopes.clear();
      for (const { provider, scope } of listed) scopes.set(provider, scope);
      for (const o of own.values()) if (o.keyless && !scopes.has(o.id)) scopes.set(o.id, o.keyless);
    } catch (e) {
      report(`The stored keys could not be read: ${messageOf(e)}`);
    }
    loaded = true;
    publish();
  };

  const adapters = new Map<string, ModelAdapter>();
  const resolveApiKey = (provider: string) => async () => {
    const key = await keys.get(provider);
    if (key) return key;
    // The adapters refuse an empty key; a local server given none ignores what is sent.
    if (own.get(provider)?.keyless) return KeylessPlaceholder;
    throw new Error(`No key for ${provider}.`);
  };

  const labelOf = (provider: string) =>
    own.get(provider)?.label ?? Builtin.find((b) => b.id === provider)?.label ?? provider;
  /** An own endpoint no longer without a key (one given, or the endpoint's removed). */
  const keyed = (provider: string) => {
    const o = own.get(provider);
    if (!o?.keyless) return;
    own.set(provider, { id: o.id, label: o.label, baseURL: o.baseURL, api: o.api, models: o.models });
    prefs.set('providers', JSON.stringify([...own.values()]));
  };
  const setKey = (provider: string, key: string, scope: KeyScope) => {
    keyed(provider);
    scopes.set(provider, scope);
    adapters.delete(provider);
    publish();
    keys.set(provider, key, scope).catch((e) => {
      report(`The key for ${labelOf(provider)} could not be kept: ${messageOf(e)}`);
      void resync();
    });
  };

  return {
    providers: providers.asReadonly,
    setup: setup.asReadonly,
    load: resync,
    adapterFor(provider) {
      const given = provided.get(provider);
      if (given) return given.adapter;
      if (!scopes.has(provider)) return undefined;
      let adapter = adapters.get(provider);
      if (!adapter) {
        const o = own.get(provider);
        const key = resolveApiKey(provider);
        adapter =
          provider === 'anthropic'
            ? lazyAdapter(async () =>
                (await import('@molstar/ai-model-adapters/anthropic/adapter')).createAnthropicModelAdapter({
                  resolveApiKey: key,
                }),
              )
            : provider === 'google'
              ? lazyAdapter(async () =>
                  (await import('@molstar/ai-model-adapters/google/adapter')).createGoogleModelAdapter({
                    resolveApiKey: key,
                  }),
                )
              : lazyAdapter(async () =>
                  (await import('@molstar/ai-model-adapters/openai/adapter')).createOpenAIModelAdapter({
                    resolveApiKey: key,
                    ...(o ? { baseURL: o.baseURL, api: o.api, provider: o.id, label: o.label } : {}),
                  }),
                );
        adapters.set(provider, adapter);
      }
      return adapter;
    },
    setKey,
    removeKey(provider) {
      keyed(provider);
      scopes.delete(provider);
      adapters.delete(provider);
      publish();
      keys.remove(provider).catch((e) => {
        report(`The key for ${labelOf(provider)} could not be removed: ${messageOf(e)}`);
        void resync();
      });
    },
    choose(next) {
      chosen = next;
      prefs.set('setup', JSON.stringify(next));
      publish();
    },
    async refreshModels() {
      // Not built yet: the catalog's models stand (asking each provider's `/models` comes with their adapters).
      publish();
    },
    async listModels() {
      // Not built yet: an own endpoint's models are typed in.
      return [];
    },
    addProvider(p) {
      const problem = checkEndpoint(p.baseURL);
      if (problem) {
        report(`${p.label || 'The endpoint'} was not added: ${problem}`);
        return;
      }
      // Its id from its name; another endpoint whose name reads alike keeps its own (the same name replaces it).
      const base = `own-${p.label.toLowerCase().replace(/[^a-z0-9]+/g, '-')}`;
      let id = base;
      for (let n = 2; own.has(id) && own.get(id)!.label !== p.label; n++) id = `${base}-${n}`;
      // A local server takes no key: none is kept (an empty one would show as kept, and fail every turn).
      const keyless = p.key.trim().length === 0;
      own.set(id, {
        id,
        label: p.label,
        baseURL: p.baseURL,
        api: p.api,
        models: p.models,
        ...(keyless ? { keyless: p.scope } : {}),
      });
      prefs.set('providers', JSON.stringify([...own.values()]));
      if (!keyless) return setKey(id, p.key, p.scope);
      scopes.set(id, p.scope);
      adapters.delete(id);
      publish();
      // A key this endpoint had before is dropped.
      keys.remove(id).catch((e) => {
        report(`The key for ${labelOf(id)} could not be removed: ${messageOf(e)}`);
        void resync();
      });
    },
  };
}

const KeylessPlaceholder = 'no-key';

const messageOf = (e: unknown) => (e instanceof Error ? e.message : String(e));

function parse<T>(text: string | undefined): T | undefined {
  if (!text) return undefined;
  try {
    return JSON.parse(text) as T;
  } catch {
    return undefined;
  }
}

/**
 * A provider's adapter whose SDK is imported on its first use: a person with one provider's key does not download the
 * others' (each its own chunk). A failed import is tried again on the next use.
 */
export function lazyAdapter(load: () => Promise<ModelAdapter>): ModelAdapter {
  let loading: Promise<ModelAdapter> | undefined;
  const adapter = () => {
    loading ??= load().catch((e) => {
      loading = undefined;
      throw e;
    });
    return loading;
  };
  return {
    capabilitiesFor: async (configuration) => (await adapter()).capabilitiesFor(configuration),
    project: async (input) => (await adapter()).project(input),
    async *invoke(input) {
      yield* (await adapter()).invoke(input);
    },
  };
}
