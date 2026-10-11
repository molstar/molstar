/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map } from 'rxjs';
import { Select } from '../../components/select/select.js';
import { TextInput } from '../../components/text-input/text-input.js';
import { For, Show } from '../../dom/control-flow.js';
import type { AIProvider, AISettingsModel, KeyScope } from '../../model/agent.js';
import { type ObservableState, observableState } from '../../reactive/observable-state.js';
import { classes } from './_generated/ai-setup.module.css.js';

const Scopes: readonly (readonly [KeyScope, string])[] = [
  ['reload', 'Until reload'],
  ['session', 'This session'],
  ['browser', 'This browser'],
];

const scopeText: Record<KeyScope, string> = {
  reload: 'key until reload',
  session: 'key for this session',
  browser: 'key kept in this browser',
};

/**
 * Setting up the assistant's model (as Stories untold does, over any-model-adapters): every provider with its key (one
 * for each, as many as one likes), how long a key is kept, the model in use, and one's own OpenAI-compatible endpoint
 * (OpenRouter, Ollama, vLLM, …). With no key yet, a word on what it is first, and the first provider's key open.
 */
export function AISetup(props: { settings: AISettingsModel }) {
  const { settings } = props;
  const anyKey = settings.providers.pipe(
    map((ps) => ps.some((p) => p.key)),
    distinctUntilChanged(),
  );
  const options = settings.providers.pipe(
    map((ps) =>
      ps.filter((p) => p.key).flatMap((p) => p.models.map((m) => [`${p.id}|${m}`, `${p.label} · ${m}`] as const)),
    ),
  );
  const own = observableState(false);
  // Which provider's key is being entered (here: rows are made anew when their provider changes). No key yet: the
  // first provider's, so the way in is obvious.
  const editing = observableState<string | undefined>(
    anyKeyNow(settings) ? undefined : settings.providers.value[0]?.id,
  );
  return (
    <section class={classes.setup} aria-label="AI assistant">
      <Show when={anyKey.pipe(map((k) => !k))}>
        {() => (
          <div class={classes.intro}>
            <p class={classes.lead}>Add a model key to ask the assistant</p>
            <p class={classes.copy}>
              It runs here with your own key for Claude, Gemini, GPT, or a server of your own; add as many as you like
              and pick the model. Everything else works without it.
            </p>
          </div>
        )}
      </Show>
      <Show when={anyKey}>
        {() => (
          <div class={classes.row}>
            <span class={classes.label}>Model</span>
            <Show when={options}>
              {(opts) => (
                <Select
                  size="sm"
                  class={classes.grow}
                  aria-label="Model"
                  options={opts}
                  value={settings.setup.pipe(map((s) => (s ? `${s.provider}|${s.model}` : '')))}
                  onChange={(v) => {
                    const [provider, model] = v.split('|');
                    settings.choose({ provider, model });
                  }}
                />
              )}
            </Show>
          </div>
        )}
      </Show>
      <p class={classes.caption}>Keys</p>
      <div class={classes.providers}>
        <For each={settings.providers} key={(p) => p.id}>
          {(p) => <ProviderRow settings={settings} provider={p} editing={editing} />}
        </For>
      </div>
      <button type="button" class={classes.fold} aria-expanded={own.asReadonly} onClick={() => own.set(!own.value)}>
        <span>Your own endpoint</span>
      </button>
      <Show when={own.asReadonly}>{() => <OwnEndpoint settings={settings} onDone={() => own.set(false)} />}</Show>
      <p class={classes.note}>
        Keys go only to their provider. Kept ones are stored scrambled, so reading the browser’s storage as it is gives
        no working key (it does not guard against code running in this page). The model applies from the next request.
      </p>
    </section>
  );
}

const anyKeyNow = (settings: AISettingsModel) => settings.providers.value.some((p) => p.key);

/** A provider: its name, its key's state (how long it is kept, or from the server), and adding, changing, removing it. */
function ProviderRow(props: {
  settings: AISettingsModel;
  provider: AIProvider;
  editing: ObservableState<string | undefined>;
}) {
  const { settings, provider } = props;
  const live = settings.providers.pipe(
    map((ps) => ps.find((x) => x.id === props.provider.id)!),
    distinctUntilChanged(),
  );
  const editing = props.editing.asReadonly.pipe(
    map((id) => id === provider.id),
    distinctUntilChanged(),
  );
  const toggle = () => props.editing.set(props.editing.value === provider.id ? undefined : provider.id);
  return (
    <div class={classes.provider}>
      <div class={classes.providerHead}>
        <span class={classes.providerName}>{props.provider.label}</span>
        <span class={classes.providerState} data-has={live.pipe(map((p) => (p.key ? '' : undefined)))}>
          {live.pipe(
            map((p) => (!p.key ? 'no key' : p.key.from === 'server' ? 'key from the server' : scopeText[p.key.scope])),
          )}
        </span>
        <Show when={live.pipe(map((p) => p.key?.from !== 'server'))}>
          {() => (
            <button type="button" class={classes.link} onClick={toggle}>
              {live.pipe(map((p) => (p.key ? 'Change' : 'Add key')))}
            </button>
          )}
        </Show>
        <Show when={live.pipe(map((p) => p.key?.from === 'here'))}>
          {() => (
            <button type="button" class={classes.link} onClick={() => settings.removeKey(props.provider.id)}>
              Remove
            </button>
          )}
        </Show>
      </div>
      <Show when={editing}>
        {() => (
          <KeyForm
            label={`${props.provider.label} API key`}
            onSave={(key, scope) => {
              settings.setKey(props.provider.id, key, scope);
              props.editing.set(undefined);
            }}
          />
        )}
      </Show>
    </div>
  );
}

/** A key and how long to keep it. */
function KeyForm(props: {
  label: string;
  action?: string;
  optional?: boolean;
  onSave(key: string, scope: KeyScope): void;
}) {
  let key = '';
  const scope = observableState<KeyScope>('browser');
  const save = () => {
    if (key.trim() || props.optional) props.onSave(key.trim(), scope.value);
  };
  return (
    <div class={classes.keyForm}>
      <div class={classes.keyRow}>
        <TextInput
          size="sm"
          type="password"
          class={classes.grow}
          aria-label={props.label}
          placeholder={props.label}
          onInput={(v) => {
            key = v;
          }}
          onSubmit={save}
        />
        <button type="button" class={classes.action} onClick={save}>
          {props.action ?? 'Save'}
        </button>
      </div>
      <div class={classes.scopes} role="radiogroup" aria-label="Keep the key">
        <span class={classes.scopesLabel}>Keep</span>
        {Scopes.map(([id, label]) => (
          <button
            type="button"
            role="radio"
            class={classes.scope}
            aria-checked={scope.asReadonly.pipe(map((s) => s === id))}
            onClick={() => scope.set(id)}
          >
            {label}
          </button>
        ))}
      </div>
    </div>
  );
}

/** One's own OpenAI-compatible endpoint (Stories' "your own"): a name, the URL, which API, the key, a model. */
function OwnEndpoint(props: { settings: AISettingsModel; onDone(): void }) {
  const values = { label: '', baseURL: '', model: '' };
  const api = observableState<'chat' | 'responses'>('chat');
  const field = (name: keyof typeof values, label: string, placeholder: string) => (
    <div class={classes.row}>
      <span class={classes.label}>{label}</span>
      <TextInput
        size="sm"
        class={classes.grow}
        aria-label={label}
        placeholder={placeholder}
        onInput={(v) => {
          values[name] = v;
        }}
      />
    </div>
  );
  const add = (key: string, scope: KeyScope) => {
    if (!values.label || !values.baseURL || !values.model) return;
    props.settings.addProvider({ ...values, key, scope, api: api.value });
    props.onDone();
  };
  return (
    <div class={classes.own}>
      {field('label', 'Name', 'e.g. Lab')}
      {field('baseURL', 'Endpoint', 'https://…/v1, or http://localhost:11434/v1')}
      <div class={classes.row}>
        <span class={classes.label}>API</span>
        <Select
          size="sm"
          class={classes.grow}
          aria-label="API"
          options={[
            ['chat', 'Chat completions'],
            ['responses', 'Responses'],
          ]}
          value={api.asReadonly}
          onChange={(v) => api.set(v as 'chat' | 'responses')}
        />
      </div>
      {field('model', 'Model', 'e.g. llama-4-scout')}
      <p class={classes.note}>
        The key is optional for a local server; the server must allow this page to call it (CORS).
      </p>
      <KeyForm label="API key (optional)" action="Add" optional onSave={add} />
    </div>
  );
}
