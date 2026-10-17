/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Icon, Select, TextInput, For, Show, type ObservableState, observableState } from '@molstar/prism-ui';
import { combineLatest, map } from 'rxjs';
import { type AIProvider, type AISettingsModel, checkLimits, type KeyScope } from '../model.js';
import { classes } from './_generated/ai-setup.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
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
 * Setting up the assistant's model (as Stories untold does): every provider with its key (one
 * for each, as many as one likes), how long a key is kept, the model in use, and one's own OpenAI-compatible endpoint
 * (OpenRouter, Ollama, vLLM, …). With no key yet, a word on what it is first, and the first provider's key open.
 */
export function AISetup(props: { settings: AISettingsModel }) {
  const { settings } = props;
  const anyKey = derive(settings.providers, (ps) => ps.some((p) => p.key));
  const options = derive(settings.providers, (ps) =>
    ps.filter((p) => p.key).flatMap((p) => p.models.map((m) => [`${p.id}|${m}`, `${p.label} · ${m}`] as const)),
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
            <p class={classes.lead}>Add a model key to ask AI assistant</p>
            <p class={classes.copy}>
              It runs here with your own key for Claude, Gemini, GPT, or a server of your own; add as many as you like
              and pick the model. Everything else works without it.
            </p>
          </div>
        )}
      </Show>
      <Show when={anyKey}>
        {() => (
          <div class={`${classes.row} ${classes.modelRow}`}>
            <span class={classes.label}>Model</span>
            <Show when={options}>
              {(opts) => (
                <Select
                  over
                  size="sm"
                  class={classes.grow}
                  aria-label="Model"
                  options={opts}
                  value={derive(settings.setup, (s) => (s ? `${s.provider}|${s.model}` : ''))}
                  onChange={(v) => {
                    const [provider, model] = v.split('|');
                    settings.choose({ provider, model });
                  }}
                />
              )}
            </Show>
            <RefreshModels settings={settings} />
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
      <Show when={own.asReadonly} reveal>
        {() => <OwnEndpoint settings={settings} onDone={() => own.set(false)} />}
      </Show>
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
  const live = derive(settings.providers, (ps) => ps.find((x) => x.id === props.provider.id)!);
  const editing = derive(props.editing.asReadonly, (id) => id === provider.id);
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
              {combineLatest([live, editing]).pipe(map(([p, e]) => (e ? 'Cancel' : p.key ? 'Change' : 'Add key')))}
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
      <Show when={editing} reveal>
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
      <ScopeChoice scope={scope} />
    </div>
  );
}

/** How long a key is kept: until reload, this session, this browser. */
function ScopeChoice(props: { scope: ObservableState<KeyScope> }) {
  return (
    <div class={classes.scopes} role="radiogroup" aria-label="Keep the key">
      <span class={classes.scopesLabel}>Keep</span>
      {Scopes.map(([id, label]) => (
        <button
          type="button"
          role="radio"
          class={classes.scope}
          aria-checked={derive(props.scope.asReadonly, (s) => s === id)}
          onClick={() => props.scope.set(id)}
        >
          {label}
        </button>
      ))}
    </div>
  );
}

/** Asks the providers with a key for the models they offer now (their `/models`). */
function RefreshModels(props: { settings: AISettingsModel }) {
  const busy = observableState(false);
  return (
    <button
      type="button"
      class={classes.refresh}
      aria-label="Ask the providers for their models"
      data-tip="Ask the providers for their models"
      disabled={busy.asReadonly}
      data-busy={derive(busy.asReadonly, (b) => (b ? '' : undefined))}
      onClick={async () => {
        busy.set(true);
        try {
          await props.settings.refreshModels();
        } finally {
          busy.set(false);
        }
      }}
    >
      <Icon name="refresh-cw" size={14} />
    </button>
  );
}

/**
 * One's own OpenAI-compatible endpoint (Stories' "your own"): a name, the URL, which OpenAI API it speaks, its key and
 * how long to keep it, its models (typed, or fetched from its `/models`), and optionally their context window and max
 * output (32,768 and 4,096 tokens when left empty).
 */
function OwnEndpoint(props: { settings: AISettingsModel; onDone(): void }) {
  const values = { label: '', baseURL: '', key: '', contextWindow: '', maxOutputTokens: '' };
  const problem = observableState<string | undefined>(undefined);
  const api = observableState<'chat' | 'responses'>('chat');
  const scope = observableState<KeyScope>('browser');
  const models = observableState('');
  const fetching = observableState<'idle' | 'busy' | 'failed'>('idle');
  const field = (name: keyof typeof values, label: string, placeholder: string, type: 'text' | 'password' = 'text') => (
    <div class={classes.row}>
      <span class={classes.label}>{label}</span>
      <TextInput
        size="sm"
        type={type}
        class={classes.grow}
        aria-label={label}
        placeholder={placeholder}
        onInput={(v) => {
          values[name] = v.trim();
        }}
      />
    </div>
  );
  const fetchModels = async () => {
    fetching.set('busy');
    try {
      const ids = await props.settings.listModels({ baseURL: values.baseURL, key: values.key });
      models.set(ids.join(', '));
      fetching.set('idle');
    } catch {
      fetching.set('failed');
    }
  };
  const add = () => {
    const list = models.value
      .split(/[,\s]+/)
      .map((m) => m.trim())
      .filter(Boolean);
    if (!values.label || !values.baseURL || !list.length) return;
    // Empty is the default; anything else is a whole number of tokens.
    const tokens = (text: string) => (text ? Number(text.replaceAll(/[,_\s]/g, '')) : undefined);
    const limits = { contextWindow: tokens(values.contextWindow), maxOutputTokens: tokens(values.maxOutputTokens) };
    const wrong = checkLimits(limits);
    problem.set(wrong);
    if (wrong) return;
    const { label, baseURL, key } = values;
    props.settings.addProvider({
      label,
      baseURL,
      key,
      scope: scope.value,
      api: api.value,
      models: list,
      ...(limits.contextWindow !== undefined ? { contextWindow: limits.contextWindow } : {}),
      ...(limits.maxOutputTokens !== undefined ? { maxOutputTokens: limits.maxOutputTokens } : {}),
    });
    props.onDone();
  };
  return (
    <div class={classes.own}>
      {field('label', 'Name', 'e.g. Lab')}
      {field('baseURL', 'Endpoint', 'https://…/v1, or http://localhost:11434/v1')}
      <div class={classes.row}>
        <span class={classes.label}>API</span>
        <Select
          over
          size="sm"
          class={classes.grow}
          aria-label="API"
          options={[
            ['chat', 'OpenAI chat completions'],
            ['responses', 'OpenAI Responses'],
          ]}
          value={api.asReadonly}
          onChange={(v) => api.set(v as 'chat' | 'responses')}
        />
      </div>
      {field('key', 'Key', 'Optional for a local server', 'password')}
      <ScopeChoice scope={scope} />
      <div class={classes.row}>
        <span class={classes.label}>Models</span>
        <div class={classes.keyRow}>
          <TextInput
            size="sm"
            class={classes.grow}
            aria-label="Models"
            placeholder="e.g. llama-4-scout, qwen3-coder"
            value={models.asReadonly}
            onInput={(v) => models.set(v)}
          />
          <button
            type="button"
            class={classes.link}
            disabled={derive(fetching.asReadonly, (f) => f === 'busy')}
            onClick={fetchModels}
          >
            {derive(fetching.asReadonly, (f) => (f === 'busy' ? 'Fetching…' : 'Fetch'))}
          </button>
        </div>
      </div>
      <Show when={derive(fetching.asReadonly, (f) => f === 'failed')}>
        {() => (
          <p class={classes.error}>Could not list its models: check the endpoint, the key, and that it allows CORS.</p>
        )}
      </Show>
      {field('contextWindow', 'Context window', 'Optional: 32,768 tokens if empty')}
      {field('maxOutputTokens', 'Max output', 'Optional: 4,096 tokens if empty')}
      <Show when={derive(problem.asReadonly, (p) => p !== undefined)}>
        {() => <p class={classes.error}>{derive(problem.asReadonly, (p) => p ?? '')}</p>}
      </Show>
      <p class={classes.note}>It must speak an OpenAI API, and allow this page to call it (CORS).</p>
      <div class={classes.ownFoot}>
        <button type="button" class={classes.action} onClick={add}>
          Add
        </button>
      </div>
    </div>
  );
}
