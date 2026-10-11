/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map } from 'rxjs';
import { Choices } from '../../components/choices/choices.js';
import { Select } from '../../components/select/select.js';
import { TextInput } from '../../components/text-input/text-input.js';
import { For, Show } from '../../dom/control-flow.js';
import type { AIProvider, AISettingsModel } from '../../model/agent.js';
import { observableState } from '../../reactive/observable-state.js';
import { classes } from './_generated/ai-setup.module.css.js';

const KeyCopy = 'Kept only in this browser; it goes only to that provider.';

/**
 * Setting up the assistant's model (as Stories untold does, over any-model-adapters): with no key yet, what it is and
 * a key to add; then the providers with their keys, the model it uses, and one's own OpenAI-compatible endpoint
 * (OpenRouter, Ollama, vLLM, …). The Settings tab's first section; asking without a model opens it.
 */
export function AISetup(props: { settings: AISettingsModel }) {
  const { settings } = props;
  const anyKey = settings.providers.pipe(
    map((ps) => ps.some((p) => p.key)),
    distinctUntilChanged(),
  );
  return (
    <section class={classes.setup} aria-label="AI assistant">
      <Show when={anyKey} fallback={() => <KeyNeeded settings={settings} />}>
        {() => <Configured settings={settings} />}
      </Show>
    </section>
  );
}

/** No key yet: what the assistant is, a provider, its key. */
function KeyNeeded(props: { settings: AISettingsModel }) {
  const { settings } = props;
  const provider = observableState('anthropic');
  let key = '';
  const save = () => settings.setKey(provider.value, key);
  return (
    <div class={classes.needed}>
      <p class={classes.lead}>Add a model key to ask the assistant</p>
      <p class={classes.copy}>
        It runs here with your own key for Claude, Gemini or GPT. The key stays in this browser and goes only to that
        provider. Everything else works without it.
      </p>
      <Choices
        choices={[
          { id: 'anthropic', label: 'Anthropic' },
          { id: 'google', label: 'Gemini' },
          { id: 'openai', label: 'OpenAI' },
        ]}
        active={provider.asReadonly}
        onPick={(id) => provider.set(id)}
      />
      <div class={classes.keyRow}>
        <TextInput
          size="sm"
          type="password"
          class={classes.grow}
          aria-label="API key"
          placeholder={provider.asReadonly.pipe(map((p) => `${labelOf(p)} API key`))}
          onInput={(v) => {
            key = v;
          }}
          onSubmit={save}
        />
        <button type="button" class={classes.action} onClick={save}>
          Save
        </button>
      </div>
      <p class={classes.note}>{KeyCopy}</p>
    </div>
  );
}

const labelOf = (id: string) => ({ anthropic: 'Anthropic', google: 'Gemini', openai: 'OpenAI' })[id] ?? id;

/** Set up: the model in use, each provider's key, one's own endpoint. */
function Configured(props: { settings: AISettingsModel }) {
  const { settings } = props;
  const options = settings.providers.pipe(
    map((ps) =>
      ps.filter((p) => p.key).flatMap((p) => p.models.map((m) => [`${p.id}|${m}`, `${p.label} · ${m}`] as const)),
    ),
  );
  const own = observableState(false);
  return (
    <div class={classes.configured}>
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
      <p class={classes.caption}>Keys</p>
      <div class={classes.providers}>
        <For each={settings.providers} key={(p) => p.id}>
          {(p) => <ProviderRow settings={settings} provider={p} />}
        </For>
      </div>
      <button type="button" class={classes.fold} aria-expanded={own.asReadonly} onClick={() => own.set(!own.value)}>
        <span>Your own endpoint</span>
      </button>
      <Show when={own.asReadonly}>{() => <OwnEndpoint settings={settings} onDone={() => own.set(false)} />}</Show>
      <p class={classes.note}>{KeyCopy} The choice applies from the next request.</p>
    </div>
  );
}

/** A provider: its name, whether it has a key (here, or from the server), and adding or removing one. */
function ProviderRow(props: { settings: AISettingsModel; provider: AIProvider }) {
  const { settings } = props;
  const live = settings.providers.pipe(
    map((ps) => ps.find((x) => x.id === props.provider.id)!),
    distinctUntilChanged(),
  );
  const editing = observableState(false);
  let key = '';
  return (
    <div class={classes.provider}>
      <div class={classes.providerHead}>
        <span class={classes.providerName}>{props.provider.label}</span>
        <span class={classes.providerState} data-has={live.pipe(map((p) => (p.key ? '' : undefined)))}>
          {live.pipe(
            map((p) => (p.key === 'browser' ? 'key kept here' : p.key === 'server' ? 'from the server' : 'no key')),
          )}
        </span>
        <Show
          when={live.pipe(map((p) => p.key === 'browser'))}
          fallback={() => (
            <button type="button" class={classes.link} onClick={() => editing.set(!editing.value)}>
              Add key
            </button>
          )}
        >
          {() => (
            <button type="button" class={classes.link} onClick={() => settings.removeKey(props.provider.id)}>
              Remove
            </button>
          )}
        </Show>
      </div>
      <Show when={editing.asReadonly}>
        {() => (
          <div class={classes.keyRow}>
            <TextInput
              size="sm"
              type="password"
              class={classes.grow}
              aria-label={`${props.provider.label} API key`}
              placeholder={`${props.provider.label} API key`}
              onInput={(v) => {
                key = v;
              }}
              onSubmit={() => {
                settings.setKey(props.provider.id, key);
                editing.set(false);
              }}
            />
            <button
              type="button"
              class={classes.action}
              onClick={() => {
                settings.setKey(props.provider.id, key);
                editing.set(false);
              }}
            >
              Save
            </button>
          </div>
        )}
      </Show>
    </div>
  );
}

/** One's own OpenAI-compatible endpoint (Stories' "your own"): a name, the URL, which API, the key, a model. */
function OwnEndpoint(props: { settings: AISettingsModel; onDone(): void }) {
  const values = { label: '', baseURL: '', key: '', model: '' };
  const api = observableState<'chat' | 'responses'>('chat');
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
          values[name] = v;
        }}
      />
    </div>
  );
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
      {field('key', 'Key', 'Optional for a local server', 'password')}
      {field('model', 'Model', 'e.g. llama-4-scout')}
      <p class={classes.note}>The server must allow this page to call it (CORS).</p>
      <div class={classes.ownFoot}>
        <button
          type="button"
          class={classes.action}
          onClick={() => {
            if (!values.label || !values.baseURL || !values.model) return;
            props.settings.addProvider({ ...values, api: api.value });
            props.onDone();
          }}
        >
          Add
        </button>
      </div>
    </div>
  );
}
