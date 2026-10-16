/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Message, Provider } from '../../messages/index.js';
import { type ModelStep, type ScriptedModel, type ScriptedModelOptions, scriptedModel } from '../../testing/index.js';
import type { Clock } from '../../harness/index.js';
import { createToolkit, type DisclosureModule, type Toolkit, toolFor } from '../../tools/index.js';
import { AgentSession, type AgentSessionOptions, type ModelSource } from '../index.js';

/** A toy host: a counter with a version, a log, and what a conversation holds. */
export interface Toy {
  state: { version: number; value: number; bump(by?: number): void };
  log: { lines: string[] };
  /** How the flaky tools behave. */
  flaky: { calls: number; failures: number; error: 'none' | 'applied' | 'throw' };
  disclosure?: DisclosureModule;
}

export function toy(): Toy {
  const state = {
    version: 0,
    value: 0,
    bump(by = 1) {
      state.value += by;
      state.version += 1;
    },
  };
  return { state, log: { lines: [] }, flaky: { calls: 0, failures: 0, error: 'none' } };
}

export const tool = toolFor<Toy>();
const none = { kind: 'object', fields: {} } as const;

export const guidance = 'You operate a counter.';

/** The toy host's tools, as a toolkit: a core, a skill with a tool, a single tool, a tool that needs approval. */
export function toyToolkit(extra: { tools?: Parameters<typeof createToolkit<Toy>>[0]['tools'] } = {}): Toolkit<Toy> {
  const get = tool({
    id: 'get_state',
    description: 'Reads the value.',
    input: { kind: 'object', fields: { look: { kind: 'integer', optional: true } } },
    modules: ['state'],
    read: true,
    run: ({ api }) => ({ ok: true, did: `Value ${api.state.value}.` }),
  });
  const set = tool({
    id: 'set_value',
    description: 'Bumps the value.',
    input: { kind: 'object', fields: { by: { kind: 'integer', optional: true } } },
    modules: ['state'],
    run: ({ args, api }) => {
      api.state.bump(args.by ?? 1);
      return { ok: true, did: 'Bumped.', changed: ['value'] };
    },
  });
  // Fails as the host's `flaky` says: `failures` times, then works.
  const behave = (api: Pick<Toy, 'flaky' | 'log'>, name: string) => {
    api.flaky.calls += 1;
    api.log.lines.push(`${name} ${api.flaky.calls}`);
    if (api.flaky.calls > api.flaky.failures) return undefined;
    if (api.flaky.error === 'throw')
      throw Object.assign(new Error('The service is unavailable (503).'), { transient: true });
    return {
      ok: false as const,
      code: 'service.unavailable',
      message: 'The service is unavailable (503).',
      transient: true,
      ...(api.flaky.error === 'applied' ? { effect: 'applied' as const } : {}),
    };
  };
  const fetchRead = tool({
    id: 'fetch_read',
    description: 'Reads from a service.',
    input: none,
    modules: ['flaky', 'log'],
    read: true,
    run: ({ api }) => behave(api, 'fetch_read') ?? { ok: true, did: 'Fetched.' },
  });
  const fetchWrite = tool({
    id: 'fetch_write',
    description: 'Writes to a service.',
    input: none,
    modules: ['flaky', 'log'],
    run: ({ api }) => behave(api, 'fetch_write') ?? { ok: true, did: 'Written.' },
  });
  const fetchSafe = tool({
    id: 'fetch_safe',
    description: 'Writes to a service, safe to repeat.',
    input: none,
    modules: ['flaky', 'log'],
    idempotent: true,
    run: ({ api }) => behave(api, 'fetch_safe') ?? { ok: true, did: 'Written.' },
  });
  const archive = tool({
    id: 'archive',
    description: 'Archives the value. Loaded by name.',
    input: none,
    modules: ['state'],
    run: () => ({ ok: true, did: 'Archived.' }),
  });
  const loadMap = tool({
    id: 'load_map',
    description: 'Loads a map.',
    input: none,
    modules: ['log'],
    run: ({ api }) => {
      api.log.lines.push('load_map');
      return { ok: true, did: 'Loaded.' };
    },
  });
  const ask = tool({
    id: 'ask_service',
    description: 'Calls a service, with the person’s approval.',
    input: none,
    modules: ['log'],
    approval: true,
    run: ({ api }) => {
      api.log.lines.push('ask_service');
      return { ok: true, did: 'Asked.' };
    },
  });
  return createToolkit<Toy>({
    tools: [get, set, fetchRead, fetchWrite, fetchSafe, ask, archive, loadMap, ...(extra.tools ?? [])],
    core: ['get_state', 'set_value', 'fetch_read', 'fetch_write', 'fetch_safe', 'ask_service'],
    groups: [{ id: 'maps', title: 'Maps', summary: 'maps', guidance: 'Skill maps: load a map.', tools: ['load_map'] }],
  });
}

/** Real time, for the few tests that wait. */
export const realClock: Clock = {
  now: () => Date.now(),
  sleep: (ms, signal) =>
    new Promise<void>((resolve, reject) => {
      if (signal?.aborted) return reject(signal.reason);
      const timer = setTimeout(resolve, ms);
      signal?.addEventListener(
        'abort',
        () => {
          clearTimeout(timer);
          reject(signal.reason);
        },
        { once: true },
      );
    }),
};

export const DEFAULT_MODEL = { provider: 'scripted', model: 'scripted-1' } as const;

/** Models for a session: the providers by id, the scripted one under `scripted`. */
export function modelsOf(scripted: Provider, others: Record<string, Provider> = {}): ModelSource {
  return {
    providerFor: (provider) => (provider === 'scripted' ? scripted : others[provider]),
    default: { ...DEFAULT_MODEL },
  };
}

export interface Fixture {
  session: AgentSession<Toy>;
  model: ScriptedModel;
  api: Toy;
  toolkit: Toolkit<Toy>;
}

/** A session over the toy host and a scripted model (the script is sequential: one step per model call). */
export function fixture(
  script: ModelStep[] = [],
  options: Partial<AgentSessionOptions<Toy>> = {},
  model: ScriptedModelOptions = {},
  api: Toy = toy(),
): Fixture {
  const scripted = scriptedModel(script, model);
  const toolkit = options.toolkit ?? toyToolkit();
  const session = new AgentSession<Toy>({
    models: modelsOf(scripted.provider),
    toolkit,
    api,
    guidance,
    toolRetryDelayMs: 1,
    retryDelayMs: 1,
    ...options,
  });
  return { session, model: scripted, api, toolkit };
}

export const textOf = (m: Message): string => m.parts.map((p) => (p.kind === 'text' ? p.text : '')).join('');

/** The messages of the request the model saw (the default payload holds them). */
export const messagesOf = (model: ScriptedModel, n: number): Message[] =>
  (model.invocations[n]?.payload as unknown as { messages: Message[] }).messages;

/** The tool results the conversation holds, as text, in order. */
export function resultsIn(messages: readonly Message[]): string[] {
  const out: string[] = [];
  for (const m of messages)
    if (m.role === 'user')
      for (const p of m.parts)
        if (p.kind === 'tool-result') out.push(p.parts.map((x) => (x.kind === 'text' ? x.text : '')).join(''));
  return out;
}

/** The committed messages of a session's conversation. */
export async function headOf(session: AgentSession<Toy>): Promise<Message[]> {
  const client = await session.client();
  const [head] = await client.heads([(await client.execution()).main]);
  return head === undefined ? [] : await client.entries(head.entries);
}

export const callOf = (tool: string, args: Record<string, unknown> = {}) => ({
  tool,
  arguments: args as never,
});
