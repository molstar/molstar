/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, ToolDefinition } from '@molstar/ai-agent';
import type { RunContext, ToolFrame } from '@molstar/ai-harness';
import { type DisclosureModule, type HandlerOptions, type ToolRegistration, toolFor } from '../index.js';

/** A toy host: a counter with a version, a log, and what a conversation holds. */
export interface Toy {
  state: { version: number; value: number; bump(by?: number): void };
  log: { lines: string[] };
  disclosure: DisclosureModule;
  /** A module a group may need, absent in some hosts. */
  extra?: { name: string };
}

export function toy(overrides: Partial<Toy> = {}): Toy & { loaded: string[]; keys: string[] } {
  const state = {
    version: 0,
    value: 0,
    bump(by = 1) {
      state.value += by;
      state.version += 1;
    },
  };
  const held = { loaded: [] as string[], keys: [] as string[] };
  return Object.assign(
    {
      state,
      log: { lines: [] as string[] },
      disclosure: { loaded: async () => [...held.loaded], skills: async () => [...held.keys] },
    },
    held,
    overrides,
  );
}

export const tool = toolFor<Toy>();

/** A call's context for calling a handler directly: no harness (journal empty, no nested runs unless given). */
export function runContext(definition: ToolDefinition, args: Json, overrides: Partial<RunContext> = {}): RunContext {
  return {
    harnessId: 'test',
    turnId: 'test-turn',
    callId: 'call-1',
    path: [],
    key: 'key-1',
    prior: 'none',
    attempt: 1,
    definition,
    arguments: args,
    signal: new AbortController().signal,
    journal: { read: async () => ({ messages: [] }), append() {}, edit() {} },
    run: async () => {
      throw new Error('no nested runs when a tool is called directly');
    },
    ...overrides,
  };
}

export type ResultFrame = Extract<ToolFrame, { kind: 'result' }>;

export interface Called {
  frames: ToolFrame[];
  progress: Json[];
  result: ResultFrame;
  status: 'succeeded' | 'failed';
  text: string;
}

/** Runs a registration's handler once and gathers its frames. */
export async function call<Api>(
  registration: ToolRegistration<Api>,
  api: Api,
  args: Json,
  options: HandlerOptions<Api> = {},
  context: Partial<RunContext> = {},
): Promise<Called> {
  const frames: ToolFrame[] = [];
  for await (const frame of registration
    .handler(api, options)
    .invoke(runContext(registration.definition, args, context)))
    frames.push(frame);
  const result = frames.find((f): f is ResultFrame => f.kind === 'result');
  if (!result) throw new Error('no result');
  const text = result.parts.map((p) => (p.kind === 'text' ? p.text : '')).join('');
  return {
    frames,
    progress: frames.flatMap((f) => (f.kind === 'progress' ? [f.value] : [])),
    result,
    status: result.status,
    text,
  };
}

/** The ids of the tools a contribution's tool-definition parts declare. */
export const definedIn = (contributions: ResultFrame['contributions']): string[] =>
  (contributions ?? []).flatMap((m) =>
    m.role === 'system' ? m.parts.flatMap((p) => (p.kind === 'tool-definition' ? [p.definition.tool.id] : [])) : [],
  );

/** The text parts of a contribution's user messages. */
export const textsIn = (contributions: ResultFrame['contributions']): string[] =>
  (contributions ?? []).flatMap((m) =>
    m.role === 'user' ? m.parts.flatMap((p) => (p.kind === 'text' ? [p.text] : [])) : [],
  );
