/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  canonicalJson,
  type Json,
  type Message,
  type ProviderBindings,
  type ProjectInput,
} from '../../messages/index.js';
import { createSession as session, type Session, type SessionSetup } from '../../testing/index.js';
import { declareTools, type ScriptedModel, settle, system, user } from '../../testing/index.js';
import { defaultDrive, drive, type DriveHandlers, type TurnConfiguration } from '../../harness/index.js';

/** A short stable digest of any value (the trace compares them, nothing more). */
function short(value: Json): string {
  let h = 5381;
  for (const c of canonicalJson(value)) h = (h * 33 + c.charCodeAt(0)) >>> 0;
  return h.toString(16);
}

/** What a test varies between turns, and a provider would report as the request's parts. */
export const adapterState = { systemVersion: 0 };

/**
 * A request's parts as a provider caches them: the tools, the system text, then each message (what a real provider
 * reports as `prefix`; the scripted model reports none).
 */
export function prefixOf(input: ProjectInput): string[] {
  const systemText = input.messages.flatMap((m) =>
    m.role === 'system' ? m.parts.map((p) => (p.kind === 'text' ? p.text : '')) : [],
  );
  return [
    short(input.tools.map((t) => [t.tool.id, t.inputSchema]) as unknown as Json),
    short([systemText, adapterState.systemVersion]),
    ...input.messages.filter((m) => m.role !== 'system').map((m) => short(m as unknown as Json)),
  ];
}

/** The session's scripted model, reporting the request's `prefix`. */
export function prefixing(model: ScriptedModel): ProviderBindings {
  const bindings: ProviderBindings = {
    resolve: async () => ({
      ...model.provider,
      project: async (input: ProjectInput) => ({
        ...(await model.provider.project(input)),
        prefix: prefixOf(input),
      }),
    }),
  };
  return bindings;
}

export const with_ = (handlers: DriveHandlers) => drive(handlers, defaultDrive());

/** Whole turns over the harness: the scripted model (reporting a prefix) and tools, a simulated clock. */
export const createSession = (setup: SessionSetup = {}): Session => session({ models: prefixing, ...setup });

/**
 * A request as a session writes it: the text, labelled with the person's words and focus (the trace reads the label,
 * never the text). A text without `Request: ` stays a plain message.
 */
function requestMessage(text: string): Message {
  const at = text.startsWith('Request: ') ? 0 : text.indexOf('\n\nRequest: ') + 2;
  if (at <= 1 && !text.startsWith('Request: ')) return user(text);
  const focus = text
    .slice(0, at)
    .split('\n')
    .find((line) => line.startsWith('Focus: '))
    ?.slice('Focus: '.length);
  return {
    ...user(text),
    label: {
      kind: 'session:request',
      data: { prompt: text.slice(at + 'Request: '.length), ...(focus !== undefined ? { focus } : {}) },
    },
  };
}

/** The first prompt: guidance, the tools (all loaded unless named as deferred), then the request. */
export function first(
  s: Session,
  request: string,
  options: { deferred?: readonly string[]; guidance?: string } = {},
): Message[] {
  const deferred = new Set(options.deferred ?? []);
  const defs = s.tools.definitions;
  return [
    system(options.guidance ?? 'Be brief.'),
    declareTools(defs.filter((d) => !deferred.has(d.tool.id))),
    ...(deferred.size
      ? [
          declareTools(
            defs.filter((d) => deferred.has(d.tool.id)),
            true,
          ),
        ]
      : []),
    requestMessage(request),
  ];
}

/** Sends a prompt and runs its turn, with a configuration of its own if given. */
export const turn = (s: Session, prompt: string | Message[], options: { configuration?: TurnConfiguration } = {}) =>
  s.send(typeof prompt === 'string' ? [requestMessage(prompt)] : prompt, options);

/** A request as the session composes it. */
export const request = (prompt: string, parts: { state?: string; focus?: string; notes?: string[] } = {}): string =>
  [
    parts.state,
    parts.focus !== undefined ? `Focus: ${parts.focus}` : undefined,
    parts.notes?.length ? `Notes since your last turn:\n${parts.notes.map((n) => `- ${n}`).join('\n')}` : undefined,
    `Request: ${prompt}`,
  ]
    .filter((x) => x !== undefined)
    .join('\n\n');

/** Lets real work proceed until `condition` holds. */
export async function until(condition: () => boolean, limit = 500): Promise<void> {
  for (let i = 0; i < limit && !condition(); i += 1) await settle(2);
  if (!condition()) throw new Error('until: the condition never held');
}

export const configured = (
  s: Session,
  patch: Partial<TurnConfiguration['model']>,
): { configuration: TurnConfiguration } => {
  const ids = s.tools.definitions.map((d) => d.tool.id);
  return {
    configuration: {
      model: { target: { provider: 'scripted', model: 'scripted-1' }, parameters: {}, ...patch },
      tools: { enabled: ids, preauthorized: ids },
      limits: {},
    },
  };
};

export type { Session, SessionSetup };
