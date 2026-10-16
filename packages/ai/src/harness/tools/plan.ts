/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type AssistantMessage,
  canonicalJson,
  digest,
  type Json,
  type Message,
  type ToolCallPart,
  type ToolDefinition,
} from '../../messages/index.js';
import type { CallView } from '../contracts/drive.js';
import type { EntryId, PlannedCall } from '../contracts/records.js';
import { messageOf, type Run } from '../loop/run.js';

/** A planned call, with the definition it runs under when it is bound. */
export interface Planned extends PlannedCall {
  readonly bound?: ToolDefinition;
}

const resultIds = (message: Message | undefined) =>
  message?.role === 'user' ? message.parts.flatMap((p) => (p.kind === 'tool-result' ? [p.requestId] : [])) : [];
const callIds = (message: Message | undefined) =>
  message?.role === 'assistant' ? message.parts.flatMap((p) => (p.kind === 'tool-call' ? [p.requestId] : [])) : [];

/** Every call has its result in the next user message, and every result has its call (harness.md §2). */
export function projectable(messages: readonly Message[]): boolean {
  return messages.every((message, i) => {
    const answered = resultIds(messages[i + 1]);
    const asked = callIds(messages[i - 1]);
    return (
      callIds(message).every((id) => answered.includes(id)) && resultIds(message).every((id) => asked.includes(id))
    );
  });
}

export const toolCalls = (message: AssistantMessage): ToolCallPart[] =>
  message.parts.filter((part): part is ToolCallPart => part.kind === 'tool-call');

interface Loaded {
  readonly definition: ToolDefinition;
  readonly entry?: EntryId;
}

/** The tools a head declares and loads, in load order; a nested run's own definitions first (harness.md §8). */
export function scan(run: Run): { loaded: Map<string, Loaded>; declared: Set<string> } {
  const loaded = new Map<string, Loaded>();
  const declared = new Set<string>();
  for (const definition of run.definitions ?? []) loaded.set(definition.tool.id, { definition });
  for (const entry of run.context.entries) {
    const message = messageOf(run, entry);
    if (message?.role !== 'system') continue;
    for (const part of message.parts) {
      if (part.kind !== 'tool-definition') continue;
      declared.add(part.definition.tool.id);
      if (!part.deferred) loaded.set(part.definition.tool.id, { definition: part.definition, entry });
    }
  }
  return { loaded, declared };
}

/** The available tools: enabled and loaded (harness.md §8). */
export const available = (run: Run, scanned = scan(run)): ToolDefinition[] =>
  [...scanned.loaded.values()]
    .filter(({ definition }) => run.tools.enabled.includes(definition.tool.id))
    .map(({ definition }) => definition);

const isObject = (value: Json) => typeof value === 'object' && value !== null && !Array.isArray(value);

/** Plans each call of an assistant message (harness.md §8); keys per harness.md §8. */
export function planCalls(run: Run, message: AssistantMessage, step: string): Planned[] {
  const scanned = scan(run);
  const planned: Planned[] = [];
  for (const call of toolCalls(message)) {
    const toolId = call.tool.id;
    const callId = `${step}:${run.core.ids.next('call')}`;
    const loaded = scanned.loaded.get(toolId);
    const unbound =
      loaded === undefined
        ? scanned.declared.has(toolId)
          ? 'not-loaded'
          : 'unknown-tool'
        : !run.tools.enabled.includes(toolId)
          ? 'not-enabled'
          : run.core.tools.resolve(loaded.definition) === undefined
            ? 'no-handler'
            : !isObject(call.arguments)
              ? 'invalid-arguments'
              : undefined;
    let key = callId;
    if (run.depth > 0) {
      const same = canonicalJson({ toolId, arguments: call.arguments });
      const n = run.seen.get(same) ?? 0;
      run.seen.set(same, n + 1);
      key = digest({ parentKey: run.parentKey ?? '', toolId, arguments: call.arguments, n });
    }
    planned.push({
      callId,
      key,
      requestId: call.requestId,
      toolId,
      arguments: call.arguments,
      ...(loaded?.entry !== undefined ? { definition: loaded.entry } : {}),
      ...(unbound !== undefined ? { unbound } : {}),
      ...(unbound === undefined && loaded !== undefined ? { bound: loaded.definition } : {}),
    });
  }
  return planned;
}

/** A stored plan's definitions, read back from the entries it names (the tools step, and its re-run on resume). */
export function bind(run: Run, plan: readonly PlannedCall[]): Planned[] {
  return plan.map((call) => {
    if (call.unbound !== undefined || call.definition === undefined) return call;
    const message = messageOf(run, call.definition);
    const part = message?.parts.findLast(
      (p) => p.kind === 'tool-definition' && !p.deferred && p.definition.tool.id === call.toolId,
    );
    return part?.kind === 'tool-definition' ? { ...call, bound: part.definition } : { ...call, unbound: 'no-handler' };
  });
}

export const viewOf = (call: Planned): CallView => ({
  callId: call.callId,
  toolId: call.toolId,
  arguments: call.arguments,
  ...(call.bound !== undefined ? { definition: call.bound } : {}),
});

/** A plan as stored with the pending step: no definitions, only their entries. */
export const stored = (plan: readonly Planned[]): PlannedCall[] => plan.map(({ bound: _, ...call }) => call);
