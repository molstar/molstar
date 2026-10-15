/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { EntryLabel, JsonObject, Message, Reference, ToolDefinition } from '@molstar/ai-agent';
import { fingerprint, LABEL_DEFINITIONS } from '@molstar/ai-tools';

/**
 * The messages of a turn (tools-and-session.md §4): the first turn sends the guidance, then the tool definitions
 * (the core loaded, the rest declared), then the request; later turns send what changed since (new or changed
 * definitions, changed guidance as a user message), then the request.
 */

/** Opens the person's words in a request: compaction and the trace find them by it. */
export const REQUEST_MARKER = 'Request: ';
/** Opens the Focus line of a request. */
export const FOCUS_MARKER = 'Focus: ';
/** Heads guidance sent again because it changed since the conversation began. */
export const GUIDANCE_UPDATE = 'Updated guidance (it replaces the earlier version):';
/** Heads the host's notes in a request (the default; a host sets its own). */
export const NOTES_HEADING = 'Notes since your last turn:';

/**
 * The labels the session gives the entries it writes (an `EntryLabel`'s `kind`): what compaction, `catchUp` and a
 * trace find them by. The markers above stay in the text, for the model.
 */
/** A turn's request; its data is `{ prompt, focus? }`. */
export const LABEL_REQUEST = 'session:request';
/** The guidance; its key is the text's fingerprint. */
export const LABEL_GUIDANCE = 'session:guidance';
/** Guidance sent again because it changed (`GUIDANCE_UPDATE` and the text); its key is the text's fingerprint. */
export const LABEL_GUIDANCE_UPDATE = 'session:guidance-update';
/** What compaction carried over from the earlier conversation. */
export const LABEL_CARRY_OVER = 'session:carry-over';
/** The host's memory, as compaction placed it. */
export const LABEL_MEMORY = 'session:memory';

/** What a conversation that began earlier has not seen: tools added or changed since, guidance that changed. */
export interface CatchUp {
  tools: readonly ToolDefinition[];
  deferred?: readonly ToolDefinition[];
  guidance: boolean;
}

export interface TurnInput {
  prompt: string;
  state?: string;
  focus?: string;
  notes?: readonly string[];
  /** Heads `notes` (default `NOTES_HEADING`). */
  notesHeading?: string;
  preamble?: string;
  /** Files sent with the request: references that carry their bytes. */
  attachments?: readonly Reference[];
  /** First turn of this conversation: add guidance and tool definitions. */
  first: boolean;
  tools: readonly ToolDefinition[];
  /**
   * First turn: tools declared but not loaded: a provider that loads tools in place lists them from the start, so
   * loading one later appends to the context instead of changing it.
   */
  deferred?: readonly ToolDefinition[];
  guidance: string;
  /** Later turns: what the conversation has not seen yet. Undefined: nothing to add. */
  catchUp?: CatchUp;
}

/**
 * The request as the model reads it: `[preamble][state]\n\nFocus: [focus]\n\n[notes]Request: [prompt]`, a part left out
 * (with its separator) when the request has none.
 */
export function requestOf(
  r: Pick<TurnInput, 'prompt' | 'state' | 'focus' | 'notes' | 'notesHeading' | 'preamble'>,
): string {
  const notes = r.notes?.length ? `${r.notesHeading ?? NOTES_HEADING}\n${r.notes.map((n) => `- ${n}`).join('\n')}` : '';
  return [
    r.preamble?.trimEnd(),
    r.state,
    r.focus !== undefined ? `${FOCUS_MARKER}${r.focus}` : undefined,
    notes,
    `${REQUEST_MARKER}${r.prompt}`,
  ]
    .filter((part) => part !== undefined && part !== '')
    .join('\n\n');
}

/** Loaded definitions, then declared ones, in one system message (never text). */
const definitions = (tools: readonly ToolDefinition[], deferred: readonly ToolDefinition[]): Message[] => {
  const parts = [
    ...tools.map((definition) => ({ kind: 'tool-definition' as const, definition, deferred: false })),
    ...deferred.map((definition) => ({ kind: 'tool-definition' as const, definition, deferred: true })),
  ];
  return parts.length ? [{ role: 'system', parts, label: { kind: LABEL_DEFINITIONS } }] : [];
};

/** The request's label: the person's words and focus as data, so no reader parses them out of the text. */
const requestLabel = (r: Pick<TurnInput, 'prompt' | 'focus'>): JsonObject => ({
  prompt: r.prompt,
  ...(r.focus !== undefined ? { focus: r.focus } : {}),
});

/** The guidance's label, by its text. */
export const guidanceLabel = (guidance: string, update = false): EntryLabel => ({
  kind: update ? LABEL_GUIDANCE_UPDATE : LABEL_GUIDANCE,
  key: fingerprint(guidance),
});

export function turnInput(options: TurnInput): Message[] {
  const attached = (options.attachments ?? []).map((reference) => ({ kind: 'reference' as const, reference }));
  const user: Message = {
    role: 'user',
    parts: [{ kind: 'text', text: requestOf(options) }, ...attached],
    label: { kind: LABEL_REQUEST, data: requestLabel(options) },
  };
  const { guidance } = options;
  // The first turn: the guidance, then the tools (loaded, then declared), each its own system message.
  if (options.first)
    return [
      { role: 'system', parts: [{ kind: 'text', text: guidance }], label: guidanceLabel(guidance) },
      ...definitions(options.tools, options.deferred ?? []),
      user,
    ];
  const c = options.catchUp;
  if (!c) return [user];
  // Later: definitions in a system message (never text); changed guidance as a user message, since a provider
  // without mid-conversation system text could not place it.
  return [
    ...definitions(c.tools, c.deferred ?? []),
    ...(c.guidance
      ? [
          {
            role: 'user' as const,
            parts: [{ kind: 'text' as const, text: `${GUIDANCE_UPDATE}\n\n${guidance}` }],
            label: guidanceLabel(guidance, true),
          },
        ]
      : []),
    user,
  ];
}
