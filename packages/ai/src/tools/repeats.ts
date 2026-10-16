/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, ToolDefinition } from '../messages/index.js';
import { DISCLOSE_SKILLS, DISCLOSE_TOOLS } from './names.js';
import type { ToolOutcome } from './result.js';

export interface RepeatedReadsOptions {
  /** The host's state version: a read is repeated only while it is the same value (compared with `===`). */
  version(): unknown;
  /** Reads in a row before the note to act (default 6). */
  streak?: number;
  /** Tools that count as neither a read nor an action (by default the disclosure tools: setting up is not acting). */
  exempt?: readonly string[];
}

export const DEFAULT_STREAK = 6;

/**
 * Guards against reading instead of acting (Gemini, thinking more, gathers far more than it needs: neighboring scenes,
 * the state again in full, one method at a time; tools.md §4):
 * - a read called again in the same turn with the same arguments, while the host's state version is unchanged, answers
 *   in one line instead of repeating its result (a changed state reads fresh);
 * - after `streak` reads in a row with no other call, each read's result says to act now.
 * Cleared at the start of each turn. Setting up a skill does not count as either.
 */
export class RepeatedReads {
  private readonly seen = new Map<string, unknown>();
  private streak = 0;

  constructor(private readonly options: RepeatedReadsOptions) {}

  reset(): void {
    this.seen.clear();
    this.streak = 0;
  }

  /** A call that is not a read: the streak starts again (not for the exempt ones). */
  acted(definition: ToolDefinition): void {
    if (!(this.options.exempt ?? [DISCLOSE_SKILLS, DISCLOSE_TOOLS]).includes(definition.tool.id)) this.streak = 0;
  }

  /** Counts a read; past the streak, the note to add to its result. */
  nudge(): string | undefined {
    this.streak += 1;
    if (this.streak < (this.options.streak ?? DEFAULT_STREAK)) return undefined;
    return `(${this.streak} reads in a row without a change: you have what you need. Act now with what you know; read more only for one specific missing fact.)`;
  }

  check(definition: ToolDefinition, args: Record<string, Json>): { outcome?: ToolOutcome; remember(): void } {
    let version: unknown;
    try {
      version = this.options.version();
    } catch {
      return { remember() {} };
    }
    const key = `${definition.tool.id} ${stableJson(args)}`;
    const remember = () => void this.seen.set(key, version);
    if (!this.seen.has(key) || this.seen.get(key) !== version) return { remember };
    return {
      outcome: {
        ok: true,
        text: `Unchanged: nothing changed since your earlier ${definition.name} call with these arguments, so its result above is current. Do not read it again; act on it.`,
      },
      remember,
    };
  }
}

const isRecord = (v: Json): v is Record<string, Json> => typeof v === 'object' && v !== null && !Array.isArray(v);

function stableJson(v: Json): string {
  if (Array.isArray(v)) return `[${v.map(stableJson).join(',')}]`;
  if (isRecord(v))
    return `{${Object.keys(v)
      .sort()
      .map((k) => `${JSON.stringify(k)}:${stableJson(v[k] as Json)}`)
      .join(',')}}`;
  return JSON.stringify(v);
}
