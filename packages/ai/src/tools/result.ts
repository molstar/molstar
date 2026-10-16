/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Message } from '../messages/index.js';
import type { ToolFrame } from '../harness/index.js';

/**
 * What a call did to the world (harness.md §8): `none` (a read, or a failure before any change), `applied`, or
 * `unknown` (it threw and it is not a read). A tool may say it; otherwise it is derived (tools.md §3).
 */
export type Effect = NonNullable<Extract<ToolFrame, { kind: 'result' }>['effect']>;

/** What a tool returns: what it did, what changed, and what else the model should know (tools.md §3). */
export type ToolResult =
  | {
      ok: true;
      /** What it did, one sentence ("Added scene pocket after §2."); "…: already so." when held. */
      did: string;
      /** What changed, by name or address ("§4 pocket", "chain A"): rendered `Changed: a, b.`. */
      changed?: readonly string[];
      /** What the person now sees, for tools that move the view; rendered `View: …`. */
      view?: string;
      warnings?: readonly string[];
      /** A read's content (a scene's detail, a table's rows, a query's elements). */
      body?: string;
      /** The obvious next step, when there is one. */
      next?: string;
      contributions?: readonly Message[];
      effect?: Effect;
    }
  | { ok: false; code: string; message: string; hint?: string; transient?: boolean; effect?: Effect };

/** A tool's result before it becomes model-visible content. */
export type ToolOutcome =
  | { ok: true; text: string; contributions?: readonly Message[]; effect?: Effect }
  | { ok: false; code: string; message: string; hint?: string; transient?: boolean; effect?: Effect };

const sentence = (s: string) => (/[.!?]$/.test(s) ? s : `${s}.`);

/** The result as the model reads it: always in this order. */
export function render(r: ToolResult): ToolOutcome {
  if (!r.ok) return r;
  const lines = [r.did];
  if (r.changed?.length) lines.push(`Changed: ${r.changed.join(', ')}.`);
  if (r.view) lines.push(`View: ${sentence(r.view)}`);
  // Labelled with what to do (guidance rule 11): small models act on errors and pass over notes.
  if (r.warnings?.length)
    lines.push(`Warnings (fix if it affects the request, else tell the user): ${r.warnings.join(' ')}`);
  if (r.body) lines.push(r.body);
  if (r.next) lines.push(`Next: ${r.next}`);
  return {
    ok: true,
    text: lines.join('\n'),
    ...(r.contributions ? { contributions: r.contributions } : {}),
    ...(r.effect ? { effect: r.effect } : {}),
  };
}

/** A failure as a result. */
export const failed = (code: string, message: string, hint?: string): ToolResult => ({
  ok: false,
  code,
  message,
  ...(hint ? { hint } : {}),
});

/** A tool's outcome as the harness's result frame. */
export function toFrame(o: ToolOutcome, effect: Effect): Extract<ToolFrame, { kind: 'result' }> {
  if (o.ok)
    return {
      kind: 'result',
      status: 'succeeded',
      parts: [{ kind: 'text', text: o.text }],
      effect,
      ...(o.contributions?.length ? { contributions: o.contributions } : {}),
    };
  return {
    kind: 'result',
    status: 'failed',
    parts: [{ kind: 'text', text: o.hint ? `${o.message} ${o.hint}` : o.message }],
    effect,
    ...(o.transient ? { transient: true } : {}),
  };
}
