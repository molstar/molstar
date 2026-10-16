/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '../messages/index.js';
import type { PlannedCall } from '../harness/index.js';
import type { CompactionReport } from './compaction/index.js';
import type { ToolEvent, TurnHooks } from './types.js';

/**
 * Everything the session keeps for the turn that is running, in one object made afresh for each turn (and each
 * resume, which keeps the events shown before the park and seeds `planned` from the parked step): the drive and the
 * observers write it, the result reads it.
 */
export interface TurnState {
  hooks: TurnHooks;
  /** The turn id the host has been told of. */
  announced?: string;
  /** The tool events shown to the person, in order. */
  events: ToolEvent[];
  /** When each call started (its first attempt), by call id. */
  starts: Map<string, number>;
  /** A transiently failed attempt per call, reported only if its retries run out. */
  retrying: Map<string, ToolEvent>;
  /** The streamed text per model attempt. */
  previews: Map<string, string>;
  lastText: string;
  /** The text of each top-level model call, by call id, in order: its latest attempt's. */
  replies: Map<string, string>;
  /** The latest top-level model call (the one a `model-result` replace answers). */
  lastModelCall: string | undefined;
  /** A reply cut short was asked for again with more room (once per turn). */
  widened: boolean;
  /** What compactions kept during the turn (reported in its result). */
  compactions: CompactionReport[];
  /** The size the last compaction in this turn left (a step's own is not in the event log yet). */
  compactedSize: number | undefined;
  /** The calls of the pending tools step, so the ones that never ran show as skipped. */
  planned: { callId: string; toolId: string; arguments: Json }[];
  /** What the person sees for each call skipped ("Not run: add_scene failed before it."), if the step was cut. */
  skipShown: string | undefined;
}

export function newTurnState(hooks: TurnHooks): TurnState {
  return {
    hooks,
    events: [],
    starts: new Map(),
    retrying: new Map(),
    previews: new Map(),
    lastText: '',
    replies: new Map(),
    lastModelCall: undefined,
    widened: false,
    compactions: [],
    compactedSize: undefined,
    planned: [],
    skipShown: undefined,
  };
}

export type Planned = Pick<PlannedCall, 'callId' | 'toolId' | 'arguments'>;
