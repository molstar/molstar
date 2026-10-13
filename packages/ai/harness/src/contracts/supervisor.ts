/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '@molstar/ai-agent';
import type { ItemView } from './drive.js';
import type { HarnessObservable } from './observe.js';
import type { HarnessId, LaunchKey, Progress, RunnerState, TurnId } from './records.js';

/** Liveness (harness.md §4): implementations are components, one per runtime. */
export interface Supervisor {
  /** Starts a runner for the key, idempotently; the report is lazy and unsubscribing never cancels. */
  ensure(harnessId: HarnessId, launch: LaunchKey): HarnessObservable<LaunchReport>;
}

export type LaunchReport =
  | { readonly phase: 'launching'; readonly key: LaunchKey; readonly detail?: Json }
  /** Then the report completes. */
  | { readonly phase: 'picked-up'; readonly runnerId: string; readonly turnId?: TurnId }
  | { readonly phase: 'failed'; readonly error: Json; readonly relaunch: boolean };

export interface SupervisorPolicy {
  /** Default: true. */
  shouldStart?(waiting: readonly ItemView[]): boolean;
  /** Undefined: never stalled. */
  stalledAfterMs?(phase: Progress['phase']): number | undefined;
  /** Once per stall; no default action. */
  onStalled?(state: Extract<RunnerState, { phase: 'live' }>): void;
}
