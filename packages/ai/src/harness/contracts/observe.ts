/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '../../messages/index.js';
import type { DecisionKind } from './drive.js';
import type { AttemptId, AttemptRecord, CallId, HarnessEvent, Seq, StepId } from './records.js';

/** The harness's own minimal observable; apps wrap it in rxjs if they want to. */
export interface HarnessObservable<T> {
  /** Returns the teardown. Whether the current value comes first is each source's to state. */
  subscribe(onNext: (value: T) => void, onError?: (error: unknown) => void, onComplete?: () => void): () => void;
}

/** Ids only, never payloads: observations are meant to be serialized (harness.md §13). */
export interface DecisionIds {
  readonly calls?: readonly CallId[];
  readonly step?: StepId;
  readonly seq?: Seq;
}

interface AttemptScope {
  readonly path: readonly CallId[];
  readonly callId: CallId;
  readonly attemptId: AttemptId;
}

export type Observation =
  | (AttemptScope & {
      readonly kind: 'model.delta';
      readonly delta: { readonly text: string } | { readonly reasoning: string };
    })
  | (AttemptScope & { readonly kind: 'tool.progress'; readonly value: Json })
  | (AttemptScope & { readonly kind: 'attempt.started' })
  | (AttemptScope & { readonly kind: 'attempt.ended'; readonly outcome: AttemptRecord['outcome'] })
  | {
      readonly kind: 'decision.opened';
      readonly decisionKind: DecisionKind;
      readonly path: readonly CallId[];
      readonly ids: DecisionIds;
    }
  | {
      readonly kind: 'decision.answered';
      readonly decisionKind: DecisionKind;
      readonly path: readonly CallId[];
      readonly answer: string;
      readonly ids: DecisionIds;
    }
  | { readonly kind: 'event'; readonly event: HarnessEvent };
