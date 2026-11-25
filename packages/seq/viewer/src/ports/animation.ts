/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Seq*'s animation port (_spec/README.md "Ports", Animation): the history port's twin. Seq* owns no clock; the
// application does. Two faces:
//
// - the participant (the application drives Seq*): `progress(plan, t)` or `show(plan, frame)` when it decides, `clear()`
//   back to the state; `drawn` moves per frame drawn, `durationOf(plan)` its natural length;
// - the animation (Seq* asks the application): Seq*'s own motions (a jump, a step's view) are plans handed over with
//   `request(plan, { duration })`, so a host that animates a structure's camera runs both together, or cuts both.
//
// Adapters live in the composition roots (`@molstar/seqx-viewer/animation`): the page's frames by default, a fake clock
// headless.

import type { ReadonlyState } from '@molstar/core/reactive';
import type { Frame, TransitionPlan } from '@molstar/seqx-core/story';

/** The participant face: what the viewer gives the application's clock. */
export interface SeqAnimationParticipant {
  /** Shows the plan at `t` in [0, 1] (sampled into the plan's frame) over the shown state: no revision, no history. */
  progress(plan: TransitionPlan, t: number): void;
  /** Shows a frame the application sampled itself. */
  show(plan: TransitionPlan, frame: Frame): void;
  /** Leaves the frame: the state is shown again. */
  clear(): void;
  /** Draws what is shown now (a headless clock's tick; a mounted viewer draws on its own frames). */
  draw(): void;
  /** Moves once per frame drawn while a frame is shown. */
  readonly drawn: ReadonlyState<number>;
  /** A plan's natural duration, ms (from the steps' timing). */
  durationOf(plan: TransitionPlan): number;
}

/** A plan handed to the application: what it may stop. */
export interface AnimationRequest {
  /** Stops it where it is (the participant is not cleared: the caller settles what is shown). */
  cancel(): void;
}

export interface AnimationRequestOptions {
  readonly duration: number;
  /** Told once the plan ran to its end (not when cancelled). */
  readonly done?: () => void;
}

/** The animation face: the application's clock as Seq* asks it. */
export interface SeqAnimation {
  /** Runs a plan from 0 to 1 over `duration` ms on the application's clock, through the participant's `progress`. */
  request(plan: TransitionPlan, options: AnimationRequestOptions): AnimationRequest;
  /** False: the host owns it (one clock for several views), and the viewer leaves it. */
  readonly disposeWithViewer?: boolean;
  dispose?(): void;
}

/** The port in a viewer's spec: a factory given the participant, called once per viewer. */
export type SeqAnimationPort = (participant: SeqAnimationParticipant) => SeqAnimation;
