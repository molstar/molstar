/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Motion is the application's (dom §5): `Show`'s `reveal` marks opening content with `data-reveal` (its style animates
// it) and asks `Motion.close` to play a closing back; `swap` asks `Motion.swap` to run a replacement as a transition.
// Unset, content closes at once and a swap is a plain update.

export type RevealKind = 'grow' | 'side' | 'fade';

export const Motion: {
  /** Plays `el` closing; undefined when nothing plays (reduced motion, no support): it goes at once. */
  close?(el: Element, kind: RevealKind): Animation | undefined;
  /** Runs `update` as a transition scoped to `scope`, where the platform has them; else runs it now. */
  swap?(scope: Element | null | undefined, update: () => void): void;
} = {};

export function playClose(el: Element, kind: RevealKind): Animation | undefined {
  return Motion.close?.(el, kind);
}

export function swapIn(scope: Element | null | undefined, update: () => void): void {
  if (Motion.swap) Motion.swap(scope, update);
  else update();
}
