/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// How far ahead of the scroll a virtual list keeps rows made (rendering.md §2.7, the margin): the browser scrolls on the
// compositor before the scroll event reaches the list, so a fast wheel or flick uncovers what was not made yet (the
// skeleton shows there). Ahead of the scroll, at least `atLeast` rows, and where the scroll will be two frames on at its
// speed now: a wheel's large steps (1,200 px a frame) are covered from the second step. Capped at three screens. Only the
// wheel's and touch's speed counts: a thumb's drag and `scrollToIndex` set the position from script, and the list
// follows in the same frame.

/** A scroll this long after a wheel or touch event is the compositor's (its speed counts). */
export const NativeMs = 200;

/** Screens ahead at most. */
export const MaxScreensAhead = 3;

/** Rows to keep made ahead of the scroll, moving `speed` px since the last frame. */
export function rowsAhead(speed: number, rowHeight: number, screen: number, atLeast: number): number {
  const bySpeed = Math.ceil((2 * speed) / rowHeight);
  return Math.max(atLeast, Math.min(bySpeed, MaxScreensAhead * Math.max(1, screen)));
}

/** How long an opened run's rows take to come in (`--msp-motion-base`), and a frame's margin for its end. */
export const RevealMs = 240;

/**
 * The run of rows `next` inserted into `prev` and nothing else (a node opened, its children loaded): `[at, at + count)`
 * in `next`, or undefined for any other change (rows removed, moved, or inserted in more than one place). A list's
 * rows animate in only for that change (a scroll never animates). One pass over the keys, per change of the rows.
 */
export function insertedRun<T>(
  prev: readonly T[],
  next: readonly T[],
  key: (item: T) => unknown,
): { at: number; count: number } | undefined {
  const count = next.length - prev.length;
  if (count <= 0 || prev.length === 0) return undefined;
  let at = 0;
  while (at < prev.length && key(prev[at]) === key(next[at])) at++;
  for (let i = at; i < prev.length; i++) if (key(prev[i]) !== key(next[i + count])) return undefined;
  return { at, count };
}
