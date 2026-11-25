/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a layer over a source has read, as its track's sizing knows it (_spec/README.md "Render", Track heights): the
// lanes a ranges layer packed, by the layer's prep (another file or region is another prep and starts again), noted by
// the layout once the view settles. A count moves each time one grew: the sizing's key, a section's reuse and the
// layout's second pass read it.

const Lanes = new WeakMap<object, number>();
const Empty = new WeakSet<object>();
let Count = 0;

/** Notes that a source's layer packed `lanes` lanes: the track grows to them, never shrinks by a pan. */
export function noteLanes(prep: object, lanes: number): void {
  if (lanes > (Lanes.get(prep) ?? 1)) {
    Lanes.set(prep, lanes);
    Empty.delete(prep);
    Count++;
  }
}

/**
 * Notes that a source's layer, settled over its reachable region, holds no item: its track is empty, not pending. Not
 * when lanes are known (a block's layer shares its key with its other blocks': one empty block empties no track).
 */
export function noteEmpty(prep: object): void {
  if (!Empty.has(prep) && !Lanes.has(prep)) {
    Empty.add(prep);
    Count++;
  }
}

/** Whether a source's layer settled with no item (`noteEmpty`). */
export function knownEmpty(prep: object): boolean {
  return Empty.has(prep);
}

/** The lanes a layer has packed (one while it is pending: the loading lane). */
export function knownLanes(prep: object): number {
  return Lanes.get(prep) ?? 1;
}

/** A count that moves each time some layer's known lanes grew (a ranges layer type's `known`). */
export function knownLanesVersion(): number {
  return Count;
}

/** Forgets what a prep's source had packed (its data was replaced): the track sizes again from what is read. */
export function resetLanes(prep: object): void {
  if (Lanes.delete(prep) || Empty.delete(prep)) Count++;
}
