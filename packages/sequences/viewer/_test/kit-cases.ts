/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The conformance kit's `heap` hook for the guards that run renderers through the kit (conformance-cases.ts,
// draws-cases.ts, interactions/fingerprint-cases.ts) and their quick twins: the kit hands each frame it draws (a layer
// at a zoom, `tiny` too) to the hook, here told the geometry it opened; a quick twin measures each as a quick case
// (`KIT` counts), its id from `frameIds`. No viewer module is imported here, so a cases module that loads the viewer
// afresh (`vi.resetModules`) can import it.

import type { QuickOptions } from './heap.js';

/** The kit's `heap` hook, told the geometry of the frame (the last one the kit opened). */
export type KitHeap = (frame: () => void, geometry: string) => void;

/**
 * A kit frame's quick counts: a frame is 0.05 to 0.3 ms (M4) and a kit hands over up to 30 a run, so windows of 250
 * steps (exact accounting reads a byte a step in 250 bytes); tier-up is done by 1,000 steps under the quick flags.
 */
export const KIT: QuickOptions = { warmSteps: 1_000, steps: 250, windows: 3 };

/**
 * The quick cases' ids for a kit run's frames: `<prefix>, <geometry> <k>`, `k` the frame's order in its geometry (the
 * kit's order is fixed). Each call of the returned function names the next frame.
 */
export function frameIds(prefix: string): (geometry: string) => string {
  const seen: Record<string, number> = {};
  return (geometry) => `${prefix}, ${geometry} ${(seen[geometry] = (seen[geometry] ?? -1) + 1)}`;
}
