/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The chevron of a group's or split track's header: a thin stroked chevron after the label's text, which the composer
// draws and the headers behavior hits, from these numbers (CSS px). The label starts where a track's does (10 px, 14
// per group level); the text is cut `CHEVRON_ROOM` short, so the chevron always has its place.

import { AxisGeometry } from "./axis.ts";
import type { ViewModel } from "./graph.ts";
import type { SectionLayout, SplitFeatures } from "./layout.ts";

/** The gap between the label's text and the chevron's box. */
export const CHEVRON_GAP = 4;
/** Half the chevron's box width (the chevron is 7 px wide). */
export const CHEVRON_HALF = 4;
/** What the label gives up to the chevron: the gap and its box. */
export const CHEVRON_ROOM = CHEVRON_GAP + 2 * CHEVRON_HALF + 2;
/** The hit box reaches this far beyond the chevron's box on each side. */
export const CHEVRON_SLOP = 3;

/**
 * Whether entry `e` is a split track's label entry whose chevron is drawn: a scrolled entry of a linear section whose
 * split has values. A group's header always has one.
 */
export function splitChevron(m: ViewModel, l: SectionLayout, e: number): boolean {
  const t = l.entryTrack[e];
  const split = t >= 0 ? m.resolved.tracks.split[t] : -1;
  if (split < 0 || e < l.pinnedCount || l.geometry !== AxisGeometry.linear) return false;
  return ((m.prepared[split].state as SplitFeatures).split?.labels.length ?? 0) > 0;
}
