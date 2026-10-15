/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Reveal (seqspec:interaction#focus, reveal and focus): a target brought into view in every section that shows it (or
// one), changing no selection, focus or zoom. Wrapped text scrolls its lines in (`Layouts.reveal`); a linear section
// pans at its cell (`Viewport.revealBlocks`), each only where the target is not shown already. Edge-rate: a target goes
// through the projector to extents, then to each section's hull.

import { Extents } from '@molstar/seq-core/model';
import type { LayoutManager } from './layout.js';
import type { SeqLociAny } from './lineage.js';
import { BlockHulls, type Coordinates } from './projection.js';
import type { Origin } from './state.js';
import type { MoveOpts, Viewports } from './viewport.js';

export interface RevealOptions {
  /** One section; null or left out: every section that shows the target. */
  section?: string | null;
  /** Default true: over the viewports' animation (`Viewport.AnimationMs`, or `durationMs`). */
  animate?: boolean;
  durationMs?: number;
  /** Default `api`. */
  origin?: Origin;
  /** Px clear of the edge it is brought to (`REVEAL_PAD`). */
  pad?: number;
}
/** A reveal's default room past the target, px (the `focusViewport` behavior's `revealPad`). */
export const REVEAL_PAD = 12;

const hulls = BlockHulls.create();

/** Loci, `Coordinates` or a list of them (their projections' hull) revealed; false: no section shows them. */
export function revealTarget(
  viewports: Viewports,
  layout: LayoutManager,
  target: SeqLociAny | Coordinates | readonly (SeqLociAny | Coordinates)[],
  opts: RevealOptions = {},
): boolean {
  const ext = Extents.create(4);
  for (const part of Array.isArray(target) ? target : [target]) viewports.projector.fromLoci(part, ext);
  return revealExtents(viewports, layout, ext, opts);
}

/** Extents (the projector's form) revealed; false: no section shows them. */
export function revealExtents(viewports: Viewports, layout: LayoutManager, ext: Extents, opts: RevealOptions = {}) {
  const move: MoveOpts = { animate: opts.animate !== false, durationMs: opts.durationMs, origin: opts.origin ?? 'api' };
  const pad = opts.pad ?? REVEAL_PAD;
  let any = false;
  for (const ref of opts.section ? [opts.section] : viewports.current) {
    const n = viewports.hull(viewports.sectionIndex(ref), ext, hulls);
    if (n === 0) continue;
    const ranges = hulls.data.subarray(0, 3 * n);
    if (layout.isWrapped(ref)) {
      if (layout.reveal(ref, ranges, n, pad, move)) any = true;
    } else if (viewports.get(ref)?.revealBlocks(ranges, pad, move)) any = true;
  }
  return any;
}
