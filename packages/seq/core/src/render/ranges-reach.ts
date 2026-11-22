/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { type Issue, issueOf } from '@molstar/core/repr/tree';
import type { QueryItem } from '@molstar/seqx-core/query';
import { projectRanges } from '@molstar/seqx-core/resolve';
import type { Want } from '@molstar/seqx-core/data';
import { RenderIssues } from './issues.js';
import { type LayerSource, boundedLevel, mergeNear, readTiles } from './layer-data.js';
import { collectItems, newItemRows } from './items.js';
import { noteLanes, resetLanes } from './known.js';
import type { LayoutContext, StoreReader } from './layer-type.js';
import { itemHulls, packedLanes } from './ranges-content.js';
import type { RangesPrep } from './ranges-lay.js';

// The lanes of a source's layer over its reachable region (_spec/README.md "Render", Track heights): the items of the
// section's whole extent, its tiles at the level that fits it, packed once they are all read; a track of a source is
// as tall as they need from the first settle, and a pan or a zoom keeps it. The layer's own lanes in view are noted
// beside (`ranges.ts`): a lane the packing here does not model (a label's room) still grows the track.

/** What a prep keeps of its reachable read: the verdict for the fit it was made for, and what was packed. */
interface Reach {
  // the fit the verdict is for
  pieces: ArrayLike<number> | undefined;
  width: number;
  dpr: number;
  coverage: object | undefined;
  /** The fit level holds items and the layer's bounds hold them: the region is read. */
  read: boolean;
  /** The region is read at the fit level but past the layer's bounds: said (an info) on every layout. */
  issue: Issue | undefined;
  regions: number[];
  level: number;
  /** The entries the lanes were packed from (their versions, in the order read): another set is other data. */
  versions: number[];
  packed: boolean;
}

const Reaches = new WeakMap<RangesPrep, Reach>();

/**
 * Wants the tiles of the section's whole extent at the fit's level (a margin want: the view's come first; every frame,
 * a want set is the whole set) and, once every one is read as items and the view has settled, notes the lanes they
 * pack into (`noteLanes`). Nothing is read where the fit level holds no items (a density), past the layer's bounds
 * (`layerLoad` over the same wants as the view's: said, `seq.render.reach`), or where a tile failed. A verdict is kept
 * for its fit and allocates nothing a frame.
 */
export function noteReachable(prep: RangesPrep, src: LayerSource, ctx: LayoutContext, data: StoreReader): void {
  const coverage = data.coverage(src.source);
  if (coverage === undefined) return;
  let r = Reaches.get(prep);
  if (!r) {
    r = {
      pieces: undefined,
      width: -1,
      dpr: -1,
      coverage: undefined,
      read: false,
      issue: undefined,
      regions: [],
      level: 1,
      versions: [],
      packed: false,
    };
    Reaches.set(prep, r);
  }
  const axis = ctx.axis;
  if (r.coverage !== coverage || r.pieces !== axis.pieces || r.width !== ctx.box.width || r.dpr !== axis.dpr) {
    // the data's identity changed (another reading of its source): the lanes it packed are another data's
    if (r.coverage !== undefined && r.coverage !== coverage) {
      resetLanes(prep);
      r.packed = false;
    }
    r.coverage = coverage;
    r.pieces = axis.pieces;
    r.width = ctx.box.width;
    r.dpr = axis.dpr;
    fit(r, src, ctx, data, coverage.levels ?? [1]);
  }
  if (r.issue && !data.layerRefused?.(ctx.key)) ctx.report(r.issue);
  if (!r.read) return;
  const read = readTiles(ctx, data, src, r.regions, r.level, { priority: 'margin' });
  if (read.tooMany || read.failed.length || read.pending.length) return;
  const entries = read.entries;
  if (r.packed && sameEntries(r.versions, entries)) return;
  if (r.packed) resetLanes(prep);
  if (!ctx.settled) return;
  const rows = newItemRows();
  const densities = collectItems(rows, entries, [], src.filter, src.base);
  if (densities.length) return;
  const items: QueryItem[] = [];
  for (let j = 0; j < rows.n; j++) items.push({ key: rows.keys[j], start: rows.starts[j], end: rows.ends[j] });
  noteLanes(prep, Math.min(packedLanes(prep.line === 'gene' ? itemHulls(items) : items), ctx.config.ranges.maxLanes));
  r.versions.length = 0;
  for (const e of entries) r.versions.push(e.version);
  r.packed = true;
}

function sameEntries(versions: readonly number[], entries: readonly { readonly version: number }[]): boolean {
  if (versions.length !== entries.length) return false;
  for (let i = 0; i < entries.length; i++) if (versions[i] !== entries[i].version) return false;
  return true;
}

/** The verdict for a fit: the regions, the level, whether the layer's bounds hold them (`Reach.read`, `issue`). */
function fit(r: Reach, src: LayerSource, ctx: LayoutContext, data: StoreReader, levels: readonly number[]) {
  r.read = false;
  r.issue = undefined;
  // the section's own pieces: its window of the system (a range's), never the whole system's length; the stretches
  // near one another read as one
  const whole = Array.from(r.pieces ?? []);
  r.regions = mergeNear(HalfOpenRanges.normalize(src.back ? projectRanges(src.back, whole) : whole));
  let span = 0;
  for (let i = 0; i + 1 < r.regions.length; i += 2) span += r.regions[i + 1] - r.regions[i];
  if (!(span > 0)) return;
  const level = boundedLevel(levels, span / Math.max(1, r.width * r.dpr), span, data);
  // the finest level holds the items: a coarser one is a density, whose lanes nobody knows
  if (level !== levels[0]) return;
  r.level = level;
  // within the layer's bounds as the view's reads are (the same layer, its tiles by key: the view's are among them)
  const wants: Want[] = [];
  for (let i = 0; i + 1 < r.regions.length; i += 2)
    wants.push({ source: src.source, system: src.system, start: r.regions[i], end: r.regions[i + 1], level, layer: ctx.key });
  if (data.layerLoad && !data.layerLoad(wants).fits) {
    r.issue = issueOf(
      RenderIssues,
      'seq.render.reach',
      'The track’s section is too large to read whole for its height: it is as tall as what is in view',
      { address: ctx.key, hint: RenderIssues.codes['seq.render.reach'].hint },
    );
    return;
  }
  r.read = true;
}
