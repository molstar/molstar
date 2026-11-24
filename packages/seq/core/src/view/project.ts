/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A Location carried onto a view's sections (_spec/README.md "Resolved view", Projection): what is hovered, selected
// or highlighted on one section's system shows on every section it links to. Each part is carried along the view's
// link from a section on its system (the paths composed when the view was resolved, so nothing composes here), else,
// for a system no section shows (a structure's chain, an item system), along the model's static path. One Location keyed
// by system: two sections on one system show the same.

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { EmptyLocation, type Location, type LocationPart, location } from '@molstar/seqx-core/query';
import { type Path, projectRanges, type StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedView } from './view.js';

export interface ProjectOptions {
  /**
   * Close the Location first among the view's systems (a nucleotide gains its residue and lifts to its codon, EC-30,
   * EC-144): what a hover lights. False by default.
   */
  readonly close?: boolean;
  /**
   * Keep the parts already on a section's system (`true`: the closed Location as each section shows it, a hover's), or
   * only those carried from other systems onto systems the Location has no part on (`false`, the default: what a
   * selection adds to the sections it is not on).
   */
  readonly own?: boolean;
  /** Carry parts along links at all (true by default); false keeps the own parts alone. */
  readonly linked?: boolean;
  /**
   * What carries a part through a map read by region (a liftover's, EC-212), where no static link does: from the runs
   * the store holds now (`regionProjector`, the render's). None: such parts are not carried.
   */
  readonly regions?: RegionProjector;
}

/**
 * A part carried through a map read by region (_spec/README.md "Resolved view", Projection: a liftover), from the runs
 * held now: what a viewer gives the projections so hover and selection reach across a liftover both ways (EC-20).
 */
export interface RegionProjector {
  /**
   * Its value moves when the runs it reads may have changed (the store's version): a projection kept is made again, when
   * the model has a map read by region (`reads`).
   */
  readonly version: { readonly value: number };
  /** Whether the model has a map read by region: else nothing it holds changes a projection (and none is made again). */
  reads(model: StaticModel): boolean;
  /**
   * `ranges` on `from` carried onto `to` through a map read by region joining them, every image appended to `out`
   * (unsorted); false when no such map joins them.
   */
  project(model: StaticModel, from: string, ranges: ArrayLike<number>, to: string, out: number[]): boolean;
  /**
   * `[start, end)` on `system` carried onto every block `view` shows: its places named when there are at least two (a
   * card's `Places`, a duplicated stretch's), else undefined.
   */
  places?(model: StaticModel, view: ResolvedView, system: string, start: number, end: number): string | undefined;
}

/**
 * `loc` on the view's sections' systems: per system a section shows, the parts carried to it, a rect's two sides with
 * them (a grid's rectangle lights both ranges on the sections of its systems, and where they link). Empty when nothing
 * reaches a section. Its work is in the parts and the links from their systems, never in every system the view shows
 * times every part (a genome of 3,000 blocks shows 3,000): the view's systems and links are indexed once a resolved
 * view (`indexOf`).
 */
export function projectOntoView(
  model: StaticModel,
  view: ResolvedView,
  loc: Location,
  options: ProjectOptions = {},
): Location {
  if (!loc.parts.length && !loc.rects.length) return EmptyLocation;
  const index = indexOf(view);
  const own = new Set(loc.parts.map((p) => p.system));
  // a rect's sides as positions, carried as the parts are; a side is no own part (its system's section shows it)
  const flat = loc.rects.length ? location({ parts: [...loc.parts, ...loc.rects.flatMap((r) => [r.a, r.b])] }) : loc;
  const src = options.close ? model.close(flat, index.systems) : flat;
  // onto a system holding a part of `loc` only with `own`
  const takes = (system: string) => !!options.own || !own.has(system);
  // a map read by region carries nothing on a model without one (`RegionProjector.reads`)
  const regions = options.regions?.reads(model) ? options.regions : undefined;
  const parts: { system: string; ranges: ArrayLike<number> }[] = [];
  const buffer: number[] = [];
  const carry = (part: LocationPart, system: string, path: Path | undefined) => {
    if (path) parts.push({ system, ranges: projectRanges(path, part.ranges, buffer).slice() });
    else if (regions) {
      // through a map read by region (a liftover): every image of the runs held now (EC-20)
      buffer.length = 0;
      if (regions.project(model, part.system, part.ranges, system, buffer) && buffer.length)
        parts.push({ system, ranges: HalfOpenRanges.normalize(buffer) });
    }
  };
  for (const part of src.parts) {
    const shown = index.shown.has(part.system);
    if (shown && takes(part.system)) {
      parts.push(part);
      // the item's other images on its own system, through a map read by region and back (EC-20)
      buffer.length = 0;
      if (options.linked !== false && regions?.project(model, part.system, part.ranges, part.system, buffer))
        if (buffer.length) parts.push({ system: part.system, ranges: HalfOpenRanges.normalize(buffer) });
    }
    if (options.linked === false) continue;
    const links = index.links.get(part.system);
    if (links || (shown && !regions)) {
      // from a section showing the part's system: along the view's links (a section of blocks shows each block's);
      // the other systems only through a map read by region
      if (regions) {
        for (const system of index.systems)
          if (system !== part.system && takes(system)) carry(part, system, links?.get(system));
      } else if (links)
        for (const [system, path] of links) if (system !== part.system && takes(system)) carry(part, system, path);
      continue;
    }
    // a system no section shows (a structure's chain, an item system): the model's static path to each system shown
    for (const system of index.systems) {
      if (system === part.system || !takes(system)) continue;
      carry(part, system, shown ? undefined : staticPath(model, part.system, system));
    }
  }
  return parts.length ? location({ parts }) : EmptyLocation;
}

/** The model's static path between two systems, when it has both. */
function staticPath(model: StaticModel, from: string, to: string): Path | undefined {
  if (!model.systems.has(from) || !model.systems.has(to)) return undefined;
  return model.path(from, to).path;
}

/**
 * A resolved view's systems in its sections' order (a section of blocks shows its blocks'), as a list and a set, and
 * its links by their systems, the first link from one system to another taken: made once a view (a view is immutable).
 */
interface ViewIndex {
  readonly systems: readonly string[];
  readonly shown: ReadonlySet<string>;
  readonly links: ReadonlyMap<string, ReadonlyMap<string, Path>>;
}

const Indexes = new WeakMap<ResolvedView, ViewIndex>();

function indexOf(view: ResolvedView): ViewIndex {
  let index = Indexes.get(view);
  if (index) return index;
  const shown = new Set<string>();
  for (const s of view.sections)
    if (s.blocks) for (const b of s.blocks) shown.add(b.axis.system);
    else shown.add(s.axis.system);
  const links = new Map<string, Map<string, Path>>();
  for (const l of view.links) {
    let to = links.get(l.fromSystem);
    if (!to) links.set(l.fromSystem, (to = new Map()));
    // only onto a system a section shows (as the projection walks them), the first link kept
    if (shown.has(l.toSystem) && !to.has(l.toSystem)) to.set(l.toSystem, l.path);
  }
  index = { systems: [...shown], shown, links };
  Indexes.set(view, index);
  return index;
}
