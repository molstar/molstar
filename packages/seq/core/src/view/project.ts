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
}

/**
 * `loc` on the view's sections' systems: per system a section shows, the parts carried to it, a rect's two sides with
 * them (a grid's rectangle lights both ranges on the sections of its systems, and where they link). Empty when nothing
 * reaches a section.
 */
export function projectOntoView(
  model: StaticModel,
  view: ResolvedView,
  loc: Location,
  options: ProjectOptions = {},
): Location {
  if (!loc.parts.length && !loc.rects.length) return EmptyLocation;
  const systems: string[] = [];
  // a section of blocks shows its blocks' systems (its own axis is their units, on no system)
  for (const s of view.sections)
    for (const system of s.blocks ? s.blocks.map((b) => b.axis.system) : [s.axis.system])
      if (!systems.includes(system)) systems.push(system);
  const own = new Set(loc.parts.map((p) => p.system));
  // a rect's sides as positions, carried as the parts are; a side is no own part (its system's section shows it)
  const flat = loc.rects.length ? location({ parts: [...loc.parts, ...loc.rects.flatMap((r) => [r.a, r.b])] }) : loc;
  const src = options.close ? model.close(flat, systems) : flat;
  const parts: { system: string; ranges: ArrayLike<number> }[] = [];
  const buffer: number[] = [];
  for (const system of systems) {
    if (!options.own && own.has(system)) continue;
    for (const part of src.parts) {
      if (part.system === system) {
        parts.push(part);
        continue;
      }
      if (options.linked === false) continue;
      const path = pathOf(model, view, part, system);
      if (path) parts.push({ system, ranges: projectRanges(path, part.ranges, buffer).slice() });
    }
  }
  return parts.length ? location({ parts }) : EmptyLocation;
}

/**
 * The view's link from a section showing the part's system to one showing `system` (a section of blocks shows each
 * block's); else the model's static path (a system no section shows).
 */
function pathOf(model: StaticModel, view: ResolvedView, part: LocationPart, system: string): Path | undefined {
  let shown = false;
  for (const l of view.links) {
    if (l.fromSystem !== part.system) continue;
    shown = true;
    if (l.toSystem === system) return l.path;
  }
  if (shown || view.sections.some((s) => shows(s, part.system))) return undefined;
  if (!model.systems.has(part.system) || !model.systems.has(system)) return undefined;
  return model.path(part.system, system).path;
}

/** Whether a section shows a system: its axis's, or one of its blocks'. */
const shows = (s: ResolvedView['sections'][number], system: string) =>
  s.blocks ? s.blocks.some((b) => b.axis.system === system) : s.axis.system === system;
