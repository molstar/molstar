/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Locations carried onto a structure's units and back (_spec/README.md "Structures", Translation): along the model's
// paths (static or the overlay's) under the copy rule (EC-25) and a binding's `via`, a molecule reaching every unit of
// it (EC-48); unobserved positions cut (nothing to light, EC-118); a molecule unit whole unless atoms are asked; the
// size gate past its bounds (whole units or nothing, EC-167). Pure: a hover a frame reads kept paths and indices.

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { EmptyLocation, type Location, type LocationPart, location, type PartInput } from '@molstar/seqx-core/query';
import { atomsSystemOf, collapseAtoms, expandAtoms, projectRanges, type UnitEntry } from '@molstar/seqx-core/resolve';
import { observedOf } from './facts.js';
import { StructureSyncIssues } from './issues.js';
import {
  bindingsOf,
  type StructureBinding,
  type StructureModel,
  type StructureReads,
  structureReads,
} from './model.js';

/** What crosses to a structure viewer: a Location whose parts are on units (or a unit's atoms). Structured-cloneable. */
export type StructureTarget = Location;

/** The size gate (EC-167): past `maxRuns` runs or `maxLength` positions, whole units (`units`) or nothing (`skip`). */
export interface SyncLimits {
  readonly maxRuns: number;
  readonly maxLength: number;
  readonly over: 'units' | 'skip';
}

/** The link's defaults: a hover a frame is smaller than a selection sent once. */
export const DefaultSyncLimits: { readonly hover: SyncLimits; readonly select: SyncLimits } = {
  hover: { maxRuns: 2_000, maxLength: 20_000, over: 'units' },
  select: { maxRuns: 20_000, maxLength: 200_000, over: 'units' },
};

export interface ToUnitsOptions {
  /** The gate; absent: `DefaultSyncLimits.select`; `null`: none (colors, rows). */
  readonly limits?: SyncLimits | null;
  /** A molecule unit (not a sequence's) lit whole (`whole`, the default) or on its atoms where it has them (`atoms`). */
  readonly ligand?: 'whole' | 'atoms';
}

export interface UnitsResult {
  readonly target: StructureTarget;
  /** The unobserved positions cut from it, on units. */
  readonly cut: Location;
  /** Whether the gate held it: whole units were sent, or nothing. */
  readonly gated?: 'units' | 'skip';
  readonly issues: readonly Issue[];
}

const ForwardedPathCodes = new Set(['seq.resolve.path.ambiguous', 'seq.resolve.path.limit']);

/** A Location carried onto the units of the bound complexes (every complex when `complexes` is absent). */
export function toUnits(
  model: StructureModel,
  loc: Location,
  complexes?: readonly StructureBinding[],
  options: ToUnitsOptions = {},
): UnitsResult {
  const r = structureReads(model);
  const issues: Issue[] = [];
  const said = new Set<string>();
  const say = (i: Issue) => {
    const k = `${i.code}\u0000${i.address ?? ''}`;
    if (said.has(k)) return;
    said.add(k);
    issues.push(i);
  };
  const atomsMode = options.ligand === 'atoms';
  // a rect's sides are positions as any part's
  const src: readonly LocationPart[] = loc.rects.length
    ? [...loc.parts, ...loc.rects.flatMap((x) => [x.a, x.b])]
    : loc.parts;
  const parts: PartInput[] = [];
  const cut: PartInput[] = [];
  const touched: UnitEntry[] = [];
  let runs = 0;
  let positions = 0;
  if (src.length)
    for (const b of bindingsOf(r, complexes)) {
      const c = r.complex(b.complex);
      if (!c) {
        say(
          issueOf(StructureSyncIssues, 'seq.structure.complex', `${b.complex} is no complex of the model`, {
            address: b.complex,
            hint: StructureSyncIssues.codes['seq.structure.complex'].hint,
          }),
        );
        continue;
      }
      const carry = carrier(r, b.via, say);
      for (const u of c.units) {
        const im = imageOn(r, u, src, carry, atomsMode);
        if (!im) continue;
        touched.push(u);
        const observed = observedOf(u, u.end! - u.start);
        const kept = HalfOpenRanges.intersect(im.residues, observed);
        const lost = HalfOpenRanges.subtract(im.residues, observed);
        if (lost.length) cut.push({ system: u.address, ranges: lost });
        if (im.atoms) {
          const a = r.atoms(atomsSystemOf(u.address))!;
          const ranges = HalfOpenRanges.subtract(im.atoms, HalfOpenRanges.normalize(expandAtoms(a, lost)));
          if (ranges.length) parts.push({ system: a.address, ranges });
          runs += ranges.length / 2;
          positions += HalfOpenRanges.size(ranges);
        } else if (kept.length) {
          parts.push({ system: u.address, ranges: kept });
          runs += kept.length / 2;
          positions += HalfOpenRanges.size(kept);
        }
      }
    }
  const limits = options.limits === undefined ? DefaultSyncLimits.select : options.limits;
  const cutLoc = cut.length ? location({ parts: cut }) : EmptyLocation;
  if (limits && (runs > limits.maxRuns || positions > limits.maxLength)) {
    say(
      issueOf(
        StructureSyncIssues,
        'seq.structure.gated',
        `${runs} runs over ${positions} positions are past the gate (${limits.maxRuns} runs, ${limits.maxLength} positions): ${limits.over === 'units' ? 'whole units sent' : 'nothing sent'}`,
        { hint: StructureSyncIssues.codes['seq.structure.gated'].hint },
      ),
    );
    if (limits.over === 'skip') return { target: EmptyLocation, cut: cutLoc, gated: 'skip', issues };
    const whole = touched.map((u) => ({ system: u.address, ranges: observedOf(u, u.end! - u.start) }));
    return { target: location({ parts: whole }), cut: cutLoc, gated: 'units', issues };
  }
  return { target: parts.length ? location({ parts }) : EmptyLocation, cut: cutLoc, issues };
}

/** A part's image on a unit's positions (unclipped), or undefined when no path reaches it. */
type Carry = (part: number, from: string, ranges: ArrayLike<number>, u: UnitEntry) => ArrayLike<number> | undefined;

/**
 * Parts carried onto units through their molecule: one search a part and molecule, then the unit's free step down
 * (its range's start subtracted), so an assembly of thousands of copies searches once, not once a copy (a search from
 * a sequence reads every unit's step). The copy rule (EC-25) is kept as the search keeps it: a part on a unit of the
 * molecule, or a path climbing out of one onto it, never steps down to another; a unit with maps of its own is
 * searched on its own (its map may be the shorter path). The round trip tests it against the model's projection.
 */
function carrier(r: StructureReads, via: readonly string[] | undefined, say: (i: Issue) => void): Carry {
  const onMolecule = new Map<string, ArrayLike<number> | null>();
  const viaOf = (p: ReturnType<StructureReads['path']>) => {
    if (p.issue && ForwardedPathCodes.has(p.issue.code)) say(p.issue);
    return p.path;
  };
  return (k, from, ranges, u) => {
    if (r.model.edgesFrom(u.address).some((e) => !e.unit)) {
      const p = viaOf(r.path(from, u.address, via));
      return p && projectRanges(p, ranges);
    }
    if (r.unit(from)?.of === u.of) return undefined;
    let image: ArrayLike<number> | null | undefined;
    if (from === u.of) image = ranges;
    else {
      const key = `${k}\u0000${u.of}`;
      image = onMolecule.get(key);
      if (image === undefined) {
        const p = viaOf(r.path(from, u.of, via));
        image = p && !p.edges.some((e) => e.unit === 'up' && e.to === u.of) ? projectRanges(p, ranges) : null;
        onMolecule.set(key, image);
      }
    }
    if (!image) return undefined;
    if (!u.start) return image;
    const out = new Array<number>(image.length);
    for (let i = 0; i < image.length; i++) out[i] = image[i] - u.start;
    return out;
  };
}

/** The image of the parts on one unit: its residues (normalized, within it) and, in atoms mode, its atoms. */
function imageOn(
  r: StructureReads,
  u: UnitEntry,
  src: readonly LocationPart[],
  carry: Carry,
  atomsMode: boolean,
): { residues: Int32Array; atoms?: Int32Array } | undefined {
  if (u.end === undefined) return undefined;
  const length = u.end - u.start;
  const unitAtoms = atomsMode && u.class !== 'polymer' ? r.atoms(atomsSystemOf(u.address)) : undefined;
  const residues: number[] = [];
  const atoms: number[] = [];
  for (let k = 0; k < src.length; k++) {
    const part = src[k];
    if (part.system === u.address) {
      push(residues, part.ranges);
      continue;
    }
    if (unitAtoms && part.system === unitAtoms.address) {
      push(atoms, part.ranges);
      continue;
    }
    const a = r.atoms(part.system);
    const from = a ? a.of : part.system;
    const ranges = a ? collapseAtoms(a, part.ranges) : part.ranges;
    if (from === u.address) {
      push(residues, ranges);
      continue;
    }
    const image = carry(k, from, ranges, u);
    if (image) push(residues, image);
  }
  if (!residues.length && !atoms.length) return undefined;
  let res = HalfOpenRanges.intersect(HalfOpenRanges.normalize(residues), Int32Array.of(0, length));
  if (atoms.length && unitAtoms) {
    // an atom part lights its residues as any part does (for the cut and the whole-unit rule)
    res = HalfOpenRanges.union(res, HalfOpenRanges.normalize(collapseAtoms(unitAtoms, atoms)));
  }
  if (!res.length) return undefined;
  if (u.class !== 'polymer' && !unitAtoms) return { residues: Int32Array.of(0, length) };
  if (!unitAtoms) return { residues: res };
  const n = unitAtoms.offsets[unitAtoms.offsets.length - 1];
  const own = HalfOpenRanges.intersect(HalfOpenRanges.normalize(atoms), Int32Array.of(0, n));
  // residues reached by a path light their atoms; an atom part lights its own atoms only
  const fromPaths = HalfOpenRanges.intersect(HalfOpenRanges.normalize(residues), Int32Array.of(0, length));
  return {
    residues: res,
    atoms: HalfOpenRanges.union(own, HalfOpenRanges.normalize(expandAtoms(unitAtoms, fromPaths))),
  };
}

function push(out: number[], ranges: ArrayLike<number>) {
  for (let i = 0; i < ranges.length; i++) out.push(ranges[i]);
}

export interface FromUnitsOptions {
  /** A unit's atoms kept as atoms (a section shows the `atoms` form); false: collapsed onto its residues. */
  readonly atoms?: boolean;
}

export interface FromUnitsResult {
  readonly location: Location;
  /** The positions on nothing the bound complexes hold (counted, said once). */
  readonly unplaced: number;
  readonly issues: readonly Issue[];
}

/**
 * A structure viewer's target as a Seq\* Location: parts on units kept, on a unit's atoms collapsed onto its residues
 * (unless `atoms`), on a molecule spread over every bound unit of it; the rest counted in `seq.structure.unplaced`.
 */
export function fromUnits(
  model: StructureModel,
  target: StructureTarget,
  complexes?: readonly StructureBinding[],
  options: FromUnitsOptions = {},
): FromUnitsResult {
  const r = structureReads(model);
  const bound = bindingsOf(r, complexes)
    .map((b) => r.complex(b.complex))
    .filter((c) => c !== undefined);
  const inBound = new Set(bound.map((c) => c.address));
  const parts: PartInput[] = [];
  let unplaced = 0;
  const clipped = (ranges: Int32Array, length: number) => {
    const kept = HalfOpenRanges.intersect(ranges, Int32Array.of(0, length));
    unplaced += HalfOpenRanges.size(ranges) - HalfOpenRanges.size(kept);
    return kept;
  };
  for (const part of target.parts) {
    const ranges = part.ranges;
    const u = r.unit(part.system);
    if (u && inBound.has(u.complex) && u.end !== undefined) {
      parts.push({ system: u.address, ranges: clipped(ranges, u.end - u.start) });
      continue;
    }
    const a = r.atoms(part.system);
    const au = a && r.unit(a.of);
    if (a && au && inBound.has(au.complex)) {
      const kept = clipped(ranges, a.offsets[a.offsets.length - 1]);
      parts.push(
        options.atoms ? { system: a.address, ranges: kept } : { system: au.address, ranges: collapseAtoms(a, kept) },
      );
      continue;
    }
    let placed = false;
    for (const c of bound)
      for (const m of r.unitsOf(c, part.system)) {
        if (m.end === undefined) continue;
        placed = true;
        const shifted = new Int32Array(ranges.length);
        for (let i = 0; i < ranges.length; i++) shifted[i] = ranges[i] - m.start;
        parts.push({ system: m.address, ranges: HalfOpenRanges.intersect(shifted, Int32Array.of(0, m.end - m.start)) });
      }
    if (!placed) unplaced += HalfOpenRanges.size(ranges);
  }
  for (const x of target.rects) unplaced += HalfOpenRanges.size(x.a.ranges) + HalfOpenRanges.size(x.b.ranges);
  const issues: Issue[] = unplaced
    ? [
        issueOf(
          StructureSyncIssues,
          'seq.structure.unplaced',
          `${unplaced} positions are on no unit, unit atoms or molecule of the bound complexes`,
        ),
      ]
    : [];
  return { location: parts.length ? location({ parts }) : EmptyLocation, unplaced, issues };
}
