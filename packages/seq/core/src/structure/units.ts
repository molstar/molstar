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
import {
  atomsSystemOf,
  type ComplexEntry,
  collapseAtoms,
  expandAtoms,
  projectRanges,
  type UnitEntry,
} from '@molstar/seqx-core/resolve';
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
  let said: Set<string> | undefined;
  const say = (i: Issue) => {
    const k = `${i.code}\u0000${i.address ?? ''}`;
    if (said?.has(k)) return;
    (said ??= new Set()).add(k);
    issues.push(i);
  };
  const limits = options.limits === undefined ? DefaultSyncLimits.select : options.limits;
  const sources = sourcesOf(r, loc);
  const out: Out = { parts: [], cut: [], touched: [], runs: 0, positions: 0, over: false };
  if (sources.length)
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
      const images = imagesOn(r, c, sources, b.via, say);
      for (const u of c.units) {
        const im = images.get(u.address);
        if (im) place(r, u, im, options.ligand === 'atoms', limits, out);
      }
      if (out.over && limits?.over === 'skip') break;
    }
  const cut = out.cut.length ? location({ parts: out.cut }) : EmptyLocation;
  if (limits && out.over) {
    say(
      issueOf(
        StructureSyncIssues,
        'seq.structure.gated',
        `More than the gate's ${limits.maxRuns} runs or ${limits.maxLength} positions: ${limits.over === 'units' ? 'whole units sent' : 'nothing sent'}`,
        { hint: StructureSyncIssues.codes['seq.structure.gated'].hint },
      ),
    );
    if (limits.over === 'skip') return { target: EmptyLocation, cut, gated: 'skip', issues };
    const whole = out.touched.map((u) => ({ system: u.address, ranges: observedOf(u) }));
    return { target: location({ parts: whole }), cut, gated: 'units', issues };
  }
  return { target: out.parts.length ? location({ parts: out.parts }) : EmptyLocation, cut, issues };
}

/** A part of the Location, read once: its system, its ranges, and its unit or atoms when it is on one. */
interface Source {
  readonly system: string;
  readonly ranges: ArrayLike<number>;
  /** The unit it is on (its residues, or its atoms collapsed onto them). */
  readonly unit?: UnitEntry;
  /** Its atoms' ranges when it is on a unit's atoms. */
  readonly atoms?: ArrayLike<number>;
}

/** The parts (a rect's sides as parts), each read once: an atoms part collapsed onto its residues system. */
function sourcesOf(r: StructureReads, loc: Location): Source[] {
  const all: LocationPart[] = loc.rects.length
    ? [...loc.parts, ...loc.rects.flatMap((x) => [x.a, x.b])]
    : [...loc.parts];
  const out: Source[] = [];
  for (const part of all) {
    const a = r.atoms(part.system);
    if (!a) {
      const unit = r.unit(part.system);
      out.push({ system: part.system, ranges: part.ranges, ...(unit ? { unit } : {}) });
      continue;
    }
    const unit = r.unit(a.of);
    out.push({ system: a.of, ranges: collapseAtoms(a, part.ranges), ...(unit ? { unit, atoms: part.ranges } : {}) });
  }
  return out;
}

/** One unit's image, unclipped: residues reached (paths and parts on it) and atoms parts on its own atoms. */
interface Image {
  residues: number[];
  atoms?: number[];
}

/**
 * Each part's image on a complex's units, linear in the parts and the units reached (_spec/README.md "Structures", To
 * units): a part on a unit lights it; carried through each molecule once (one search a part and molecule, then each
 * unit's free step down, its start subtracted), so an assembly of thousands of copies searches once, not once a copy.
 * The copy rule (EC-25) as the search keeps it: a part on a unit of the molecule, or a path climbing out of one onto
 * it, never steps down to another; a unit with maps of its own (the spec's or the overlay's) is searched on its own.
 * A unit no part reaches makes nothing.
 */
function imagesOn(
  r: StructureReads,
  c: ComplexEntry,
  sources: readonly Source[],
  via: readonly string[] | undefined,
  say: (i: Issue) => void,
): Map<string, Image> {
  const images = new Map<string, Image>();
  const imageOf = (u: UnitEntry) => {
    let im = images.get(u.address);
    if (!im) images.set(u.address, (im = { residues: [] }));
    return im;
  };
  const searched = (from: string, to: string) => {
    const p = r.path(from, to, via);
    if (p.issue && ForwardedPathCodes.has(p.issue.code)) say(p.issue);
    return p.path;
  };
  const own = r.ownMapped(c);
  for (const s of sources) {
    if (s.unit && s.unit.complex === c.address) {
      const im = imageOf(s.unit);
      // an atoms part is its atoms (its residues, collapsed, are `place`'s to light)
      if (s.atoms) {
        im.atoms ??= [];
        push(im.atoms, s.atoms, 0);
      } else push(im.residues, s.ranges, 0);
    }
    for (const m of r.moleculesOf(c)) {
      // a part on a unit of this molecule reaches no sibling
      if (s.unit?.of === m) continue;
      let image: ArrayLike<number> | undefined;
      if (s.system === m) image = s.ranges;
      else {
        const p = searched(s.system, m);
        // the scratch: an image is pushed onto the units it reaches, never kept
        if (p && !p.edges.some((e) => e.unit === 'up' && e.to === m)) image = projectRanges(p, s.ranges, Scratch);
      }
      if (!image?.length) continue;
      let lo = image[0];
      let hi = image[1];
      for (let i = 2; i + 1 < image.length; i += 2) {
        if (image[i] < lo) lo = image[i];
        if (image[i + 1] > hi) hi = image[i + 1];
      }
      for (const u of r.unitsOf(c, m)) {
        if (own.has(u.address) || u.end === undefined || hi <= u.start || lo >= u.end) continue;
        push(imageOf(u).residues, image, u.start);
      }
    }
    for (const a of own) {
      const u = r.unit(a)!;
      if (u === s.unit || s.unit?.of === u.of) continue;
      const p = searched(s.system, a);
      if (p) push(imageOf(u).residues, projectRanges(p, s.ranges, Scratch), 0);
    }
  }
  return images;
}

interface Out {
  readonly parts: PartInput[];
  readonly cut: PartInput[];
  /** The units lit (an observed position or atom reached), in order: what `over: units` sends whole. */
  readonly touched: UnitEntry[];
  runs: number;
  positions: number;
  /** Past the gate: nothing more is counted or kept, only which units are lit. */
  over: boolean;
}

/** A unit's image clipped to it, unobserved positions cut, a molecule unit whole or on its atoms; counted for the gate. */
function place(r: StructureReads, u: UnitEntry, im: Image, atomsMode: boolean, limits: SyncLimits | null, out: Out) {
  if (u.end === undefined) return;
  const length = u.end - u.start;
  const atomsOf = im.atoms || (atomsMode && u.class !== 'polymer') ? r.atoms(atomsSystemOf(u.address)) : undefined;
  const unitAtoms = atomsMode && u.class !== 'polymer' ? atomsOf : undefined;
  const reached = clip(normalized(im.residues), length);
  const ownAtoms = atomsOf && im.atoms ? HalfOpenRanges.normalize(im.atoms) : undefined;
  let residues =
    ownAtoms && atomsOf
      ? HalfOpenRanges.union(reached, HalfOpenRanges.normalize(collapseAtoms(atomsOf, ownAtoms)))
      : reached;
  if (!residues.length) return;
  if (u.class !== 'polymer' && !unitAtoms) residues = Int32Array.of(0, length);
  const observed = observedOf(u);
  // all reached observed (a hover inside the modeled residues, the common case): nothing to cut, nothing made
  const all = inside(residues, observed);
  const kept = all ? residues : HalfOpenRanges.intersect(residues, observed);
  const lost = all ? HalfOpenRanges.Empty : HalfOpenRanges.subtract(residues, observed);
  if (lost.length) out.cut.push({ system: u.address, ranges: lost });
  if (!kept.length) return;
  let system = u.address;
  let ranges = kept;
  if (unitAtoms) {
    // residues reached by a path light their atoms; an atom part its own atoms only; unobserved residues' none
    const n = Int32Array.of(0, unitAtoms.offsets[unitAtoms.offsets.length - 1]);
    const fromPaths = HalfOpenRanges.intersect(HalfOpenRanges.normalize(expandAtoms(unitAtoms, reached)), n);
    const atoms = ownAtoms ? HalfOpenRanges.union(HalfOpenRanges.intersect(ownAtoms, n), fromPaths) : fromPaths;
    ranges = HalfOpenRanges.subtract(atoms, HalfOpenRanges.normalize(expandAtoms(unitAtoms, lost)));
    system = unitAtoms.address;
    if (!ranges.length) return;
  }
  out.touched.push(u);
  if (out.over) return;
  out.runs += ranges.length / 2;
  out.positions += HalfOpenRanges.size(ranges);
  if (limits && (out.runs > limits.maxRuns || out.positions > limits.maxLength)) {
    out.over = true;
    out.parts.length = 0;
    return;
  }
  out.parts.push({ system, ranges });
}

/** The projections' scratch (an image is read at once, pushed onto units, never kept). */
const Scratch: number[] = [];

/** Ranges normalized: copied once when they already are (sorted, disjoint, non-empty runs), else `normalize`d. */
function normalized(r: ArrayLike<number>): Int32Array {
  if (r.length % 2) return HalfOpenRanges.normalize(r);
  for (let i = 0; i + 1 < r.length; i += 2)
    if (r[i + 1] <= r[i] || (i > 0 && r[i] <= r[i - 1])) return HalfOpenRanges.normalize(r);
  return Int32Array.from(r);
}

/** Normalized ranges clipped to `[0, length)`: the same array when they lie within it. */
function clip(r: Int32Array, length: number): Int32Array {
  return !r.length || (r[0] >= 0 && r[r.length - 1] <= length)
    ? r
    : HalfOpenRanges.intersect(r, Int32Array.of(0, length));
}

/** Whether every run of `a` lies in a run of `b` (both normalized). */
function inside(a: ArrayLike<number>, b: ArrayLike<number>): boolean {
  let j = 0;
  for (let i = 0; i + 1 < a.length; i += 2) {
    while (j + 1 < b.length && b[j + 1] <= a[i]) j += 2;
    if (j + 1 >= b.length || b[j] > a[i] || b[j + 1] < a[i + 1]) return false;
  }
  return true;
}

/** `ranges` shifted down by `start` onto `out`. */
function push(out: number[], ranges: ArrayLike<number>, start: number) {
  for (let i = 0; i < ranges.length; i++) out.push(ranges[i] - start);
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
