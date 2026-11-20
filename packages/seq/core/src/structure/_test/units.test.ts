/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Locations onto units and back (_spec/README.md "Structures", Translation) on a 2OCJ-shaped structure: the copy rule
// (EC-25), a molecule reaching every unit of it (EC-48), copies (EC-65), unobserved positions cut (EC-118), a molecule
// unit whole or on its atoms, the size gate (EC-167) monotone in its bounds, the order of bindings and units
// irrelevant, `via` and `priority` choosing between two maps of one pair; both round trips as properties.

import { PCG } from '@molstar/core/data/util/hash-functions';
import {
  intersectLocations,
  type Location,
  location,
  locationsEqual,
  partOf,
  positionCount,
  rangesLocation,
} from '@molstar/seqx-core/query';
import { resolveSpec, type StaticModel } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { fromUnits, toUnits } from '../units.js';
import { chainId, Complex, Ent, EntLength, EntStart, Observed, p53Spec, Up, UpLength, Zn } from './p53.js';

const unit = (k: number) => `${Complex}/${chainId(k)}`;
const zinc = (k: number) => `${Complex}/Z${chainId(k)}`;
const ranges = (l: Location, system: string) => Array.from(partOf(l, system)?.ranges ?? []);
const systems = (l: Location) => l.parts.map((p) => p.system);

const model = resolveSpec(p53Spec());

describe('toUnits', () => {
  it('a hover on P04637 lights the residue on every copy (EC-48, EC-65), shifted onto each unit', () => {
    expect(model.issues.filter((i) => i.severity !== 'info')).toEqual([]);
    const r = toUnits(model, rangesLocation(Up, [EntStart + 10, EntStart + 13]));
    expect(systems(r.target)).toEqual([0, 1, 2, 3].map(unit));
    for (let k = 0; k < 4; k++) expect(ranges(r.target, unit(k))).toEqual([10, 13]);
    expect(r.issues).toEqual([]);
  });

  it('unobserved positions are cut and said in `cut` (EC-118); a Location wholly unobserved sends nothing', () => {
    const r = toUnits(model, rangesLocation(Up, [EntStart, EntStart + 5]));
    expect(ranges(r.target, unit(0))).toEqual([Observed[0], 5]);
    expect(ranges(r.cut, unit(0))).toEqual([0, Observed[0]]);
    expect(positionCount(toUnits(model, rangesLocation(Up, [EntStart, EntStart + 2])).target)).toBe(0);
    // positions of P04637 the entity does not hold reach no unit
    expect(positionCount(toUnits(model, rangesLocation(Up, [0, 50])).target)).toBe(0);
  });

  it('a Location on one copy stays on it: no path climbs to the molecule and down to a sibling (EC-25)', () => {
    const r = toUnits(model, rangesLocation(unit(1), [20, 30]));
    expect(systems(r.target)).toEqual([unit(1)]);
    expect(ranges(r.target, unit(1))).toEqual([20, 30]);
  });

  it('a Location on the entity reaches every unit of it (EC-48)', () => {
    const r = toUnits(model, rangesLocation(Ent, [50, 51]));
    expect(systems(r.target)).toEqual([0, 1, 2, 3].map(unit));
  });

  it('a molecule unit is whole by default; with `ligand: atoms` it is its atoms', () => {
    const on = rangesLocation(Zn, [0, 1]);
    expect(systems(toUnits(model, on).target)).toEqual([0, 1, 2, 3].map(zinc));
    const atoms = toUnits(model, on, undefined, { ligand: 'atoms' }).target;
    expect(systems(atoms)).toEqual([0, 1, 2, 3].map((k) => `${zinc(k)}/(atoms)`));
    expect(ranges(atoms, `${zinc(0)}/(atoms)`)).toEqual([0, 1]);
  });

  it('only the bound complexes, each with its `via`; a binding naming no complex is said once', () => {
    const r = toUnits(model, rangesLocation(Up, [100, 101]), [{ complex: '/data/none' }, { complex: '/data/none' }]);
    expect(positionCount(r.target)).toBe(0);
    expect(r.issues.map((i) => i.code)).toEqual(['seq.structure.complex']);
  });

  it('two maps of one pair: ambiguous until ranked or named (EC-23, EC-24), the issue forwarded once', () => {
    const tied = resolveSpec(p53Spec({ aln: 0 }));
    const loc = rangesLocation(Up, [EntStart + 20, EntStart + 21]);
    const r = toUnits(tied, loc);
    expect(positionCount(r.target)).toBe(0);
    expect(r.issues.map((i) => i.code)).toEqual(['seq.resolve.path.ambiguous']);
    const viaSifts = toUnits(tied, loc, [{ complex: Complex, via: ['/data/sifts'] }]).target;
    expect(ranges(viaSifts, unit(0))).toEqual([20, 21]);
    const viaAln = toUnits(tied, loc, [{ complex: Complex, via: ['/data/aln'] }]).target;
    expect(ranges(viaAln, unit(0))).toEqual([19, 20]);
    const ranked = toUnits(resolveSpec(p53Spec({ aln: -1 })), loc).target;
    expect(ranges(ranked, unit(0))).toEqual([20, 21]);
  });

  it('the order of bindings, units and parts is irrelevant', () => {
    const reversed = resolveSpec(p53Spec({ reversed: true }));
    const rng = new PCG(7);
    for (let n = 0; n < 30; n++) {
      const loc = randomLocation(rng);
      const a = toUnits(model, loc).target;
      expect(locationsEqual(a, toUnits(reversed, loc).target)).toBe(true);
      const two = [{ complex: Complex }, { complex: Complex, via: ['/data/sifts'] }];
      expect(locationsEqual(toUnits(model, loc, two).target, toUnits(model, loc, [...two].reverse()).target)).toBe(
        true,
      );
    }
  });
});

describe('the size gate (EC-167)', () => {
  const big = location({
    parts: [
      { system: Up, ranges: Array.from({ length: 40 }, (_, i) => [EntStart + 4 + 4 * i, EntStart + 6 + 4 * i]).flat() },
    ],
  });

  it('past `maxRuns`: whole units (`over: units`, their observed runs) or nothing (`skip`), said', () => {
    const units = toUnits(model, big, undefined, { limits: { maxRuns: 10, maxLength: 1e9, over: 'units' } });
    expect(units.gated).toBe('units');
    expect(ranges(units.target, unit(0))).toEqual([...Observed]);
    expect(units.issues.map((i) => i.code)).toEqual(['seq.structure.gated']);
    const skip = toUnits(model, big, undefined, { limits: { maxRuns: 1e9, maxLength: 10, over: 'skip' } });
    expect([skip.gated, positionCount(skip.target)]).toEqual(['skip', 0]);
  });

  it('monotone: a bound raised never gates what a lower one passed; a gated target holds the exact one', () => {
    const rng = new PCG(11);
    for (let n = 0; n < 40; n++) {
      const loc = randomLocation(rng);
      const exact = toUnits(model, loc, undefined, { limits: null }).target;
      let gatedBefore = true;
      for (const maxRuns of [0, 1, 2, 4, 8, 16, 64, 1e9]) {
        const r = toUnits(model, loc, undefined, { limits: { maxRuns, maxLength: 1e9, over: 'units' } });
        const gated = r.gated !== undefined;
        if (!gatedBefore) expect(gated).toBe(false);
        gatedBefore = gated;
        expect(locationsEqual(intersectLocations(r.target, exact), exact)).toBe(true);
        expect(systems(r.target)).toEqual(systems(exact));
      }
    }
  });
});

describe('round trips', () => {
  it('fromUnits(toUnits(loc)) is the Location’s observed part on the units, by the model’s own projection (EC-118)', () => {
    const rng = new PCG(3);
    for (let n = 0; n < 60; n++) {
      const loc = randomLocation(rng);
      const there = toUnits(model, loc, undefined, { limits: null });
      const back = fromUnits(model, there.target);
      expect(back.unplaced).toBe(0);
      const expected = location({
        parts: [0, 1, 2, 3].flatMap((k) => expectedOn(model, loc, unit(k))),
      });
      expect(locationsEqual(back.location, expected), JSON.stringify(ranges(loc, Up))).toBe(true);
    }
  });

  it('toUnits(fromUnits(target)) is the target’s observed part; a molecule part spreads over its units', () => {
    const rng = new PCG(5);
    for (let n = 0; n < 60; n++) {
      const k = rng.int() % 4;
      const s = rng.int() % EntLength;
      const t = rangesLocation(unit(k), [s, Math.min(EntLength, s + 1 + (rng.int() % 30))]);
      const back = toUnits(model, fromUnits(model, t).location, undefined, { limits: null }).target;
      const observed = rangesLocation(unit(k), [...Observed]);
      expect(locationsEqual(back, intersectLocations(t, observed))).toBe(true);
    }
    const spread = fromUnits(model, rangesLocation(Ent, [5, 7])).location;
    expect(systems(spread)).toEqual([0, 1, 2, 3].map(unit));
  });

  it('a unit’s atoms collapse onto its residues unless kept; positions on nothing bound are counted', () => {
    const atoms = rangesLocation(`${zinc(2)}/(atoms)`, [0, 1]);
    expect(systems(fromUnits(model, atoms).location)).toEqual([zinc(2)]);
    expect(systems(fromUnits(model, atoms, undefined, { atoms: true }).location)).toEqual([`${zinc(2)}/(atoms)`]);
    const r = fromUnits(
      model,
      location({
        parts: [
          { system: '/data/elsewhere', ranges: [0, 4] },
          { system: unit(0), ranges: [EntLength - 1, EntLength + 2] },
        ],
      }),
    );
    expect([r.unplaced, r.issues.map((i) => i.code)]).toEqual([6, ['seq.structure.unplaced']]);
  });
});

/** The observed part of the model's own projection of `loc` onto a unit (EC-118). */
function expectedOn(m: StaticModel, loc: Location, u: string) {
  const parts = loc.parts.flatMap((p) => {
    if (p.system === u) return [p];
    const r = m.path(p.system, u);
    return r.path ? [m.project(location({ parts: [p] }), u).location.parts[0]].filter((x) => x !== undefined) : [];
  });
  const image = location({ parts: parts.map((p) => ({ system: u, ranges: p.ranges })) });
  const observed = rangesLocation(u, [...Observed]);
  return intersectLocations(image, observed).parts;
}

/** A Location on P04637, the entity or a unit: a few ranges anywhere. */
function randomLocation(rng: PCG): Location {
  const parts = [];
  for (let i = 1 + (rng.int() % 3); i > 0; i--) {
    const which = rng.int() % 3;
    const system = which === 0 ? Up : which === 1 ? Ent : unit(rng.int() % 4);
    const length = which === 0 ? UpLength : EntLength;
    const r: number[] = [];
    for (let j = 1 + (rng.int() % 4); j > 0; j--) {
      const s = rng.int() % length;
      r.push(s, Math.min(length, s + 1 + (rng.int() % 40)));
    }
    parts.push({ system, ranges: r });
  }
  return location({ parts });
}
