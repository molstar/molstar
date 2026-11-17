/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The path search over small hand-made graphs: shortest, bounded, at most one codon step, ambiguity, the copy rule
// for units, a path ending at a ranged unit, sparse composition (EC-19 to EC-27).

import { describe, expect, it } from 'vitest';
import { codonLinear } from '../codon.js';
import { ImplicitLinear, type Linear } from '../linear.js';
import { type Edge, findPath, type PathGraph, type PathSearch, projectRanges, type Step } from '../paths.js';

interface MapDef {
  readonly from: string;
  readonly to: string;
  readonly linear: Linear;
  readonly step?: Step;
  readonly unit?: boolean;
}

/** A graph of maps, each read both ways; a unit map is the unit's place on its molecule (`up` forward). */
function graph(maps: Record<string, MapDef>): PathGraph {
  const out = new Map<string, Edge[]>();
  const add = (e: Edge) => {
    const list = out.get(e.from) ?? [];
    list.push(e);
    out.set(e.from, list);
  };
  for (const [map, m] of Object.entries(maps)) {
    const step = m.step ?? 0;
    add({ map, from: m.from, to: m.to, backward: false, step, linear: m.linear, ...(m.unit ? { unit: 'up' } : {}) });
    add({
      map,
      from: m.to,
      to: m.from,
      backward: true,
      step: -step as Step,
      linear: m.linear.inverse(),
      ...(m.unit ? { unit: 'down' } : {}),
    });
  }
  return { edgesFrom: (s) => out.get(s) ?? [] };
}

const offset = (from: string, to: string, by: number, length = Infinity): MapDef => ({
  from,
  to,
  linear: new ImplicitLinear(by, 1, Math.max(0, -by), length),
});
/** A unit `[start, end)` of its molecule. */
const unit = (u: string, m: string, start: number, end: number): MapDef => ({
  from: u,
  to: m,
  linear: new ImplicitLinear(start, 1, 0, end - start),
  unit: true,
});

const via = (r: PathSearch) => (r.ok ? r.path.edges.map((e) => e.map) : r.why);

describe('paths: search (EC-19 to EC-27)', () => {
  it('takes the shortest path; equal compositions count once (EC-22, EC-23)', () => {
    // a → b → c (+2, +3) and a → d → c (+1, +4): both +5, one map
    const g = graph({
      ab: offset('/a', '/b', 2),
      bc: offset('/b', '/c', 3),
      ad: offset('/a', '/d', 1),
      dc: offset('/d', '/c', 4),
      ae: offset('/a', '/e', 7),
    });
    const r = findPath(g, '/a', '/c');
    expect(r.ok && r.path.linear.forward(0, 1)[0].to).toBe(5);
    expect(via(findPath(g, '/a', '/e'))).toEqual(['ae']);
    // a third way of two steps that maps otherwise makes it ambiguous
    const g2 = graph({
      ab: offset('/a', '/b', 2),
      bc: offset('/b', '/c', 3),
      ae: offset('/a', '/e', 7),
      ec: offset('/e', '/c', 1),
    });
    expect(via(findPath(g2, '/a', '/c'))).toBe('ambiguous');
    expect(via(findPath(g, '/a', '/a'))).toEqual([]);
  });

  it('two maps joining the same systems differently are ambiguous (EC-24)', () => {
    const g = graph({ m1: offset('/a', '/b', 2), m2: offset('/a', '/b', 3), bc: offset('/b', '/c', 0) });
    expect(findPath(g, '/a', '/b')).toEqual({ ok: false, why: 'ambiguous', candidates: 2 });
    expect(via(findPath(g, '/a', '/c'))).toBe('ambiguous');
  });

  it('a candidate that pairs nothing is dropped first', () => {
    const g = graph({
      m1: offset('/a', '/b', 2),
      empty: { from: '/a', to: '/b', linear: new ImplicitLinear(5, 1, 0, 0) },
    });
    expect(via(findPath(g, '/a', '/b'))).toEqual(['m1']);
  });

  it('at most four counted steps, and a bounded number of candidates (EC-22)', () => {
    const chain = (n: number) =>
      graph(Object.fromEntries(Array.from({ length: n }, (_, i) => [`m${i}`, offset(`/s${i}`, `/s${i + 1}`, 1)])));
    expect(via(findPath(chain(4), '/s0', '/s4'))).toHaveLength(4);
    expect(via(findPath(chain(5), '/s0', '/s5'))).toBe('none');
    // three equal maps a → b, three b → c: nine candidates over a bound of four
    const many = graph({
      ab1: offset('/a', '/b', 1),
      ab2: offset('/a', '/b', 1),
      ab3: offset('/a', '/b', 1),
      bc1: offset('/b', '/c', 1),
      bc2: offset('/b', '/c', 1),
      bc3: offset('/b', '/c', 1),
    });
    expect(via(findPath(many, '/a', '/c', { maxCandidates: 4 }))).toBe('limit');
    expect(via(findPath(many, '/a', '/c'))).toHaveLength(2);
  });

  it('at most one codon step: two proteins on one gene have no path (EC-19)', () => {
    const g = graph({
      a: { from: '/pa', to: '/g', linear: codonLinear([0], [30], '+'), step: 1 },
      b: { from: '/pb', to: '/g', linear: codonLinear([30], [60], '+'), step: 1 },
    });
    expect(via(findPath(g, '/pa', '/g'))).toEqual(['a']);
    expect(via(findPath(g, '/pa', '/pb'))).toBe('none');
  });

  it('a protein through its gene onto a chromosome is as sparse as its coding sequence (EC-27)', () => {
    const g = graph({
      codon: { from: '/p', to: '/g', linear: codonLinear([100, 160], [150, 200], '-'), step: 1 },
      chr: { from: '/g', to: '/chr', linear: new ImplicitLinear(7_000_000, 1, 0, 1_000) },
    });
    const r = findPath(g, '/p', '/chr');
    expect(r.ok && r.path.step).toBe(1);
    expect(r.ok && r.path.linear.runs()).toHaveLength(2);
    expect(r.ok && projectRanges(r.path, [0, 1])).toEqual([7_000_197, 7_000_200]);
  });
});

describe('paths: units and the copy rule (EC-25, EC-26)', () => {
  // a molecule of 100 with two units, [0, 40) and [40, 80); a protein mapped onto the molecule
  const maps = {
    u1: unit('/u1', '/m', 0, 40),
    u2: unit('/u2', '/m', 40, 80),
    pm: offset('/p', '/m', 0),
  };

  it('no path climbs from a unit to its molecule and down to a sibling (EC-25)', () => {
    expect(via(findPath(graph(maps), '/u1', '/u2'))).toBe('none');
    // unit steps are free: a unit reaches its molecule at no cost
    expect(via(findPath(graph(maps), '/u1', '/m'))).toEqual(['u1']);
  });

  it('a longer legal path beats a shorter illegal one: the rule is part of the search (EC-25)', () => {
    const g = graph({ ...maps, ux: offset('/u1', '/x', 0), xu: offset('/x', '/u2', 0) });
    expect(via(findPath(g, '/u1', '/u2'))).toEqual(['ux', 'xu']);
  });

  it('a path ending at a ranged unit is clipped to the unit (EC-26)', () => {
    const r = findPath(graph(maps), '/p', '/u1');
    expect(r.ok && projectRanges(r.path, [30, 50])).toEqual([30, 40]);
    const r2 = findPath(graph(maps), '/p', '/u2');
    expect(r2.ok && projectRanges(r2.path, [30, 50])).toEqual([0, 10]);
  });
});
