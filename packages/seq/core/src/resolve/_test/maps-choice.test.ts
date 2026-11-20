/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Two maps joining one pair of systems (structure-syncing.md, Choosing between two maps; _spec/README.md "Resolve",
// Paths): a map's `priority` ranks the distinct shortest candidates by their maps' priorities sorted weakest first
// and compared in turn, a full tie stays ambiguous; a
// use's `via` leaves out the named maps' rivals. EC-22 to EC-27 on every encoding (implicit, runs, expanded) with a
// ranked pair, a tie and a `via`; EC-24 as amended: two maps of one pair make paths through them ambiguous unless
// ranked or named.

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { partOf, rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { ImplicitLinear, isEmptyLinear, linearsEqual } from '../linear.js';
import { viaRivals, withoutRivals } from '../maps-via.js';
import { createResolver, type ResolveConfig, resolveSpec, type StaticModel } from '../model.js';
import { composeEdges, type Edge, findPath, type Path, type PathGraph, type PathSearch } from '../paths.js';

const Encodings: ResolveConfig['encoding'][] = ['cheapest', 'runs', 'expanded'];
const UP = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN';
const ENT = 'PQSDPSVEPPLSQETFSDLWKLLP';

interface ChoiceOptions {
  /** `sifts`' priority (absent: none written). */
  readonly sifts?: number;
  /** `aln`'s priority; `null` leaves the map out. */
  readonly aln?: number | null;
  /** `aln` the same runs as `sifts` (equal compositions). */
  readonly same?: boolean;
  /** A per-unit map `pa` from the sequence onto unit A, with this priority. */
  readonly perUnit?: number;
  /** A longer route through `/data/x` (two maps of this priority) beside a direct map `direct` of priority -10. */
  readonly longer?: number;
  readonly title?: string;
}

/**
 * UniProt-like `up` (30 residues) and an entity `ent` (24) with units A and B (copies) and C (`[10, 20)`); `sifts`
 * a structure map up `[4, 24)` onto ent from 2, `aln` a segments map from 3 (one off): two maps of one pair.
 */
function choiceSpec(o: ChoiceOptions = {}): Node {
  const b = createSeqSpec({ title: o.title ?? 'two maps of one pair' });
  const up = b.data.sequence('up', { alphabet: 'protein', letters: UP });
  const ent = b.data.sequence('ent', { alphabet: 'protein', letters: ENT });
  let unitA = '';
  b.data.complex('c', {}, (c) => {
    unitA = c.unit('A', { of: ent, chain: 'A' });
    c.unit('B', { of: ent, chain: 'B' });
    c.unit('C', { of: ent, range: [10, 20], chain: 'C' });
  });
  b.data.map('sifts', {
    from: up,
    to: ent,
    method: 'structure',
    segments: { from: [4], to: [2], length: [20], provenance: ['sifts'] },
    status: { start: [8], end: [9], status: ['conflict'] },
    ...(o.sifts !== undefined ? { priority: o.sifts } : {}),
  });
  if (o.aln !== null)
    b.data.map('aln', {
      from: up,
      to: ent,
      method: 'segments',
      segments: { from: [4], to: [o.same ? 2 : 3], length: [20] },
      ...(o.aln !== undefined ? { priority: o.aln } : {}),
    });
  if (o.perUnit !== undefined)
    b.data.map('pa', {
      from: up,
      to: unitA,
      method: 'segments',
      segments: { from: [4], to: [4], length: [20] },
      priority: o.perUnit,
    });
  if (o.longer !== undefined) {
    const x = b.data.sequence('x', { alphabet: 'protein', letters: UP });
    b.data.map('ux', { from: up, to: x, method: 'offset', by: 0, priority: o.longer });
    b.data.map('xe', { from: x, to: ent, method: 'offset', by: -5, priority: o.longer });
    b.data.map('direct', { from: up, to: ent, method: 'offset', by: -4, priority: -10 });
  }
  return b.build();
}

const project = (m: StaticModel, from: string, r: number[], to: string, via?: string[]) =>
  Array.from(partOf(m.project(rangesLocation(from, r), to, { via }).location, to)?.ranges ?? []);
const maps = (m: StaticModel, from: string, to: string, via?: string[]) => {
  const r = m.path(from, to, via);
  return r.path ? r.path.edges.filter((e) => !e.unit).map((e) => e.map) : r.issue?.code;
};

describe('two maps of one pair: priority (EC-22 to EC-27, every encoding)', () => {
  for (const encoding of Encodings) {
    const resolve = (o: ChoiceOptions) => resolveSpec(choiceSpec(o), { config: { encoding } });

    it(`a tie stays ambiguous, and its hint names priority and via (EC-23, EC-24; ${encoding})`, () => {
      const m = resolve({});
      expect(m.issues).toEqual([]);
      const r = m.path('/data/up', '/data/ent');
      expect(r.issue?.code).toBe('seq.resolve.path.ambiguous');
      expect(r.issue?.hint).toMatch(/priority/);
      expect(r.issue?.hint).toMatch(/via/);
      // a tie at 2 against 2 is still a tie
      expect(maps(resolve({ sifts: 2, aln: 2 }), '/data/up', '/data/c/A')).toBe('seq.resolve.path.ambiguous');
    });

    it(`a ranked pair: the higher-ranked map is taken, onto the molecule and down to each copy (EC-24, EC-25; ${encoding})`, () => {
      const m = resolve({ sifts: 1 });
      expect(maps(m, '/data/up', '/data/ent')).toEqual(['/data/sifts']);
      expect(project(m, '/data/up', [10, 11], '/data/ent')).toEqual([8, 9]);
      // each copy apart (the copy rule): A and B alike, each from the sequence
      expect(project(m, '/data/up', [10, 11], '/data/c/A')).toEqual([8, 9]);
      expect(project(m, '/data/up', [10, 11], '/data/c/B')).toEqual([8, 9]);
      expect(m.path('/data/c/A', '/data/c/B').issue?.code).toBe('seq.resolve.path.none');
      // below the default: the other one
      const n = resolve({ aln: -1 });
      expect(project(n, '/data/up', [10, 11], '/data/c/A')).toEqual([8, 9]);
      const o = resolve({ sifts: -1 });
      expect(project(o, '/data/up', [10, 11], '/data/c/A')).toEqual([9, 10]);
    });

    it(`equal compositions count once, whatever their ranks (EC-23; ${encoding})`, () => {
      const m = resolve({ same: true, aln: 5 });
      expect(project(m, '/data/up', [4, 24], '/data/ent')).toEqual([2, 22]);
      expect(resolve({ same: true }).path('/data/up', '/data/ent').path).toBeDefined();
    });

    it(`a priority never makes a longer path beat a shorter one (EC-22; ${encoding})`, () => {
      const m = resolve({ aln: null, sifts: 0, longer: 10 });
      // up → ent: sifts and direct (one map each) beat the route through x, whatever its rank
      expect(maps(m, '/data/up', '/data/ent')).toEqual(['/data/sifts']);
      const tied = resolve({ aln: null, sifts: -10, longer: 10 });
      expect(maps(tied, '/data/up', '/data/ent')).toBe('seq.resolve.path.ambiguous');
    });

    it(`a per-unit map against a map onto the molecule is ranked, never chosen by kind (EC-24, EC-25; ${encoding})`, () => {
      expect(maps(resolve({ aln: null, perUnit: 0 }), '/data/up', '/data/c/A')).toBe('seq.resolve.path.ambiguous');
      const general = resolve({ aln: null, sifts: 1, perUnit: 0 });
      expect(project(general, '/data/up', [10, 11], '/data/c/A')).toEqual([8, 9]);
      const specific = resolve({ aln: null, perUnit: 1 });
      expect(project(specific, '/data/up', [10, 11], '/data/c/A')).toEqual([10, 11]);
      // the per-unit map never reaches B: A's molecule is not climbed out of and down to a copy
      expect(project(specific, '/data/up', [10, 11], '/data/c/B')).toEqual([8, 9]);
    });

    it(`a path ending at a ranged unit is clipped, through the ranked map (EC-26; ${encoding})`, () => {
      const m = resolve({ sifts: 1 });
      expect(project(m, '/data/up', [0, 30], '/data/c/C')).toEqual([0, 10]);
      expect(project(m, '/data/up', [0, 30], '/data/c/C', ['/data/aln'])).toEqual([0, 10]);
      expect(project(m, '/data/c/C', [0, 10], '/data/up', ['/data/aln'])).toEqual([11, 21]);
    });

    it(`a composed path is as sparse as its maps (EC-27; ${encoding})`, () => {
      const m = resolve({ sifts: 1 });
      const p = m.path('/data/up', '/data/c/C').path!;
      expect(p.linear.runs()?.length).toBe(1);
      if (encoding === 'cheapest') expect(p.linear.encoding).toBe('implicit');
      expect(p.linear.encoding).not.toBe('expanded');
    });
  }
});

describe('ranks compared as sorted priorities (EC-24)', () => {
  /** A gene `g` onto `up` by an unranked offset, then up onto `ent` by `sifts` or `aln`. */
  function throughOffset(sifts: number, aln: number): StaticModel {
    const b = createSeqSpec();
    const g = b.data.sequence('g', { alphabet: 'protein', letters: UP });
    const up = b.data.sequence('up', { alphabet: 'protein', letters: UP });
    const ent = b.data.sequence('ent', { alphabet: 'protein', letters: ENT });
    b.data.map('gu', { from: g, to: up, method: 'offset', by: 0 });
    b.data.map('sifts', {
      from: up,
      to: ent,
      method: 'segments',
      segments: { from: [4], to: [2], length: [20] },
      priority: sifts,
    });
    b.data.map('aln', {
      from: up,
      to: ent,
      method: 'segments',
      segments: { from: [4], to: [3], length: [20] },
      priority: aln,
    });
    return resolveSpec(b.build());
  }

  it('a shared unranked map ties, and the rival maps decide; a full tie stays ambiguous', () => {
    // [0, 1] against [0, 0]: the weakest maps tie, the next weakest decide
    expect(maps(throughOffset(1, 0), '/data/g', '/data/ent')).toEqual(['/data/gu', '/data/sifts']);
    expect(maps(throughOffset(-1, 0), '/data/g', '/data/ent')).toEqual(['/data/gu', '/data/aln']);
    expect(maps(throughOffset(2, 2), '/data/g', '/data/ent')).toBe('seq.resolve.path.ambiguous');
  });
});

describe('a use’s via (EC-24)', () => {
  const tie = resolveSpec(choiceSpec());
  const ranked = resolveSpec(choiceSpec({ sifts: 1 }));

  it('names one map of a tied pair: its rivals are left out for that use alone', () => {
    expect(maps(tie, '/data/up', '/data/c/A', ['/data/aln'])).toEqual(['/data/aln']);
    expect(maps(tie, '/data/up', '/data/c/A', ['/data/sifts'])).toEqual(['/data/sifts']);
    expect(maps(tie, '/data/up', '/data/c/A')).toBe('seq.resolve.path.ambiguous');
    // a via beats a priority for its use; the rest keep the ranked one
    expect(project(ranked, '/data/up', [10, 11], '/data/c/A', ['/data/aln'])).toEqual([9, 10]);
    expect(project(ranked, '/data/up', [10, 11], '/data/c/A')).toEqual([8, 9]);
    // naming both leaves the tie
    expect(maps(tie, '/data/up', '/data/ent', ['/data/sifts', '/data/aln'])).toBe('seq.resolve.path.ambiguous');
  });

  it('keeps one entry a pair and via set: the same path object, its inverse read the other way', () => {
    const a = tie.path('/data/up', '/data/ent', ['/data/aln']).path!;
    expect(tie.path('/data/up', '/data/ent', ['/data/aln']).path).toBe(a);
    expect(tie.path('/data/up', '/data/ent', ['/data/aln', '/data/aln']).path).toBe(a);
    const back = tie.path('/data/ent', '/data/up', ['/data/aln']).path!;
    expect(linearsEqual(back.linear, a.linear.inverse())).toBe(true);
  });

  it('leaves alone a path that reaches no named pair, and names nothing when it names no map', () => {
    const m = resolveSpec(choiceSpec({ longer: 0 }));
    // up → x reaches no pair `aln` joins (its rivals `sifts` and `direct` join up and ent)
    expect(maps(m, '/data/up', '/data/x', ['/data/aln'])).toEqual(maps(m, '/data/up', '/data/x'));
    expect(maps(tie, '/data/up', '/data/ent', ['/data/nothing'])).toBe('seq.resolve.path.ambiguous');
  });
});

describe('a carried search and priority', () => {
  it('is dropped when a priority changes (the param is in the map entry’s hash), kept on an edit elsewhere', () => {
    const r = createResolver();
    const a = r.resolve(choiceSpec({ sifts: 1 }));
    const pa = a.path('/data/up', '/data/c/A').path!;
    // an edit elsewhere: the maps and units the search read are the same objects
    const b = r.resolve(choiceSpec({ sifts: 1, title: 'retitled' }));
    expect(b).not.toBe(a);
    expect(b.path('/data/up', '/data/c/A').path).toBe(pa);
    const d = r.resolve(choiceSpec({ sifts: -1 }));
    expect(d.maps.get('/data/sifts')).not.toBe(a.maps.get('/data/sifts'));
    expect(project(d, '/data/up', [10, 11], '/data/c/A')).toEqual([9, 10]);
  });
});

// ---- properties over random graphs ----

/**
 * A random graph of offset maps among `n` systems, a few joining one pair twice, each with a random priority; the last
 * system a molecule with two units (`/u0` its `[0, 50)`, `/u1` its `[50, 100)`), their steps free under the copy rule,
 * and maps reaching the units too.
 */
function randomGraph(rng: PCG, n: number, ranked: boolean): { graph: PathGraph; maps: Edge[]; systems: string[] } {
  const out = new Map<string, Edge[]>();
  const add = (e: Edge) => out.set(e.from, [...(out.get(e.from) ?? []), e]);
  const molecule = `/s${n - 1}`;
  for (const [u, start] of [
    ['/u0', 0],
    ['/u1', 50],
  ] as const) {
    const linear = new ImplicitLinear(start, 1, 0, 50);
    add({ map: u, from: u, to: molecule, backward: false, step: 0, linear, unit: 'up' });
    add({ map: u, from: molecule, to: u, backward: true, step: 0, linear: linear.inverse(), unit: 'down' });
  }
  const systems = [...Array.from({ length: n }, (_, i) => `/s${i}`), '/u0', '/u1'];
  const maps: Edge[] = [];
  const count = 4 + (rng.int() % 6);
  for (let i = 0; i < count; i++) {
    const a = rng.int() % systems.length;
    let b = rng.int() % systems.length;
    if (a === b) b = (b + 1) % systems.length;
    const by = rng.int() % 3;
    // drawn either way, so the ranked graph is the unranked one with priorities
    const drawn = (rng.int() % 5) - 2;
    const priority = ranked ? drawn : 0;
    const linear = new ImplicitLinear(by, 1, 0, 100);
    const fwd: Edge = { map: `/m${i}`, from: systems[a], to: systems[b], backward: false, step: 0, linear, priority };
    maps.push(fwd);
    add(fwd);
    add({ ...fwd, from: fwd.to, to: fwd.from, backward: true, linear: linear.inverse() });
  }
  return { graph: { edgesFrom: (s) => out.get(s) ?? [] }, maps, systems };
}

/**
 * The oracle: every candidate of the shortest length that has any, by the search's rules written out plainly (each
 * system once, unit steps free, no step down from a molecule the path climbed out of, an empty composition dropped).
 */
function shortestCandidates(graph: PathGraph, from: string, to: string, max = 4): { length: number; paths: Path[] } {
  for (let length = 0; length <= max; length++) {
    const paths: Path[] = [];
    const walk = (at: string, edges: Edge[], seen: Set<string>, climbed: Set<string>, cost: number) => {
      if (at === to) {
        const p = cost === length ? composeEdges(from, edges) : undefined;
        if (p && !isEmptyLinear(p.linear)) paths.push(p);
        return;
      }
      for (const e of graph.edgesFrom(at)) {
        const c = cost + (e.unit ? 0 : 1);
        if (seen.has(e.to) || c > length || (e.unit === 'down' && climbed.has(e.from))) continue;
        walk(e.to, [...edges, e], new Set([...seen, e.to]), e.unit === 'up' ? new Set([...climbed, e.to]) : climbed, c);
      }
    };
    walk(from, [], new Set([from]), new Set(), 0);
    if (paths.length) return { length, paths };
  }
  return { length: -1, paths: [] };
}

const counted = (r: PathSearch) => (r.ok ? r.path.edges.filter((e) => !e.unit).length : undefined);

/** Whether some simple path of at most `max` counted steps from `from` to `to` steps through one of `maps`. */
function touches(graph: PathGraph, from: string, to: string, maps: ReadonlySet<string>, max = 4): boolean {
  const go = (at: string, seen: Set<string>, depth: number, hit: boolean): boolean => {
    if (at === to) return hit;
    for (const e of graph.edgesFrom(at)) {
      if (seen.has(e.to) || depth + (e.unit ? 0 : 1) > max) continue;
      seen.add(e.to);
      // a unit's step is free, as in the search
      const r = go(e.to, seen, depth + (e.unit ? 0 : 1), hit || maps.has(e.map));
      seen.delete(e.to);
      if (r) return true;
    }
    return false;
  };
  return go(from, new Set([from]), 0, false);
}

describe('paths: properties at small seeds (EC-22, EC-24)', () => {
  it('a priority never lengthens a path: the ranked path is one of the unranked search’s shortest candidates (EC-25 too)', () => {
    const rng = new PCG(23);
    let settled = 0;
    let crossed = 0;
    for (let c = 0; c < 60; c++) {
      const seed = rng.int();
      const plain = randomGraph(new PCG(seed), 6, false).graph;
      const { graph: ranked, systems } = randomGraph(new PCG(seed), 6, true);
      for (const a of systems)
        for (const b of systems) {
          if (a === b) continue;
          const p = findPath(plain, a, b);
          const q = findPath(ranked, a, b);
          const oracle = shortestCandidates(plain, a, b);
          // ranking only settles ties: never a path where none is, never another length
          if (!p.ok && p.why === 'none') expect(q.ok).toBe(false);
          if (q.ok) {
            expect(counted(q)).toBe(oracle.length);
            expect(oracle.paths.some((o) => linearsEqual(o.linear, q.path.linear))).toBe(true);
            if (q.path.edges.some((e) => e.unit)) crossed++;
          }
          if (q.ok && !p.ok) {
            expect(p.why).not.toBe('none');
            settled++;
          }
          // an unambiguous unranked path stays the same map
          if (p.ok) expect(q.ok && linearsEqual(q.path.linear, p.path.linear)).toBe(true);
        }
    }
    // the seeds hold ties that a rank settles, and ranked paths through units
    expect(settled).toBeGreaterThan(0);
    expect(crossed).toBeGreaterThan(0);
  });

  it('a via leaves alone every path that no rival lies on', () => {
    const rng = new PCG(29);
    let compared = 0;
    for (let c = 0; c < 60; c++) {
      const { graph, maps, systems } = randomGraph(rng, 6, true);
      const via = [maps[rng.int() % maps.length].map];
      const rivals = viaRivals(
        via,
        maps.map((m) => ({ address: m.map, from: m.from, to: m.to })),
      );
      const filtered = withoutRivals(graph, rivals);
      for (const a of systems)
        for (const b of systems) {
          if (a === b || touches(graph, a, b, rivals)) continue;
          const p = findPath(graph, a, b);
          const q = findPath(filtered, a, b);
          expect(q.ok).toBe(p.ok);
          if (p.ok && q.ok) expect(linearsEqual(q.path.linear, p.path.linear)).toBe(true);
          else if (!p.ok && !q.ok) expect(q.why).toBe(p.why);
          if (p.ok) compared++;
        }
    }
    expect(compared).toBeGreaterThan(0);
  });
});
