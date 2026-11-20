/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The translation past the 2OCJ shape (_spec/README.md "Structures", To units, From units): units that start past their
// molecule's 0 (the step down and the spread back shift by the start, EC-26), atoms of a molecule of several residues,
// the gate never lighting a unit reached only where nothing was modeled (EC-167, EC-118), and the two consumers in one
// test: a host's description through `structureNodes` into an overlay, and the same units through toUnits and fromUnits,
// with a map onto one of the overlay's units read as its own.

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node } from '@molstar/core/repr/tree';
import { type NodeEntry, nodeKey } from '@molstar/seqx-core/data';
import { structureNodes } from '@molstar/seqx-core/io/structure';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import {
  intersectLocations,
  type Location,
  location,
  locationsEqual,
  partOf,
  rangesLocation,
} from '@molstar/seqx-core/query';
import {
  buildOverlay,
  composeEdges,
  type Edge,
  ImplicitLinear,
  nodeAddress,
  type Overlay,
  resolveSpec,
  type StaticModel,
} from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { Source, zincFinger } from '../../io/_test/fixtures.js';
import { unitFacts } from '../facts.js';
import { structureReads } from '../model.js';
import { fromUnits, toUnits } from '../units.js';

const ranges = (l: Location, system: string) => Array.from(partOf(l, system)?.ranges ?? []);
const systems = (l: Location) => l.parts.map((p) => p.system);

/**
 * A sequence `s` (40) reached from `up` (50, by an offset of -5) with three units: A and B on `[10, 30)` (copies, A
 * modeled on `[2, 15)`, B on `[0, 20)`), C on all of it; a glycan of three residues with six atoms (unit G).
 */
function rangedSpec(): Node {
  const b = createSeqSpec();
  const up = b.data.sequence('up', { alphabet: 'protein', letters: 'A'.repeat(50) });
  const s = b.data.sequence('s', { alphabet: 'protein', letters: 'C'.repeat(40) });
  const gly = b.data.molecule('gly', {
    class: 'branched',
    residues: ['NAG', 'NAG', 'MAN'],
    atoms: {
      residue: [0, 0, 1, 1, 1, 2],
      name: ['C1', 'O4', 'C1', 'C2', 'O4', 'C1'],
      element: ['C', 'O', 'C', 'C', 'O', 'C'],
    },
  });
  b.data.complex('c', {}, (c) => {
    c.unit('A', { of: s, range: [10, 30], observed: { start: [2], end: [15] } });
    c.unit('B', { of: s, range: [10, 30], observed: { start: [0], end: [20] } });
    c.unit('C', { of: s });
    c.unit('G', { of: gly });
  });
  b.data.map('m', { from: up, to: s, method: 'offset', by: -5 });
  return b.build();
}

const ranged = resolveSpec(rangedSpec());

/** The model's own projection of `loc` onto a unit, its observed part (EC-118). */
function expectedOn(m: StaticModel, loc: Location, u: string): Location {
  const parts = loc.parts.flatMap((p) => {
    if (p.system === u) return [p];
    return m.path(p.system, u).path ? m.project(location({ parts: [p] }), u).location.parts : [];
  });
  const image = location({ parts: parts.map((p) => ({ system: u, ranges: p.ranges })) });
  const observed = unitFacts(structureReads(m), structureReads(m).unit(u)!).observed;
  return intersectLocations(image, rangesLocation(u, observed));
}

describe('units past their molecule’s start (EC-26)', () => {
  it('the step down shifts by the start; the model’s projection agrees, for every unit', () => {
    const rng = new PCG(21);
    for (let n = 0; n < 60; n++) {
      const system = rng.int() % 2 ? '/data/up' : '/data/s';
      const len = system === '/data/up' ? 50 : 40;
      const s = rng.int() % len;
      const loc = rangesLocation(system, [s, Math.min(len, s + 1 + (rng.int() % 20))]);
      const got = toUnits(ranged, loc, undefined, { limits: null }).target;
      const want = location({
        parts: ['A', 'B', 'C'].flatMap((u) => expectedOn(ranged, loc, `/data/c/${u}`).parts),
      });
      expect(locationsEqual(got, want), `${system} ${ranges(loc, system)}`).toBe(true);
    }
    // `s` 12..14 is A's and B's 2..4 and C's 12..14
    const r = toUnits(ranged, rangesLocation('/data/s', [12, 14])).target;
    expect([ranges(r, '/data/c/A'), ranges(r, '/data/c/B'), ranges(r, '/data/c/C')]).toEqual([
      [2, 4],
      [2, 4],
      [12, 14],
    ]);
  });

  it('fromUnits spreads a molecule’s part over its units, shifted by each start and clipped to it', () => {
    const r = fromUnits(ranged, rangesLocation('/data/s', [5, 12])).location;
    expect([ranges(r, '/data/c/A'), ranges(r, '/data/c/B'), ranges(r, '/data/c/C')]).toEqual([
      [0, 2],
      [0, 2],
      [5, 12],
    ]);
  });
});

describe('atoms of a molecule of several residues', () => {
  const g = '/data/c/G';
  const ga = `${g}/(atoms)`;

  it('whole: a residue reached lights the unit whole; atoms: its residue’s atoms, an atom part its own atoms', () => {
    expect(ranges(toUnits(ranged, rangesLocation('/data/gly', [1, 2])).target, g)).toEqual([0, 3]);
    const atoms = (loc: Location) => ranges(toUnits(ranged, loc, undefined, { ligand: 'atoms' }).target, ga);
    expect(atoms(rangesLocation('/data/gly', [1, 2]))).toEqual([2, 5]);
    expect(atoms(rangesLocation(ga, [3, 4]))).toEqual([3, 4]);
    // the molecule's own atoms collapse onto its residues, then light the unit's atoms of them
    expect(atoms(rangesLocation('/data/gly/(atoms)', [5, 6]))).toEqual([5, 6]);
    expect(
      atoms(
        location({
          parts: [
            { system: ga, ranges: [0, 1] },
            { system: '/data/gly', ranges: [2, 3] },
          ],
        }),
      ),
    ).toEqual([0, 1, 5, 6]);
    expect(systems(fromUnits(ranged, rangesLocation(ga, [3, 5])).location)).toEqual([g]);
    expect(ranges(fromUnits(ranged, rangesLocation(ga, [3, 5])).location, g)).toEqual([1, 2]);
  });
});

describe('the gate over units reached only where nothing was modeled (EC-167, EC-118)', () => {
  it('a unit reached only in unobserved positions is never sent whole, at any bound, units or skip', () => {
    // A is modeled on 2..15: its part 0..2 reaches nothing modeled; C's runs are lit (a sibling's part stays on it)
    const loc = location({
      parts: [
        { system: '/data/c/C', ranges: [10, 12, 14, 16, 18, 20, 22, 24] },
        { system: '/data/c/A', ranges: [0, 2] },
      ],
    });
    const exact = toUnits(ranged, loc, undefined, { limits: null });
    expect(systems(exact.target)).toEqual(['/data/c/C']);
    const rng = new PCG(5);
    for (let n = 0; n < 40; n++) {
      const maxRuns = rng.int() % 12;
      const maxLength = rng.int() % 30;
      const over = rng.int() % 2 ? 'units' : 'skip';
      const r = toUnits(ranged, loc, undefined, { limits: { maxRuns, maxLength, over } });
      if (!r.gated) {
        expect(locationsEqual(r.target, exact.target)).toBe(true);
        continue;
      }
      if (over === 'skip') expect(r.target.parts).toEqual([]);
      else expect(systems(r.target)).toEqual(systems(exact.target));
    }
    // A only through its unmodeled 0..2: no unit lit, nothing gated
    const only = toUnits(ranged, rangesLocation('/data/c/A', [0, 2]), undefined, {
      limits: { maxRuns: 0, maxLength: 0, over: 'units' },
    });
    expect([only.target.parts, only.gated, ranges(only.cut, '/data/c/A')]).toEqual([[], undefined, [0, 2]]);
  });
});

/** The zinc finger described by a host, read by `structureNodes` into a source's overlay over a spec naming it. */
function hostOverlay(): { model: StaticModel; overlay: Overlay } {
  const b = createSeqSpec();
  b.data.source('zf', { format: 'mmcif', uri: 'zf.cif' });
  const model = resolveSpec(b.build());
  const nodes = new Map<string, NodeEntry>();
  for (const n of structureNodes(zincFinger(), { source: Source }).nodes)
    nodes.set(nodeKey(Source, n.key), { source: Source, key: n.key, node: n.node, version: 1 });
  return { model, overlay: buildOverlay(model, { nodes, nodeVersion: { value: 1 } }) };
}

describe('two consumers: a host’s description and the translation, in one test', () => {
  const complex = nodeAddress(Source, 'complex');
  const A = `${complex}/A-1_555`;
  const B = `${complex}/B-1_555`;
  const seq = nodeAddress(Source, '1');

  it('structureNodes’ units through toUnits and fromUnits over the overlay: copies, the unmodeled residue cut', () => {
    const host = hostOverlay();
    const r = toUnits(host, rangesLocation(seq, [0, 3]));
    // A's residue 0 was not modeled (the description lists 1..4 for chain A)
    expect([ranges(r.target, A), ranges(r.target, B), ranges(r.cut, A)]).toEqual([
      [1, 3],
      [0, 3],
      [0, 1],
    ]);
    expect(locationsEqual(fromUnits(host, r.target).location, r.target)).toBe(true);
    const back = fromUnits(host, rangesLocation(seq, [2, 4])).location;
    expect([ranges(back, A), ranges(back, B)]).toEqual([
      [2, 4],
      [2, 4],
    ]);
    // a copy's part stays on it (EC-25)
    expect(systems(toUnits(host, rangesLocation(B, [0, 2])).target)).toEqual([B]);
  });

  it('a map onto one of the overlay’s units is its own: searched on its own, the copies reached from the molecule', () => {
    const host = hostOverlay();
    const x = '/data/x';
    // a source's map onto unit A (SIFTS per chain): x i ↦ A i + 1, as the overlay would hold it
    const edge: Edge = {
      map: '/data/xa',
      from: A,
      to: x,
      backward: true,
      step: 0,
      linear: new ImplicitLinear(-1, 1, 1, 5),
    };
    const forward = composeEdges(x, [
      { ...edge, from: x, to: A, backward: false, linear: new ImplicitLinear(1, 1, 0, 4) },
    ])!;
    const overlay: Overlay = {
      ...host.overlay,
      molecular: host.overlay.molecular,
      system: (a) => host.overlay.system(a),
      annotation: (a) => host.overlay.annotation(a),
      node: (a) => host.overlay.node(a),
      nested: (n) => host.overlay.nested(n),
      systemNamed: (n, near) => host.overlay.systemNamed(n, near),
      holds: (a) => host.overlay.holds(a),
      edgesFrom: (s) => [...(host.overlay.edgesFrom?.(s) ?? []), ...(s === A ? [edge] : [])],
      path: (from, to, via) => (from === x && to === A ? { path: forward } : host.overlay.path(from, to, via)),
    };
    const r = toUnits({ model: host.model, overlay }, rangesLocation(x, [0, 2]));
    expect(systems(r.target)).toEqual([A]);
    expect(ranges(r.target, A)).toEqual([1, 3]);
    // without the overlay's edges A is no unit of its own: its map is never asked
    const plain = toUnits(host, rangesLocation(x, [0, 2]));
    expect(plain.target.parts).toEqual([]);
  });
});
