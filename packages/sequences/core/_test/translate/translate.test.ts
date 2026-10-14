/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// TR1 (seqspec:libraries#translate): the registry, translation through one and composed maps, copies, chimeras, complex
// runs across units, keys kept, the builders and the index caches.

import {
  type CoordinateMap,
  type CoordinateSystem,
  CoordinateSystemKind,
  Coords,
  Extents,
  NO_LABEL,
  Numbering,
  StringTable,
} from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { IndexKinds, type SystemDescriptor, TranslateMaps, Translator } from '../../src/translate/index.js';

const K = CoordinateSystemKind;
const unit = (id: number): CoordinateSystem => ({ kind: K.unit, id, row: -1 });
const molecule = (id: number): CoordinateSystem => ({ kind: K.molecule, id, row: -1 });
const desc = (key: string, stop: number, more: Partial<SystemDescriptor> = {}): SystemDescriptor => ({
  key,
  kind: 'external',
  start: 0,
  stop,
  refs: {},
  members: null,
  structure: null,
  numberings: [],
  owner: 'test',
  ...more,
});
const shift = (from: CoordinateSystem, to: CoordinateSystem, at: number, length: number, toLo: number) =>
  Coords.create(from, to, 1, [at, length, toLo, 1]);

/** Runs in one system, as (kind, id, start, stop) tuples; `runs` the input as [system, start, stop]. */
function translate(
  t: Translator,
  from: CoordinateSystem,
  to: CoordinateSystem,
  ...runs: [CoordinateSystem, number, number][]
) {
  const input = Extents.create(1);
  for (const [s, a, b] of runs) Extents.add(input, s, a, b);
  const out = Extents.create(1);
  t.translate(input, from, to, out);
  return Array.from({ length: out.count }, (_, i) => [out.kind[i], out.id[i], out.start[i], out.stop[i]]);
}
function at(t: Translator, from: CoordinateSystem, index: number, to: CoordinateSystem) {
  const seen: string[] = [];
  const n = t.forEachAt(from, index, to, (i, c, s) => c.push(`${s.kind}:${s.id}:${i}`), seen);
  expect(n).toBe(seen.length);
  return seen;
}

describe('translate', () => {
  it('goes through one map, both ways, and allocates its route once', () => {
    const t = new Translator();
    const u = unit(7);
    const x = t.register(desc('host:x', 20), (s) => [shift(s, u, 0, 10, 100)]);
    expect(x).toEqual({ kind: K.external, id: 0, row: -1 });
    expect(translate(t, x, u, [x, 2, 5])).toEqual([[K.unit, 7, 102, 105]]);
    expect(translate(t, u, x, [u, 98, 103])).toEqual([[K.external, 0, 0, 3]]);
    expect(translate(t, x, u, [x, 12, 15])).toEqual([]); // outside the support
    expect(at(t, x, 4, u)).toEqual(['1:7:104']);
    expect(at(t, x, x.id, x)).toEqual(['5:0:0']); // the identity route
  });

  it('composes maps, a codon step either way, split by an intron', () => {
    const t = new Translator();
    const [p, g] = [molecule(1), molecule(2)];
    // Protein p over gene g: residues 0–4 on [6, 21), 5–9 on [30, 45).
    t.addMap(Coords.create(p, g, 3, [0, 15, 6, 1, 15, 15, 30, 1]));
    const x = t.register(desc('host:x', 10), (s) => [shift(s, p, 0, 10, 0)]);
    expect(translate(t, x, g, [x, 4, 6])).toEqual([
      [K.molecule, 2, 18, 21],
      [K.molecule, 2, 30, 33],
    ]);
    expect(translate(t, g, x, [g, 19, 31])).toEqual([[K.external, 0, 4, 6]]); // backward through the span map
    expect(at(t, g, 31, x)).toEqual(['5:0:5']);
    expect(at(t, x, 4, g)).toEqual(['0:2:18', '0:2:19', '0:2:20']);
  });

  it('visits every copy of a source used twice', () => {
    const t = new Translator();
    const [ub, diub] = [molecule(1), molecule(2)];
    t.addMap(Coords.create(diub, ub, 1, [0, 76, 0, 1, 76, 76, 0, 1]));
    expect(translate(t, ub, diub, [ub, 47, 48])).toEqual([
      [K.molecule, 2, 47, 48],
      [K.molecule, 2, 123, 124],
    ]);
    expect(at(t, ub, 47, diub)).toEqual(['0:2:47', '0:2:123']);
    expect(at(t, diub, 123, ub)).toEqual(['0:1:47']);
  });

  it("takes a chimera's partial maps, each part to its source", () => {
    const t = new Translator();
    const [gfp, p53] = [molecule(1), molecule(2)];
    const c = t.register(desc('host:chimera', 300), (s) => [shift(s, gfp, 0, 238, 0), shift(s, p53, 250, 50, 0)]);
    expect(translate(t, c, gfp, [c, 200, 260])).toEqual([[K.molecule, 1, 200, 238]]);
    expect(translate(t, c, p53, [c, 200, 260])).toEqual([[K.molecule, 2, 0, 10]]);
    expect(translate(t, gfp, p53, [gfp, 0, 238])).toEqual([]); // their parts do not overlap
    expect(at(t, c, 245, gfp)).toEqual([]); // the linker
  });

  it('translates complex runs across units in one call, member to member', () => {
    const t = new Translator();
    const [ua, ub, mol] = [unit(1), unit(2), molecule(9)];
    // The copy rule's trap: unit A → molecule → unit B is a longer route, never taken while A's own map is shorter.
    t.addMap(shift(ua, mol, 0, 100, 0));
    t.addMap(shift(ub, mol, 0, 100, 0));
    t.register(desc('seqstar:A', 100, { kind: 'unit' }), [], ua);
    t.register(desc('seqstar:B', 100, { kind: 'unit' }), [], ub);
    const seq = t.register(desc('seqstar:1brs', 0, { kind: 'complex', members: ['seqstar:A', 'seqstar:B'] }), [], {
      kind: K.complex,
      id: 0,
      row: -1,
    });
    const labels = Int32Array.from({ length: 100 }, (_, i) => i + 1);
    const chain = (c: string, u: CoordinateSystem) =>
      t.register(desc(`molstar:3d/${c}@1`, 98, { structure: { labelAsymId: c } }), (s) => [
        TranslateMaps.byLabels(s, labels.subarray(2), u, labels),
      ]);
    const [ca, cb] = [chain('A', ua), chain('B', ub)];
    const structure = t.register(
      desc('molstar:3d', 0, { kind: 'complex', members: ['molstar:3d/A@1', 'molstar:3d/B@1'] }),
    );
    expect(translate(t, structure, seq, [ca, 0, 5], [cb, 10, 12])).toEqual([
      [K.unit, 1, 2, 7],
      [K.unit, 2, 12, 14],
    ]);
    expect(translate(t, seq, structure, [ua, 2, 7], [ub, 0, 3])).toEqual([
      [K.external, ca.id, 0, 5],
      [K.external, cb.id, 0, 1],
    ]);
    expect(at(t, ua, 50, structure)).toEqual([`5:${ca.id}:48`]);
    expect(t.index(structure, IndexKinds.chainUnits)).toEqual(
      new Map([
        ['A', [ca]],
        ['B', [cb]],
      ]),
    );
  });

  it("keeps a key's id through unregister and register, and replaces the maps it brought", () => {
    const t = new Translator();
    const u = unit(3);
    const maps = (s: CoordinateSystem) => [shift(s, u, 0, 10, 0)];
    const x = t.register(desc('host:x', 10), maps);
    t.register(desc('host:y', 10));
    expect(translate(t, x, u, [x, 0, 2])).toEqual([[K.unit, 3, 0, 2]]);
    t.unregister('host:x');
    expect(t.system('host:x')).toBe(x);
    expect(t.descriptor('host:x')).toBeNull();
    expect(translate(t, x, u, [x, 0, 2])).toEqual([]);
    expect(at(t, x, 0, u)).toEqual([]);
    expect(t.register(desc('host:x', 10), maps)).toBe(x);
    t.register(desc('host:x', 10), maps); // again: its maps replaced, not added
    expect(translate(t, x, u, [x, 0, 2])).toEqual([[K.unit, 3, 0, 2]]);
    expect(t.system('host:none')).toBeNull();
  });

  it("names a system's neighbors, one map away either way, each once, until a change (IN1)", () => {
    const t = new Translator();
    const [u, v] = [unit(3), unit(4)];
    const x = t.register(desc('host:x', 10), (s) => [
      shift(s, u, 0, 10, 0),
      shift(s, u, 0, 5, 20),
      shift(s, v, 0, 4, 0),
    ]);
    expect(t.neighbors(x)).toEqual([u, v]);
    expect(t.neighbors(x)).toBe(t.neighbors(x));
    expect(t.neighbors(u)).toEqual([x]);
    t.unregister('host:x');
    expect([t.neighbors(x), t.neighbors(u), t.neighbors(molecule(9))]).toEqual([[], [], []]);
  });
});

describe('TranslateMaps', () => {
  const [a, b] = [unit(1), unit(2)];
  const segments = (m: CoordinateMap) =>
    Array.from({ length: m.count }, (_, k) => [m.at[k], m.length[k], m.toLo[k], m.dir[k]]);

  it('byLabels pairs equal labels, the first to position of each, and skips the unlabeled', () => {
    const from = Int32Array.of(5, 6, 7, NO_LABEL, 9, 10, 42);
    const to = Int32Array.of(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10);
    expect(segments(TranslateMaps.byLabels(a, from, b, to))).toEqual([
      [0, 3, 4, 1],
      [4, 2, 8, 1],
    ]);
  });

  it('byNumbering reads auth numbers and insertion codes', () => {
    const strings = new StringTable();
    const n = Numbering.build({ name: 'auth', method: 'labels', labels: [10, 11, '11A', 12, 13] }, 5, strings, [], '');
    const labels = Int32Array.of(11, 11, 12, 13, 99, 11);
    const map = TranslateMaps.byNumbering(a, n as Numbering, b, labels, [null, 'A', null, null, null, 'Z']);
    expect(segments(map)).toEqual([[1, 4, 0, 1]]);
    const offset = Numbering.build({ name: 'o', method: 'offset', first: 100 }, 50, strings, [], '') as Numbering;
    expect(segments(TranslateMaps.byNumbering(a, offset, b, Int32Array.of(100, 101, 200, 120)))).toEqual([
      [0, 2, 0, 1],
      [20, 1, 3, 1],
    ]);
  });

  it('bySegments takes a SIFTS-style table with holes', () => {
    const map = TranslateMaps.bySegments(a, b, Int32Array.of(30, 0, 10, 0, 100, 20));
    expect(segments(map)).toEqual([
      [0, 20, 100, 1],
      [30, 10, 0, 1],
    ]);
  });
});

describe('Translator.index', () => {
  it('builds each kind once per system, again after the system is registered again', () => {
    const t = new Translator();
    let reads = 0;
    const atoms = IndexKinds.atomRanges('atoms', () => {
      reads++;
      return Int32Array.of(0, 0, 0, 1, 2, 2);
    });
    const x = t.register(desc('host:x', 4));
    expect(Array.from(t.index(x, atoms))).toEqual([0, 3, 4, 6, 6]);
    t.index(x, atoms);
    expect(reads).toBe(1);
    t.register(desc('host:x', 4));
    t.index(x, atoms);
    expect(reads).toBe(2);
    const labels = IndexKinds.labelIndex('labels', () => Int32Array.of(5, 6, NO_LABEL, 6));
    expect(t.index(x, labels)).toEqual(
      new Map([
        [5, 0],
        [6, 1],
      ]),
    );
    const sifts = IndexKinds.siftsTable('sifts', () => Int32Array.of(30, 0, 10, 0, 100, 20));
    expect(Array.from(t.index(x, sifts))).toEqual([0, 100, 20, 30, 0, 10]);
    const own = { name: 'own', build: (s: CoordinateSystem) => ({ id: s.id }) };
    expect(t.index(x, own)).toBe(t.index(x, own));
  });
});
