/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Molecules, complexes and units in the static model: a unit a system on its molecule's positions from 0, its steps
// free under the copy rule, copies derived, links held on per-row units, the atoms group map at a path's ends only,
// the modified index; `maxUnits`.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { type Node, printAddress, parseAddress } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { location, partOf, rangesLocation, subtractLocations, isEmptyLocation } from '@molstar/seqx-core/query';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { moleculesSpec } from '../../lang/_test/molecules.js';
import { linearsEqual } from '../linear.js';
import { resolveSpec } from '../model.js';
import { atomsSystemOf, collapseAtoms, componentsAt, expandAtoms, modifiedAt } from '../molecules.js';

const ranges = (l: { parts: readonly { system: string; ranges: Int32Array }[] }, system: string) =>
  Array.from(partOf(l as never, system)?.ranges ?? []);

describe('molecules and units as systems', () => {
  const model = resolveSpec(moleculesSpec());

  it('a molecule is its residues and atoms systems; a unit its molecule’s positions over its range, from 0 (EC-35)', () => {
    expect(model.issues.filter((i) => i.severity === 'error')).toEqual([]);
    expect(model.systems.get('/data/glycan')).toMatchObject({ kind: 'molecule', length: 3, cyclic: false });
    expect(model.systems.get('/data/atp/(atoms)')).toMatchObject({ kind: 'atoms', length: 3 });
    expect(model.molecular.molecules.get('/data/atp')).toMatchObject({ class: 'ligand', name: 'ATP' });
    expect(model.molecular.molecules.get('/data/glycan')?.identifiers[0]).toEqual({
      namespace: 'glytoucan',
      value: 'G15407YE',
    });
    const a = model.systems.get('/data/c/A')!;
    expect(a).toMatchObject({
      kind: 'unit',
      alphabet: 'protein',
      length: 20,
      parent: { system: '/data/fc', start: 0 },
    });
    expect(a.letters).toBe('MKTAYIAKQRQISFVKSHFS');
    expect(a.numberings).toEqual(['index', 'ordinal', 'auth']);
    expect(model.systems.get('/data/c/G')).toMatchObject({ kind: 'unit', length: 3 });
    expect(model.systems.get('/data/c/L/(atoms)')).toMatchObject({ kind: 'atoms', length: 3 });
    expect(model.molecular.units.get('/data/c/A')).toMatchObject({
      class: 'polymer',
      chain: 'A',
      authChain: 'H',
      operator: '1_555',
    });
    expect(model.molecular.units.get('/data/c/G')?.class).toBe('branched');
    expect(model.molecular.units.get('/data/c/L')?.group).toBe('PRD_000001');
    expect(Array.from(model.molecular.units.get('/data/c/A')!.observed!)).toEqual([2, 18]);
  });

  it('an atoms system’s address is its system’s, then (atoms), as addresses print it', () => {
    const a = atomsSystemOf('/data/c/L');
    expect(printAddress(parseAddress(a).pattern!.segments)).toBe(a);
  });

  it('a ranges annotation targets a unit; a unit’s own numbering labels its positions', () => {
    expect(model.annotations.get('/data/sites')?.target).toBe('/data/c/A');
    expect(model.numberingOf('/data/c/A', 'auth')?.numbering?.labelAt(10)).toBe('240');
    expect(model.system('/data/c/A')?.parent?.system.address).toBe('/data/fc');
  });

  it('a complex’s links are held with each end’s unit a row’s; a molecule’s on its residues', () => {
    const links = model.annotations.get('/data/c/links')!;
    expect(links.units).toEqual(['/data/c/A', '/data/c/B', '/data/c/G', '/data/c/L']);
    expect(Array.from(links.pairs!.aUnit!)).toEqual([0, 0]);
    expect(Array.from(links.pairs!.bUnit!)).toEqual([2, 1]);
    expect(links.fields).toEqual(['type', 'a_atom', 'b_atom']);
    const q = model.annotation('/data/c/links')!;
    expect(q.items().items.map((it) => [it.key, it.system, it.start])).toEqual([
      [0, '/data/c/A', 10],
      [0, '/data/c/G', 0],
      [1, '/data/c/A', 4],
      [1, '/data/c/B', 4],
    ]);
    // an end on unit B lights its pair alone
    expect(q.items(rangesLocation('/data/c/B', [4, 5])).items.map((it) => [it.key, it.system])).toEqual([
      [1, '/data/c/B'],
    ]);
    const bonds = model.annotations.get('/data/glycan/bonds')!;
    expect([bonds.target, bonds.ySystem, bonds.pairs?.n]).toEqual(['/data/glycan', '/data/glycan', 2]);
  });

  it('copies are derived: units of one complex with equal `of` and range (EC-65)', () => {
    const c = model.molecular.complexes.get('/data/c')!;
    expect([...c.copies.values()]).toEqual([['/data/c/A', '/data/c/B'], ['/data/c/G'], ['/data/c/L']]);
    expect(c.links).toEqual(['/data/c/links']);
  });

  it('the modified index: rows in order, alternatives at one position, a component over several', () => {
    const fc = model.modified('/data/fc')!;
    expect(modifiedAt(fc, 0)).toEqual(['MSE']);
    expect([modifiedAt(fc, 3), modifiedAt(fc, 5), modifiedAt(fc, 6)]).toEqual([['CRO'], ['CRO'], []]);
    const trna = model.modified('/data/trna')!;
    expect(modifiedAt(trna, 12)).toEqual(['H2U', 'U']);
    expect(componentsAt(fc, 'protein', 'MKTAYIAK', 1)).toEqual(['LYS']);
    // a unit's is its sequence's, shifted and clipped
    expect(Array.from(model.modified('/data/c/A')!.position)).toEqual([0, 3]);
  });
});

/** A sequence `s` and `t` (t = s by an offset map), units of each in one complex. */
function unitsSpec(units: readonly { id: string; of: 's' | 't'; range?: [number, number] }[], map = false): Node {
  const b = createSeqSpec();
  const s = b.data.sequence('s', { alphabet: 'protein', letters: 'ABCDEFGHIK' });
  const t = b.data.sequence('t', { alphabet: 'protein', letters: 'ABCDEFGHIK' });
  if (map) b.data.map('m', { from: s, to: t, method: 'offset', by: 0 });
  b.data.complex('c', {}, (c) => {
    for (const u of units) c.unit(u.id, { of: u.of === 's' ? s : t, ...(u.range ? { range: u.range } : {}) });
  });
  return b.build();
}

describe('paths through units', () => {
  it('a unit steps up and down freely; a path to a ranged unit is clipped by its step (EC-26)', () => {
    const m = resolveSpec(unitsSpec([{ id: 'U', of: 's', range: [2, 8] }]));
    expect(ranges(m.project(rangesLocation('/data/s', [0, 10]), '/data/c/U').location, '/data/c/U')).toEqual([0, 6]);
    expect(ranges(m.project(rangesLocation('/data/c/U', [0, 2]), '/data/s').location, '/data/s')).toEqual([2, 4]);
    expect(m.path('/data/c/U', '/data/s').path?.edges.map((e) => e.unit)).toEqual(['up']);
  });

  it('the copy rule: no path climbs to a molecule and down to a sibling; a longer legal path wins (EC-25)', () => {
    const m = resolveSpec(
      unitsSpec(
        [
          { id: 'U', of: 's', range: [2, 8] },
          { id: 'V', of: 's', range: [2, 8] },
          { id: 'X', of: 't' },
        ],
        true,
      ),
    );
    expect(m.path('/data/c/U', '/data/c/V').issue?.code).toBe('seq.resolve.path.none');
    const p = m.path('/data/c/U', '/data/c/X').path!;
    expect(p.edges.map((e) => e.unit ?? e.map)).toEqual(['up', '/data/m', 'down']);
    expect(ranges(m.project(rangesLocation('/data/c/U', [0, 1]), '/data/c/X').location, '/data/c/X')).toEqual([2, 3]);
  });

  it('a target that is a molecule reaches every unit over it (EC-48)', () => {
    const m = resolveSpec(moleculesSpec());
    const loc = rangesLocation('/data/fc', [3, 6]);
    expect(ranges(m.project(loc, '/data/c/A').location, '/data/c/A')).toEqual([3, 6]);
    expect(ranges(m.project(loc, '/data/c/B').location, '/data/c/B')).toEqual([3, 6]);
    expect(m.project(rangesLocation('/data/c/A', [3, 6]), '/data/c/B').issues[0]?.code).toBe('seq.resolve.path.none');
  });

  it('a carried search is dropped when a unit it read changes, kept otherwise', async () => {
    const { createResolver } = await import('../model.js');
    const r = createResolver();
    const a = r.resolve(unitsSpec([{ id: 'U', of: 's', range: [2, 8] }]));
    const pa = a.path('/data/c/U', '/data/s').path!;
    const b = r.resolve(
      unitsSpec([
        { id: 'U', of: 's', range: [2, 8] },
        { id: 'W', of: 't' },
      ]),
    );
    expect(b.path('/data/c/U', '/data/s').path).toBe(pa);
    const c = r.resolve(unitsSpec([{ id: 'U', of: 's', range: [1, 8] }]));
    expect(ranges(c.project(rangesLocation('/data/c/U', [0, 1]), '/data/s').location, '/data/s')).toEqual([1, 2]);
  });

  it('a complex past maxUnits is refused with an issue naming a subset', () => {
    const spec = unitsSpec(Array.from({ length: 5 }, (_, i) => ({ id: `U${i}`, of: 's' as const })));
    const m = resolveSpec(spec, { config: { maxUnits: 4 } });
    expect(m.molecular.complexes.size).toBe(0);
    expect(m.systems.has('/data/c/U0')).toBe(false);
    expect(m.issues.map((i) => i.code)).toEqual(['seq.resolve.complex.units']);
    expect(resolveSpec(spec).molecular.complexes.get('/data/c')?.units.length).toBe(5);
  });
});

describe('a reference to a unit that made no system', () => {
  /** A complex of two units of `s` (or of a source's node), and a ranges, a map and pairs naming its unit A. */
  function namingUnit(of?: string): Node {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'ABCDEFGHIK' });
    b.data.source('src', { format: 'mmcif', uri: 'x.cif' });
    b.data.complex('c', {}, (c) => {
      c.unit('A', { of: (of ?? s) as never });
      c.unit('B', { of: s });
    });
    b.data.ranges('r', { target: '/data/c/A' as never, rows: { start: [0], end: [1] } });
    b.data.map('m', { from: '/data/c/A' as never, to: s, method: 'offset', by: 0 });
    b.data.pairs('p', { a: '/data/c/A' as never, rows: { a_start: [0], a_end: [1], b_start: [2], b_end: [3] } });
    return b.build();
  }

  const codes = (m: ReturnType<typeof resolveSpec>) =>
    m.issues.filter((i) => i.code !== 'seq.resolve.source.runtime').map((i) => [i.code, i.address]);

  it('a unit of a complex refused past maxUnits: what names it is left out with an issue, never a throw', () => {
    const m = resolveSpec(namingUnit(), { config: { maxUnits: 1 } });
    expect([m.annotations.has('/data/r'), m.maps.has('/data/m'), m.annotations.has('/data/p')]).toEqual([
      false,
      false,
      false,
    ]);
    expect(codes(m)).toEqual([
      ['seq.resolve.complex.units', '/data/c'],
      ['seq.resolve.ref.system', '/data/m'],
      ['seq.resolve.ref.system', '/data/p'],
      ['seq.resolve.ref.system', '/data/r'],
    ]);
  });

  it('a unit whose `of` names nothing held (a source’s dynamic node) is no system either', () => {
    const m = resolveSpec(namingUnit('/data/src/(node :key "1")'));
    expect(m.systems.has('/data/c/A')).toBe(false);
    expect(m.annotations.has('/data/r')).toBe(false);
    expect(codes(m).filter(([c]) => c === 'seq.resolve.ref.system')).toEqual([
      ['seq.resolve.ref.system', '/data/m'],
      ['seq.resolve.ref.system', '/data/p'],
      ['seq.resolve.ref.system', '/data/r'],
    ]);
  });
});

describe('the atoms group map, at a path’s ends only', () => {
  const model = resolveSpec(moleculesSpec());

  it('atoms collapse onto their residues before a path, and expand onto atoms after it (EC-18)', () => {
    const fromAtoms = model.project(rangesLocation('/data/atp/(atoms)', [1, 2]), '/data/c/L');
    expect(ranges(fromAtoms.location, '/data/c/L')).toEqual([0, 1]);
    const toAtoms = model.project(rangesLocation('/data/atp', [0, 1]), '/data/c/L/(atoms)');
    expect(ranges(toAtoms.location, '/data/c/L/(atoms)')).toEqual([0, 3]);
  });

  it('a follow needs a fixed scale: no path enters an atoms system (EC-113)', () => {
    expect(model.path('/data/atp/(atoms)', '/data/atp').issue?.code).toBe('seq.resolve.path.none');
  });
});

// ---- properties ----

const Cases = 40;

/** Units of `s` and `t` with ranges from a small set, so copies are common. */
function randomUnits(rng: PCG) {
  const n = 1 + (rng.int() % 8);
  const choices: [number, number][] = [
    [0, 10],
    [2, 8],
    [2, 6],
  ];
  return Array.from({ length: n }, (_, i) => ({
    id: `U${i}`,
    of: rng.int() % 2 ? ('s' as const) : ('t' as const),
    range: choices[rng.int() % choices.length],
  }));
}

function shuffle<T>(rng: PCG, xs: readonly T[]): T[] {
  const out = xs.slice();
  for (let i = out.length - 1; i > 0; i--) {
    const j = rng.int() % (i + 1);
    [out[i], out[j]] = [out[j], out[i]];
  }
  return out;
}

/** Copy groups as sets of unit ids, order-free. */
function groups(spec: Node): string[] {
  const c = resolveSpec(spec).molecular.complexes.get('/data/c')!;
  return [...c.copies.values()]
    .map((g) =>
      g
        .map((a) => a.split('/').pop()!)
        .sort()
        .join(','),
    )
    .sort();
}

describe('properties of molecules and units', () => {
  it('copies are derived equal for any order of the units (EC-65)', () => {
    const rng = new PCG(65);
    for (let k = 0; k < Cases; k++) {
      const units = randomUnits(rng);
      expect(groups(unitsSpec(shuffle(rng, units)))).toEqual(groups(unitsSpec(units)));
    }
  });

  /** A molecule of 1–8 residues with 0–4 heavy atoms each. */
  function randomMolecule(rng: PCG, withAtoms: boolean) {
    const n = 1 + (rng.int() % 8);
    const residue: number[] = [];
    for (let r = 0; r < n; r++) for (let a = rng.int() % 5; a > 0; a--) residue.push(r);
    const b = createSeqSpec();
    const m = b.data.molecule('m', {
      residues: Array.from({ length: n }, () => 'NAG'),
      ...(withAtoms && residue.length
        ? { atoms: { residue, name: residue.map((_, i) => `C${i}`), element: residue.map(() => 'C') } }
        : {}),
    });
    const start = rng.int() % n;
    const end = start + 1 + (rng.int() % (n - start));
    b.data.complex('c', {}, (c) => c.unit('U', { of: m, range: [start, end] }));
    return { spec: b.build(), n, atoms: residue.length };
  }

  it('atoms collapsed then expanded cover what they started from', () => {
    const rng = new PCG(18);
    for (let k = 0; k < Cases; k++) {
      const { spec, atoms } = randomMolecule(rng, true);
      if (!atoms) continue;
      const a = resolveSpec(spec).molecular.atoms.get('/data/m/(atoms)')!;
      const s = rng.int() % atoms;
      const e = s + 1 + (rng.int() % (atoms - s));
      const back = expandAtoms(a, collapseAtoms(a, [s, e]));
      const left = subtractLocations(
        rangesLocation('/data/m/(atoms)', [s, e]),
        location({ parts: [{ system: '/data/m/(atoms)', ranges: back }] }),
      );
      expect(isEmptyLocation(left)).toBe(true);
    }
  });

  it('a path’s result is unchanged by a molecule’s atoms form', () => {
    const rng = new PCG(64);
    for (let k = 0; k < Cases; k++) {
      const seed = rng.int();
      const plain = resolveSpec(randomMolecule(new PCG(seed), false).spec);
      const withAtoms = resolveSpec(randomMolecule(new PCG(seed), true).spec);
      const a = plain.path('/data/c/U', '/data/m').path!;
      const b = withAtoms.path('/data/c/U', '/data/m').path!;
      expect(linearsEqual(a.linear, b.linear)).toBe(true);
      const n = plain.systems.get('/data/m')!.length!;
      const loc = rangesLocation('/data/m', [0, n]);
      expect(ranges(withAtoms.project(loc, '/data/c/U').location, '/data/c/U')).toEqual(
        ranges(plain.project(loc, '/data/c/U').location, '/data/c/U'),
      );
    }
  });
});

describe('200 merged copies (EC-67)', () => {
  /** A chain in `copies` units over one range, one copy group, and a glycan on each. */
  function copiesSpec(copies: number): Node {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MKTAYIAKQRQISFVKSHFSRQ'.repeat(10) });
    const g = b.data.molecule('g', { class: 'branched', residues: ['NAG', 'NAG', 'BMA'] });
    b.data.complex('c', {}, (c) => {
      for (let i = 0; i < copies; i++) c.unit(`A${i}`, { of: s, range: [10, 200], chain: `A${i}` });
      for (let i = 0; i < copies; i++) c.unit(`G${i}`, { of: g, chain: `G${i}` });
    });
    return b.build();
  }

  /** Resolves the copies and carries one hover on the sequence onto every unit of it. */
  function hoverAll(spec: Node) {
    const m = resolveSpec(spec);
    const hover = rangesLocation('/data/s', [40, 41]);
    let lit = 0;
    for (const u of m.molecular.complexes.get('/data/c')!.units)
      if (u.class === 'polymer') lit += ranges(m.project(hover, u.address).location, u.address).length / 2;
    return { m, lit };
  }

  it('two copy groups of 200; a hover on the sequence lights each copy once, on its own positions', () => {
    const { m, lit } = hoverAll(copiesSpec(200));
    expect([...m.molecular.complexes.get('/data/c')!.copies.values()].map((g) => g.length)).toEqual([200, 200]);
    expect(lit).toBe(200);
    expect(ranges(m.project(rangesLocation('/data/s', [40, 41]), '/data/c/A7').location, '/data/c/A7')).toEqual([
      30, 31,
    ]);
  });

  itWithin('resolving 200 copies and carrying a hover onto each', (within) => {
    const spec = copiesSpec(200);
    hoverAll(spec); // warm
    // measured: a few ms (an M-series laptop); 40 ms leaves room for a collection and a parallel heavy run
    expect(within(40, () => hoverAll(spec)).lit).toBe(200);
  });
});
