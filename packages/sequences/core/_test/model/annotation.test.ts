/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  type AlignmentInput,
  type Annotation,
  type AnnotationId,
  type AnnotationInput,
  Annotations,
  Conversion,
  ConversionStatus,
  CoordinateSystemKind,
  CustomData,
  type CustomProperty,
  ItemLocation,
  Loci,
  OrderedSet,
  Pairs,
  Providers,
  type SeqIndex,
  SeqLocation,
  type SeqSet,
  SeqSets,
  type SetInput,
  type UnitId,
} from '../../src/model/index.js';
import { contacts, IEP, iep1, LIGAND_UNITS, POSE_UNITS } from './annotation-fixtures.js';
import { brs1 } from './fixtures.js';

const codes = (issues: { code: string; path: string }[]) => issues.map((i) => `${i.code} ${i.path}`);
const unitId = (set: SeqSet, key: string) => set.units.key.indexOf(set.strings.find(key));
const byKey = (set: SeqSet, key: string) =>
  set.annotations.find((a) => a && set.strings.get(a.key) === key) as Annotation;

function withAnnotations(base: SetInput, annotations: AnnotationInput[]) {
  const { set, issues } = SeqSets.build({ ...base, annotations });
  return { set, issues, ann: set.annotations[annotations.length - 1] as Annotation };
}

/** Each row's [slot, start, stop] for an annotation's rows. */
const rowsOf = (a: Annotation) => Array.from({ length: a.count }, (_, r) => [a.slot[r], a.start[r], a.stop[r]]);

/** The positions of a ranges annotation's rows per target key. */
function positions(input: AnnotationInput | null): Record<string, number[]> {
  const out: Record<string, number[]> = {};
  const c = input?.columns as Record<string, unknown[]>;
  c.target.forEach((t, k) => {
    out[t as string] ??= [];
    for (let i = c.start[k] as number; i < (c.stop[k] as number); i++) out[t as string].push(i);
  });
  return out;
}

describe('annotations from columns', () => {
  it('convert labels once per target, keep row order, and group rows by id into items', () => {
    const { set, issues, ann } = withAnnotations(brs1(), [
      {
        key: 'features',
        geometry: 'ranges',
        numbering: 'auth',
        columns: {
          target: ['1brs/A', '1brs/D', '1brs/A', '1brs/A', '1brs/entity:2'],
          beg: [10, 5, 30, 12, '1'],
          end: [20, null, 40, 999, 3],
          id: ['helix', null, 'helix', null, null],
          value: [1, 2, null, 4, 5],
        },
      },
    ]);
    // Row 3 ends at an absent label (999): placed at the last position at or before it.
    expect(codes(issues)).toEqual(['model.position.form /annotations/0/columns/beg/4']);
    expect(rowsOf(ann)).toEqual([
      [0, 9, 20],
      [1, 4, 5],
      [0, 29, 40],
      [0, 11, 110],
    ]);
    expect([ann.itemCount, Array.from(ann.item), set.strings.get(ann.itemKey[0])]).toEqual([3, [0, 1, 0, 2], 'helix']);
    expect([Array.from(ann.slots.unit), Array.from(ann.slots.molecule)]).toEqual([
      [0, 3],
      [0, 1],
    ]);
    // Missing values stay missing: NaN, never zero.
    expect(Array.from(ann.sparse ?? [])).toEqual([1, 2, Number.NaN, 4]);
    const seen: number[] = [];
    Annotations.itemsOver(ann, 0, 15, 16, (row) => seen.push(row), null);
    expect(seen.sort()).toEqual([0, 3]);
    // The last covering row wins.
    expect([Annotations.valueAt(ann, 0, 15), Annotations.valueAt(ann, 0, 35), Annotations.valueAt(ann, 0, 5)]).toEqual([
      4,
      4,
      Number.NaN,
    ]);
    expect(Annotations.valueAt(ann, 1, 4)).toBe(2);
    expect(Annotations.slotOf(ann, { kind: CoordinateSystemKind.unit, id: 3, row: -1 })).toBe(1);
    expect(Annotations.slotOf(ann, { kind: CoordinateSystemKind.molecule, id: 0, row: -1 })).toBe(-1);
  });

  it('read fields: numbers, strings, dictionary-coded and JSON, missing values kept missing', () => {
    const json = [{ score: 1, tags: ['a'] }, [1, 2], null, 'x'];
    const { set, issues, ann } = withAnnotations(brs1(), [
      {
        key: 'f',
        geometry: 'ranges',
        target: '1brs/entity:1',
        links: [
          { label: 'UniProt', url: 'https://example.org/{acc}/{start}' },
          { label: 'bad', url: 'http://x' },
        ],
        columns: {
          start: [0, 5, 10, 20],
          stop: [3, 6, 12, 21],
          category: { dict: ['helix', 'strand'], codes: [0, 1, null, 0] },
          acc: ['P1', null, 'P 3', 'P4'],
          weight: [0.5, null, Number.NaN, 2],
          extra: { json },
          mixed: [1, 'a', true, null],
        },
      },
    ]);
    expect(codes(issues)).toEqual(['model.annotation.link /annotations/0/links/1']);
    expect(Array.from(ann.category ?? [], (c) => set.strings.get(c))).toEqual(['helix', 'strand', '', 'helix']);
    expect(ann.category?.[2]).toBe(-1);
    expect(ann.fields.map((f) => [set.strings.get(f.name), f.kind])).toEqual([
      ['acc', 'string'],
      ['weight', 'number'],
      ['extra', 'json'],
      ['mixed', 'json'],
    ]);
    const row = (r: number, name: string) => Annotations.rowValue(ann, set, r, name);
    expect([0, 1, 2, 3].map((r) => row(r, 'extra'))).toEqual([json[0], json[1], undefined, 'x']);
    expect(JSON.parse(JSON.stringify([0, 1, 3].map((r) => row(r, 'extra'))))).toEqual([json[0], json[1], 'x']);
    expect([0, 1, 2, 3].map((r) => [row(r, 'acc'), row(r, 'weight')])).toEqual([
      ['P1', 0.5],
      [undefined, undefined],
      ['P 3', undefined],
      ['P4', 2],
    ]);
    expect([0, 1, 2].map((r) => Annotations.linkOf(ann, set, 0, r))).toEqual([
      'https://example.org/P1/0',
      null,
      'https://example.org/P%203/10',
    ]);
  });

  it('hold dense values as one array per target, no rows; more values than positions reported', () => {
    const n = 110 + 89 + 5;
    const target = Array.from({ length: n }, (_, i) => (i < 110 ? '1brs/A' : i < 199 ? '1brs/D' : '1brs/A'));
    const value = Float64Array.from({ length: n }, (_, i) => i / 10);
    const { set, issues } = SeqSets.build({
      ...brs1(),
      annotations: [{ key: 'v', geometry: 'values', columns: { target, value } }],
    });
    const v = set.annotations[0] as Annotation;
    expect([v.count, v.dense?.map((d) => d.length), v.slots.unit.length]).toEqual([0, [110, 89], 2]);
    expect([Annotations.valueAt(v, 0, 3), Annotations.valueAt(v, 1, 0)]).toEqual([0.3, 11]);
    expect(codes(issues)).toEqual(
      Array.from({ length: 5 }, (_, k) => `model.position.range /annotations/0/columns/value/${199 + k}`),
    );
  });

  it('hold boundaries, dense values per target and a matrix with its maximum', () => {
    const { set, issues } = SeqSets.build({
      ...brs1(),
      annotations: [
        {
          key: 'b',
          geometry: 'ranges',
          target: '1brs/A',
          columns: { boundary: [4, null], start: [null, 2], stop: [null, 9] },
        },
        { key: 'v', geometry: 'values', columns: { target: ['1brs/A', '1brs/D', '1brs/A'], value: [0.5, null, 0.25] } },
        {
          key: 'm',
          geometry: 'matrix',
          target: '1brs/D',
          states: ['x', 'y'],
          columns: {
            pos: [3, 5],
            value: [
              [0.1, 0.9],
              [null, 0.2],
            ],
          },
          numbering: 'auth',
        },
      ],
    });
    expect(codes(issues)).toEqual([]);
    const [b, v, m] = set.annotations as Annotation[];
    expect(rowsOf(b)).toEqual([
      [0, 4, 4],
      [0, 2, 9],
    ]);
    const at = (i: number) => {
      const rows: number[] = [];
      Annotations.itemsOver(b, 0, i, i + 1, (r) => rows.push(r), null);
      return rows;
    };
    expect([at(3), at(4), at(5)]).toEqual([[1], [1, 0], [1]]);
    expect(v.dense?.map((d) => d.length)).toEqual([110, 89]);
    expect([0, 1, 2].map((i) => Annotations.valueAt(v, 0, i))).toEqual([0.5, 0.25, Number.NaN]);
    expect(Annotations.valueAt(v, 1, 0)).toBeNaN();
    expect(rowsOf(m)).toEqual([
      [0, 2, 3],
      [0, 4, 5],
    ]);
    expect(Array.from(m.matrix?.values ?? [], (x) => Math.round(x * 10) / 10)).toEqual([0.1, 0.9, Number.NaN, 0.2]);
    expect(Array.from(m.matrix?.max ?? [], (x) => Math.round(x * 10) / 10)).toEqual([0.9, 0.2]);
  });

  it('leave out rows they cannot place, with instructive issues', () => {
    const { issues, ann } = withAnnotations(brs1(), [
      {
        key: 'f',
        geometry: 'ranges',
        numbering: 'auth',
        columns: { target: ['1brs/A', 'nope', '1brs/A', '1brs/A'], beg: [1, 1, 500, null], end: [2, 2, 501, null] },
      },
    ]);
    expect(codes(issues)).toEqual([
      'model.annotation.target /annotations/0/columns/target/1',
      'model.position.unknown /annotations/0/columns/beg/2',
      'model.position.missing /annotations/0/columns/beg/3',
    ]);
    expect(ann.count).toBe(1);
    const mixed = SeqSets.build({
      ...brs1(),
      annotations: [{ key: 'g', geometry: 'ranges', target: '1brs/A', columns: { start: [1], end: [3] } }],
    });
    expect(codes(mixed.issues)).toEqual(['model.position.form /annotations/0/columns']);
    expect(mixed.set.annotations).toEqual([]);
  });

  it("accept an empty column set in every geometry: no rows, no issues (the model's half)", () => {
    for (const geometry of ['ranges', 'values', 'pairs', 'matrix', 'rows'] as const) {
      const states = geometry === 'matrix' ? ['a', 'b'] : undefined;
      const { issues, ann } = withAnnotations(brs1(), [{ key: 'e', geometry, target: '1brs/A', columns: {}, states }]);
      expect([geometry, issues, ann.count]).toEqual([geometry, [], 0]);
    }
  });

  it('keep ids across versions, rebuild over a changed target, and attach runtime properties', () => {
    const input: AnnotationInput = {
      key: 'f',
      geometry: 'ranges',
      target: '1brs/entity:2',
      columns: { start: [3] },
      props: { 'ex:note': 1 },
    };
    const v1 = SeqSets.build({ ...brs1(), annotations: [input] }).set;
    const a1 = v1.annotations[0] as Annotation;
    expect(a1.custom.props.get('ex:note')).toBe(1);
    const prop: CustomProperty<Annotation, number> = { name: 'ex:rows', compute: (a) => a.count * 10 };
    expect(CustomData.get(a1.custom, prop, a1, v1)).toBe(10);
    CustomData.attach(a1.custom, prop, 7);
    expect(CustomData.get(a1.custom, prop, a1, v1)).toBe(7);
    const v2 = SeqSets.withChanges(v1, { upsert: { molecules: [{ key: 'x', residues: 'ACD' }] } }).set;
    expect(v2.annotations[0]).toBe(a1);
    const mol = brs1().molecules[1];
    const v3 = SeqSets.withChanges(v2, { upsert: { molecules: [{ ...mol, label: 'Barstar (renamed)' }] } }).set;
    const a3 = v3.annotations[0] as Annotation;
    expect([a3 === a1, a3.id, a3.custom.props.has('ex:rows'), a3.custom.props.get('ex:note')]).toEqual([
      false,
      0,
      false,
      1,
    ]);
    const r4 = SeqSets.withChanges(v3, { remove: ['1brs/entity:2'] });
    expect(codes(r4.issues).filter((c) => c.startsWith('model.annotation'))).toEqual([
      'model.annotation.orphaned /remove/0',
    ]);
    expect(r4.set.annotations).toEqual([null]);
  });
});

describe('pairs', () => {
  const { set, issues } = SeqSets.build({ ...iep1(), annotations: [contacts()] });
  const pairs = set.annotations[0] as Annotation;
  const kinase = unitId(set, '1iep/A') as UnitId;

  it('read both sides in their own numberings', () => {
    expect(codes(issues)).toEqual([]);
    const total = IEP.poses.reduce((n, p) => n + p.contacts.length, 0);
    expect([pairs.count, pairs.slots.unit.length, pairs.pairs?.bTarget.unit[0]]).toEqual([total, 1, 1]);
    // Crystal T315 hbond: residue index 90, atom 12 ("N13") of the crystal STI.
    const r = IEP.poses[0].contacts.findIndex((c) => c[0] === 315);
    expect([pairs.start[r], pairs.pairs?.bStart[r], set.strings.get(pairs.category?.[r] ?? -1)]).toEqual([
      90,
      12,
      'hbond',
    ]);
  });

  it('find every pair touching a position, on either side', () => {
    const rows: [number, number][] = [];
    Pairs.touching(pairs, kinase, 90 as SeqIndex, (row, side) => rows.push([row, side]), null);
    expect(rows.length).toBe(IEP.poses.filter((p) => p.contacts.some((c) => c[0] === 315)).length);
    expect(rows.every(([, side]) => side === 0)).toBe(true);
    const crystal = unitId(set, '1iep/C') as UnitId;
    const atoms: number[] = [];
    Pairs.touching(
      pairs,
      crystal,
      31 as SeqIndex,
      (row, side) => atoms.push(pairs.start[row] + 225 + side * 1000),
      null,
    );
    expect(atoms.sort()).toEqual([1360, 1361]);
  });

  it('give the 1IEP pocket from `paired`: every residue in contact with any pose', () => {
    const pocket = Providers.paired(set, { key: 'pocket', pairs: ['contacts'], partner: LIGAND_UNITS });
    expect(pocket.issues).toEqual([]);
    const expected = [...new Set(IEP.poses.flatMap((p) => p.contacts.map((c) => c[0] - IEP.first)))].sort(
      (a, b) => a - b,
    );
    expect(positions(pocket.annotation)).toEqual({ '1iep/A': expected });
    const built = SeqSets.withChanges(set, { upsert: { annotations: [pocket.annotation as AnnotationInput] } });
    expect(codes(built.issues)).toEqual([]);
    const hbonds = Providers.paired(set, {
      key: 'hb',
      pairs: ['contacts'],
      partner: ['1iep/C'],
      categories: ['hbond'],
    });
    expect(positions(hbonds.annotation)['1iep/A'].map((i) => i + IEP.first)).toEqual([286, 315, 318, 360, 361, 381]);
  });

  it("fold contacts by copy group: the poses' group gives the union of their pockets", () => {
    const group = set.units.group[unitId(set, POSE_UNITS[0])];
    const members = Array.from(set.groups.members[group], (u) => set.strings.get(set.units.key[u]));
    expect(members).toEqual(POSE_UNITS);
    const byGroup = Providers.paired(set, { key: 'g', pairs: ['contacts'], partner: members });
    const expected = [...new Set(IEP.poses.slice(1).flatMap((p) => p.contacts.map((c) => c[0] - IEP.first)))];
    expect(positions(byGroup.annotation)['1iep/A']).toEqual(expected.sort((a, b) => a - b));
    // partnersOf aggregates per position: three of four poses touch T315, all four F382.
    const perPosition = new Map<number, number>();
    for (const m of members) {
      Pairs.partnersOf(
        pairs,
        unitId(set, m) as UnitId,
        (_u, i) => perPosition.set(i, (perPosition.get(i) ?? 0) + 1),
        null,
      );
    }
    expect([perPosition.get(90), perPosition.get(382 - 225)]).toEqual([3, 4]);
  });

  it('reach every copy from a molecule-level side, and every endpoint as an interface', () => {
    const { set: s, issues: is } = SeqSets.build({
      ...brs1(),
      annotations: [
        {
          key: 'c',
          geometry: 'pairs',
          target: '1brs/entity:1',
          columns: {
            start: [34, 58],
            b_target: ['1brs/entity:2', '1brs/entity:2'],
            b_start: [28, 37],
            direction: ['ab', null],
          },
        },
      ],
    });
    expect(codes(is)).toEqual([]);
    const c = s.annotations[0] as Annotation;
    expect(Array.from(c.pairs?.direction ?? [])).toEqual([1, 0]);
    const hits: number[] = [];
    Pairs.touching(c, 1 as UnitId, 58 as SeqIndex, (row) => hits.push(row), null);
    const partners: string[] = [];
    Pairs.partnersOf(c, 4 as UnitId, (u, i) => partners.push(`${s.strings.get(s.units.key[u])}:${i}`), null);
    expect([hits, partners]).toEqual([
      [1],
      ['1brs/A:34', '1brs/B:34', '1brs/C:34', '1brs/A:58', '1brs/B:58', '1brs/C:58'],
    ]);
    const iface = Providers.paired(s, { key: 'i', pairs: ['c'], partner: null });
    expect(positions(iface.annotation)).toEqual({ '1brs/entity:1': [34, 58], '1brs/entity:2': [28, 37] });
  });
});

describe('items at locations', () => {
  it('visit each item once at a unit position, and convert loci to items', () => {
    const { set } = SeqSets.build({
      ...brs1(),
      annotations: [
        {
          key: 'f',
          geometry: 'ranges',
          columns: {
            target: ['1brs/A', '1brs/entity:1', '1brs/A', '1brs/D'],
            start: [0, 5, 8, 3],
            stop: [10, 20, 9, 4],
            id: ['a', 'b', 'a', 'd'],
          },
        },
      ],
    });
    const items: number[] = [];
    const out = ItemLocation.create(set);
    const loc = SeqLocation.set(SeqLocation.create(set), 0, 8);
    expect(SeqLocation.itemsAt(loc, 0 as AnnotationId, out, (l) => items.push(l.item), null)).toBe(2);
    expect(items).toEqual([0, 1]);
    expect(SeqLocation.itemsAt(SeqLocation.set(loc, 1, 8), 0 as AnnotationId, out, () => {}, null)).toBe(1);
    const c = Conversion.create();
    const loci = {
      kind: 'seq-loci' as const,
      set,
      elements: [{ unit: 0 as UnitId, indices: OrderedSet.ofRange(0 as SeqIndex, 30 as SeqIndex) }],
    };
    const il = Loci.toItemLoci(loci, 0 as AnnotationId, c);
    expect([
      il.elements.map((e) => Array.from({ length: OrderedSet.size(e.items) }, (_, i) => OrderedSet.getAt(e.items, i))),
      c.status,
      c.dropped,
      c.multiplied,
    ]).toEqual([[[0, 1]], ConversionStatus.partial, 10, 5]);
  });
});

const ALIGNMENT: AlignmentInput = {
  key: 'msa',
  flavor: 'aligned',
  alphabet: 'protein',
  gapped: ['AL-W', 'AK-W', 'CKAW', '-L-W', 'X--B'],
};

describe('column statistics, one fixture column per rule', () => {
  const { set } = SeqSets.build({
    molecules: [],
    alignments: [
      ALIGNMENT,
      { key: 'rna', flavor: 'aligned', alphabet: 'rna', gapped: ['U', 'U', 'A', 'N'] },
      { key: 'dna', flavor: 'aligned', alphabet: 'dna', gapped: ['T', 'T', 'A', 'N'] },
    ],
  });
  const stats = (measure: 'conservation' | 'consensus' | 'identity' | 'gap_fraction', more = {}) => {
    const r = Providers.columnStats(set, { key: 's', alignment: 'msa', measure, ...more });
    expect(r.issues).toEqual([]);
    return r.annotation as AnnotationInput;
  };
  const values = (a: AnnotationInput) => Array.from(a.columns.value as ArrayLike<number>);
  const H = (...p: number[]) => -p.reduce((h, x) => h + x * Math.log2(x), 0);

  it('conservation: k + 1 symbols with gaps, X and B left out of every count', () => {
    const [c0, , , c3] = values(stats('conservation'));
    expect(c0).toBeCloseTo(1 - H(0.5, 0.25, 0.25) / Math.log2(21), 12);
    expect(c3).toBe(1);
    const ignore = values(stats('conservation', { gaps: 'ignore' }));
    expect(ignore[0]).toBeCloseTo(1 - H(2 / 3, 1 / 3) / Math.log2(20), 12);
    expect(ignore[2]).toBe(1);
  });

  it('identity: NaN at a reference gap; consensus: ties to the lower code, `value` its share', () => {
    expect(values(stats('identity'))).toEqual([0.5, 0.4, Number.NaN, 1]);
    expect(values(stats('identity', { gaps: 'ignore', referenceRow: 2 }))[2]).toBe(1);
    const consensus = stats('consensus');
    expect([consensus.geometry, consensus.columns.label, values(consensus)]).toEqual([
      'ranges',
      ['A', 'K', 'A', 'W'],
      [0.5, 0.4, 0.2, 1],
    ]);
    expect(values(stats('consensus', { gaps: 'ignore' }))).toEqual([2 / 3, 0.5, 1, 1]);
  });

  it('gap fraction: over all rows or the subset, whatever `gaps` says', () => {
    expect(values(stats('gap_fraction', { gaps: 'ignore' }))).toEqual([0.2, 0.2, 0.8, 0]);
    expect(values(stats('gap_fraction', { rowSubset: [2, 3] }))).toEqual([0.5, 0, 0.5, 0]);
  });

  it('nucleotides: T counted as U (one symbol, DNA writes T), N left out; states in alphabet order', () => {
    const one = (alignment: string) => {
      const c = Providers.columnStats(set, { key: 'c', alignment, measure: 'consensus' }).annotation;
      const p = Providers.profile(set, { key: 'p', alignment }).annotation as AnnotationInput;
      return [c?.columns.label, values(c as AnnotationInput), p.states, p.columns.value];
    };
    expect(one('rna')).toEqual([['U'], [2 / 3], ['A', 'C', 'G', 'U'], [[1 / 3, 0, 0, 2 / 3]]]);
    expect(one('dna')).toEqual([['T'], [2 / 3], ['A', 'C', 'G', 'T'], [[1 / 3, 0, 0, 2 / 3]]]);
    const p = Providers.profile(set, { key: 'p', alignment: 'rna' }).annotation as AnnotationInput;
    const built = SeqSets.withChanges(set, { upsert: { annotations: [p, stats('conservation')] } });
    expect(codes(built.issues)).toEqual([]);
    const [profile, cons] = built.set.annotations as Annotation[];
    expect([profile.geometry, profile.count, cons.dense?.[0].length]).toEqual(['matrix', 1, 4]);
    expect(Annotations.valueAt(cons, 0, 3)).toBe(1);
  });

  it('records per row: `rows` annotations, rows out of range left out', () => {
    const input: AnnotationInput = {
      key: 'species',
      geometry: 'rows',
      target: 'msa',
      columns: { row: [3, 0, 7], species: ['mouse', 'human', 'fish'] },
    };
    const { set: s, issues } = SeqSets.withChanges(set, { upsert: { annotations: [input] } });
    expect(codes(issues)).toEqual(['model.annotation.target /upsert/annotations/0/columns/target/2']);
    const a = byKey(s, 'species');
    expect([Array.from(a.target.row), Array.from(a.slots.alignment), Annotations.rowValue(a, s, 1, 'species')]).toEqual(
      [[3, 0], [0, 0], 'human'],
    );
  });

  it('key tables: a `rows` record named by `state` or `node` has no target, and the name is a string field', () => {
    const input: AnnotationInput = {
      key: 'samples',
      geometry: 'rows',
      columns: { state: ['S1', 'S2', 'S3'], population: ['x', 'y', 'x'], age: [30, 20, 10] },
    };
    const { set: s, issues } = SeqSets.withChanges(set, { upsert: { annotations: [input] } });
    expect(issues).toEqual([]);
    const a = byKey(s, 'samples');
    expect([a.count, Array.from(a.target.row), Array.from(a.target.alignment)]).toEqual([
      3,
      [-1, -1, -1],
      [-1, -1, -1],
    ]);
    expect([Annotations.rowValue(a, s, 2, 'state'), Annotations.rowValue(a, s, 0, 'population')]).toEqual(['S3', 'x']);
    expect(Annotations.rowValue(a, s, 1, 'age')).toBe(20);
    // Without a row, a state or a node, a `rows` record still needs its row.
    const bad = SeqSets.withChanges(set, {
      upsert: { annotations: [{ key: 'x', geometry: 'rows', target: 'msa', columns: { age: [1] } }] },
    });
    expect(codes(bad.issues)).toEqual(['model.annotation.row /upsert/annotations/0/columns/row/0']);
  });

  it('are unavailable for a generic alphabet and reject unknown rows', () => {
    const g = SeqSets.build({ molecules: [], alignments: [{ ...ALIGNMENT, alphabet: 'generic' }] }).set;
    expect(codes(Providers.profile(g, { key: 'p', alignment: 'msa' }).issues)).toEqual([
      'model.provider.alphabet /alignment',
    ]);
    expect(
      codes(Providers.columnStats(set, { key: 's', alignment: 'msa', measure: 'identity', rowSubset: [9] }).issues),
    ).toEqual(['model.provider.input /row_subset/0']);
  });
});

describe('coverage', () => {
  it('gives observed, unobserved and outside per unit, and conflicts against a reference', () => {
    const BARSTAR = brs1().molecules[1].residues as string;
    const reference = `M${BARSTAR.slice(0, 10)}G${BARSTAR.slice(11)}`;
    const { set, issues } = SeqSets.build({
      ...brs1(),
      molecules: [...brs1().molecules, { key: 'P11540', residues: reference }],
      units: [...(brs1().units ?? []), { key: 'part', molecule: '1brs/entity:2', range: [5, 50] }],
      mappings: [{ key: 'm', from: 'P11540', to: '1brs/entity:2', method: 'segments', segments: [[1, 0, 89]] }],
    });
    expect(codes(issues)).toEqual([]);
    const own = Providers.coverage(set, { key: 'cov', of: ['1brs/D', 'part'] }).annotation as AnnotationInput;
    expect(own.columns).toEqual({
      target: ['1brs/D', '1brs/D', '1brs/D', 'part', 'part', 'part'],
      start: [0, 63, 65, 0, 5, 50],
      stop: [63, 65, 89, 5, 50, 89],
      category: ['observed', 'unobserved', 'observed', 'outside', 'observed', 'outside'],
      id: ['1brs/D', '1brs/D', '1brs/D', 'part', 'part', 'part'],
    });
    const ref = Providers.coverage(set, { key: 'cov', of: ['1brs/D'], reference: 'P11540' })
      .annotation as AnnotationInput;
    expect([ref.columns.start, ref.columns.category]).toEqual([
      [0, 1, 11, 12, 64, 66],
      ['outside', 'observed', 'conflict', 'observed', 'unobserved', 'observed'],
    ]);
    const built = SeqSets.withChanges(set, { upsert: { annotations: [ref] } });
    expect([codes(built.issues), byKey(built.set, 'cov').slots.molecule[0]]).toEqual([[], 2]);
  });
});
