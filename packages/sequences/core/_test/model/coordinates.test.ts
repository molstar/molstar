/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Coordinate maps (seqspec:model#coordinate-maps, seqspec:model#codon-maps; FR1): what a codon mapping holds
// over a long molecule (G1-1), a source used twice (linear di-ubiquitin), every index of the codon fixtures through
// `Coords` against the residues' own spans, boundaries on both strands and composition clipped to a domain.

import { describe, expect, it } from 'vitest';
import {
  type CoordinateMap,
  CoordinateSystemKind,
  Coords,
  type MapHit,
  type Mapping,
  type MappingInput,
  Mappings,
  type PathStep,
  Paths,
  type SeqSet,
  SeqSets,
} from '../../src/model/index.js';
import { CDS, GENE, PROTEIN } from './codon-fixtures.js';
import { gene as dmd } from './support/generators.js';

const COMPLEMENT: Record<string, string> = { A: 'T', T: 'A', C: 'G', G: 'C' };
const revcomp = (s: string) =>
  [...s]
    .reverse()
    .map((c) => COMPLEMENT[c] ?? c)
    .join('');
const system = (set: SeqSet, key: string) => {
  const s = { kind: -1, id: -1, row: -1 };
  if (!Mappings.systemOf(set, key, s)) throw new Error(`no system ${key}`);
  return s;
};
const runsTo = (m: CoordinateMap, s: number, t: number) => {
  const out: number[][] = [];
  Coords.forEachTo(m, s, t, (a, b, acc) => acc.push([a, b]), out);
  return out;
};
const runsFrom = (m: CoordinateMap, s: number, t: number) => {
  const out: number[][] = [];
  Coords.forEachFrom(m, s, t, (a, b, acc) => acc.push([a, b]), out);
  return out;
};
const holders = (m: CoordinateMap, to: number) => {
  const out: number[][] = [];
  const h: MapHit = { element: -1, phase: -1 };
  const n = Coords.forEachFromAt(m, to, h, (x, acc) => acc.push([x.element, x.phase]), out);
  expect(n).toBe(out.length);
  return out;
};
/** Every typed array a mapping holds, in bytes. */
const bytes = (m: Mapping) =>
  [m.forward, m.status, m.spans, m.map.at, m.map.length, m.map.toLo, m.map.dir, m.map.byTo, m.map.byToMax].reduce(
    (n, a) => n + (a?.byteLength ?? 0),
    0,
  );

function coding(cds: [number, number][], strand: '+' | '-', dna: string, residues = PROTEIN) {
  const { set, issues } = SeqSets.build({
    molecules: [
      { key: 'gene', alphabet: 'dna', residues: dna, ownUnit: true },
      { key: 'protein', alphabet: 'protein', residues, ownUnit: true },
    ],
    mappings: [{ key: 'c', from: 'protein', to: 'gene', method: 'codon', cds, strand }],
  });
  expect(issues).toEqual([]);
  return { set, m: set.mappings[0] as Mapping };
}
const mirrored = CDS.map(([s, t]) => [GENE.length - t, GENE.length - s] as [number, number]);

describe('coordinate maps', () => {
  it('holds O(CDS) bytes for a codon mapping over a 1M-nucleotide molecule (G1-1)', () => {
    // A 3 kb CDS in three exons in the middle of 1M nucleotides: nothing is held per nucleotide.
    const protein = 'M'.repeat(999);
    const codon = 'ATG'.repeat(1000);
    const dna = `${'C'.repeat(400_000)}${codon.slice(0, 1000)}${'G'.repeat(100_000)}${codon.slice(1000, 2000)}${'G'.repeat(100_000)}${codon.slice(2000)}`;
    const full = dna + 'C'.repeat(1_000_000 - dna.length);
    const cds: [number, number][] = [
      [400_000, 401_000],
      [501_000, 502_000],
      [602_000, 603_000],
    ];
    const { m } = coding(cds, '+', full, protein);
    expect([m.map.count, m.map.toStart, m.map.toStop]).toEqual([3, 400_000, 603_000 - 3]);
    // forward (4 per residue), status (1), spans (12) and the segment tables: about 17 per residue, under 20 KB.
    expect(bytes(m)).toBeLessThan(20 * 1024);
    expect(holders(m.map, 501_000)).toEqual([[333, 1]]);
    expect(Coords.toAt(m.map, 998, 2)).toBe(602_996);
  });

  it("gives every index of the codon fixtures through Coords as the residues' spans say, on both strands", () => {
    // Case 14 (a codon split by the intron), its mirror on strand -, and CDSs at either end of the molecule.
    const ends = `${GENE.slice(12, 46)}${GENE.slice(76, 105)}`;
    const fixtures = [
      coding(CDS, '+', GENE),
      coding(mirrored, '-', revcomp(GENE)),
      coding([[0, 63]], '+', ends),
      coding([[0, 63]], '-', revcomp(ends)),
      coding([[57, 120]], '+', `${'A'.repeat(57)}${ends}`),
    ];
    for (const { set, m } of fixtures) {
      const n = set.molecules[0]?.length ?? 0;
      const spans = m.spans as Int32Array;
      const expected = Array.from({ length: n }, () => [] as number[][]);
      for (let q = 0; q < spans.length; q++) if (spans[q] >= 0) expected[spans[q]].push([Math.floor(q / 3), q % 3]);
      for (let j = 0; j < n; j++) expect(holders(m.map, j), `nucleotide ${j}`).toEqual(expected[j]);
      for (let q = 0; q < spans.length; q++) expect(Coords.toAt(m.map, Math.floor(q / 3), q % 3)).toBe(spans[q]);
      // Each residue's runs are its codon's positions, in to order.
      for (let r = 0; r < PROTEIN.length; r++) {
        const own = Array.from(spans.subarray(3 * r, 3 * r + 3)).sort((a, b) => a - b);
        expect(runsTo(m.map, r, r + 1).flatMap(([s, t]) => Array.from({ length: t - s }, (_, k) => s + k))).toEqual(
          own,
        );
      }
      expect(m.map.count).toBeLessThanOrEqual(2);
    }
  });

  it('walks a residue range as runs in to order, an intron splitting them, on both strands', () => {
    const plus = coding(CDS, '+', GENE).m.map;
    expect(runsTo(plus, 4, 14)).toEqual([
      [24, 46],
      [76, 84],
    ]);
    const minus = coding(mirrored, '-', revcomp(GENE)).m.map;
    expect(runsTo(minus, 4, 14)).toEqual([
      [36, 44],
      [74, 96],
    ]);
    // Nucleotides back to the residues holding them: one run when they touch across the intron.
    expect(runsFrom(plus, 40, 80)).toEqual([[9, 13]]);
    expect(runsFrom(minus, 40, 80)).toEqual([[9, 13]]);
    expect([runsFrom(plus, 0, 12), runsFrom(plus, 102, 120)]).toEqual([[], []]);
  });

  it('maps boundaries on both strands, at the end and across an intron', () => {
    const at = (m: CoordinateMap, k: number) => runsTo(m, k, k);
    // Exons end on a codon: residue 11 starts the second exon; the stop codon (103–105) pairs nothing.
    const whole: [number, number][] = [
      [12, 45],
      [76, 106],
    ];
    const plus = coding(whole, '+', GENE).m.map;
    const minus = coding(
      whole.map(([s, t]) => [120 - t, 120 - s] as [number, number]),
      '-',
      revcomp(GENE),
    ).m.map;
    // Before residue 0, before residue 11 (the next exon's first nucleotide), at the end (after residue 19's last).
    expect([at(plus, 0), at(plus, 11), at(plus, 20)]).toEqual([[[12, 12]], [[76, 76]], [[103, 103]]]);
    // Mirrored: the boundary after the first nucleotide by index, and the end before the last.
    expect([at(minus, 0), at(minus, 11), at(minus, 20)]).toEqual([[[108, 108]], [[44, 44]], [[17, 17]]]);
    // A span-1 map: the end of the domain after the last position of the last element.
    const row = Coords.create(sys(0), sys(1), 1, [0, 10, 5, 1]);
    expect([at(row, 0), at(row, 4), at(row, 10), at(row, 11)]).toEqual([[[5, 5]], [[9, 9]], [[15, 15]], []]);
  });

  it('holds linear di-ubiquitin as one mapping with two segments onto ubiquitin: not ambiguous, both copies found', () => {
    const UB = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG';
    const { set, issues } = SeqSets.build({
      molecules: [
        { key: 'ub', residues: UB, ownUnit: true },
        { key: 'diub', residues: UB + UB, ownUnit: true },
      ],
      mappings: [
        {
          key: 'diub.ub',
          from: 'diub',
          to: 'ub',
          method: 'segments',
          segments: [
            [0, 0, 76],
            [76, 0, 76],
          ],
        },
      ],
    });
    expect(issues).toEqual([]);
    const m = set.mappings[0] as Mapping;
    // K48 (index 47) is held by diub 47 and 123.
    expect(holders(m.map, 47)).toEqual([
      [47, 0],
      [123, 0],
    ]);
    expect(runsFrom(m.map, 40, 50)).toEqual([
      [40, 50],
      [116, 126],
    ]);
    expect([m.forward[47], m.forward[123], Array.from(m.status).every((s) => s === 0)]).toEqual([47, 47, true]);
    // Ub's own unit to diub's: one path, the copies' side as the map's `from`, ub 47 lighting both copies.
    const steps: PathStep[] = [];
    expect(Mappings.path(set, system(set, 'ub/unit'), system(set, 'diub/unit'), steps)).toBe(3);
    const p = Mappings.compose(set, system(set, 'ub/unit'), steps);
    expect([Coords.sameSystem(p.map.from, system(set, 'diub/unit')), Paths.at(p, 47), p.map.order]).toEqual([
      true,
      47,
      0,
    ]);
    expect(runsFrom(p.map, 47, 48)).toEqual([
      [47, 48],
      [123, 124],
    ]);
    // Written from the source to the construct, the from ranges would overlap: the builder refuses it.
    const back: MappingInput = {
      key: 'x',
      from: 'ub',
      to: 'diub',
      method: 'segments',
      segments: [
        [0, 0, 76],
        [0, 76, 76],
      ],
    };
    expect(SeqSets.withChanges(set, { upsert: { mappings: [back] } }).issues.map((i) => i.path)).toEqual([
      '/upsert/mappings/0/segments/1',
    ]);
  });

  it('composes through a ranged unit, clipped to its domain at either end', () => {
    const P = 'ACDEFGHIKLMNPQRSTVWY'.repeat(6);
    const { set } = SeqSets.build({
      molecules: [
        { key: 'P', residues: P, ownUnit: true },
        { key: 'entity', residues: P },
      ],
      units: [{ key: 'A', molecule: 'entity', range: [10, 100] }],
      mappings: [{ key: 'sifts', from: 'P', to: 'A', method: 'identity' }],
    });
    const steps: PathStep[] = [];
    expect(Mappings.path(set, system(set, 'P/unit'), system(set, 'A'), steps)).toBe(2);
    const p = Mappings.compose(set, system(set, 'P/unit'), steps);
    expect([Paths.at(p, 9), Paths.at(p, 10), Paths.at(p, 99), Paths.at(p, 100), p.map.toStart, p.map.toStop]).toEqual([
      -1, 10, 99, -1, 10, 100,
    ]);
    expect(runsTo(p.map, 0, 120)).toEqual([[10, 100]]);
    // And from the unit: its domain clips the start.
    const out: PathStep[] = [];
    Mappings.path(set, system(set, 'A'), system(set, 'P/unit'), out);
    const q = Mappings.compose(set, system(set, 'A'), out);
    expect([Paths.at(q, 5), Paths.at(q, 50), q.map.fromStart, q.map.fromStop]).toEqual([-1, 50, 10, 100]);
  });

  it("builds maps from each method's own runs, never per position", () => {
    const g = dmd();
    const { set, issues } = SeqSets.build({
      molecules: [
        { key: 'dmd', alphabet: 'dna', residues: g.dna },
        { key: 'p', alphabet: 'protein', residues: g.protein, ownUnit: true },
        { key: 'q', alphabet: 'protein', residues: g.protein, ownUnit: true },
      ],
      mappings: [
        { key: 'c', from: 'p', to: 'dmd', method: 'codon', cds: g.exons.map(([s, t]) => [s, t]) },
        { key: 'aln', from: 'p', to: 'q', method: 'pairwise', cigar: `10M5I5D${g.protein.length - 15}M` },
      ],
    });
    expect(issues).toEqual([]);
    const [c, aln] = set.mappings as Mapping[];
    expect([c.map.count, aln.map.count, bytes(c) < 30 * g.protein.length]).toEqual([79, 2, true]);
    expect(Array.from(aln.map.at)).toEqual([0, 15]);
    expect(c.status.every((s) => s === 0)).toBe(true);
  });
});

function sys(id: number) {
  return { kind: CoordinateSystemKind.molecule, id, row: -1 };
}

describe('forEachTo over a map whose segments overlap in to (CM-08)', () => {
  it('gives runs in to order, merged, as its contract says', () => {
    const A = { kind: 0, id: 0, row: -1 };
    const B = { kind: 0, id: 1, row: -1 };
    const at = (m: ReturnType<typeof Coords.create>, s: number, t: number) => {
      const runs: number[][] = [];
      Coords.forEachTo(m, s, t, (x, y) => runs.push([x, y]), null);
      return runs;
    };
    // A tandem: construct 0..100 and 100..200 both onto GFP 0..100.
    const tandem = Coords.create(A, B, 1, [0, 100, 0, 1, 100, 100, 0, 1]);
    expect(at(tandem, 90, 110)).toEqual([
      [0, 10],
      [90, 100],
    ]);
    // Segments (0,100→0), (100,100→50), (200,10→10): [90, 210) clips to [90,100), [50,150) and [10,20); the
    // first lies inside the second (it came first, unsorted, before).
    const three = Coords.create(A, B, 1, [0, 100, 0, 1, 100, 100, 50, 1, 200, 10, 10, 1]);
    expect(at(three, 90, 210)).toEqual([
      [10, 20],
      [50, 150],
    ]);
  });
});
