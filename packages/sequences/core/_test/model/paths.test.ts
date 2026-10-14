/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Sparse composed paths (ST7, seqspec:model#paths): a path holds only what it pairs on its `from` side, so a protein's
// path through a chromosome holds its CDS; `Paths` reads it per element and by ranges either way, as the dense arrays
// it replaced did (checked against the map itself).

import { describe, expect, it } from 'vitest';
import {
  type ComposedPath,
  type CoordinateSystem,
  Coords,
  type MappingInput,
  Mappings,
  type PathStep,
  Paths,
  type SeqSet,
  SeqSets,
} from '../../src/model/index.js';
import { CDS, GENE, PROTEIN } from './codon-fixtures.js';

const CHR = 200_000_000;
/** A protein coded by three exons of a lazy 200 Mb chromosome, on either strand. */
function genome(strand: '+' | '-' = '+') {
  const cds: [number, number][] = [
    [120_000_000, 120_000_030],
    [120_004_000, 120_004_031],
    [120_090_000, 120_090_002],
  ];
  const input: MappingInput = {
    key: 'p.codon',
    from: 'p',
    to: 'chr',
    method: 'codon',
    cds: strand === '+' ? cds : cds.slice().reverse(),
    strand,
  };
  return SeqSets.build({
    molecules: [
      { key: 'chr', alphabet: 'dna', lazy: { length: CHR }, ownUnit: true },
      { key: 'p', alphabet: 'protein', residues: 'M'.repeat(21), ownUnit: true },
    ],
    mappings: [input],
  }).set;
}
const space = (set: SeqSet, key: string): CoordinateSystem => {
  const s = { kind: -1, id: -1, row: -1 };
  if (!Mappings.systemOf(set, key, s)) throw new Error(`no system ${key}`);
  return s;
};
function compose(set: SeqSet, from: string, to: string): ComposedPath {
  const steps: PathStep[] = [];
  if (Mappings.path(set, space(set, from), space(set, to), steps) < 0) throw new Error(`no path ${from} → ${to}`);
  return Mappings.compose(set, space(set, from), steps);
}
const runsTo = (p: ComposedPath, a: number, e: number) => {
  const out: number[][] = [];
  Paths.forEachTo(p, a, e, (s, t, acc) => acc.push([s, t]), out);
  return out;
};
const runsFrom = (p: ComposedPath, a: number, e: number) => {
  const out: number[][] = [];
  Paths.forEachFrom(p, a, e, (s, t, acc) => acc.push([s, t]), out);
  return out;
};
/** What the dense arrays held per from element, read from the map: its first to position, or its span's j-th. */
function reference(p: ComposedPath, i: number, j: number): number {
  const m = p.map;
  if (Coords.sameSystem(m.from, p.from)) return Coords.toAt(m, i, j);
  if (j > 0) return -1;
  let first = -1;
  const hit = { element: -1, phase: -1 };
  Coords.forEachFromAt(
    m,
    i,
    hit,
    (h) => {
      if (first < 0 || h.element < first) first = h.element;
    },
    null,
  );
  return first;
}

describe('sparse composed paths', () => {
  it("a protein's path through a 200 Mb chromosome holds its CDS, not the chromosome", () => {
    const set = genome();
    const down = compose(set, 'chr/unit', 'p/unit');
    expect([down.length, down.heldTo.length, down.heldSpans, Array.from(down.runs)]).toEqual([
      CHR,
      63,
      null,
      [120_000_000, 120_000_030, 120_004_000, 120_004_031, 120_090_000, 120_090_002],
    ]);
    // Codon 10 is split by the first intron: 120,000,030 is intronic.
    const at = [120_000_000, 120_000_029, 120_000_030, 120_004_000, 120_004_002, 119_999_999, CHR - 1, -1];
    expect(at.map((i) => Paths.at(down, i))).toEqual([0, 9, -1, 10, 10, -1, -1, -1]);
    const up = compose(set, 'p/unit', 'chr/unit');
    expect([up.length, up.heldTo.length, up.heldSpans?.length, up.span]).toEqual([21, 21, 63, 3]);
    expect([0, 1, 2].map((j) => Paths.spanAt(up, 9, j))).toEqual([120_000_027, 120_000_028, 120_000_029]);
    expect([0, 1, 2].map((j) => Paths.spanAt(up, 10, j))).toEqual([120_004_000, 120_004_001, 120_004_002]);
  });

  it("reads ranges either way: a residue's codon split across an intron, the exons of a stretch of residues", () => {
    const down = compose(genome(), 'chr/unit', 'p/unit');
    // To (residues) → from (nucleotides): the request's direction.
    expect(runsFrom(down, 9, 11)).toEqual([
      [120_000_027, 120_000_030],
      [120_004_000, 120_004_003],
    ]);
    expect(runsFrom(down, 0, 21)).toEqual([
      [120_000_000, 120_000_030],
      [120_004_000, 120_004_031],
      [120_090_000, 120_090_002],
    ]);
    // From (nucleotides) → to (residues): an item's image; an intronic one reaches nothing.
    expect([runsTo(down, 120_000_020, 120_004_010), runsTo(down, 120_001_000, 120_002_000)]).toEqual([[[6, 14]], []]);
    expect([Paths.ratio(down), Paths.ratio(compose(genome(), 'p/unit', 'chr/unit'))]).toEqual([3, 1 / 3]);
    // On the minus strand residue 0 is at the top of the last exon: runs still come in to order.
    const minus = compose(genome('-'), 'chr/unit', 'p/unit');
    expect([Paths.at(minus, 120_090_001), Paths.at(minus, 120_000_000)]).toEqual([0, 20]);
    expect([runsFrom(minus, 0, 1), runsFrom(minus, 10, 11)]).toEqual([
      [
        [120_004_030, 120_004_031],
        [120_090_000, 120_090_002],
      ],
      [[120_004_000, 120_004_003]],
    ]);
  });

  it("reads a residue's bases back through the path, in reading order (fromAt)", () => {
    const down = compose(genome(), 'chr/unit', 'p/unit');
    expect([0, 1, 2].map((k) => Paths.fromAt(down, 10, k))).toEqual([120_004_000, 120_004_001, 120_004_002]);
    const minus = compose(genome('-'), 'chr/unit', 'p/unit');
    // Residue 0 on the minus strand: two bases at the top of the last exon, the third at the top of the middle one.
    expect([0, 1, 2, 3].map((k) => Paths.fromAt(minus, 0, k))).toEqual([120_090_001, 120_090_000, 120_004_030, -1]);
    expect([Paths.fromAt(down, -1, 0), Paths.fromAt(down, 21, 0)]).toEqual([-1, -1]);
    // From the protein side (the map read from its own from): a base's residue, k 0 only.
    const up = compose(genome(), 'p/unit', 'chr/unit');
    expect([
      Paths.fromAt(up, 120_004_001, 0),
      Paths.fromAt(up, 120_004_001, 1),
      Paths.fromAt(up, 120_000_030, 0),
    ]).toEqual([10, -1, -1]);
  });

  it('agrees with the dense arrays it replaced, element by element, both ways and on both strands', () => {
    const gene = SeqSets.build({
      molecules: [
        { key: 'gene', alphabet: 'dna', residues: GENE, ownUnit: true },
        { key: 'protein', alphabet: 'protein', residues: PROTEIN, ownUnit: true },
      ],
      mappings: [{ key: 'c', from: 'protein', to: 'gene', method: 'codon', cds: CDS }],
    }).set;
    const cases: [SeqSet, string, string][] = [
      [gene, 'gene/unit', 'protein/unit'],
      [gene, 'protein/unit', 'gene/unit'],
      [gene, 'gene', 'protein'],
      [genome('-'), 'p/unit', 'chr/unit'],
    ];
    for (const [set, a, b] of cases) {
      const p = compose(set, a, b);
      const lo = p.runs.length ? p.runs[0] - 2 : 0;
      const hi = p.runs.length ? p.runs[p.runs.length - 1] + 2 : p.length;
      for (let i = Math.max(-1, lo); i <= Math.min(p.length, hi); i++)
        for (let j = 0; j < p.span; j++)
          expect([a, i, j, Paths.spanAt(p, i, j)]).toEqual([a, i, j, reference(p, i, j)]);
    }
  });
});
