/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqQL's runtime fixture (QR1): 1BRS with a molecule-level domain, an antibody chain with insertion codes and
// unlabeled positions (a ranged unit and a full one), a molecule with no unit, a small alignment, and an annotation of
// every geometry a query reads.

import { readFileSync } from 'node:fs';
import type { SetInput } from '../../src/model/index.js';
import { CDS, GENE, PROTEIN } from './codon-fixtures.js';
import { brs1 } from './fixtures.js';
import { prng } from './support/generators.js';

/** Kabat-like labels: insertion codes (3A, 3B), decimals (9.1, 9.2) and gaps (null), on 20 positions. */
export const KABAT = [1, 2, 3, '3A', '3B', 4, null, null, 7, 8, 9, '9.1', '9.2', 10, 11, 12, null, 14, 15, 16];

export function queryFixture(): SetInput {
  const b = brs1();
  return {
    molecules: [
      ...b.molecules,
      {
        key: 'ab',
        residues: 'QVQLVESGGGLVQPGGSLRL',
        numberings: [
          { name: 'kabat', method: 'labels', labels: KABAT, isDefault: true },
          { name: 'seq', method: 'offset', first: -5, skipZero: true },
        ],
      },
      { key: 'lonely', residues: 'ACDEFGHIKL' },
    ],
    units: [...(b.units ?? []), { key: 'ab/H', molecule: 'ab', range: [2, 18] }, { key: 'ab/L', molecule: 'ab' }],
    complexes: b.complexes,
    alignments: [
      { key: 'msa', gapped: ['AC-DE-FG', 'A--DEKFG', 'ACCD--FG'], flavor: 'aligned', labels: ['x/11-16', 'y', 'z'] },
    ],
    annotations: [
      {
        key: 'sites',
        geometry: 'ranges',
        target: '1brs/A',
        columns: {
          start: [5, 30, 60, 8],
          stop: [10, 40, 70, 12],
          label: ['a', 'b', 'c', 'd'],
          id: ['s1', 's2', 's3', 's4'],
        },
      },
      { key: 'dom', geometry: 'ranges', target: '1brs/entity:1', columns: { start: [20, 90], stop: [50, 105] } },
      { key: 'lone', geometry: 'ranges', target: 'lonely', columns: { start: [2], stop: [6] } },
      { key: 'cols', geometry: 'ranges', target: 'msa', columns: { start: [1, 5], stop: [3, 8] } },
      { key: 'rowsite', geometry: 'ranges', target: 'msa/1', columns: { start: [0], stop: [2] } },
      {
        key: 'multi',
        geometry: 'ranges',
        columns: { target: ['ab/L', '1brs/D', 'ab/H'], start: [15, 0, 0], stop: [19, 4, 3], id: ['m', 'm', 'n'] },
      },
      {
        key: 'cons',
        geometry: 'values',
        target: 'msa',
        columns: { value: [0.9, 0.2, 0.85, 0.95, Number.NaN, 0.5, 0.81, 1] },
      },
      {
        key: 'plddt',
        geometry: 'values',
        target: '1brs/A',
        columns: { value: Float64Array.from({ length: 110 }, (_, i) => 50 + 45 * Math.sin(i / 7)) },
      },
      {
        key: 'sv',
        geometry: 'values',
        target: '1brs/D',
        columns: { start: [0, 10, 12, 40], stop: [5, 14, 20, 41], value: [0.9, 0.7, 0.95, 0.1] },
      },
      {
        key: 'contacts',
        geometry: 'pairs',
        target: '1brs/A',
        columns: {
          start: [33, 58, 101],
          b_target: ['1brs/D', '1brs/D', '1brs/E'],
          b_start: [28, 37, null],
          b_stop: [30, 38, null],
        },
      },
      { key: 'species', geometry: 'rows', target: 'msa', columns: { row: [2, 0], species: ['mouse', 'human'] } },
      {
        key: 'm',
        geometry: 'matrix',
        target: '1brs/D',
        states: ['x', 'y'],
        columns: { index: [3], value: [[0.1, 0.9]] },
      },
    ],
    mappings: [],
  };
}

/** Refs a random query draws from. */
export const REFS = [
  '1brs/A',
  '1brs/D',
  '1brs',
  '1brs/entity:1',
  'ab/H',
  'ab',
  'lonely',
  'msa',
  'msa/1',
  'sites',
  'dom',
  'lone',
  'cols',
  'rowsite',
  'multi',
  'contacts',
  'species',
];

/**
 * QR2a's fixture for `onto` and `partners`: 1BRS with its 43 contacts within 4 Å (lib/seqstar/core/_test/fixtures/1brs.json,
 * chain A to chain D, those under 3.2 Å as `hbond`) and a pairs annotation whose ends are a whole unit, a molecule over
 * three units, a molecule with no unit, a row and the columns; case 14's protein over its spliced CDS; a family
 * alignment whose row 0 is bound to barnase's molecule (so chains A, B and C), row 1 to chain D, row 2 to case 14's
 * protein and chain E, row 3 to nothing; an A3M alignment with inserts; a protein mapped onto case 14's (a path through
 * the codon step), one mapped twice (an ambiguous path), and a molecule with no unit.
 */
export function ontoFixture(): SetInput {
  const b = brs1();
  const rand = prng(20261012);
  const COLUMNS = 130;
  /** A sequence's residues at sorted random columns, gaps elsewhere. */
  const spread = (seq: string) => {
    const cols = new Set<number>();
    while (cols.size < seq.length) cols.add(Math.floor(rand() * COLUMNS));
    const row = Array(COLUMNS).fill('-');
    for (const [i, c] of [...cols].sort((x, y) => x - y).entries()) row[c] = seq[i];
    return row.join('');
  };
  const barnase = b.molecules[0].residues as string;
  const barstar = b.molecules[1].residues as string;
  const json = JSON.parse(readFileSync(new URL('../fixtures/1brs.json', import.meta.url), 'utf8'));
  const contacts = json.contacts as [number, number, number][];
  return {
    molecules: [
      ...b.molecules,
      { key: 'gene', residues: GENE, alphabet: 'dna', ownUnit: true },
      { key: 'protein', residues: PROTEIN, ownUnit: true },
      { key: 'p2', residues: PROTEIN.slice(2, 18), ownUnit: true },
      { key: 'amb', residues: 'SKELTPRDG', ownUnit: true },
      { key: 'lonely', residues: 'ACDEFGHIKL' },
    ],
    units: b.units,
    complexes: b.complexes,
    alignments: [
      {
        key: 'fam',
        gapped: [spread(barnase), spread(barstar.slice(0, 40)), spread(PROTEIN), spread(barstar.slice(50, 80))],
        flavor: 'aligned',
      },
      { key: 'ins', gapped: ['ACdeFG-H', 'A-CDeeFG'], flavor: 'a3m' },
    ],
    mappings: [
      { key: 'r0', from: 'fam/0', to: '1brs/entity:1', method: 'segments', segments: [[0, 0, 110]] },
      { key: 'r1', from: '1brs/D', to: 'fam/1', method: 'segments', segments: [[20, 0, 40]] },
      { key: 'r2a', from: 'fam/2', to: 'protein', method: 'segments', segments: [[0, 0, 20]] },
      { key: 'r2b', from: 'fam/2', to: '1brs/E', method: 'segments', segments: [[0, 5, 20]] },
      { key: 'codon', from: 'protein', to: 'gene', method: 'codon', cds: CDS },
      { key: 'm-p2', from: 'p2', to: 'protein', method: 'segments', segments: [[0, 2, 16]] },
      { key: 'amb1', from: 'amb', to: 'protein', method: 'segments', segments: [[0, 3, 9]] },
      { key: 'amb2', from: 'amb', to: 'protein', method: 'segments', segments: [[0, 5, 9]] },
      { key: 'm-lone', from: 'lonely', to: 'p2', method: 'segments', segments: [[0, 1, 10]] },
    ],
    annotations: [
      {
        key: 'contacts',
        geometry: 'pairs',
        columns: {
          target: contacts.map(() => '1brs/A'),
          start: contacts.map((c) => c[0] - 48),
          b_target: contacts.map(() => '1brs/D'),
          b_start: contacts.map((c) => c[1] - 2),
          category: contacts.map((c) => (c[2] < 3.2 ? 'hbond' : 'contact')),
          distance: contacts.map((c) => c[2]),
          label: contacts.map((c) => `${c[0] - 47}-${c[1] - 1}`),
        },
      },
      {
        key: 'links',
        geometry: 'pairs',
        columns: {
          target: ['1brs/B', '1brs/entity:2', 'fam/3', 'lonely', '1brs/C'],
          start: [5, 10, 2, 1, 40],
          stop: [6, 12, 4, 3, 45],
          b_target: ['1brs/F', 'fam', 'lonely', '1brs/A', '1brs/C'],
          b_start: [null, 3, 1, 60, 50],
          b_stop: [null, 5, 3, 62, 52],
          category: ['whole', 'molecule', 'row', 'lonely', 'self'],
          id: ['w', 'm', 'r', 'l', 's'],
        },
      },
      { key: 'dom', geometry: 'ranges', target: '1brs/entity:1', columns: { start: [20, 90], stop: [50, 105] } },
      {
        key: 'sites',
        geometry: 'ranges',
        columns: {
          target: ['1brs/A', '1brs/A', '1brs/D', 'protein', 'gene', 'fam'],
          start: [5, 30, 60, 4, 30, 10],
          stop: [10, 40, 70, 9, 60, 20],
          label: ['a', 'b', 'c', 'kinase-like', 'exon site', 'Kinase motif'],
          id: ['s1', 's2', 's3', 's4', 's5', 's6'],
          category: ['x', 'y', 'x', 'y', 'x', 'z'],
          value: [0.5, 2, 3.5, 1, 4, 0.1],
        },
      },
    ],
  };
}

/** Refs a random QR2a query draws from. */
export const ONTO_REFS = [
  '1brs/A',
  '1brs/D',
  '1brs/E',
  '1brs',
  '1brs/entity:1',
  'gene',
  'protein',
  'p2',
  'amb',
  'lonely',
  'fam',
  'fam/0',
  'fam/1',
  'fam/2',
  'fam/3',
  'ins',
  'ins/0',
  'ins/1',
  'dom',
  'sites',
  'contacts',
  'links',
];
