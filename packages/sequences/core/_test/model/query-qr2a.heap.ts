/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for QR2a's runtimes (seqspec:seqql#performance): `onto` through a row's binding and a codon mapping,
// `partners` over 100,000 pairs, `expand`, `group`, `pick`, `first`, `ids` and an item test each run over 100,000 runs
// or 50,000 groups and allow 1 MB of heap growth (10 bytes a run). The positive control is the template's: an object
// a run, through `(leak Q)` (query-leak.ts), caught by the same measurement. `each`
// is not guarded: it builds a selection per group by definition (its `(current)` and E's result).

import { checkText } from '@molstar/seq-core/query';
import { describe, expect, it } from 'vitest';
import { SeqSets } from '../../src/model/index.js';
import { compileQuery, createQueryEnv } from '../../src/model/query/index.js';
import { controlCaught, heapGrowth, LIMIT } from './heap.js';
import { leakExtension } from './query-leak.js';

const N = 1_000_000;
const P = 300_000;
const PAIRS = 100_000;
const ACGT = 'ACGT'.repeat(N / 4);
/** The bound row: the molecule's residues with a gap column after every tenth. */
const row = ACGT.replace(/(.{10})/g, '$1-');
const { set } = SeqSets.build({
  molecules: [
    { key: 'big', residues: ACGT, ownUnit: true, alphabet: 'dna' },
    { key: 'prot', residues: 'MA'.repeat(P / 2), ownUnit: true },
    { key: 'gene', residues: 'ATG'.repeat(P), ownUnit: true, alphabet: 'dna' },
    { key: 'x', residues: 'A'.repeat(2 * PAIRS), ownUnit: true },
    { key: 'y', residues: 'A'.repeat(2 * PAIRS), ownUnit: true },
  ],
  alignments: [{ key: 'aln', gapped: [row], flavor: 'aligned' }],
  mappings: [
    { key: 'bind', from: 'aln/0', to: 'big', method: 'identity' },
    { key: 'codon', from: 'prot', to: 'gene', method: 'codon', start: 0 },
  ],
  annotations: [
    // passes in runs of 5 every 10: 100,000 runs on big, 30,000 on prot
    {
      key: 'v',
      geometry: 'values',
      target: 'big',
      columns: { value: Float64Array.from({ length: N }, (_, i) => (i % 10 < 5 ? 1 : 0)) },
    },
    {
      key: 'pv',
      geometry: 'values',
      target: 'prot',
      columns: { value: Float64Array.from({ length: P }, (_, i) => (i % 10 < 5 ? 1 : 0)) },
    },
    // 50,000 items of 4 positions, every 20
    {
      key: 'reps',
      geometry: 'ranges',
      target: 'big',
      columns: {
        start: Int32Array.from({ length: N / 20 }, (_, i) => 20 * i + 3),
        stop: Int32Array.from({ length: N / 20 }, (_, i) => 20 * i + 7),
        label: Array.from({ length: N / 20 }, (_, i) => (i % 2 ? 'even' : 'odd')),
        id: Array.from({ length: N / 20 }, (_, i) => `r${i}`),
      },
    },
    // 100,000 pairs between x and y, every other position
    {
      key: 'pp',
      geometry: 'pairs',
      columns: {
        target: Array(PAIRS).fill('x'),
        start: Int32Array.from({ length: PAIRS }, (_, i) => 2 * i),
        b_target: Array(PAIRS).fill('y'),
        b_start: Int32Array.from({ length: PAIRS }, (_, i) => 2 * i + 1),
      },
    },
  ],
});

const env = createQueryEnv([leakExtension()]);

const QUERIES: [string, number][] = [
  ['(onto (values :of v :above 0.5) aln)', 1],
  ['(onto (values :of pv :above 0.5) gene)', 1],
  ['(onto (sel "aln@0:1000000") big)', 1],
  ['(partners y :via pp)', 90_000],
  ['(expand (values :of v :above 0.5) 2)', 90_000],
  ['(group (values :of v :above 0.5) :by run)', 90_000],
  ['(group (items :of reps) :by system)', 1],
  ['(pick (items :of reps) :where (matches (label) "^even$"))', 20_000],
  ['(items :of reps :where (> (length) 3))', 40_000],
  ['(first (items :of reps) 40000)', 40_000],
  ['(ids (items :of reps) r1 r2 r3)', 3],
  ['(items :of reps :where (in (id) r1 r5 r9))', 3],
  ['(partners y :via pp :where (not (= (category) x)))', 90_000],
];

/** The least of three measured runs after ten warm-ups (a loaded machine's late compiles can box in one run). */
function growth(text: string): { groups: number; growth: number } {
  const q = compileQuery(checkText(text, env.language).expr as never, env);
  let groups = 0;
  const run = () => {
    groups = q.run({ set }).selection.groupCount;
  };
  for (let k = 0; k < 10; k++) run();
  return { groups, growth: Math.min(heapGrowth(run), heapGrowth(run), heapGrowth(run)) };
}

describe("QR2a's runtimes allocate nothing per run", () => {
  for (const [text, groups] of QUERIES) {
    it(text, () => {
      const r = growth(text);
      expect(r.groups).toBeGreaterThanOrEqual(groups);
      expect(r.growth).toBeLessThan(LIMIT);
    });
  }

  it("an object a run is caught by the same measurement (the template's control)", () => {
    for (const text of ['(onto (values :of v :above 0.5) aln)', '(partners y :via pp)']) {
      expect([text, controlCaught(() => growth(`(leak ${text})`).growth)]).toEqual([text, true]);
    }
  });
});
