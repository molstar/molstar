/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QR2b's runtimes' cases, for their heap guard and its quick twin (query-qr2b.heap.ts, query-qr2b.quick.ts)
// (seqspec:seqql#performance): `letters`, `positions` with a test per position, `search` with each matcher but the
// regex one (which allocates its match per hit, seqspec:libraries#search), the group tests over a selection unioned
// once, and system tests over a 100,000-row alignment; a step is one run of a compiled query over 100,000 runs, groups
// or rows (a million positions scanned or tested). Label tests (`=`, `!=`, `in` on `(label)`, ordinal and through a
// labels and an offset numbering) compare as the numbering holds labels, no string per position (CM-04). The
// environment has the `(leak Q)` extension of the guard's positive control (query-leak.ts).

import { checkText } from '@molstar/seq-core/query';
import { SeqSets } from '../../src/model/index.js';
import { compileQuery, createQueryEnv } from '../../src/model/query/index.js';
import { leakExtension } from './query-leak.js';

const N = 1_000_000;
const ROWS = 100_000;
const { set } = SeqSets.build({
  // ACGT: 250,000 runs of each letter, 250,000 hits of ACG
  // Numbered two ways (LA-15: the numbered path of `labelEquals`, where `Numbering.format` made a string a position):
  // author labels with a letter or a multi-digit `.n` suffix every tenth, and an offset.
  molecules: [
    {
      key: 'big',
      residues: 'ACGT'.repeat(N / 4),
      ownUnit: true,
      alphabet: 'dna',
      numberings: [
        {
          name: 'auth',
          method: 'labels',
          labels: Array.from({ length: N }, (_, i) => (i % 10 === 3 ? `${i}A` : i % 10 === 7 ? `${i}.12` : i + 1)),
        },
        { name: 'off', method: 'offset', first: -47 },
      ],
    },
  ],
  alignments: [
    { key: 'msa', gapped: Array.from({ length: ROWS }, (_, r) => (r % 2 ? 'AC-GT' : 'ACGGT')), flavor: 'aligned' },
  ],
  annotations: [
    // passes in runs of 5 every 10, values that are no small integers (a boxed double would show)
    {
      key: 'v',
      geometry: 'values',
      target: 'big',
      columns: { value: Float64Array.from({ length: N }, (_, i) => (i % 10 < 5 ? 0.75 : 0.25)) },
    },
    // 50,000 items of 4 positions, every 20
    {
      key: 'reps',
      geometry: 'ranges',
      target: 'big',
      columns: {
        start: Int32Array.from({ length: N / 20 }, (_, i) => 20 * i + 3),
        stop: Int32Array.from({ length: N / 20 }, (_, i) => 20 * i + 7),
      },
    },
    {
      key: 'stats',
      geometry: 'rows',
      target: 'msa',
      columns: {
        row: Int32Array.from({ length: ROWS }, (_, r) => r),
        identity: Float64Array.from({ length: ROWS }, (_, r) => (r % 4) * 0.3 + 0.05),
      },
    },
  ],
});

const env = createQueryEnv([leakExtension({ dna: { N: 'ACGT', R: 'AG', Y: 'CT' } })]);

export const QUERIES: [string, number][] = [
  ['(letters A :in big)', 250_000],
  ['(and big (letters CG))', 1],
  ['(positions :in big :where (< (index) 500000))', 1],
  ['(positions :in big :where (= (letter) T))', 250_000],
  ['(positions :in big :where (> (value-of v) 0.5))', 100_000],
  // CM-04: a label compared as the numbering holds it (an ordinal label here), never a string a position.
  ['(positions :in big :where (= (label) "500000"))', 1],
  ['(positions :in big :where (!= (label) "1"))', 1],
  ['(positions :in big :where (in (label) "10" "20" "999999"))', 3],
  // LA-15: the same through a numbering with suffixes, and an offset one.
  ['(positions :in big :where (= (label :numbering auth) "500001"))', 1],
  ['(positions :in big :where (= (label :numbering auth) "13A"))', 1],
  ['(positions :in big :where (in (label :numbering auth) "10" "23A" "57.12"))', 3],
  ['(positions :in big :where (= (label :numbering off) "100"))', 1],
  ['(positions :in big :where (in (label :numbering off) "-40" "20" "999000"))', 3],
  ['(search ACG :in big)', 250_000],
  ['(search NRY :in big :classes dna)', 250_000],
  ['(search ACGA :in big :edits 1)', 200_000],
  ['(pick (values :of v :above 0.5) :where (overlaps reps))', 50_000],
  ['(pick (values :of v :above 0.5) :where (within reps 2))', 50_000],
  ['(pick (values :of v :above 0.5) :where (>= (count reps) 2))', 50_000],
  ['(pick (values :of v :above 0.5) :where (and (= (runs) 1) (> (length) 4)))', 100_000],
  ['(rows :of msa :where (and (> (length) 3) (= (kind) row) (< (score) 1)))', 0],
];

/** A query compiled, and its run: the step, and the groups the last run selected. */
export function prepare(text: string): { run: () => void; groups: () => number } {
  const q = compileQuery(checkText(text, env.language).expr as never, env);
  let groups = 0;
  const run = () => {
    groups = q.run({ set }).selection.groupCount;
  };
  return { run, groups: () => groups };
}
