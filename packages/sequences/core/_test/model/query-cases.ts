/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqQL's runtime cases, for their heap guard and its quick twin (query.heap.ts, query.quick.ts)
// (seqspec:seqql#performance, QR1): a query allocates its typed result arrays (off the JS heap), its scratch and one closure tree at compile, never
// an object per run or per position; a step is one run of a compiled query over 100,000 runs or groups (a million
// positions scanned). The environment has the `(leak Q)` extension of the guard's positive control (query-leak.ts).

import { checkText } from '@molstar/seq-core/query';
import { SeqSets } from '../../src/model/index.js';
import { compileQuery, createQueryEnv, SeqSelection } from '../../src/model/query/index.js';
import { leakExtension } from './query-leak.js';

const N = 1_000_000;
const { set } = SeqSets.build({
  molecules: [{ key: 'big', residues: 'ACGT'.repeat(N / 4), ownUnit: true, alphabet: 'dna' }],
  annotations: [
    // passes in runs of 5 every 10: 100,000 runs
    {
      key: 'v',
      geometry: 'values',
      target: 'big',
      columns: { value: Float64Array.from({ length: N }, (_, i) => (i % 10 < 5 ? 1 : 0)) },
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
  ],
});
const env = createQueryEnv([leakExtension()]);
const SELECTION = SeqSelection.fromExtents(set, {
  count: 1,
  kind: Uint8Array.of(1),
  id: Int32Array.of(0),
  row: Int32Array.of(-1),
  start: Int32Array.of(1000),
  stop: Int32Array.of(900_000),
});

export const QUERIES = [
  '(values :of v :above 0.5)',
  '(and (values :of v :above 0.5) (not reps))',
  '(and (values :of v :above 0.5) reps (input selection))',
  '(range (values :of v :above 0.5) 1-500000 600000-900000)',
  '(slice (values :of v :above 0.5) 0 500000 600000 900000)',
  '(or reps (values :of v :above 0.5) reps)',
  '(not reps :in big)',
  '(items :of reps)',
];

/** A query compiled, and its run: the step, and the groups the last run selected. */
export function prepare(text: string): { run: () => void; groups: () => number } {
  const q = compileQuery(checkText(text, env.language).expr as never, env);
  let groups = 0;
  const run = () => {
    groups = q.run({ set, inputs: { selection: SELECTION } }).selection.groupCount;
  };
  return { run, groups: () => groups };
}
