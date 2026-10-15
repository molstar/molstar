/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST7's cases, for its heap guard and its quick twin (paths-st7.heap.ts, paths-st7.quick.ts): reading a sparse composed
// path (seqspec:model#paths) allocates nothing per step. A DMD-sized protein (79 exons over 2.2 Mb) coded on a lazy
// 200 Mb chromosome, its path read per element both ways (`Paths.at`, `Paths.spanAt`) and by ranges either way
// (`forEachFrom`: a streamed request's runs; `forEachTo`: an item's image), as the stream manager's pass and its draws
// read it.

import { type CoordinateSystem, Mappings, type PathStep, Paths, SeqSets } from '../../src/model/index.js';
import { gene as dmd } from './support/generators.js';

const MASK = 0xfffff;
const AT = 100_000_000;
const g = dmd();
const { set } = SeqSets.build({
  molecules: [
    { key: 'chr', alphabet: 'dna', lazy: { length: 200_000_000 }, ownUnit: true },
    { key: 'p', alphabet: 'protein', residues: g.protein, ownUnit: true },
  ],
  mappings: [{ key: 'c', from: 'p', to: 'chr', method: 'codon', cds: g.exons.map(([s, t]) => [AT + s, AT + t]) }],
});
const sys = (key: string): CoordinateSystem => {
  const s = { kind: -1, id: -1, row: -1 };
  Mappings.systemOf(set, key, s);
  return s;
};
const steps: PathStep[] = [];
Mappings.path(set, sys('chr/unit'), sys('p/unit'), steps);
export const down = Mappings.compose(set, sys('chr/unit'), steps);
Mappings.path(set, sys('p/unit'), sys('chr/unit'), steps);
export const up = Mappings.compose(set, sys('p/unit'), steps);
const span = g.exons[g.exons.length - 1][1] - g.exons[0][0];

const ctx = { n: 0 };
let acc = 0;
function runSum(s: number, t: number, c: { n: number }) {
  c.n = (c.n + s + t) & MASK;
}
/** A step: one element read both ways, and ranges either way. */
export const read = (i: number) => {
  const r = (i * 7919) % 3685;
  acc = (acc + Paths.at(down, AT + g.exons[0][0] + ((i * 104_729) % span))) & MASK;
  acc = (acc + Paths.spanAt(up, r, i % 3)) & MASK;
  Paths.forEachFrom(down, r, r + (i % 40) + 1, runSum, ctx);
  if ((i & 3) === 0) Paths.forEachTo(down, AT + ((i * 31) % span), AT + ((i * 31) % span) + 5000, runSum, ctx);
};

/** The steps' results: finite, or a step that reads nothing would pass blind. */
export const sums = () => acc + ctx.n;
