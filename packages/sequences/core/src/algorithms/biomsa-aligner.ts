/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `biomsa` aligner (seqspec:format#mappings): biomsa's pairwise alignment (vendored in ./biomsa, VENDORED.md),
// complete Needleman–Wunsch for short sequences and a diagonal heuristic for long ones, its gapped rows read as a
// CIGAR.

import {
  type Aligner,
  AlignerOptionsError,
  type AlignerParams,
  alignable,
  cigarFromRows,
  detectType,
  type ResolvedOptions,
  readOptions,
  type SequenceType,
} from './aligner.js';
import BioMSA from './biomsa/index.js';

/** The matrices biomsa ships, each for one sequence type (biomsa/align/params.ts). */
export const BiomsaMatrices = { VTML240: 'amino', BLASTZ: 'nucleic' } as const;
export type BiomsaMatrix = keyof typeof BiomsaMatrices;

/** biomsa's default penalties per type, as positive costs (its own values are the negated scores). */
const DEFAULTS = { amino: { gap_open: 300, gap_extend: 11 }, nucleic: { gap_open: 400, gap_extend: 60 } } as const;
const GAP = '-';

const params: AlignerParams = {
  matrix: {
    type: 'choice',
    choices: Object.keys(BiomsaMatrices),
    nullable: true,
    default: null,
    description: "The substitution matrix: VTML240 for amino acids, BLASTZ for nucleotides; null: the type's matrix.",
  },
  gap_open: {
    type: 'float',
    nullable: true,
    min: 0,
    default: null,
    description: "The cost of opening a gap, on the matrix's scale; null: 300 for amino acids, 400 for nucleotides.",
  },
  gap_extend: {
    type: 'float',
    nullable: true,
    min: 1,
    default: null,
    description:
      'The cost of extending a gap (biomsa adds twice it to every substitution score, as MUSCLE does); ' +
      'null: 11 for amino acids, 60 for nucleotides.',
  },
  method: {
    type: 'choice',
    choices: ['auto', 'complete', 'diag'],
    nullable: false,
    default: 'auto',
    description:
      'complete: global Needleman–Wunsch; diag: diagonals of shared k-mers first, global alignment between them; ' +
      'auto: diag when a sequence is longer than 1,600.',
  },
  type: {
    type: 'choice',
    choices: ['auto', 'amino', 'nucleic'],
    nullable: false,
    default: 'auto',
    description:
      "What the sequences are; auto: the matrix's type when a matrix is given, otherwise detected (both must agree). Ambiguous letters score as a fixed standard one.",
  },
};

/** The type `seqs` are aligned as under `o` (the matrix names it, else `type`, else detected: all must agree); null: not alignable. */
function typeOf(seqs: readonly string[], o: ResolvedOptions): SequenceType | null {
  const matrix = o.matrix as BiomsaMatrix | null;
  // A matrix scores one type of sequence: with `type: auto` it names the type; an explicit other type is an error.
  if (matrix !== null) {
    if (o.type === 'auto' || o.type === BiomsaMatrices[matrix]) return BiomsaMatrices[matrix];
    throw new AlignerOptionsError(
      `biomsa: matrix ${matrix} scores ${BiomsaMatrices[matrix]} sequences; type is ${o.type}.`,
    );
  }
  if (o.type === 'amino' || o.type === 'nucleic') return o.type;
  const first = detectType(seqs[0]);
  return first !== null && seqs.every((s) => detectType(s) === first) ? first : null;
}

/** Biomsa's gapped rows of `seqs` (upper case), awaiting `step` between the alignments of a multiple alignment. */
async function rowsOf(
  seqs: readonly string[],
  type: SequenceType,
  o: ResolvedOptions,
  step?: () => void | Promise<void>,
): Promise<string[]> {
  try {
    return await BioMSA.align(
      [...seqs],
      {
        type,
        method: o.method as 'auto' | 'complete' | 'diag',
        gapopen: -(o.gap_open ?? DEFAULTS[type].gap_open),
        gapextend: -(o.gap_extend ?? DEFAULTS[type].gap_extend),
        gapchar: GAP,
      },
      step,
    );
  } catch (e) {
    // biomsa rejects with plain strings on its argument checks.
    throw typeof e === 'string' ? new Error(`biomsa: ${e}`) : e;
  }
}

export const biomsa: Aligner = {
  name: 'biomsa',
  description: 'biomsa (Paul Pillot, MIT): progressive alignment as MUSCLE does it, pairwise here.',
  params,
  async align(a, b, options) {
    const o = readOptions('biomsa', params, options);
    if (!alignable(a, b)) return null;
    const A = a.toUpperCase();
    const B = b.toUpperCase();
    const type = typeOf([A, B], o);
    if (type === null) return null;
    const rows = await rowsOf([A, B], type, o);
    return cigarFromRows(rows[0], rows[1], GAP);
  },
  async alignMany(seqs, options, signal) {
    const o = readOptions('biomsa', params, options);
    if (seqs.length < 2) throw new AlignerOptionsError('biomsa: a multiple alignment needs at least two sequences.');
    if (!seqs.every((s) => alignable(s, s))) return null;
    const S = seqs.map((s) => s.toUpperCase());
    const type = typeOf(S, o);
    if (type === null) return null;
    signal?.throwIfAborted();
    let last = Date.now();
    return rowsOf(S, type, o, async () => {
      signal?.throwIfAborted();
      // Yield to the event loop about every 20 ms of work, so a newer request can abort the run.
      if (Date.now() - last < 20) return;
      await new Promise<void>((r) => setTimeout(r, 0));
      signal?.throwIfAborted();
      last = Date.now();
    });
  },
};
