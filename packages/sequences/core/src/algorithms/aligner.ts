/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The aligner interface (A1; seqspec:format#mappings): an aligner takes two sequences and gives a CIGAR from the first
// to the second (seqspec:model#coordinate-maps: M a pair, I a letter only in `a`, D a letter only in `b`). Its params
// are plain data, so the format (V2) builds its registered params table from them without this folder importing
// the format.

import { SEQUENCE_TYPE } from './biomsa/align/params.js';
import { getSequenceType } from './biomsa/sequence/sequence.js';

/** One option of an aligner, described as data. */
export interface AlignerParam {
  /** `choice`: one of `choices`; `str`: a string, matching `pattern` when given. */
  readonly type: 'int' | 'float' | 'choice' | 'str';
  readonly choices?: readonly string[];
  /** A `str`'s format: a regular expression (its source) the whole string must match. */
  readonly pattern?: string;
  /** Whether `null` is allowed; `null` means "the default for the alphabet" where the description says so. */
  readonly nullable: boolean;
  /** The inclusive lower bound of a number. */
  readonly min?: number;
  readonly default: string | number | null;
  readonly description: string;
}
export type AlignerParams = { readonly [name: string]: AlignerParam };
/** Options as a document gives them: missing (or `undefined`) means the param's default. */
export type AlignerOptions = { readonly [name: string]: string | number | null | undefined };
/** Options after `readOptions`: every param present. */
export type ResolvedOptions = { readonly [name: string]: string | number | null };

export interface Aligner {
  readonly name: string;
  readonly description: string;
  readonly params: AlignerParams;
  /**
   * A CIGAR from `a` to `b` covering both sequences exactly; null when these sequences cannot be aligned (empty, a
   * character other than a letter, or alphabets the aligner cannot pair). Deterministic: the same sequences and
   * options give the same CIGAR. Letters are case-insensitive. Throws an `AlignerOptionsError` on invalid options and an
   * `AlignerSizeError` (limits.ts) when a full-matrix DP would be too large.
   */
  align(a: string, b: string, options?: AlignerOptions): Promise<string | null>;
  /**
   * A multiple alignment, for an aligner that has one: the gapped rows of `seqs` (at least two), in their order and
   * upper case, `-` a gap, all of one length; null when they cannot be aligned. Deterministic. Honors `signal`: a long
   * run yields to the event loop between its steps and throws the signal's reason when aborted. Throws as `align` does.
   */
  alignMany?(seqs: readonly string[], options?: AlignerOptions, signal?: AbortSignal): Promise<string[] | null>;
}

/** Thrown for an invalid option or option combination: the only error a caller reports as the options' fault. */
export class AlignerOptionsError extends RangeError {
  constructor(message: string) {
    super(message);
    this.name = 'AlignerOptionsError';
  }
}

/** What an aligner treats a sequence as. */
export type SequenceType = 'amino' | 'nucleic';

const LETTERS = /^[A-Za-z]+$/;

/** Whether both sequences are non-empty and letters only. */
export function alignable(a: string, b: string): boolean {
  return LETTERS.test(a) && LETTERS.test(b);
}

/** The sequence's type as biomsa detects it; null when it recognizes neither (an upper-case sequence). */
export function detectType(seq: string): SequenceType | null {
  try {
    return getSequenceType(seq) === SEQUENCE_TYPE.PROTEIN ? 'amino' : 'nucleic';
  } catch {
    return null;
  }
}

/**
 * The type both sequences share: the `type` option when it is not `auto`, otherwise detected; null when the
 * detected types differ or either is not recognized.
 */
export function pairType(a: string, b: string, type: string | number | null): SequenceType | null {
  if (type === 'amino' || type === 'nucleic') return type;
  const ta = detectType(a);
  return ta !== null && ta === detectType(b) ? ta : null;
}

/** Checks `options` against `params` and fills in defaults; throws an `AlignerOptionsError` naming the aligner and option. */
export function readOptions(aligner: string, params: AlignerParams, options: AlignerOptions = {}): ResolvedOptions {
  const out: { [name: string]: string | number | null } = {};
  for (const name of Object.keys(options)) {
    if (!Object.hasOwn(params, name)) {
      throw new AlignerOptionsError(
        `${aligner}: unknown option "${name}"; the options are ${Object.keys(params).join(', ')}.`,
      );
    }
  }
  for (const [name, p] of Object.entries(params)) {
    const v = options[name] === undefined ? p.default : options[name];
    const fail = (expected: string) => {
      throw new AlignerOptionsError(`${aligner}: option "${name}" must be ${expected}; got ${JSON.stringify(v)}.`);
    };
    const str = p.pattern === undefined ? 'a string' : `a string matching /${p.pattern}/`;
    if (v === null) {
      if (!p.nullable)
        fail(p.type === 'choice' ? `one of ${p.choices?.join(', ')}` : p.type === 'str' ? str : 'a number');
    } else if (p.type === 'choice') {
      if (typeof v !== 'string' || !p.choices?.includes(v)) fail(`one of ${p.choices?.join(', ')}`);
    } else if (p.type === 'str') {
      if (typeof v !== 'string' || (p.pattern !== undefined && !new RegExp(`^(?:${p.pattern})$`).test(v))) fail(str);
    } else if (typeof v !== 'number' || !Number.isFinite(v) || (p.type === 'int' && !Number.isInteger(v))) {
      fail(p.type === 'int' ? 'an integer' : 'a finite number');
    } else if (p.min !== undefined && v < p.min) {
      fail(`at least ${p.min}`);
    }
    out[name] = v;
  }
  return out;
}

/** Run-length encodes per-column operations (`M`, `I`, `D`, one character each) as a CIGAR. */
export function cigarOf(ops: string): string {
  let cigar = '';
  for (let k = 0; k < ops.length; ) {
    let r = k + 1;
    while (r < ops.length && ops[r] === ops[k]) r++;
    cigar += `${r - k}${ops[k]}`;
    k = r;
  }
  return cigar;
}

/** A CIGAR from two gapped rows of equal length (`gap` marks a gap; a column gapped in both is skipped). */
export function cigarFromRows(a: string, b: string, gap = '-'): string {
  if (a.length !== b.length) throw new RangeError(`rows of different lengths: ${a.length} and ${b.length}`);
  let ops = '';
  for (let i = 0; i < a.length; i++) {
    const inA = a[i] !== gap;
    const inB = b[i] !== gap;
    if (inA || inB) ops += inA && inB ? 'M' : inA ? 'I' : 'D';
  }
  return cigarOf(ops);
}
