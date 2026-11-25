/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The search's contract (_spec/README.md "Search"): the matcher port an application implements, its targets, batches
// and hits, the bounds, and the kernel a matcher of kernels runs. Types and the bounds' defaults only: what the default
// viewer loads at once, the matchers themselves a lazy chunk.

import type { Observable } from '@molstar/core/reactive';
import type { Issue } from '@molstar/core/repr/tree';
import type { SequenceAlphabet } from '@molstar/core/util/sequence-constants';

/** The strands a nucleotide pattern is searched on: both (the default) or the forward one only. */
export type SearchStrands = 'both' | 'forward';

/** What a compile is given besides the pattern. */
export interface SearchOptions {
  /** Both strands of a nucleotide target (`both`, the default), or the forward one (`forward`). Proteins have one. */
  readonly strands?: SearchStrands;
}

/** Every bound of a search (_spec/README.md "Search", Bounds), each checked before the work it sizes. */
export interface SearchBounds {
  /** The longest `exact` or `degenerate` pattern, in letters (1,000). */
  readonly maxPattern: number;
  /** A regular expression's characters (256). */
  readonly maxRegexChars: number;
  /** A regular expression's letter positions once its counted repeats are expanded (128: four words). */
  readonly maxRegexPositions: number;
  /** A regular expression's nesting of groups (64). */
  readonly maxRegexDepth: number;
  /** The longest regular expression hit: the longest match from a start stops there (1,000). */
  readonly maxHitLength: number;
  /** The letters a run reads over all its targets (2^25): the target that would pass it is skipped. */
  readonly maxLetters: number;
  /** The word steps a run takes (2^28): the run ends there, incomplete. */
  readonly maxSteps: number;
  /** The hits a run holds (10,000): counting goes on past it. */
  readonly maxHits: number;
  /** The letters read at once, a batch each (2^20). */
  readonly chunkLetters: number;
}

export const DefaultSearchBounds: SearchBounds = {
  maxPattern: 1000,
  maxRegexChars: 256,
  maxRegexPositions: 128,
  maxRegexDepth: 64,
  maxHitLength: 1000,
  maxLetters: 2 ** 25,
  maxSteps: 2 ** 28,
  maxHits: 10_000,
  chunkLetters: 2 ** 20,
};

/** The bounds with the given ones over the defaults. */
export function searchBounds(over: Partial<SearchBounds> = {}): SearchBounds {
  return { ...DefaultSearchBounds, ...over };
}

/** A target's letters by region: `[start, end)` of its positions, one code a position, gaps (`-`, `.`) included. */
export interface LetterReader {
  read(start: number, end: number, signal: AbortSignal): Promise<Uint8Array>;
}

/** One system a run searches. */
export interface SearchTarget {
  /** Its address. */
  readonly system: string;
  readonly alphabet: SequenceAlphabet;
  readonly length: number;
  readonly cyclic: boolean;
  /** The ranges searched in (`[start, end)` pairs, sorted), from `within`: the whole system when absent. */
  readonly within?: Int32Array;
  readonly letters: LetterReader;
}

/**
 * Hits, columnar: no object a hit. `target` the index of its target in the run, `[start, end)` its positions (`end`
 * past the length for a hit across a circle's origin: unrolled), `strand` `1`, `-1`, or `0` (a protein, or a site the
 * same on both strands).
 */
export interface SearchHits {
  readonly length: number;
  readonly target: Int32Array;
  readonly start: Int32Array;
  readonly end: Int32Array;
  readonly strand: Int8Array;
}

/** What a run streams: the hits found since the last batch, in order, and what it read and said. */
export interface SearchBatch {
  readonly hits: SearchHits;
  /** Hits found in this batch past the held bound: counted, not held. */
  readonly counted: number;
  /** Letters read so far (for a progress said once). */
  readonly scanned: number;
  /** A target skipped (its alphabet, too long, a failed read), a bound reached. */
  readonly issues: readonly Issue[];
}

/** What a run is given: its cancellation, its bounds, and a pause between batches (the viewer's frame). */
export interface SearchRunContext {
  readonly signal: AbortSignal;
  readonly bounds: SearchBounds;
  /** Awaited after each batch: the viewer yields a frame there, so a run never holds one. */
  readonly yield?: () => Promise<void>;
}

/** A pattern checked and compiled: the matcher's own data beside the pattern and the options. */
export interface SearchCompiled {
  readonly pattern: string;
  readonly strands: SearchStrands;
}

/**
 * The matcher port (_spec/README.md "Search", Matchers): what an application registers (`search-matcher`). `compile` is
 * synchronous and reads no letters; `run` streams batches in the final order, ending when done; unsubscribing, or the
 * signal, cancels it, and a cancelled run ends with no error.
 */
export interface SearchMatcher<C extends SearchCompiled = SearchCompiled> {
  /** `exact`, `degenerate`, `regex`, or an application's `<id>:<name>`. */
  readonly name: string;
  readonly title: string;
  /** The alphabets it reads; absent: any. */
  readonly alphabets?: readonly SequenceAlphabet[];
  compile(pattern: string, options: SearchOptions, bounds: SearchBounds): C | readonly Issue[];
  run(compiled: C, targets: readonly SearchTarget[], ctx: SearchRunContext): Observable<SearchBatch>;
}

/** Whether a compile's answer is its issues. */
export function compileFailed<C>(r: C | readonly Issue[]): r is readonly Issue[] {
  return Array.isArray(r);
}

/**
 * Called once a hit: `[start, end)` in the searched array's indices. The one callback of the surface, on purpose: the
 * hot loop makes no object a hit (an exception to returning observables).
 */
export type HitVisitor = (start: number, end: number) => void;

/** The word steps a run may still take: a kernel adds what it does past the letters it was charged. */
export interface StepBudget {
  steps: number;
  readonly limit: number;
}

/** A kernel's buffer edges: its first letter is the sequence's first (`^`), its last the sequence's last (`$`). */
export const EdgeStart = 1;
export const EdgeEnd = 2;

/**
 * A kernel (_spec/README.md "Search", Kernels): the hot loop, synchronous, over bytes it knows nothing of. Calls `visit`
 * for every hit lying in `[0, n)` of `seq` that starts in `[from, to)`; nothing allocated a hit. The letters it reads
 * are charged before it runs (`words` a letter); what it does past them it adds to `budget`, and it returns false when
 * the budget ran out. Not reentrant.
 */
export interface SearchKernel {
  /** The longest hit it reports. */
  readonly maxLength: number;
  /** The word steps a letter costs. */
  readonly words: number;
  find(
    seq: Uint8Array,
    n: number,
    from: number,
    to: number,
    edges: number,
    budget: StepBudget,
    visit: HitVisitor,
  ): boolean;
}
