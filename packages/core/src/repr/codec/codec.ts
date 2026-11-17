/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The codec interface (_spec/README.md "Codecs"; tree-core.md): a tree to text or bytes and back. No codec lives here:
// JSON and the S-expression printing are `repr/codecs/*`, MessagePack is `@molstar/io`'s, beside its library. A codec
// is schema-free: it needs no language, so it may give a column back as a plain array, and the read pipeline's
// normalization (`schema/normalize.ts`) types it. What every codec guarantees is checked by the conformance suite
// (conformance.ts), run on every codec, built-in or not.

import { defineIssues, type Issue, type Node } from '../tree/index.js';

/** A decoded tree: `node` when there is no issue. */
export interface Decoded {
  readonly node?: Node;
  readonly issues: readonly Issue[];
}

export interface DecodeOptions {
  /** Kinds whose positional items are values (`args`), for codecs whose spelling cannot tell (S-expressions). */
  readonly dataForms?: ReadonlySet<string>;
}

/**
 * A tree's encoding. Guarantees, for every tree of the value model (atoms with -0, NaN, the infinities and bigints;
 * lists, records, slots, data forms, custom data, extension kinds, typed columns):
 *
 * - `decode(encode(x))`, normalized by the language that declares `x`'s columns, is `x` (equal hashes), except where
 *   the codec says it elides (`elides`: a printing for people, whose reader refuses what it elided);
 * - `encode` is deterministic: equal trees encode to equal output;
 * - `decode` never throws: what is not a tree is issues, each with a code and a place (an offset, a path).
 */
export interface Codec<T extends string | Uint8Array = string | Uint8Array> {
  /** `json`, `sexpr`, `msgpack`. */
  readonly id: string;
  /** What it writes. */
  readonly output: T extends string ? 'text' : 'binary';
  /** Whether `decode` gives typed columns as typed arrays (else as plain arrays, typed by normalization). */
  readonly typedColumns: boolean;
  /** Whether `encode` may leave out data (long columns), which `decode` then refuses. */
  readonly elides: boolean;
  encode(node: Node): T;
  decode(data: T, options?: DecodeOptions): Decoded;
}

/** The codecs' codes: reading text or bytes that are no tree. */
export const CodecIssues = defineIssues('repr.codec', {
  'repr.codec.json.syntax': { severity: 'error', description: 'The text is not JSON' },
  'repr.codec.sexpr.syntax': {
    severity: 'error',
    description: 'The text is not S-expressions: an unbalanced bracket, an unterminated string, a malformed word',
  },
  'repr.codec.sexpr.shape': {
    severity: 'error',
    description:
      'Forms in the wrong place: not exactly one top-level node, a positional other than the id, a keyword with no value',
  },
  'repr.codec.sexpr.duplicate_key': {
    severity: 'error',
    description: 'A node or record names the same keyword twice',
  },
  'repr.codec.sexpr.elided': {
    severity: 'error',
    description: 'A column printed elided (`[float32 × 3,000: …]`): the printing is for people, its data is not there',
    hint: 'encode the tree with JSON, or print it with a longer head',
  },
});
