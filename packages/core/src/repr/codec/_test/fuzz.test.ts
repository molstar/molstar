/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A short seeded fuzz of the readers (testing.md "Parsers and readers"): random text and mutations of valid text,
// a small fixed count in the quick lane. A reader never throws: it gives a tree, or issues.

import { describe, expect, it } from 'vitest';
import { PCG } from '../../../data/util/hash-functions.js';
import { parseJson, printJson } from '../../codecs/json/index.js';
import { parseSexpr, printSexpr } from '../../codecs/sexpr/index.js';
import { parseAddress } from '../../tree/index.js';
import { conformanceTrees } from '../index.js';

const Count = 300;
const Alphabet = ['(', ')', '[', ']', '{', '}', '"', ':', ';', '/', '\\', ' ', ',', 'a', '1', '-', '+', '.', 'e', 'n'];
const Snippets = ['(range :from 1)', '/data/p53', '"x"', ':$column', ':$custom', 'nan', '7n', '×', '{"kind":"x"}', '*'];

function samples(seed: number, valid: readonly string[]): string[] {
  const rng = new PCG(seed);
  const pick = <T>(xs: readonly T[]) => xs[rng.int() % xs.length];
  const out: string[] = [];
  for (let k = 0; k < Count; k++) {
    let text = k % 2 ? pick(valid) : '';
    const edits = 1 + (rng.int() % 6);
    for (let e = 0; e < edits; e++) {
      const at = text.length ? rng.int() % (text.length + 1) : 0;
      const op = rng.int() % 3;
      if (op === 0) text = text.slice(0, at) + pick(Alphabet) + text.slice(at);
      else if (op === 1) text = text.slice(0, at) + pick(Snippets) + text.slice(at);
      else text = text.slice(0, at) + text.slice(at + 1);
    }
    out.push(text);
  }
  return out;
}

const trees = conformanceTrees().map((t) => t.node);

describe('repr: readers never throw', () => {
  it('parseAddress', () => {
    for (const text of samples(1, ['/data/p53/(range :from 94 :to [1 2] :of /a/(b :c "d"))/*', '/a/b'])) {
      const r = parseAddress(text);
      expect(r.pattern !== undefined || r.issues.length > 0, text).toBe(true);
    }
  });

  it('the S-expression reader', () => {
    const valid = trees.map((t) => printSexpr(t));
    for (const text of samples(2, valid)) {
      const r = parseSexpr(text, { dataForms: new Set(['expr']) });
      expect(r.node !== undefined || r.issues.length > 0, text).toBe(true);
    }
  });

  it('the JSON reader', () => {
    const valid = trees.map((t) => printJson(t));
    for (const text of samples(3, valid)) {
      const r = parseJson(text);
      expect(r.node !== undefined || r.issues.length > 0, text).toBe(true);
    }
  });
});
