/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The spelling of atoms and keywords, one for addresses (address.ts) and the S-expression codec (codecs/sexpr), so an
// address reads the same inside a spec and alone (_spec/README.md "Atoms"; locations-and-queries.md):
//
//   true false nil            the booleans and null (JSON's `null`)
//   nan +inf -inf             the numbers JSON lacks
//   123n -5n                  bigints
//   -1 1e5 0.25 .5 +2         a word that parses as a JavaScript decimal number is a number
//   94A HBA_HUMAN 2nm         any other word is a symbol: a string, written bare
//   "a b" "HBA_HUMAN/1-141"   a string in JSON's grammar, quoted when it would not read back bare
//
// A word is printable, non-space code points, no delimiter (`()[]";{},`) and not opening with `:` (a keyword) or `$`
// (reserved). Changed from meso*'s reader: a word that looks like a number but is not JSON's (`+1`, `.5`) was an
// issue there; here it is the number JavaScript reads, and a word that is no number (`94A`, `0x10`) is a symbol, so a
// residue label with an insertion code needs no quotes and a quantity (`2nm`) is a string with no units table.

import type { Atom } from './node.js';

/** A JavaScript decimal number, optionally signed (no hex, no `Infinity`, no separators). */
const NumberWord = /^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$/;
const BigIntWord = /^-?[0-9]+n$/;

/** What may be a bare word. */
export const WordPattern = /^[^\p{C}\p{Z}()";[\]{},:$][^\p{C}\p{Z}()";[\]{},]*$/u;
/** A keyword's name written bare (`:name`); any other is quoted (`:"p.ALB"`). */
export const KeywordPattern = /^[A-Za-z_][A-Za-z0-9_-]*$/;
/** What ends a word. */
export const Delimiter = /[ \t\r\n()";[\]{},]/;

const Words: ReadonlyMap<string, Atom> = new Map<string, Atom>([
  ['true', true],
  ['false', false],
  ['nil', null],
  ['nan', Number.NaN],
  ['+inf', Number.POSITIVE_INFINITY],
  ['-inf', Number.NEGATIVE_INFINITY],
]);

/** The atom a bare word reads as, or undefined when it is no word. */
export function readWord(word: string): Atom | undefined {
  if (Words.has(word)) return Words.get(word);
  if (NumberWord.test(word)) return Number(word);
  if (BigIntWord.test(word)) return BigInt(word.slice(0, -1));
  return WordPattern.test(word) ? word : undefined;
}

/**
 * Whether a string prints bare: a word that reads back as the same string, with no `/` (a word opening with `/` is an
 * address, the caller's to check; one holding a `/` is quoted so that it never reads as a path: `"HBA_HUMAN/1-141"`).
 */
export function isBareString(s: string): boolean {
  return readWord(s) === s && !s.includes('/');
}

/** An atom's text: a string bare when it reads back as itself (`isBareString`), else quoted. */
export function printAtom(v: Atom): string {
  switch (typeof v) {
    case 'number':
      if (Number.isNaN(v)) return 'nan';
      if (!Number.isFinite(v)) return v > 0 ? '+inf' : '-inf';
      return Object.is(v, -0) ? '-0' : String(v);
    case 'bigint':
      return `${v}n`;
    case 'string':
      return isBareString(v) ? v : JSON.stringify(v);
  }
  return v === null ? 'nil' : String(v);
}

/** A key as the printers spell it: `:name`, or `:"…"` when it is not keyword-safe. */
export function printKey(key: string): string {
  return KeywordPattern.test(key) ? `:${key}` : `:${JSON.stringify(key)}`;
}

/**
 * A quoted string from `i` (at its `"`) in JSON's grammar, raw line breaks (CR, LF, CRLF: all `\n`) and tabs allowed:
 * the string and the offset past it, or undefined when it is unterminated or malformed.
 */
export function readQuoted(text: string, i: number, end: number): [string, number] | undefined {
  let j = i + 1;
  while (j < end && text[j] !== '"') j += text[j] === '\\' ? 2 : 1;
  if (j >= end) return undefined;
  const raw = text
    .slice(i, j + 1)
    .replace(/\r\n?/g, '\n')
    .replace(/\n/g, '\\n')
    .replace(/\t/g, '\\t');
  try {
    return [JSON.parse(raw) as string, j + 1];
  } catch {
    return undefined;
  }
}
