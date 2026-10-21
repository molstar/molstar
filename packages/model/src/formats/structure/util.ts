/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { StringLike } from '@molstar/core/util/string-like';
import { TokenBuilder, type Tokens } from '@molstar/io/reader/common/text/tokenizer';

export function guessElementSymbolTokens(tokens: Tokens, str: StringLike, start: number, end: number) {
  let s = start,
    e = end - 1;

  // trim spaces and numbers
  let c = str.charCodeAt(s);
  while ((c === 32 || (c >= 48 && c <= 57)) && s <= e) c = str.charCodeAt(++s);
  c = str.charCodeAt(e);
  while ((c === 32 || (c >= 48 && c <= 57)) && e >= s) c = str.charCodeAt(--e);

  ++e;

  if (s === e) return TokenBuilder.add(tokens, s, e); // empty
  if (s + 1 === e) return TokenBuilder.add(tokens, s, e); // one char

  c = str.charCodeAt(s);

  if (s + 2 === e) {
    // two chars
    const c2 = str.charCodeAt(s + 1);
    if (
      ((c === 78 || c === 110) && (c2 === 65 || c2 === 97)) || // NA na Na nA
      ((c === 67 || c === 99) && (c2 === 76 || c2 === 108)) || // CL
      ((c === 70 || c === 102) && (c2 === 69 || c2 === 101)) || // FE
      ((c === 83 || c === 115) && (c2 === 73 || c2 === 105)) || // SI
      ((c === 66 || c === 98) && (c2 === 82 || c2 === 114)) || // BR
      ((c === 65 || c === 97) && (c2 === 83 || c2 === 115)) || // AS
      ((c === 76 || c === 108) && (c2 === 73 || c2 === 105)) || // LI
      ((c === 71 || c === 103) && (c2 === 65 || c2 === 97)) // GA
    )
      return TokenBuilder.add(tokens, s, s + 2);
  }

  if (
    c === 67 ||
    c === 99 || // C c
    c === 72 ||
    c === 104 || // H h
    c === 78 ||
    c === 110 || // N n
    c === 79 ||
    c === 111 || // O o
    c === 80 ||
    c === 112 || // P p
    c === 83 ||
    c === 115 // S s
  )
    return TokenBuilder.add(tokens, s, s + 1);

  TokenBuilder.add(tokens, s, s); // no reasonable guess, add empty token
}

const TwoCharElementNames = new Set(['NA', 'CL', 'FE', 'SI', 'BR', 'AS', 'LI', 'GA']);
const OneCharElementNames = new Set(['C', 'H', 'N', 'O', 'P', 'S', 'F', 'B']);

const reTrimSpacesAndNumbers = /^[\s\d]+|[\s\d]+$/g;
export function guessElementSymbolString(atomId: string, compId: string) {
  // ASCII names (all in practice) by char codes, without a regex or a new string per call (perf-audit 04/P11);
  // anything else as before (Unicode spaces, upper cases longer than the letter).
  let s = 0;
  let e = atomId.length;
  for (let i = 0; i < e; i++) {
    if (atomId.charCodeAt(i) > 127) return guessElementSymbolStringUnicode(atomId, compId);
  }
  // trim spaces (what `\s` matches in ASCII) and digits
  while (s < e && isSpaceOrDigit(atomId.charCodeAt(s))) s++;
  while (e > s && isSpaceOrDigit(atomId.charCodeAt(e - 1))) e--;
  const l = e - s;

  if (l === 0) return ''; // empty
  const c0 = upperCase(atomId.charCodeAt(s));
  if (l === 1) return AsciiChars[c0]; // one char
  if (l === 2) {
    const two = TwoCharElementCodes.get(c0 * 128 + upperCase(atomId.charCodeAt(s + 1)));
    if (two !== undefined) return two; // two chars
  }
  // check for Charmm ion names where component and atom id are the same
  if (l === 3 && compId.length === 3) {
    const c1 = upperCase(atomId.charCodeAt(s + 1));
    const c2 = upperCase(atomId.charCodeAt(s + 2));
    if (compId.charCodeAt(0) === c0 && compId.charCodeAt(1) === c1 && compId.charCodeAt(2) === c2) {
      if (compId === 'SOD') return 'NA';
      if (compId === 'POT') return 'K';
      if (compId === 'CES') return 'CS';
      if (compId === 'CAL') return 'CA';
      if (compId === 'CLA') return 'CL';
    }
  }

  if (OneCharElementNames.has(AsciiChars[c0])) return AsciiChars[c0];
  return ''; // no reasonable guess, return empty string
}

/** One-character strings by ASCII code. */
const AsciiChars = Array.from({ length: 128 }, (_, c) => String.fromCharCode(c));
/** TwoCharElementNames by `first * 128 + second` char code. */
const TwoCharElementCodes = new Map<number, string>();
for (const name of TwoCharElementNames) TwoCharElementCodes.set(name.charCodeAt(0) * 128 + name.charCodeAt(1), name);

/** `toUpperCase` of an ASCII code. */
function upperCase(c: number) {
  return c >= 97 && c <= 122 ? c & ~32 : c;
}

function isSpaceOrDigit(c: number) {
  return (c >= 48 && c <= 57) || c === 32 || (c >= 9 && c <= 13);
}

function guessElementSymbolStringUnicode(atomId: string, compId: string) {
  // trim spaces and numbers, convert to upper case
  atomId = atomId.replace(reTrimSpacesAndNumbers, '').toUpperCase();
  const l = atomId.length;

  if (l === 0) return atomId; // empty
  if (l === 1) return atomId; // one char
  if (TwoCharElementNames.has(atomId)) return atomId; // two chars

  // check for Charmm ion names where component and atom id are the same
  if (l === 3 && compId === atomId) {
    if (atomId === 'SOD') return 'NA';
    if (atomId === 'POT') return 'K';
    if (atomId === 'CES') return 'CS';
    if (atomId === 'CAL') return 'CA';
    if (atomId === 'CLA') return 'CL';
  }

  if (OneCharElementNames.has(atomId[0])) return atomId[0];

  return ''; // no reasonable guess, return empty string
}
