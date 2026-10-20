/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashString } from '@molstar/core/data/util/hash-functions';
import { describe, expect, it } from 'vitest';
import { guessElementSymbolString } from '../util.js';

// Non-ASCII whitespace and letters whose upper case is longer (kept as code points: the formatter writes escapes as
// the characters themselves)
const NBSP = String.fromCharCode(0xa0);
const IDEOGRAPHIC_SPACE = String.fromCharCode(0x3000);
const LINE_SEPARATOR = String.fromCharCode(0x2028);
const SHARP_S = String.fromCharCode(0xdf);
const LIGATURE_FF = String.fromCharCode(0xfb00);
const E_ACUTE = String.fromCharCode(0xe9);

describe('guessElementSymbolString', () => {
  it('guesses from atom names', () => {
    const cases: [string, string][] = [
      ['CA', 'ALA'],
      ['CA', 'CA'],
      ['N', 'ALA'],
      ['HB1', 'ALA'],
      ['1HB', 'ALA'],
      [' OW ', 'SOL'],
      ['HW1', 'SOL'],
      ["O5'", 'DA'],
      ['na', 'NA'],
      ['Cl', 'CL'],
      ['fe2', 'HEM'],
      ['SOD', 'SOD'],
      ['POT', 'POT'],
      ['CES', 'CES'],
      ['CAL', 'CAL'],
      ['CLA', 'CLA'],
      ['SOD', 'NA'],
      ['cla', 'CLA'],
      ['ZN', 'ZN'],
      ['MG', 'MG'],
      ['', 'X'],
      ['123', 'X'],
      ['  ', 'X'],
      ['\tC1\n', 'X'],
      ['C 1', 'X'],
      ['x', 'X'],
      ['1x2', 'X'],
      ['B12', 'X'],
      ['F', 'X'],
      ['KX', 'X'],
      [`${NBSP}O${NBSP}`, 'X'],
      [`${IDEOGRAPHIC_SPACE}N1`, 'X'],
      [SHARP_S, 'X'],
      [`${LIGATURE_FF}1`, 'X'],
      [E_ACUTE, 'X'],
      [`${LINE_SEPARATOR}S2`, 'X'],
    ];
    expect(cases.map(([atom, comp]) => guessElementSymbolString(atom, comp))).toEqual([
      'C',
      'C',
      'N',
      'H',
      'H',
      'O',
      'H',
      'O',
      'NA',
      'CL',
      'FE',
      'NA',
      'K',
      'CS',
      'CA',
      'CL',
      'S',
      'CL',
      '',
      '',
      '',
      '',
      '',
      'C',
      'C',
      'X',
      'X',
      'B',
      'F',
      '',
      'O',
      'N',
      'S',
      'F',
      E_ACUTE.toUpperCase(),
      'S',
    ]);
  });

  it('keeps its guesses for seeded names', () => {
    let s = 11;
    const next = () => (s = (Math.imul(s, 1103515245) + 12345) >>> 0) >>> 16;
    const letters = ` \t0123456789ACHNOPSBFKLZacnos${NBSP}${IDEOGRAPHIC_SPACE}${SHARP_S}${E_ACUTE}`;
    const comps = ['ALA', 'SOD', 'CLA', 'POT', 'CES', 'CAL', 'HOH'];
    const out: string[] = [];
    for (let k = 0; k < 5000; k++) {
      let atom = '';
      for (let i = next() % 6; i > 0; i--) atom += letters[next() % letters.length];
      out.push(
        guessElementSymbolString(atom, next() % 3 === 0 ? atom.trim().toUpperCase() : comps[next() % comps.length]),
      );
    }
    expect(hashString(out.join('|'))).toBe(1805547253);
  });
});
