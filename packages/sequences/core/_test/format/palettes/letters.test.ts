/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  type Color,
  countColumn,
  ElementSymbolColors,
  LetterDictionaries,
  LetterDictionaryNames,
  letterColors,
} from '../../../src/format/index.js';
import { clustalXColor, createColumnComposition } from '../../../src/format/palettes/clustal-x.js';

const BLUE = 0x80a0f0;
const RED = 0xf01505;
const MAGENTA = 0xc048c0;
const GREEN = 0x15c015;
const PINK = 0xf08080;
const ORANGE = 0xf09048;
const YELLOW = 0xc0c000;
const CYAN = 0x15a4a4;

describe('letter dictionaries', () => {
  it("are the spec's, keyed by upper-case one-letter code", () => {
    expect(LetterDictionaryNames).toEqual([
      'Zappo',
      'Taylor',
      'Hydrophobicity',
      'HelixPropensity',
      'StrandPropensity',
      'TurnPropensity',
      'BuriedIndex',
      'Nucleotide',
      'PurinePyrimidine',
    ]);
    for (const name of LetterDictionaryNames) {
      for (const key of Object.keys(LetterDictionaries[name])) expect(key, name).toMatch(/^[A-Z]$/);
    }
    for (const name of LetterDictionaryNames.slice(0, 7)) {
      expect(Object.keys(LetterDictionaries[name]).length, name).toBeGreaterThanOrEqual(20);
    }
    expect(LetterDictionaries.Zappo.A).toBe(0xffafaf);
    expect(LetterDictionaries.Zappo.X).toBeUndefined();
    expect(LetterDictionaries.Nucleotide.U).toBe(LetterDictionaries.Nucleotide.T);
    expect(LetterDictionaries.PurinePyrimidine.G).toBe(LetterDictionaries.PurinePyrimidine.A);
  });

  it("letter's defaults: Zappo, Nucleotide, ElementSymbol with carbon", () => {
    expect(letterColors(null, 'protein')).toBe(LetterDictionaries.Zappo);
    expect(letterColors(null, 'dna')).toBe(LetterDictionaries.Nucleotide);
    expect(letterColors(null, 'rna')).toBe(LetterDictionaries.Nucleotide);
    expect(letterColors(null, null)).toBe(ElementSymbolColors);
    expect(letterColors(null, 'generic')).toBeUndefined();
    expect(letterColors('ElementSymbol', null)?.C).toBe(0x909090);
    expect(letterColors('Taylor', 'protein')).toBe(LetterDictionaries.Taylor);
    expect(letterColors('ResidueName', 'protein')?.ALA).toBe(0x8cff8c);
    expect(letterColors('Nope', 'protein')).toBeUndefined();
    expect(letterColors({ A: 'red', C: '#00f/#fff' }, 'dna')).toEqual({ A: 0xff0000, C: 0x0000ff });
  });
});

describe('Clustal X', () => {
  const c = createColumnComposition();
  const colors = (column: string, letters: string): (Color | undefined)[] => {
    countColumn(column, c);
    return [...letters].map((l) => clustalXColor(c, l));
  };

  it('hydrophobic columns are blue, gaps left out of the shares', () => {
    expect(colors('LLLLIV----', 'LIVA')).toEqual([BLUE, BLUE, BLUE, BLUE]);
    expect(colors('lllliv....', 'lI')).toEqual([BLUE, BLUE]);
  });

  it('thresholds are strict', () => {
    // K + R is exactly 60%: no rule holds.
    expect(colors('KKKRRRDEST', 'KRDES')).toEqual([undefined, undefined, undefined, undefined, undefined]);
    // 70%: K and R red; E and D magenta by {>60%, KR}.
    expect(colors('KKKKRRRDES', 'KRDES')).toEqual([RED, RED, MAGENTA, MAGENTA, undefined]);
  });

  it('charged and polar groups', () => {
    expect(colors('EEEEEDDKST', 'EDKST')).toEqual([MAGENTA, MAGENTA, undefined, undefined, undefined]);
    expect(colors('NNNNNNDDQQ', 'NDQ')).toEqual([GREEN, undefined, undefined]);
    expect(colors('QQQQQEEKLL', 'QEKL')).toEqual([GREEN, MAGENTA, undefined, undefined]);
    expect(colors('SSSSTTTAGG', 'STAG')).toEqual([GREEN, GREEN, undefined, ORANGE]);
    expect(colors('SSSSSSSSSY', 'SY')).toEqual([GREEN, undefined]);
  });

  it('cysteines pink above 85%, blue in a hydrophobic column', () => {
    expect(colors('CCCCCCCCCA', 'CA')).toEqual([PINK, BLUE]);
    expect(colors('CCCCCCLLLL', 'C')).toEqual([BLUE]);
    expect(colors('CCDDEEKKNN', 'C')).toEqual([undefined]);
  });

  it('glycine and proline always, aromatic in hydrophobic columns', () => {
    expect(colors('GPN-------', 'GPN')).toEqual([ORANGE, YELLOW, undefined]);
    expect(colors('HHHHHHHHHY', 'HY')).toEqual([CYAN, CYAN]);
    expect(colors('HYDDEEKKNN', 'HY')).toEqual([undefined, undefined]);
  });

  it('other letters and empty columns are not colored', () => {
    expect(colors('XXXXXXBZLL', 'XBZL')).toEqual([undefined, undefined, undefined, undefined]);
    expect(colors('----', 'A')).toEqual([undefined]);
  });
});
