/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { hashString } from '@molstar/core/data/util/hash-functions';
import { describe, expect, it } from 'vitest';
import { type AlignmentOptions, align } from '../alignment.js';

function generator(seed: number): () => number {
  let s = seed;
  return () => {
    s = (Math.imul(s, 1103515245) + 12345) >>> 0;
    return s >>> 16;
  };
}

/** Seeded pairs (lengths 0 to 24, small alphabets, so ties are common) aligned under `options`, as one string. */
function seededResults(seed: number, letters: readonly string[], options: Partial<AlignmentOptions>): string {
  const next = generator(seed);
  const seq = () => Array.from({ length: next() % 25 }, () => letters[next() % letters.length]);
  const out: string[] = [];
  for (let k = 0; k < 60; k++) {
    const { aliA, aliB, score } = align(seq(), seq(), options);
    out.push(`${aliA.join('')}/${aliB.join('')}/${score}`);
  }
  return out.join(' ');
}

describe('Alignment', () => {
  it('basic', () => {
    // 3PQR: Rhodopsin
    const seqA =
      'MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA';
    // 3SN6: Endolysin,Beta-2 adrenergic
    const seqB =
      'DYKDDDDAENLYFQGNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYAADEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGAAHILTKTWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFKYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRATHQEAINCYAEETCCDFFTNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQKIDKSEGRFHVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNLIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNTGEQSG';

    const { aliA, aliB, score } = align(seqA, seqB, {
      gapPenalty: -11,
      gapExtensionPenalty: -1,
      substMatrix: 'blosum62',
    });

    expect(aliA.join('')).toEqual(
      '------------------------------------------------------------------------------------------------------------------------------------------------MNGTEGPNFYVPFSNKTGVVRSPFEA---PQYYLAEPWQFSM--LAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLH---GYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMS-NFRFGENHAIMGVAFTWVMA-LACAAPPLVGWSRYI-PEGMQC----SCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLV----------------FTVKEAAAQQQESATTQ----------KAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA',
    );
    expect(aliB.join('')).toEqual(
      'DYKDDDDAENLYFQGNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKS-RWYNQTPNRAKRVITTFRTGTWDAYAADEVWVVGMGIVMSLIVLAIVFG---NVLVITAIAKFERLQTVTNYFITSLACADLVM---GLAVVPFGAAHILTKTWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFKYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRATHQEAINCYAEETC-CDFFT------NQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQKIDKSEGRFHVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMG-TFTLCWLPFF-IVNIVHVIQDNLIRKEVYILLNWIGYVNSGFNPLIY-CRSPDFRIAFQELLCLRRSSL--KAYGNGYSSNGNTGEQSG',
    );
    expect(score).toEqual(118);
  });

  // Pinned from the number[][] implementation that the shared DP (@molstar/core/data/util/needleman-wunsch) replaced:
  // the same alignments, ties included, and scores (perf-audit 04/P1).
  it('keeps its alignments and scores: seeded pairs under each option set', () => {
    const amino = 'ARNDCQEGHILKMFPSTWYVU'.split(''); // U: not in BLOSUM62, scores 0
    const names = ['ALA', 'GLY', 'SER', 'HOH', 'DA'];
    const cases: [string, readonly string[], Partial<AlignmentOptions>][] = [
      ['default AB', ['A', 'B'], {}],
      ['default names', names, {}],
      ['blosum62', amino, { substMatrix: 'blosum62' }],
      ['blosum62x', amino, { substMatrix: 'blosum62x' }],
      ['linear gaps', ['A', 'B', 'C'], { gapPenalty: -5, gapExtensionPenalty: -5 }],
      ['extension over open', ['A', 'B', 'C'], { gapPenalty: -2, gapExtensionPenalty: -4 }],
      ['fractional', ['A', 'C', 'G', 'T'], { gapPenalty: -10, gapExtensionPenalty: -0.5 }],
      ['free gaps', ['A', 'B'], { gapPenalty: 0, gapExtensionPenalty: 0 }],
    ];
    const digests = Object.fromEntries(
      cases.map(([name, letters, options], k) => [name, hashString(seededResults(17 + k, letters, options))]),
    );
    expect(digests).toEqual({
      'default AB': 1160935526,
      'default names': 596161677,
      blosum62: 1818237639,
      blosum62x: 1727981868,
      'linear gaps': 632748018,
      'extension over open': -1537212691,
      fractional: -1089611043,
      'free gaps': -1866580488,
    });
  });
});
