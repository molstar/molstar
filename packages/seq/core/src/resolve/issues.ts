/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The resolve pass's codes (_spec/README.md "Resolve", Issues). What no map covers maps nowhere and is no issue
// (EC-21); a map that pairs nothing is.

import { defineIssues } from '@molstar/core/repr/tree';

export const ResolveIssues = defineIssues('seq.resolve', {
  'seq.resolve.ref.missing': {
    severity: 'error',
    description: 'A reference the resolve pass follows names no node of a kind it takes (the spec was not checked)',
    hint: 'Check the spec first: checkSeqSpec reports what it names',
  },
  'seq.resolve.source.runtime': {
    severity: 'warning',
    description:
      'No format runtime the application gave reads the source’s format; the source and its data are unavailable',
    hint: 'Give a runtime for the format when composing the viewer',
  },
  'seq.resolve.map.empty': {
    severity: 'warning',
    description: 'A map pairs no position of its systems',
    hint: 'Check its offset or exons against the lengths of the sequences',
  },
  'seq.resolve.map.codon_length': {
    severity: 'warning',
    description:
      'A codon map’s exons are not whole codons, or not as many as its protein has residues (with or without a stop)',
    hint: 'A coding sequence is three nucleotides a residue; it cannot start mid-codon',
  },
  'seq.resolve.path.none': {
    severity: 'warning',
    description: 'No path of maps joins two systems within the bound (at most four maps and one codon step)',
    hint: 'Add a map between the systems, or between them and a system both reach',
  },
  'seq.resolve.path.ambiguous': {
    severity: 'warning',
    description: 'Several shortest paths join two systems and map them differently; none is taken',
    hint: 'Two maps join the same systems: drop one, or map through one system',
  },
});
