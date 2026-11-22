/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The genome chunk's `seq.resolve.assembly.*` codes (_spec/README.md "Resolve", Issues): a catalog of its own, so a
// reader that never opens a genome loads none of it (the load's own code, `genome.load`, is `ResolveIssues`').

import { defineIssues } from '@molstar/core/repr/tree';

export const GenomeIssues = defineIssues('seq.resolve.assembly', {
  'seq.resolve.assembly.lengths': {
    severity: 'error',
    description: 'An assembly’s sequence has no length',
    hint: 'Give lengths',
  },
  'seq.resolve.assembly.duplicate': {
    severity: 'warning',
    description: 'Two sequences of one name: the later dropped',
    hint: 'Rename one',
  },
  'seq.resolve.assembly.alias': {
    severity: 'warning',
    description: 'An alias names two sequences: dropped',
    hint: 'Give it to one',
  },
  'seq.resolve.assembly.override': {
    severity: 'warning',
    description: 'A sequence override names no sequence of the assembly: dropped',
    hint: 'Name a derived sequence, or an alias of one',
  },
  'seq.resolve.assembly.mismatch': {
    severity: 'warning',
    description: 'A source says another assembly',
    hint: 'Add an alias',
  },
  'seq.resolve.assembly.name': {
    severity: 'info',
    description: 'Rows name sequences the assembly lacks',
    hint: 'Add aliases',
  },
  'seq.resolve.assembly.inline': {
    severity: 'info',
    description: 'Inline pairs, grids, matrices or values over an assembly are not read yet',
    hint: 'Name a source, or target a sequence',
  },
});
