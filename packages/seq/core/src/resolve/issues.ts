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
  'seq.resolve.ref.system': {
    severity: 'warning',
    description:
      'A reference names a sequence, molecule or unit that made no system (a unit of a refused complex, or of nothing held); what reads it is left out',
    hint: 'Name a system the model holds: a unit of a complex within `maxUnits`, whose `of` is held',
  },
  'seq.resolve.source.format': {
    severity: 'error',
    description: 'A source names no format (the spec was not checked); it is unavailable',
    hint: 'Give the source a `format`',
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
  'seq.resolve.map.observed': {
    severity: 'warning',
    description:
      'A structure map onto a unit says other residues are unobserved than the unit’s `observed` does; the unit’s are taken',
    hint: 'Unobserved is the structure’s own fact: drop the map’s `unobserved` rows, or read the unit again',
  },
  'seq.resolve.map.codon_length': {
    severity: 'warning',
    description:
      'A codon map’s exons are not whole codons, or not as many as its protein has residues (with or without a stop)',
    hint: 'A coding sequence is three nucleotides a residue; it cannot start mid-codon',
  },
  'seq.resolve.map.table': {
    severity: 'warning',
    description: 'A codon map names a genetic code that is not built; its codons have no status',
    hint: 'Codes built: 1 (standard), 2 (vertebrate mitochondrial), 11 (bacterial and plastid)',
  },
  'seq.resolve.length.unknown': {
    severity: 'warning',
    description:
      'A row or an exon wraps the origin of a cyclic sequence whose length is not known yet; it is left out until it is',
    hint: 'Give the sequence its letters inline, or wait for its source',
  },
  'seq.resolve.system.unnamed': {
    severity: 'error',
    description:
      'A sequence read from a source names no system the opened source holds: its length is unknown and nothing over it is laid out',
    hint: 'Name the sequence (`name`, else its id) after a system its source holds',
  },
  'seq.resolve.grid.size': {
    severity: 'warning',
    description: 'A grid’s values are not one per pair of positions of its systems, or more than a block holds',
    hint: 'Give |x| × |y| values, x-major, at most 2^24 of them; the grid is not drawn',
  },
  'seq.resolve.path.none': {
    severity: 'warning',
    description: 'No path of maps joins two systems within the bound (at most four maps and one codon step)',
    hint: 'Add a map between the systems, or between them and a system both reach',
  },
  'seq.resolve.path.ambiguous': {
    severity: 'warning',
    description: 'Several shortest paths join two systems and map them differently; none is taken',
    hint: 'Two maps join the same systems: rank one above the other with the map’s `priority`, or name one with the section’s `via`',
  },
  'seq.resolve.via.unused': {
    severity: 'warning',
    description:
      'A section’s `via` names a map that joins no pair on any of its paths (or no map); its paths are chosen as if the via were absent',
    hint: 'Name a map between systems the section’s paths go through, or remove it from `via`',
  },
  'seq.resolve.path.limit': {
    severity: 'warning',
    description: 'More shortest paths join two systems than the search compares; none is taken',
    hint: 'Map through fewer systems, or raise the candidate bound in the configuration',
  },
  'seq.resolve.node.invalid': {
    severity: 'warning',
    description: 'A dynamic node a source sent is not a valid node of its kind; it is left out',
  },
  'seq.resolve.node.ref': {
    severity: 'warning',
    description:
      'A dynamic node names what is neither a spec node nor a held dynamic node of a kind it takes; it is left out',
  },
  'seq.resolve.item.strand': {
    severity: 'warning',
    description: 'The rows of an item with a codon map are on both strands; the first row’s strand is taken',
  },
  'seq.resolve.complex.units': {
    severity: 'warning',
    description: 'A complex has more units than the configuration allows (`maxUnits`); it is left out',
    hint: 'Read a subset of the assembly (fewer operators or chains), or raise `maxUnits` at its memory',
  },
  'seq.resolve.genome.load': {
    severity: 'error',
    description: 'The genome chunk did not load',
    hint: 'Reload',
  },
  'seq.resolve.node.depth': {
    severity: 'warning',
    description: 'A nested source is deeper than the nesting the configuration allows (EC-179); it is not requested',
    hint: 'Raise the nesting depth in the configuration, or show the nodes of a shallower source',
  },
});
