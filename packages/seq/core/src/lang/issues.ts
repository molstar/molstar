/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `seq.lang.*` codes: what the description cannot express (_spec/README.md, Issues). Each with its severity, what
// it means and the usual next step, as the tree core's catalogs are; a test fails on a code a source uses and the
// catalog lacks.

import { defineIssues } from '@molstar/core/repr/tree';

export const LangIssues = defineIssues('seq.lang', {
  'seq.lang.version.unsupported': {
    severity: 'error',
    description: 'A spec is written in a newer major version of the language than this one reads',
    hint: 'Read it with a newer version of the library, or migrate the spec down',
  },
  'seq.lang.id.required': {
    severity: 'error',
    description: 'A node of the spec has no id (an extension kind’s too), so it cannot be addressed or patched',
    hint: 'Give the node an id unique among its siblings',
  },
  'seq.lang.container.id': {
    severity: 'error',
    description: 'A `data` container has an id other than `data`, or a `views` container one other than `views`',
    hint: 'The containers have reserved ids',
  },
  'seq.lang.source.index': {
    severity: 'error',
    description: 'A source with an inline `value` has an `index`, which only a file read by region has',
    hint: 'Drop the `index`, or give a `uri` instead of the `value`',
  },
  'seq.lang.source.content': {
    severity: 'error',
    description: 'A source has both a `uri` and an inline `value`, or neither',
    hint: 'Give one of `uri` and `value`',
  },
  'seq.lang.sequence.content': {
    severity: 'error',
    description: 'A sequence has both `letters` and a `source`, or neither',
    hint: 'Give one of `letters` and `source`',
  },
  'seq.lang.sequence.letters': {
    severity: 'error',
    description: 'Inline letters of a sequence are not letters of its alphabet',
    hint: 'Protein: letters, `*` and `-`; dna: ACGTN, the ambiguity codes and `-`; rna: the same with U',
  },
  'seq.lang.ranges.content': {
    severity: 'error',
    description: 'A ranges node has both `rows` and a `source`, or neither',
    hint: 'Give one of `rows` and `source`',
  },
  'seq.lang.ranges.interval': {
    severity: 'error',
    description: 'A row of ranges starts before position 0, or ends before it starts on a sequence that is not cyclic',
    hint: 'Ranges are 0-based and half-open: 0 <= start <= end; only a cyclic sequence has ranges that wrap',
  },
  'seq.lang.ranges.strand': {
    severity: 'error',
    description: 'A strand of a row of ranges is neither `+` nor `-`',
    hint: 'The strand column holds `+` or `-`',
  },
  'seq.lang.ranges.bounds': {
    severity: 'error',
    description: 'A row of ranges ends past the inline letters of its target sequence',
    hint: 'A range ends at most at the length of the sequence',
  },
  'seq.lang.values.content': {
    severity: 'error',
    description: 'A values node has both `values` and a `source`, or neither',
    hint: 'Give one of `values` and `source`',
  },
  'seq.lang.values.length': {
    severity: 'error',
    description: 'Inline values are not one per position of the inline letters of their target sequence',
    hint: 'Give as many values as the sequence has letters; NaN is missing',
  },
  'seq.lang.map.params': {
    severity: 'error',
    description: 'A map has params its method does not take, or lacks one it needs',
    hint: '`offset` takes `by`; `codon` takes `exons` and `strand`',
  },
  'seq.lang.map.alphabet': {
    severity: 'error',
    description: 'A codon map does not go from a protein to a dna or rna sequence',
    hint: 'A codon map is a protein on its gene',
  },
  'seq.lang.map.exons': {
    severity: 'error',
    description: 'An exon of a codon map ends before it starts, or overlaps or precedes the exon before it',
    hint: 'List the exons in genomic order, each from its start to its end, none overlapping',
  },
  'seq.lang.ref.target': {
    severity: 'error',
    description: 'A reference names no node of the spec, or one of a kind the param does not take',
    hint: 'Name a node by its address, like /data/p53',
  },
  'seq.lang.layer.data': {
    severity: 'error',
    description:
      'A layer’s data does not fit its type: `ranges` draws ranges, `values` draws values, `letters` draws letters',
    hint: 'Give data of the kind the type draws',
  },
  'seq.lang.section.range': {
    severity: 'error',
    description:
      'The range of a section ends before it starts on a sequence that is not cyclic, or runs past the inline letters of its system',
    hint: 'A range is [start, end) within the system',
  },
  'seq.lang.section.follow': {
    severity: 'error',
    description: 'A section follows itself, directly or through other sections',
    hint: 'Following is a chain that ends at a section that follows none',
  },
});
