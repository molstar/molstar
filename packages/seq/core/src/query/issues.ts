/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The query forms' issue codes (_spec/README.md "Queries", Issues). Bounds, cancellation and deadlines use the
// tree's codes (`repr.tree.address.*`), as a selection does.

import { defineIssues } from '@molstar/core/repr/tree';

export const QueryIssues = defineIssues('seq.query', {
  'seq.query.form.prefix': {
    severity: 'error',
    description:
      'A form applies to what its prefix names: location forms to systems, item forms to annotations, node forms to nodes',
  },
  'seq.query.form.param': {
    severity: 'error',
    description: 'A form names a param it does not take, lacks one it needs, or gives one a value of the wrong type',
  },
  'seq.query.form.type': {
    severity: 'error',
    description:
      'A form is given another result type than it takes (a location, elements, nodes), or a combinator mixes types',
  },
  'seq.query.form.unsupported': {
    severity: 'error',
    description: 'A combination the forms do not define: the complement of a rect, or of nodes outside a node',
  },
  'seq.query.numbering.unknown': {
    severity: 'error',
    description: 'A numbering the system does not have: index, ordinal, its own or its molecule’s',
  },
  'seq.query.numbering.pending': {
    severity: 'info',
    description: 'A numbering the system declares whose labels have not arrived (a source’s): the form waits for them',
  },
  'seq.query.label.syntax': {
    severity: 'error',
    description: 'A label is not an integer or a word, or an insertion code is not letters',
    hint: 'an integer (94), an integer with a suffix (94A, 100.1), or a name',
  },
  'seq.query.label.int32': {
    severity: 'error',
    description: 'A label’s integer is outside the 32-bit range: refused, never wrapped',
  },
  'seq.query.label.unknown': {
    severity: 'error',
    description: 'A single label the numbering lacks (a range’s missing end snaps inward instead)',
  },
  'seq.query.label.insertion': {
    severity: 'error',
    description: 'An insertion code on a range, or on a label that already has a suffix',
    hint: 'write the code in the label (94A); a range takes labels with their codes',
  },
  'seq.query.label.snapped': {
    severity: 'info',
    description: 'A range end the numbering lacks snapped inward to the nearest label it has',
  },
  'seq.query.range.backwards': {
    severity: 'error',
    description: 'A range ends before it starts on a system that is not cyclic',
    hint: 'write :from at or before :to, in the numbering’s order of positions',
  },
  'seq.query.range.lap': {
    severity: 'error',
    description: 'A wrapping range on a cyclic system covers a whole lap or more',
    hint: 'a whole system is (all)',
  },
  'seq.query.range.empty': {
    severity: 'warning',
    description: 'A range holds no position once its ends are snapped to labels the numbering has',
  },
  'seq.query.range.wrap': {
    severity: 'error',
    description: 'A rect’s side wraps past a cyclic system’s origin: rect sides never wrap',
    hint: 'split the rect in two at the origin',
  },
  'seq.query.position.invalid': {
    severity: 'error',
    description:
      'Data gives a position a location cannot hold: not a 32-bit integer, or an extent ending before it starts',
  },
  'seq.query.location.unreachable': {
    severity: 'warning',
    description: 'A location is on a system the items are not on, and projection is not built yet: it overlaps nothing',
  },
  'seq.query.ref.unresolved': {
    severity: 'error',
    description: 'An operand names a ref, and refs are not resolved yet',
    hint: 'write the address the ref names',
  },
  'seq.query.row.unknown': {
    severity: 'error',
    description: 'No row of the alignment has the key',
    hint: 'a row key is its label',
  },
  'seq.query.locus.name': {
    severity: 'error',
    description: 'A locus names no sequence of its assembly',
    hint: 'chr17:7,668,402-7,687,550',
  },
  'seq.query.row.pending': {
    severity: 'info',
    description: 'The rows of an alignment read from a source have not arrived: the form waits for them',
  },
  'seq.query.property.level': {
    severity: 'error',
    description:
      'A test reads a property finer than the elements it tests (an atom’s on a residue, a residue’s in a unit test)',
    hint: 'test a finer level: (elements …) with an atom field, or groups’ :atom-test',
  },
  'seq.query.property.absent': {
    severity: 'warning',
    description:
      'Elements lack a property a test reads (a chain on a lone molecule, auth-seq with no auth numbering, atoms on a polymer): their test is false, the rest kept',
  },
  'seq.query.property.pending': {
    severity: 'info',
    description: 'A property’s values have not arrived (a values channel read by region): the form waits for them',
  },
  'seq.query.extension.missing': {
    severity: 'error',
    description: 'A namespaced form or property no registered extension provides: a query is never approximated',
    hint: 'register the extension that provides it (its query-form or query-property entry)',
  },
  'seq.query.node.missing': {
    severity: 'warning',
    description: 'No node (a dynamic node, a chain) is held under the prefix with the key',
  },
});

export type QueryIssueCode = keyof typeof QueryIssues.codes;
