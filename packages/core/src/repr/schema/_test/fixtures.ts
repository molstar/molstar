/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A small language in the shape of the first Seq* slice's (language-description.md), for the schema's and the
// builder's tests: containers, ids, refs, an enum, a default, a typed column, a table, a slot, a data form.

import { defineLanguage } from '../index.js';

export const Mini = defineLanguage({
  id: 'mini',
  version: '0.0.0',
  root: 'spec',
  doc: 'A test language.',
  kinds: {
    spec: {
      doc: 'The root.',
      parents: [],
      params: {
        version: { type: 'string', required: true, doc: 'The version.' },
        title: { type: 'string', doc: 'A title.' },
      },
    },
    data: { doc: 'The data container.', parents: ['spec'], ids: 'required' },
    sequence: {
      doc: 'A sequence.',
      parents: ['data'],
      ids: 'required',
      params: {
        alphabet: { type: 'string', enum: ['protein', 'dna', 'rna'], required: true, doc: 'The alphabet.' },
        letters: { type: 'string', doc: 'Inline letters.' },
        cyclic: { type: 'boolean', default: false, doc: 'Whether it is cyclic.' },
      },
    },
    values: {
      doc: 'A value per position.',
      parents: ['data'],
      ids: 'required',
      params: {
        target: { type: 'ref', kinds: ['sequence'], required: true, doc: 'The sequence.' },
        values: { type: 'column', of: 'float32', doc: 'The values.' },
        scale: { type: 'number', min: 0, unit: 'A', doc: 'A length.' },
      },
    },
    ranges: {
      doc: 'Ranges on a sequence.',
      parents: ['data'],
      ids: 'required',
      params: {
        target: { type: 'ref', kinds: ['sequence'], required: true, doc: 'The sequence.' },
        rows: {
          type: 'table',
          required: true,
          key: 'id',
          columns: {
            start: { type: 'column', of: 'int32', required: true, doc: 'Starts, 0-based.' },
            end: { type: 'column', of: 'int32', required: true, doc: 'Ends, exclusive.' },
            id: { type: 'column', of: 'category', doc: 'Ids.' },
            label: { type: 'list', of: { type: 'string', doc: 'A label.' }, doc: 'Labels.' },
          },
          doc: 'The rows.',
        },
        color: { type: 'node', language: 'mini', kinds: ['field'], doc: 'A theme form.' },
      },
    },
    field: {
      doc: 'A theme form: color by a field.',
      parents: [],
      data: { args: [{ name: 'name', type: 'string', required: true, doc: 'The field.' }] },
    },
  },
});
