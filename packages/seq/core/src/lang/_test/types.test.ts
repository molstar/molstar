/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Type tests: `@ts-expect-error` lines fail the type check (`tsc -p tsconfig.test.json`) when what they mark compiles.
// The functions are never called.

import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../builder.js';

const wrongKinds = [
  () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    const r = b.data.ranges('r', { target: s, rows: { start: [0], end: [1] } });
    // @ts-expect-error a ranges handle where a sequence is taken
    b.data.ranges('r2', { target: r, rows: { start: [0], end: [1] } });
    // @ts-expect-error a sequence handle where a source is taken
    b.data.sequence('s2', { alphabet: 'dna', source: s });
    // @ts-expect-error a ranges handle where a sequence is taken
    b.data.map('m', { from: r, to: s, method: 'offset', by: 1 });
    b.views.view('v', (v) => {
      const sec = v.section('a', { of: s });
      // @ts-expect-error a section where a sequence is taken
      v.section('b', { of: sec });
      // @ts-expect-error a sequence handle where a section is taken
      v.section('c', { of: s, follow: s });
      // @ts-expect-error a ranges handle where a section is taken
      v.section('d', { of: s, follow: r });
      v.section('e', { of: s, follow: sec });
      // @ts-expect-error a source's handle where the kinds of a layer's data are taken
      sec.track('t').layer('ranges', { type: 'ranges', data: b.data.source('x', { format: 'fasta' }) });
    });
  },
  () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MK' });
    const r = b.data.ranges('r', { target: s, rows: { start: [0], end: [1] } });
    const c = b.data.complex('c', {}, (cx) => {
      // @ts-expect-error a ranges handle where a sequence or molecule is taken
      cx.unit('A', { of: r });
      // @ts-expect-error links name no `a`: their parent is both ends
      cx.links('l', { a: s, rows: { a_start: [0], a_end: [1], b_start: [1], b_end: [2] } });
    });
    b.views.view('v', (v) =>
      v.section('a', {}, (sec) => {
        // @ts-expect-error a complex where a block's system is taken
        sec.block('b', { of: c });
        // @ts-expect-error not a form
        sec.block('b2', { of: s, form: 'chip' });
      }),
    );
    // @ts-expect-error not a way to draw modified residues
    b.views.view('w', (v) => v.section('a', { of: s }).track('t').layer('l', { type: 'letters', modified: 'dot' }));
  },
];

const fine = [
  () => {
    const b = createSeqSpec();
    // a handle of an allowed kind, or an address written out
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.data.ranges('r', { target: s, rows: { start: [0], end: [1] } });
    b.data.ranges('r2', { target: '/data/s', rows: { start: [0], end: [1] } });
    const address: string = `/data/${'s'}`;
    b.data.values('v', { target: address });
    const id: string = b.data.source('x', { format: 'fasta' });
    // molecules and units stand where a system does; a complex where a section's system does
    const m = b.data.molecule('m', { residues: ['ATP'] });
    const units: string[] = [];
    const c = b.data.complex('c', {}, (cx) => {
      units.push(cx.unit('A', { of: m }, (unit) => unit.numbering('auth', { method: 'labels', labels: ['1'] })));
    });
    b.data.ranges('r3', { target: m, rows: { start: [0], end: [1] } });
    b.views.view('v', (v) => v.section('a', { of: c, classes: ['polymer', 'ligand'] }));
    return id;
  },
];

const wrongParams = [
  () => {
    const b = createSeqSpec();
    // @ts-expect-error alphabet is required
    b.data.sequence('s', { letters: 'ACGT' });
    // @ts-expect-error not an alphabet
    b.data.sequence('s', { alphabet: 'cat', letters: 'ACGT' });
    // @ts-expect-error not a param of sequence
    b.data.sequence('s', { alphabet: 'dna', length: 4 });
    // @ts-expect-error a sequence has an id
    b.data.sequence({ alphabet: 'dna' });
  },
];

describe('the builder types', () => {
  it('are checked by tsc -p tsconfig.test.json', () => {
    expect(wrongKinds.length + fine.length + wrongParams.length).toBe(4);
  });
});
