/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';

/** The matrix's states system and the section's positions system. */
export const States = '/data/gt/(states)';
export const Positions = '/data/chr';
export const StateCount = 8;

/** A genotype matrix of `StateCount` states over 40 positions, a track of 160 px: its states' rows 20 px each. */
export function matrixSpec() {
  const b = createSeqSpec({ title: 'matrix' });
  const chr = b.data.sequence('chr', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
  const codes = new Uint8Array(40 * StateCount).map((_, i) => i % 3);
  b.data.matrix('gt', { target: chr, codes, states: Array.from({ length: StateCount }, (_, i) => `s${i}`) } as never);
  b.views.view('main', (v) => {
    v.section('genotypes', { of: chr }, (s) => {
      s.track('gt', { height: 160 }).layer('gt', { type: 'matrix', data: '/data/gt' } as never);
    });
  });
  return b.build();
}
