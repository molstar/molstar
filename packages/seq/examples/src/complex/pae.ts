/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { ComplexLengths, letters } from './blocks.js';

// The Complex's synthetic predicted aligned error and contacts, stand-ins until the dataset is recorded (the fixtures
// hold no coordinates, so no contact is derived): CDK2 and cyclin A at 1FIN's lengths; a synthetic PAE across the two
// chains (x CDK2's residues, y cyclin A's: an error in Å a pair, low across the interface, high elsewhere, from a seed)
// and synthetic contacts at the interface (CDK2's PSTAIRE helix and activation loop against cyclin A's first cyclin
// box, random positions there from a seed). One section a chain; the PAE a grid section whose x follows CDK2's section
// and y cyclin A's; the contacts dots on the PAE and arcs between the chains' blocks.

export const title = 'Complex PAE (synthetic)';

/** The interface's patches: CDK2's ranges against cyclin A's (0-based, half-open). */
const Interface: readonly (readonly [number, number, number, number])[] = [
  [44, 56, 60, 95],
  [150, 165, 95, 125],
];

/** How far a residue is inside a range (1 inside, falling off over `edge` residues outside, 0 past). */
function near(i: number, lo: number, hi: number, edge = 10): number {
  const d = i < lo ? lo - i : i >= hi ? i - hi + 1 : 0;
  return Math.max(0, 1 - d / edge);
}

/** The synthetic PAE, CDK2 × cyclin A, x-major, Å: about 25 to 30 across the chains, 2 to 6 at the interface. */
export function complexPae(): Float32Array {
  const nx = ComplexLengths.cdk2;
  const ny = ComplexLengths.cyclin;
  const r = new PCG(1979);
  const v = new Float32Array(nx * ny);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++) {
      let at = 0;
      for (const [a0, a1, b0, b1] of Interface) at = Math.max(at, near(i, a0, a1) * near(j, b0, b1));
      const noise = (r.int() % 1000) / 1000;
      v[i * ny + j] = 2 + 23 * (1 - at) + 4 * noise;
    }
  return v;
}

/** The synthetic contacts at the interface: random pairs within its patches. */
export function complexContacts() {
  const r = new PCG(1996);
  const rows = { a_start: [] as number[], a_end: [] as number[], b_start: [] as number[], b_end: [] as number[] };
  const add = (a0: number, a1: number, b0: number, b1: number, n: number) => {
    for (let k = 0; k < n; k++) {
      const a = a0 + (r.int() % (a1 - a0));
      const b = b0 + (r.int() % (b1 - b0));
      rows.a_start.push(a);
      rows.a_end.push(a + 1);
      rows.b_start.push(b);
      rows.b_end.push(b + 1);
    }
  };
  add(...Interface[0], 14);
  add(...Interface[1], 10);
  return { ...rows, type: rows.a_start.map(() => 'contact'), id: rows.a_start.map((_, k) => `c${k + 1}`) };
}

export interface PaeOptions {
  /** The contacts as arcs between the chains' blocks, without the PAE. */
  readonly contacts?: boolean;
}

export function build(options: PaeOptions = {}): Node {
  const b = createSeqSpec({ title });
  const cdk2 = b.data.sequence('cdk2', { alphabet: 'protein', letters: letters(2, ComplexLengths.cdk2) });
  const cyclin = b.data.sequence('cyclin', { alphabet: 'protein', letters: letters(3, ComplexLengths.cyclin) });
  const contacts = b.data.pairs('contacts', { a: cdk2, b: cyclin, rows: complexContacts() });
  const pae = b.data.grid('pae', { x: cdk2, y: cyclin, values: complexPae() });
  b.views.view('main', (v) => {
    if (options.contacts) {
      const s = v.section('complex', {});
      s.block('a', { of: cdk2, label: 'CDK2 (A)' });
      s.block('b', { of: cyclin, label: 'cyclin A (B)' });
      s.track('contacts', { label: 'Contacts (synthetic)' }).layer('arcs', { type: 'pairs', data: contacts });
      s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
      return;
    }
    const x = v.section('cdk2', { of: cdk2 });
    x.track('ruler', { label: 'CDK2' }).layer('ruler', { type: 'ruler' });
    x.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
    const y = v.section('cyclin', { of: cyclin });
    y.track('ruler', { label: 'cyclin A' }).layer('ruler', { type: 'ruler' });
    y.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
    v.section('pae', {
      layout: 'grid',
      of: cdk2,
      of_y: cyclin,
      follow: x.address,
      follow_y: y.address,
    }).track('pae', { label: 'PAE (Å, synthetic)' }, (t) => {
      // the least error a sample holds when zoomed out: a confident pair is never averaged away
      t.layer('pae', {
        type: 'grid',
        data: pae,
        aggregate: 'min',
        color: theme.field('value', { domain: [0, 30], palette: theme.continuous({ colors: 'viridis' }) }),
      });
      // the contacts in white, outlined: on viridis's dark low end (the interface) in both themes
      t.layer('contacts', { type: 'pairs', data: contacts, color: '#ffffff' });
    });
  });
  return b.build();
}

/** The contacts as arcs between the chains' blocks. */
export const contacts = () => build({ contacts: true });
