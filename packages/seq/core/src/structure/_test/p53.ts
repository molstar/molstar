/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A structure shaped as 2OCJ (p53's DNA-binding domain, not yet recorded): P04637 and the entity it maps to by a
// `structure` map (residues 94 to 312, author numbering equal to UniProt's), `copies` copies of the domain with
// unobserved ends, a zinc per chain, and optionally an `alignment`-like map of the same pair. Synthetic letters.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

export const Up = '/data/p53';
export const Ent = '/data/ent';
export const Zn = '/data/zn';
export const Complex = '/data/2ocj';
export const UpLength = 393;
/** The entity: UniProt 94 to 312. */
export const EntStart = 93;
export const EntLength = 219;
/** Each chain's modeled residues on the entity (`[2, 196)`: 96 to 289). */
export const Observed: readonly [number, number] = [2, 196];

const Letters = 'ACDEFGHIKLMNPQRSTVWY';
const letters = (n: number, seed: number) => Array.from({ length: n }, (_, i) => Letters[(i * 7 + seed) % 20]).join('');

export interface P53Options {
  readonly copies?: number;
  /** A second map of the same pair, one residue off, with this priority (absent: none). */
  readonly aln?: number;
  /** The sifts map's priority (absent: none written). */
  readonly sifts?: number;
  /** The units in reverse order (order irrelevant). */
  readonly reversed?: boolean;
  /** A numbering `auth` on the entity too (one past the units' labels), which the units' own beat (EC-8). */
  readonly entityAuth?: boolean;
  /** A sequence `/data/x` with a map onto unit A alone (its own map: EC-25's longer legal path). */
  readonly perUnit?: boolean;
  /** The units' author numbering (true by default; a scale case leaves it out). */
  readonly auth?: boolean;
  /**
   * A view of P04637: its letters, its domains (a ranges node, colored by `category` when `field`) and a values track
   * (`values[i]` at `i`, NaN at 100), each a layer `colorRuns` reads.
   */
  readonly view?: { readonly field?: boolean };
}

/** The view's values: `i / UpLength`, none at 100. */
export const ValueAt = (i: number) => (i === 100 ? Number.NaN : i / UpLength);

/** Each unit's author labels: UniProt's numbers where modeled, `104A` at entity position 10 (an insertion code). */
const AuthLabels = Array.from({ length: EntLength }, (_, i) =>
  i < Observed[0] || i >= Observed[1] ? null : i === 10 ? `${EntStart + 1 + i}A` : String(EntStart + 1 + i),
);

export const chainId = (k: number) => (k < 26 ? String.fromCharCode(65 + k) : `X${k}`);

export function p53Spec(o: P53Options = {}): Node {
  const copies = o.copies ?? 4;
  const b = createSeqSpec({ title: '2OCJ-shaped' });
  const up = b.data.sequence('p53', { alphabet: 'protein', letters: letters(UpLength, 3) });
  const entLetters = letters(UpLength, 3).slice(EntStart, EntStart + EntLength);
  const ent = b.data.sequence('ent', { alphabet: 'protein', letters: entLetters }, (s) => {
    if (o.entityAuth)
      s.numbering('auth', {
        method: 'labels',
        labels: Array.from({ length: EntLength }, (_, i) => String(EntStart + 2 + i)),
      });
  });
  const zn = b.data.molecule('zn', {
    class: 'ion',
    residues: ['ZN'],
    atoms: { residue: [0], name: ['ZN'], element: ['ZN'] },
  });
  b.data.complex('2ocj', { name: '2OCJ' }, (c) => {
    const order = Array.from({ length: copies }, (_, k) => k);
    if (o.reversed) order.reverse();
    for (const k of order) {
      const chain = chainId(k);
      c.unit(
        chain,
        { of: ent, chain, auth_chain: chain, observed: { start: [Observed[0]], end: [Observed[1]] } },
        (u) => {
          if (o.auth !== false) u.numbering('auth', { method: 'labels', labels: AuthLabels });
        },
      );
    }
    for (const k of order) c.unit(`Z${chainId(k)}`, { of: zn, chain: `Z${chainId(k)}`, auth_chain: chainId(k) });
  });
  b.data.map('sifts', {
    from: up,
    to: ent,
    method: 'structure',
    segments: { from: [EntStart], to: [0], length: [EntLength], provenance: ['sifts'] },
    ...(o.sifts !== undefined ? { priority: o.sifts } : {}),
  });
  if (o.perUnit) {
    const x = b.data.sequence('x', { alphabet: 'protein', letters: letters(50, 5) });
    b.data.map('xa', { from: x, to: `${Complex}/A` as never, method: 'offset', by: 10 });
  }
  if (o.aln !== undefined)
    b.data.map('aln', {
      from: up,
      to: ent,
      method: 'segments',
      segments: { from: [EntStart + 1], to: [0], length: [EntLength - 1] },
      priority: o.aln,
    });
  if (o.view) {
    const domains = b.data.ranges('domains', {
      target: up,
      rows: {
        start: [0, 63, 99, 320],
        end: [61, 94, 120, 360],
        id: ['TAD', 'PRR', 'DBD', 'TET'],
        category: ['activation', 'linker', 'binding', 'binding'],
      },
    });
    const values = b.data.values('scores', {
      target: up,
      values: Array.from({ length: UpLength }, (_, i) => ValueAt(i)),
    });
    b.views.view('main', (v) => {
      const s = v.section('p53', { of: up });
      s.track('letters').layer('letters', { type: 'letters' });
      s.track('domains').layer('domains', {
        type: 'ranges',
        data: domains,
        ...(o.view?.field ? { color: { kind: 'field', params: { name: 'category' } } } : { color: '#ff0000' }),
      });
      s.track('scores').layer('scores', { type: 'values', data: values });
    });
  }
  return b.build();
}
