/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { addDomains, addSecondary, addVariants, conservation, P53, seeded } from './data.js';

// p53 for the track groups' stories (seqx-core _spec "Render", Track groups): the protein section in three groups, as
// the design's Defaults draw it. Features (domains, secondary structure) open; Variants (two variant tracks, a density
// as its authored summary); Signals (conservation, a pLDDT-like and a disorder-like signal). `folds` says which groups
// open folded.

/** A signal a residue in [0, 1], smooth (a seeded walk). */
function signal(seed: number): Float32Array {
  const r = seeded(seed);
  const v = new Float32Array(P53.length);
  let x = 0.5;
  for (let i = 0; i < v.length; i++) {
    x = Math.min(1, Math.max(0, x + (r() - 0.5) * 0.2));
    v[i] = Math.round(x * 100) / 100;
  }
  return v;
}

/** The variants' density: a variant's position counted with its neighbours (a window of 9). */
function density(starts: readonly number[]): Float32Array {
  const v = new Float32Array(P53.length);
  for (const s of starts) for (let i = Math.max(0, s - 4); i <= Math.min(v.length - 1, s + 4); i++) v[i] += 1;
  return v;
}

/** The section in three groups; `folds`: the groups that open folded (`collapsed`). */
export function build(folds: readonly ('features' | 'variants' | 'signals')[] = ['variants', 'signals']): Node {
  const b = createSeqSpec({ title: 'p53, track groups' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const domains = addDomains(b, p53);
  const secondary = addSecondary(b, p53);
  const variants = addVariants(b, p53);
  const rare = b.data.ranges('rare', {
    target: p53,
    rows: {
      start: new Int32Array([47, 151, 272, 337]),
      end: new Int32Array([48, 152, 273, 338]),
      label: ['P47S', 'P152L', 'V272M', 'R337H'],
    },
  });
  const counts = b.data.values('variant-density', {
    target: p53,
    values: density([47, 151, 175, 245, 248, 272, 273, 337]),
  });
  const cons = b.data.values('conservation', { target: p53, values: conservation() });
  const plddt = b.data.values('plddt', { target: p53, values: signal(7) });
  const disorder = b.data.values('disorder', { target: p53, values: signal(11) });
  const collapsed = (id: string) => (folds.includes(id as never) ? { collapsed: true } : {});
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    const f = s.group('features', { label: 'Features', ...collapsed('features') });
    f.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
    f.track('secondary', { label: 'Secondary' }).layer('secondary', { type: 'secondary', data: secondary });
    const va = s.group('variants', { label: 'Variants', summary: 'density', ...collapsed('variants') });
    va.track('clinvar', { label: 'ClinVar' }).layer('variants', { type: 'variants', data: variants });
    va.track('rare', { label: 'Rare' }).layer('variants', { type: 'variants', data: rare });
    va.track('density', { label: 'Density' }).layer('values', { type: 'values', data: counts });
    const sg = s.group('signals', { label: 'Signals', ...collapsed('signals') });
    sg.track('conservation', { label: 'Conservation' }).layer('values', { type: 'values', data: cons });
    sg.track('plddt', { label: 'pLDDT' }).layer('values', { type: 'values', data: plddt });
    sg.track('disorder', { label: 'Disorder' }).layer('values', { type: 'values', data: disorder });
  });
  return b.build();
}

/** Every group open. */
export const open = (): Node => build([]);

/** Every group folded. */
export const folded = (): Node => build(['features', 'variants', 'signals']);
