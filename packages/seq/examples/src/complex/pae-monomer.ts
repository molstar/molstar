/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { Cdk2PaeDigits, Cdk2PaeSize } from './cdk2-pae-data.js';
import { prediction } from './data.js';
import { PaeColor } from './pae.js';

// CDK2's own PAE as AlphaFold DB publishes it (AF-P24941-F1, model version 6): real values, the whole chain, in the
// colors of AlphaFold DB's plot, so it compares with https://alphafold.ebi.ac.uk/entry/P24941 side by side. The cross-chain
// PAE of the complex is not in AlphaFold DB (no permissive prediction of the complex, none in ModelArchive was found
// 2026-11-21): the complex's grids stay synthetic and say so; nothing here mixes the two.

export const title = 'CDK2: AlphaFold DB’s PAE';

/** The PAE's values as a grid holds them (the file's order: a row the aligned residue). */
export function cdk2Pae(): Float32Array {
  const n = Cdk2PaeSize;
  const out = new Float32Array(n * n);
  for (let i = 0; i < out.length; i++) out[i] = Number.parseInt(Cdk2PaeDigits[i], 32);
  return out;
}

export function build(): Node {
  const b = createSeqSpec({ title });
  const p = prediction(b).cdk2;
  const pae = b.data.grid('pae-cdk2', { x: p.protein, y: p.protein, values: cdk2Pae() });
  b.views.view('main', (v) => {
    const x = v.section('cdk2', { of: p.protein });
    x.track('ruler', { label: 'CDK2' }).layer('ruler', { type: 'ruler' });
    x.track('confidence', { label: 'pLDDT' }).layer('plddt', { type: 'values', data: p.plddt });
    x.track('letters', { label: 'P24941' }).layer('letters', { type: 'letters' });
    v.section('pae', { layout: 'grid', of: p.protein, follow: x.address }).track(
      'pae',
      { label: 'PAE (Å, AlphaFold DB v6)' },
      (t) => {
        t.layer('pae', { type: 'grid', data: pae, aggregate: 'min', color: PaeColor() });
      },
    );
  });
  return b.build();
}
