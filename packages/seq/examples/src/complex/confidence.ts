/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { ofChain } from './blocks.js';
import type { Derived } from './derived.js';
import { constructOffset, crystal, Data, Lengths, prediction } from './data.js';

// Both instances on the prediction's chains (UniProt P24941, P20248): AlphaFold DB's per-residue confidence (pLDDT) in
// AlphaFold's four bins, UniProt's secondary structure and 1FIN's (its first copy's chains, mapped through an offset
// map: CDK2's numbering is UniProt's, the cyclin A construct starts at P20248's 173), and what 1FIN holds of each chain:
// the crystal's chains as a range, gray where the crystal has no residue (cyclin A2's first 172, EC-118).

export const title = 'Complex: confidence and secondary structure';

/** AlphaFold's confidence bins: very low under 50, low to 70, confident to 90, very high above. */
export const PlddtEdges = [0, 50, 70, 90, 100] as const;
export const PlddtColors = ['#ff7d45', '#ffdb13', '#65cbf3', '#0053d6'] as const;

/** The confidence's color: AlphaFold's bins (an editing view sets other edges, one theme patch). */
export const plddtColor = (edges: readonly number[] = PlddtEdges, colors: readonly string[] = PlddtColors) =>
  theme.field('value', { palette: theme.binned({ edges: [...edges], colors: [...colors] }) });

/** The gray of a residue the crystal does not hold. */
export const Unobserved = '#9e9e9e';

export function build(d: Derived = Data): Node {
  const b = createSeqSpec({ title });
  const p = prediction(b, d);
  const x = crystal(b, d);
  const shift = constructOffset(d);
  b.data.map('cdk2-in-1fin', { from: x.cdk2, to: p.cdk2.protein, method: 'offset', by: 0 });
  b.data.map('cyclin-in-1fin', { from: x.cyclin, to: p.cyclin.protein, method: 'offset', by: shift });
  // what 1FIN holds of each chain, on UniProt's numbering: CDK2 whole, the cyclin A construct from 173
  const held = {
    cdk2: b.data.ranges('cdk2-held', {
      target: p.cdk2.protein,
      rows: { start: [0], end: [Lengths.cdk2], id: ['1FIN A, C'], label: ['1FIN chains A, C'] },
    }),
    cyclin: b.data.ranges('cyclin-held', {
      target: p.cyclin.protein,
      rows: {
        start: [shift],
        end: [shift + Lengths.construct],
        id: ['1FIN B, D'],
        label: [`1FIN chains B, D (${shift + 1}–${shift + Lengths.construct})`],
      },
    }),
  };
  // what it does not hold: cyclin A2's residues before the construct
  const unobserved = {
    cdk2: b.data.ranges('cdk2-unobserved', { target: p.cdk2.protein, rows: { start: [], end: [] } }),
    cyclin: b.data.ranges('cyclin-unobserved', {
      target: p.cyclin.protein,
      rows: { start: [0], end: [shift], id: ['not in 1FIN'], label: [`not in 1FIN (1–${shift})`] },
    }),
  };
  b.views.view('main', (v) => {
    for (const [id, chain, crystalChain, label] of [
      ['cdk2', 'cdk2', 'A', 'CDK2 (P24941)'],
      ['cyclin', 'cyclin', 'B', 'cyclin A2 (P20248)'],
    ] as const) {
      const s = v.section(id, { of: p[chain].protein });
      s.track('ruler', { label }).layer('ruler', { type: 'ruler' });
      s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
      s.track('plddt', { label: 'pLDDT (AF v6)' }).layer('plddt', {
        type: 'values',
        data: p[chain].plddt,
        color: plddtColor(),
      });
      s.track('uniprot-ss', { label: 'Secondary (UniProt)' }, (t) => {
        t.layer('helix', { type: 'secondary', data: p[chain].data('Helix') });
        t.layer('strand', { type: 'secondary', data: p[chain].data('Beta strand') });
      });
      s.track('crystal-ss', { label: `Secondary (1FIN ${crystalChain})` }).layer('ss', {
        type: 'secondary',
        data: ofChain(chain === 'cdk2' ? x.cdk2Secondary : x.cyclinSecondary, crystalChain),
      });
      s.track('held', { label: 'In 1FIN' }, (t) => {
        t.layer('held', { type: 'ranges', data: held[chain] });
        // gray where 1FIN holds no residue of the chain
        t.layer('unobserved', { type: 'ranges', data: unobserved[chain], color: Unobserved });
      });
    }
  });
  return b.build();
}
