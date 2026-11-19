/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, type SeqSpecBuilder, theme } from '@molstar/seqx-core/lang';
import type { Derived } from './derived.js';
import { complexPae, constructOffset, contactRows, Data, prediction } from './data.js';

// The Complex's prediction: CDK2 and cyclin A2 as UniProt has them (P24941, P20248), a section a chain, and the
// synthetic predicted aligned error between them (`complexPae`'s CDK2 × cyclin A block: low across 1FIN's interface,
// high where cyclin A2 is disordered) as a grid section whose x follows CDK2's section and y cyclin A's; 1FIN's contacts
// (both copies, shifted onto P20248 by the construct's offset) as dots on it.

export const title = 'Complex: PAE across the chains (synthetic)';

/** The PAE's color: viridis over 0 to 30 Å, low (confident) dark. */
export const PaeColor = () =>
  theme.field('value', { domain: [0, 30], palette: theme.continuous({ colors: 'viridis' }) });

/** The prediction's data nodes: the entries, the confidence, the PAE's four blocks and the contacts on P20248. */
export function predictionData(b: SeqSpecBuilder, d: Derived = Data) {
  const p = prediction(b, d);
  const pae = complexPae();
  const grid = (x: 'cdk2' | 'cyclin', y: 'cdk2' | 'cyclin') =>
    b.data.grid(`pae-${x}-${y}`, { x: p[x].protein, y: p[y].protein, values: pae.block(x, y) });
  const contacts = b.data.pairs('contacts', {
    a: p.cdk2.protein,
    b: p.cyclin.protein,
    rows: contactRows({ shift: constructOffset(d) }, d),
  });
  return { ...p, grid, contacts };
}

export function build(d: Derived = Data): Node {
  const b = createSeqSpec({ title });
  const p = predictionData(b, d);
  const pae = p.grid('cdk2', 'cyclin');
  b.views.view('main', (v) => {
    const x = v.section('cdk2', { of: p.cdk2.protein });
    x.track('ruler', { label: 'CDK2' }).layer('ruler', { type: 'ruler' });
    x.track('letters', { label: 'P24941' }).layer('letters', { type: 'letters' });
    const y = v.section('cyclin', { of: p.cyclin.protein });
    y.track('ruler', { label: 'cyclin A2' }).layer('ruler', { type: 'ruler' });
    y.track('letters', { label: 'P20248' }).layer('letters', { type: 'letters' });
    v.section('pae', {
      layout: 'grid',
      of: p.cdk2.protein,
      of_y: p.cyclin.protein,
      follow: x.address,
      follow_y: y.address,
    }).track('pae', { label: 'PAE (Å, synthetic)' }, (t) => {
      // the least error a sample holds when zoomed out: a confident pair is never averaged away
      t.layer('pae', { type: 'grid', data: pae, aggregate: 'min', color: PaeColor() });
      // the contacts in white, outlined: on viridis's dark low end (the interface) in both themes
      t.layer('contacts', { type: 'pairs', data: p.contacts, color: '#ffffff' });
    });
  });
  return b.build();
}
