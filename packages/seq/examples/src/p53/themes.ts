/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn, type Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { addDomains, conservation, P53, seeded } from './data.js';

// p53's themes: the protein's letters under each scheme that fits a protein; its domains colored by their category
// (the builder lists the categories, so the legend is the vocabulary); a per-residue confidence as AlphaFold's four
// pLDDT bands (binned) and as a ramp; two signals naming one shared scale (one domain, one legend); variants colored by
// consequence and sized by a log-scaled frequency; secondary structure by its scheme. Its export draws the legends
// below the data. The confidence, the second signal, the variants and the secondary structure are generated from seeds
// (stand-ins, named so): the view's point is the themes.

export const title = 'p53 themes';

/** The letters' schemes that fit a protein without column measures, in the picker's order. */
export const LetterSchemes = ['classes', 'zappo', 'taylor', 'hydrophobicity', 'helix', 'strand', 'turn', 'buried'];

/** A per-residue confidence in [20, 98], smooth (a stand-in for a model's pLDDT). */
export function confidence(): Float32Array {
  const r = seeded(4637);
  let x = 70;
  return Float32Array.from({ length: P53.length }, (_, i) => {
    // the disordered ends low, the DNA-binding domain high, as a structure model's confidence runs
    const target = i < 90 || i > 360 ? 40 : 88;
    x += (target - x) * 0.08 + (r() - 0.5) * 10;
    return Math.round(Math.max(20, Math.min(98, x)));
  });
}

/** A second signal on the scale the conservation names (its range larger: the union is the shared domain). */
export function accessibility(): Float32Array {
  const r = seeded(17);
  return Float32Array.from({ length: P53.length }, () => Math.round(r() * 150) / 100);
}

/** Variants: a position each, a consequence and a frequency (stand-ins). */
export function variants() {
  const r = seeded(273);
  const n = 40;
  const start = Int32Array.from({ length: n }, () => Math.floor(r() * P53.length)).sort();
  const kinds = ['missense', 'nonsense', 'frameshift', 'synonymous'];
  return {
    start,
    end: start.map((s) => s + 1),
    consequence: categoryColumn(Array.from({ length: n }, () => kinds[Math.floor(r() ** 2 * kinds.length)])),
    frequency: Float32Array.from({ length: n }, () => 10 ** (-6 + r() * 5)),
  };
}

/** Secondary structure: helices, strands and turns over the DNA-binding domain (stand-ins). */
export const Secondary = {
  start: new Int32Array([110, 124, 141, 156, 177, 194, 214, 230, 251, 264, 278]),
  end: new Int32Array([118, 130, 148, 164, 181, 200, 220, 236, 258, 271, 287]),
  structure: categoryColumn(['strand', 'turn', 'strand', 'strand', 'helix', 'strand', 'strand', 'turn', 'strand', 'strand', 'helix']),
};

/** The view's spec. */
export function build(): Node {
  const b = createSeqSpec({ title });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const domains = addDomains(b, p53);
  const plddt = b.data.values('confidence', { target: p53, values: confidence() });
  const cons = b.data.values('conservation', { target: p53, values: conservation() });
  const acc = b.data.values('accessibility', { target: p53, values: accessibility() });
  const v = variants();
  const vars = b.data.ranges('variants', { target: p53, rows: v });
  const ss = b.data.ranges('secondary', { target: p53, rows: Secondary });
  const shared = b.views.scale('signal', { label: 'Signal', domain: [0, null] as unknown as [number, number] });
  b.views.view('main', { export: { legend: 'below' } }, (view) => {
    const s = view.section('protein', { of: p53 });
    for (const name of LetterSchemes)
      s.track(`letters-${name}`, { label: name, height: 18 }).layer('letters', {
        type: 'letters',
        ...(name === 'classes' ? {} : { color: theme.scheme(name) }),
      });
    s.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
    s.track('confidence', { label: 'Confidence (pLDDT)' }).layer('values', {
      type: 'values',
      data: plddt,
      domain: [0, 100],
      color: theme.scheme('plddt'),
    });
    s.track('confidence-ramp', { label: 'Confidence (ramp)' }).layer('values', {
      type: 'values',
      data: plddt,
      domain: [0, 100],
      color: theme.field('value', { domain: [0, 100] }),
    });
    s.track('conservation', { label: 'Signal: conservation' }).layer('values', {
      type: 'values',
      data: cons,
      scale: shared,
      color: theme.field('value', { scale: shared }),
    });
    s.track('accessibility', { label: 'Signal: accessibility' }).layer('values', {
      type: 'values',
      data: acc,
      scale: shared,
      color: theme.field('value', { scale: shared }),
    });
    s.track('variants', { label: 'Variants' }).layer('ranges', {
      type: 'ranges',
      data: vars,
      color: theme.field('consequence'),
      size: theme.field('frequency', { type: 'log', domain: [1e-6, 0.1], range: [0.3, 1] }),
    });
    s.track('secondary', { label: 'Secondary structure' }).layer('ranges', {
      type: 'ranges',
      data: ss,
      color: theme.scheme('secondary', { field: 'structure' }),
    });
  });
  return b.build();
}
