/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch } from '@molstar/core/repr/patch';
import { type Node, printAddress } from '@molstar/core/repr/tree';
import { addRefPatch, createSeqSpec } from '@molstar/seqx-core/lang';
import { plasmidData, type PlasmidVariant, scaleData, scaleStrand } from './data.js';

// The Plasmid's views over pUC19c (INSDC L09137.2): the ring with a linear detail following its window (EC-108,
// EC-114), over the rotated variant, where bla and a primer site cross the origin on both strands (EC-9 to EC-14); the
// same with bookmarks on the two features that cross it; the rotated linear view on the published record; and the ring
// at scale. Features are arrows by strand: the plus strand's point to the end, the minus strand's to the start (the
// ranges layer's `arrow` channel, EC-232), on the piece holding the end they point to, across the origin too.

export const title = 'Plasmid ring';

/** The arrow of each strand's features: where a feature on it points. */
export const StrandArrow = { '+': 'end', '-': 'start' } as const;

/** The ring and its detail; `rotation` turns the ring (degrees clockwise: where the origin is drawn). */
export function build(options: { readonly rotation?: number; readonly variant?: PlasmidVariant } = {}): Node {
  const b = createSeqSpec({ title });
  const { seq, plus, minus, sites } = plasmidData(b, options.variant ?? 'rotated');
  b.views.view('main', (v) => {
    const ring = v.section('ring', {
      of: seq,
      layout: 'rings',
      ...(options.rotation !== undefined ? { rotation: options.rotation } : {}),
    });
    ring
      .track('plus', { label: 'Features +' })
      .layer('features', { type: 'ranges', data: plus, arrow: StrandArrow['+'] });
    ring
      .track('minus', { label: 'Features −' })
      .layer('features', { type: 'ranges', data: minus, arrow: StrandArrow['-'] });
    ring.track('sites', { label: 'Sites', height: 36 }).layer('sites', { type: 'ranges', data: sites() });
    ring.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
    const detail = v.section('detail', { of: seq, follow: ring });
    detail
      .track('plus', { label: 'Features +' })
      .layer('features', { type: 'ranges', data: plus, arrow: StrandArrow['+'] });
    detail
      .track('minus', { label: 'Features −' })
      .layer('features', { type: 'ranges', data: minus, arrow: StrandArrow['-'] });
    detail.track('letters', { label: 'pUC19c' }).layer('letters', { type: 'letters' });
    detail.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
  });
  return b.build();
}

/** The address of a feature by its id (an `item` of the features' ranges). */
export const featureAddress = (features: string, id: string) =>
  `${features}${printAddress([{ kind: 'item', params: { id } }])}`;

/** The two features that cross the rotated origin, each a bookmark that focuses it by its shortest arc. */
export const CrossingBookmarks = [
  { id: 'bla', label: 'bla (AmpR)', on: '/data/minus' },
  { id: 'primer', label: 'primer site', on: '/data/plus' },
] as const;

/** The ring of the rotated variant with bookmarks on bla and the primer site (EC-9 to EC-14). */
export function crossing(): Node {
  let root = build();
  for (const m of CrossingBookmarks) {
    const r = applyPatch(root, addRefPatch(root, { id: m.id, to: featureAddress(m.on, m.id), label: m.label }));
    if (!r.ok) throw new Error(`The bookmark ${m.id} did not apply: ${r.issues.map((i) => i.code).join(', ')}`);
    root = r.root;
  }
  return root;
}

/**
 * The rotated linear view (Benchling's "set origin"), on the published record: a tracks section whose `range` runs
 * over the origin, from inside bla to past lacZα (two pieces, the seam no break), its features cut nowhere. Its start
 * is the one param the user sets.
 */
export function rotatedLinear(range: readonly [number, number] = [2000, 700]): Node {
  const b = createSeqSpec({ title: 'Plasmid, the rotated linear view' });
  const { seq, plus, minus, sites } = plasmidData(b, 'published');
  b.views.view('main', (v) => {
    const s = v.section('linear', { of: seq, range: [...range] });
    s.track('plus', { label: 'Features +' }).layer('features', { type: 'ranges', data: plus, arrow: StrandArrow['+'] });
    s.track('minus', { label: 'Features −' }).layer('features', {
      type: 'ranges',
      data: minus,
      arrow: StrandArrow['-'],
    });
    s.track('sites', { label: 'Sites', height: 36 }).layer('sites', { type: 'ranges', data: sites() });
    s.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
  });
  return b.build();
}

/** The height of a strand's track at scale: room for the lanes of a thousand features a strand. */
const ScaleBand = 120;

/** The ring at scale: a generated 200 kb cyclic sequence with 2,000 features on both strands, a tenth across the origin. */

export function scale(): Node {
  const b = createSeqSpec({ title: 'Plasmid ring at scale' });
  const data = scaleData();
  const seq = b.data.sequence('plasmid', {
    alphabet: 'dna',
    letters: data.letters,
    cyclic: true,
    name: 'generated, 200 kb',
  });
  const plus = scaleStrand(b, seq, data, '+', 'plus');
  const minus = scaleStrand(b, seq, data, '-', 'minus');
  b.views.view('main', (v) => {
    const ring = v.section('ring', { of: seq, layout: 'rings' });
    ring
      .track('plus', { label: 'Features +', height: ScaleBand })
      .layer('features', { type: 'ranges', data: plus, arrow: StrandArrow['+'] });
    ring
      .track('minus', { label: 'Features −', height: ScaleBand })
      .layer('features', { type: 'ranges', data: minus, arrow: StrandArrow['-'] });
    ring.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
    const detail = v.section('detail', { of: seq, follow: ring });
    detail
      .track('plus', { label: 'Features +', height: ScaleBand })
      .layer('features', { type: 'ranges', data: plus, arrow: StrandArrow['+'] });
    detail
      .track('minus', { label: 'Features −', height: ScaleBand })
      .layer('features', { type: 'ranges', data: minus, arrow: StrandArrow['-'] });
    detail.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
  });
  return b.build();
}
