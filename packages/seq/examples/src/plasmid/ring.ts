/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { plasmidData, RotatedOrigin } from './data.js';

// The Plasmid as a ring with a linear detail following its window (EC-108, EC-114), over the rotated variant: bla and
// a primer site cross the origin on both strands (EC-9 to EC-14). The rings outside in: the plus strand's features, the
// minus strand's, the restriction sites, the ruler; labels outside with leaders, the name and length in the hole. The
// detail shows the window's stretch with its letters; a pan or a zoom of the detail moves the ring's window.

export const title = 'Plasmid ring';

/** The ring and its detail; `rotation` turns the ring (degrees clockwise: where the origin is drawn). */
export function build(options: { readonly rotation?: number } = {}): Node {
  const b = createSeqSpec({ title });
  const { seq, plus, minus, sites } = plasmidData(b, RotatedOrigin);
  b.views.view('main', (v) => {
    const ring = v.section('ring', {
      of: seq,
      layout: 'rings',
      ...(options.rotation !== undefined ? { rotation: options.rotation } : {}),
    });
    ring.track('plus', { label: 'Features +' }).layer('features', { type: 'ranges', data: plus });
    ring.track('minus', { label: 'Features −' }).layer('features', { type: 'ranges', data: minus });
    ring.track('sites', { label: 'Sites', height: 12 }).layer('sites', { type: 'ranges', data: sites });
    ring.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
    const detail = v.section('detail', { of: seq, follow: ring });
    detail.track('plus', { label: 'Features +' }).layer('features', { type: 'ranges', data: plus });
    detail.track('minus', { label: 'Features −' }).layer('features', { type: 'ranges', data: minus });
    detail.track('letters', { label: 'pUC19' }).layer('letters', { type: 'letters' });
    detail.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
  });
  return b.build();
}

/**
 * The rotated linear view (Benchling's "set origin"), on the published record: a tracks section whose `range` runs
 * over the origin, from inside bla to past lacZα (two pieces, the seam no break), its features cut nowhere.
 */
export function rotatedLinear(): Node {
  const b = createSeqSpec({ title: 'Plasmid, the rotated linear view' });
  const { seq, plus, minus, sites } = plasmidData(b);
  b.views.view('main', (v) => {
    const s = v.section('linear', { of: seq, range: [2000, 700] });
    s.track('plus', { label: 'Features +' }).layer('features', { type: 'ranges', data: plus });
    s.track('minus', { label: 'Features −' }).layer('features', { type: 'ranges', data: minus });
    s.track('sites', { label: 'Sites', height: 12 }).layer('sites', { type: 'ranges', data: sites });
    s.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
  });
  return b.build();
}
