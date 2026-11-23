/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewer } from '@molstar/seqx-viewer';
import type { RondoLay } from './lay.js';

// A viewer's ribbons and a stem's partners (_spec/README.md, The viewer): what the partner keys and the stem menu read.

/** The ribbons' lays of a viewer's last frame. */
export function laysOf(viewer: SeqViewer): RondoLay[] {
  const out: RondoLay[] = [];
  for (const s of viewer.frame?.layout.sections ?? [])
    for (const l of s.sectionLays) if (l.type.name === 'rondo:ribbons') out.push(l.lay as RondoLay);
  return out;
}

/** A stem's partners by count, strongest first: each the partner stem and the count between them. */
export function partnersOf(lay: RondoLay, stem: number): { stem: number; count: number; strand: number }[] {
  const out: { stem: number; count: number; strand: number }[] = [];
  for (let r = 0; r < lay.nr; r++) {
    if (lay.ribA[r] === stem) out.push({ stem: lay.ribB[r], count: lay.ribCount[r], strand: lay.ribStrand[r] });
    else if (lay.ribB[r] === stem) out.push({ stem: lay.ribA[r], count: lay.ribCount[r], strand: lay.ribStrand[r] });
  }
  return out.sort((a, b) => b.count - a.count || a.stem - b.stem);
}
