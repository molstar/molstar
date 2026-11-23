/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewer } from '@molstar/seqx-viewer';
import type { RondoLay } from './lay.js';

// A viewer's ribbons and a stem's partners (_spec/README.md, Hits, emphasis, words: The stem menu): one list, by
// chromosome, read by the stem menu and the partner keys.

/** The ribbons' lays of a viewer's last frame. */
export function laysOf(viewer: SeqViewer): RondoLay[] {
  const out: RondoLay[] = [];
  for (const s of viewer.frame?.layout.sections ?? [])
    for (const l of s.sectionLays) if (l.type.name === 'rondo:ribbons') out.push(l.lay as RondoLay);
  return out;
}

/** A partner stem of a stem: the count between them (the width's measure), its links, the ribbon's strand. */
export interface Partner {
  readonly stem: number;
  readonly count: number;
  readonly links: number;
  readonly strand: number;
}

/** A stem's partners by count, strongest first: each the partner stem and the count between them. */
export function partnersOf(lay: RondoLay, stem: number): Partner[] {
  const out: Partner[] = [];
  for (let r = 0; r < lay.nr; r++) {
    const other = lay.ribA[r] === stem ? lay.ribB[r] : lay.ribB[r] === stem ? lay.ribA[r] : -1;
    if (other >= 0) out.push({ stem: other, count: lay.ribCount[r], links: lay.ribLinks[r], strand: lay.ribStrand[r] });
  }
  return out.sort((a, b) => b.count - a.count || a.stem - b.stem);
}

/**
 * A partner chromosome of a stem: its system (an index into the lay's), its links in all, and its pieces (the ribbons
 * to its stems, by count; the stem's own links within it a piece of the stem itself on its own chromosome).
 */
export interface PartnerChromosome {
  readonly system: number;
  readonly links: number;
  /** The stem's own chromosome: links within it, never a partner. */
  readonly within: boolean;
  readonly pieces: readonly Partner[];
}

/**
 * A stem's partners by chromosome (_spec/README.md, The stem menu): a row a partner chromosome by its links, most
 * first (ties by the assembly's order), each its pieces by count; the stem's own chromosome last ("within"), its pieces
 * the ribbons to its other stems and, with links inside the stem, the stem itself. At most `max` partner rows and `max`
 * pieces a row (the within row besides): bounded however many stems a ring has. The stem menu lists them; the partner
 * keys walk their pieces in this order.
 */
export function partnerChromosomes(lay: RondoLay, stem: number, max = 8): PartnerChromosome[] {
  const own = lay.stemSystem[stem];
  const bySystem = new Map<number, Partner[]>();
  for (const p of partnersOf(lay, stem)) {
    const sys = lay.stemSystem[p.stem];
    let list = bySystem.get(sys);
    if (!list) bySystem.set(sys, (list = []));
    list.push(p);
  }
  const within = lay.withinLinks[stem];
  if (within > 0) {
    const list = bySystem.get(own) ?? [];
    list.push({ stem, count: lay.stemWithin[stem], links: within, strand: 1 });
    list.sort((a, b) => b.count - a.count || a.stem - b.stem);
    bySystem.set(own, list);
  }
  const rows: PartnerChromosome[] = [];
  let mine: PartnerChromosome | undefined;
  for (const [system, pieces] of bySystem) {
    let links = 0;
    for (const p of pieces) links += p.links;
    const row = { system, links, within: system === own, pieces: pieces.slice(0, max) };
    if (row.within) mine = row;
    else rows.push(row);
  }
  rows.sort((a, b) => b.links - a.links || a.system - b.system);
  if (rows.length > max) rows.length = max;
  if (mine) rows.push(mine);
  return rows;
}
