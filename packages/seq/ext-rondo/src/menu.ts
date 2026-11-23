/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { rangesLocation } from '@molstar/seqx-core/extension';
import type { HeaderElementAt, HeaderFacts, HeaderSubEntry, SeqViewer } from '@molstar/seqx-viewer';
import { Kind, type RondoLay } from './lay.js';
import { laysOf, partnersOf } from './partners.js';

// The stem menu's chunk (_spec/README.md, The viewer, The stem menu): what its entries (`StemMenu`, `viewer.ts`) read,
// loaded by `rondoViewer`'s behavior once composed (its `loaded` waits for it). Its entries read the lays of the viewer
// the menu is in (`HeaderFacts.viewer`): no state of their own, so two viewers of one spec each list their own.

/** How many partners the menu lists. */
const MenuPartners = 8;

/** The stem a layer's element key names in a viewer's lays (a fan is its stem), with its lay. */
export function stemOf(
  viewer: SeqViewer,
  layer: string,
  key: string | number,
): { lay: RondoLay; stem: number } | undefined {
  if (typeof key !== 'string' || !key.startsWith('stem:')) return undefined;
  for (const lay of laysOf(viewer)) {
    if (lay.key !== layer) continue;
    const i = lay.byKey.get(key);
    if (i === undefined) continue;
    const kind = lay.polyKind[i];
    if (kind === Kind.trumpet || kind === Kind.fan) return { lay, stem: lay.polyRef[i] };
  }
  return undefined;
}

/**
 * The stem at a viewer's keyboard cursor (its position on a placed stem with a trumpet, the topmost layer's), its menu's
 * target and its mouth two px inward (where the menu opens): what the menu key opens there (`HeaderEntry.atCursor`).
 */
export function stemAtCursor(viewer: SeqViewer): HeaderElementAt | undefined {
  const c = viewer.selection.cursor.value;
  if (!c || c.track !== undefined) return undefined;
  const lays = laysOf(viewer);
  // the last drawn is on top: the layer a pointer there would hit first
  for (let k = lays.length - 1; k >= 0; k--) {
    const lay = lays[k];
    if (lay.section !== c.section) continue;
    for (let s = 0; s < lay.ns; s++) {
      if (!lay.stemPlaced[s] || lay.stemPoly[s] < 0 || lay.systems[lay.stemSystem[s]] !== c.system) continue;
      if (c.position < lay.stemStart[s] || c.position >= lay.stemEnd[s]) continue;
      lay.frame.point((lay.stemT0[s] + lay.stemT1[s]) / 2, 2);
      return {
        target: {
          section: lay.section,
          track: lay.key.slice(0, lay.key.lastIndexOf('/')),
          element: { layer: lay.key, key: lay.stemKey[s], type: 'rondo:ribbons' },
        },
        x: lay.frame.px,
        y: lay.frame.py,
      };
    }
  }
  return undefined;
}

/** A menu's stem from its facts: in the viewer whose menu it is. */
export const stemOfFacts = (f: HeaderFacts) =>
  f.element?.type === 'rondo:ribbons' && f.viewer ? stemOf(f.viewer, f.element.layer, f.element.key) : undefined;

/** The names of two systems of a lay in the assembly's order (one when equal). */
function pairNames(lay: RondoLay, a: number, b: number): string[] {
  return a === b ? [lay.names[a]] : a < b ? [lay.names[a], lay.names[b]] : [lay.names[b], lay.names[a]];
}

/** The stem and its partners' stems selected (one selection); its message. */
export function selectStemPartners(viewer: SeqViewer, layer: string, key: string): { message: string } {
  const s = stemOf(viewer, layer, key);
  if (!s) return { message: 'No stem here' };
  const { lay, stem } = s;
  const loc = (k: number) => rangesLocation(lay.systems[lay.stemSystem[k]], [lay.stemStart[k], lay.stemEnd[k]]);
  viewer.selection.clear();
  viewer.selection.add(loc(stem));
  const partners = partnersOf(lay, stem);
  for (const p of partners) viewer.selection.add(loc(p.stem));
  return { message: `${partners.length} partner${partners.length === 1 ? '' : 's'} selected` };
}

/** "Show with a partner ▸": the stem's partners by count, each shown with it as a pair. */
export function partnerItems(f: HeaderFacts): HeaderSubEntry[] {
  const s = stemOfFacts(f);
  if (!s) return [];
  const { lay, stem } = s;
  const address = lay.prep?.opens ?? lay.section;
  return partnersOf(lay, stem)
    .slice(0, MenuPartners)
    .map((p) => ({
      label: `${lay.names[lay.stemSystem[p.stem]]} · ${Math.round(p.count).toLocaleString('en-US')} links${p.strand < 0 ? ', reversed' : ''}`,
      command: 'section.include',
      args: { address, include: pairNames(lay, lay.stemSystem[stem], lay.stemSystem[p.stem]) },
    }));
}

/** "Show only chr4": its stem's chromosome alone, in the section the layer opens (else its own). */
export function onlyStem(f: HeaderFacts): { label: string; args: { address: string; include: string[] } } | undefined {
  const s = stemOfFacts(f);
  if (!s) return undefined;
  const name = s.lay.names[s.lay.stemSystem[s.stem]];
  return { label: `Show only ${name}`, args: { address: s.lay.prep?.opens ?? s.lay.section, include: [name] } };
}
