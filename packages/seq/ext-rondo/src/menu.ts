/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { rangesLocation } from '@molstar/seqx-core/extension';
import type { HeaderEntry, HeaderFacts, SeqViewer } from '@molstar/seqx-viewer';
import { Kind, type RondoLay } from './lay.js';

// The stem menu (_spec/README.md, The viewer, The stem menu): a stem's element menu, the viewer's label menu's shape
// scoped to `rondo:ribbons`' stems (`HeaderFacts.element.type`): its partners by count, each shown with it as a pair
// (in the section the layer `opens`, else its own); the stem and its partners selected; its chromosome alone. Its
// entries read the lays of the viewer the menu is in (`HeaderFacts.viewer`): no state of their own, so two viewers of
// one spec each list their own.

/** How many partners the menu lists. */
const MenuPartners = 8;

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

/** A menu's stem from its facts: in the viewer whose menu it is. */
const stemOfFacts = (f: HeaderFacts) =>
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

/** The stem menu's entries (`on: element`, scoped to `rondo:ribbons`' stems). */
export const StemMenu: readonly HeaderEntry[] = [
  {
    name: 'rondo:show-partner',
    on: 'element',
    label: 'Show with a partner',
    order: 1,
    when: (f) => !!stemOfFacts(f),
    items: (f) => {
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
    },
  },
  {
    name: 'rondo:select-stem-partners',
    on: 'element',
    label: 'Select the stem and its partners',
    order: 2,
    when: (f) => !!stemOfFacts(f),
    command: 'rondo:select-partners',
    args: (f) => ({ layer: f.element!.layer, key: f.element!.key }),
  },
  {
    name: 'rondo:show-chromosome',
    order: 3,
    on: 'element',
    label: (f) => {
      const s = stemOfFacts(f);
      return s ? `Show only ${s.lay.names[s.lay.stemSystem[s.stem]]}` : 'Show only this chromosome';
    },
    when: (f) => !!stemOfFacts(f),
    command: 'section.include',
    args: (f) => {
      const s = stemOfFacts(f)!;
      return { address: s.lay.prep?.opens ?? s.lay.section, include: [s.lay.names[s.lay.stemSystem[s.stem]]] };
    },
  },
];
