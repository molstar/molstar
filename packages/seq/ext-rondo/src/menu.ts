/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { rangesLocation } from '@molstar/seqx-core/extension';
import type { HeaderElementAt, HeaderFacts, HeaderSubEntry, SeqViewer } from '@molstar/seqx-viewer';
import { Kind, type RondoLay } from './lay.js';
import { laysOf, type Partner, partnerChromosomes, partnersOf } from './partners.js';
import { links, spanText } from './words.js';

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

/**
 * A range of a system named by the cytobands a viewer draws over it on the lay's section (`q34.1`, `q21.1–q22.3`):
 * the first and the last band it meets, read from the section's `bands` layers as laid out (their items in view); its
 * span in bases when no band is drawn there (`88–104 Mb`).
 */
export function bandRange(viewer: SeqViewer, lay: RondoLay, system: string, start: number, end: number): string {
  return bandsOf(viewer, lay, system, start, end) ?? spanText(start, end);
}

/** A range's cytobands as `bandRange` names them; undefined when the section draws none there. */
export function bandsOf(
  viewer: SeqViewer,
  lay: RondoLay,
  system: string,
  start: number,
  end: number,
): string | undefined {
  const s = viewer.frame?.layout.sections.find((x) => x.section.address === lay.section);
  const st = s?.placement.strips;
  const blocks = s?.section.blocks;
  let first: { label: string; start: number } | undefined;
  let last: { label: string; start: number } | undefined;
  if (s && st && blocks)
    for (let i = 0; i < st.n; i++) {
      if (blocks[st.block[i]]?.axis.system !== system) continue;
      for (let b = 0; b < st.bands; b++)
        for (const l of s.lays[i * st.bands + b] ?? []) {
          if (l.type.name !== 'bands') continue;
          for (const it of l.type.items?.(l.lay) ?? []) {
            if (it.start >= end || it.end <= start) continue;
            if (!first || it.start < first.start) first = it;
            if (!last || it.start > last.start) last = it;
          }
        }
    }
  if (!first || !last) return undefined;
  return first.label === last.label ? first.label : `${first.label}–${last.label}`;
}

/** A piece's row: its band range, its links, reversed; the stem itself `this stem` (its links within it). */
const pieceLabel = (v: SeqViewer, lay: RondoLay, stem: number, p: Partner) =>
  `${bandRange(v, lay, lay.systems[lay.stemSystem[p.stem]], lay.stemStart[p.stem], lay.stemEnd[p.stem])} · ${links(p.links)}${p.stem === stem ? ', this stem' : ''}${p.strand < 0 ? ', reversed' : ''}`;

/**
 * "Show with a partner ▸" (_spec/README.md, The stem menu): a row a partner chromosome with its links, most first, each
 * shown with the stem's as a pair; a chromosome reached by more than one ribbon its pieces in a submenu, each named by
 * its band range and strand; the stem's own chromosome last, "Within chr22", shown alone.
 */
export function partnerItems(f: HeaderFacts): HeaderSubEntry[] {
  const s = stemOfFacts(f);
  if (!s || !f.viewer) return [];
  const { lay, stem } = s;
  const v = f.viewer;
  const address = lay.prep?.opens ?? lay.section;
  return partnerChromosomes(lay, stem, MenuPartners).map((c) => {
    const args = { address, include: pairNames(lay, lay.stemSystem[stem], c.system) };
    const name = lay.names[c.system];
    const one = c.pieces.length === 1 ? c.pieces[0] : undefined;
    const reversed = one && one.strand < 0 ? ', reversed' : '';
    return {
      label: `${c.within ? `Within ${name}` : name} · ${links(c.links)}${reversed}`,
      command: 'section.include',
      args,
      ...(one
        ? {}
        : { items: c.pieces.map((p) => ({ label: pieceLabel(v, lay, stem, p), command: 'section.include', args })) }),
    };
  });
}

/** "Show only chr4": its stem's chromosome alone, in the section the layer opens (else its own). */
export function onlyStem(f: HeaderFacts): { label: string; args: { address: string; include: string[] } } | undefined {
  const s = stemOfFacts(f);
  if (!s) return undefined;
  const name = s.lay.names[s.lay.stemSystem[s.stem]];
  return { label: `Show only ${name}`, args: { address: s.lay.prep?.opens ?? s.lay.section, include: [name] } };
}

/** Where `rondo:partner` stands: the stem it last moved to and the partner it is on, by lay. */
const Walk = new WeakMap<RondoLay, { from: number; next: number }>();

/**
 * The cursor to the next (or the previous) partner stem of its stem, in the stem menu's order (its rows' pieces in
 * turn: one list, two readers), said.
 */
export function walkPartners(v: SeqViewer, back: boolean): { message: string } {
  const c = v.selection.cursor.value;
  if (!c) return { message: 'No cursor' };
  for (const lay of laysOf(v)) {
    if (lay.section !== c.section) continue;
    let from = -1;
    for (let s = 0; s < lay.ns; s++)
      if (
        lay.stemPlaced[s] &&
        lay.systems[lay.stemSystem[s]] === c.system &&
        c.position >= lay.stemStart[s] &&
        c.position < lay.stemEnd[s]
      )
        from = s;
    const w = Walk.get(lay);
    // on a partner reached by the walk: the walk goes on from its stem
    if (w && from !== w.from && w.next >= 0) from = w.from;
    if (from < 0) return { message: 'No stem here' };
    // the stem itself, a piece of its links within, is no step
    const partners = partnerChromosomes(lay, from).flatMap((x) => x.pieces.filter((p) => p.stem !== from));
    if (!partners.length) return { message: 'No partners' };
    const k0 = w && w.from === from ? w.next : -1;
    const k = (((k0 + (back ? -1 : 1)) % partners.length) + partners.length) % partners.length;
    Walk.set(lay, { from, next: k });
    const p = partners[k];
    const position = lay.stemStart[p.stem];
    v.selection.setCursor({ section: c.section, system: lay.systems[lay.stemSystem[p.stem]], position });
    return {
      message: `${lay.names[lay.stemSystem[p.stem]]}, ${position.toLocaleString('en-US')} · ${links(p.links)}${p.strand < 0 ? ', reversed' : ''}`,
    };
  }
  return { message: 'No genome ring' };
}

/** A range's cytobands as the composed viewer whose last frame holds the lay draws them (a pair's card's ends). */
export function bandOfLay(
  viewers: ReadonlySet<SeqViewer>,
  lay: RondoLay,
  system: string,
  start: number,
  end: number,
): string | undefined {
  for (const v of viewers) if (laysOf(v).includes(lay)) return bandsOf(v, lay, system, start, end);
  return undefined;
}

/**
 * At a rest: each empty hole's words said once, politely (`viewer.notify`; `said` the viewer's, by layer), again only
 * after the hole had other words or none.
 */
export function sayHoles(v: SeqViewer, said: Map<string, string>): void {
  const now = new Set<string>();
  for (const lay of laysOf(v)) {
    if (!lay.holeText) continue;
    now.add(lay.key);
    if (said.get(lay.key) === lay.holeText) continue;
    said.set(lay.key, lay.holeText);
    v.notify(lay.holeText);
  }
  for (const k of [...said.keys()]) if (!now.has(k)) said.delete(k);
}
