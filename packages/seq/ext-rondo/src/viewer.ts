/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Hit, rangesLocation } from '@molstar/seqx-core/extension';
import {
  type BehaviorProvider,
  defineViewerExtension,
  type HeaderEntry,
  type SeqViewer,
  type SeqViewerExtension,
  type ViewerCommand,
} from '@molstar/seqx-viewer';
import { rondo, RondoParams, RondoVersion } from './index.js';
import { Kind, type RondoLay } from './lay.js';
import { laysOf, partnersOf, StemMenu, selectStemPartners, Viewers } from './menu.js';

// The genome ring extension's viewer parts (_spec/README.md, The viewer): `rondoViewer`, the headless `rondo` and
// its commands (`rondo:partner`, along the links by count from the cursor's stem; `rondo:select-partners`, the far
// ends of the selection's ribbons added; `rondo:open`, the hit's chromosomes shown in the section the layer `opens`),
// and the double-click scoped to `rondo:ribbons`' hits that runs `rondo:open` (it wins over the focus there only).

const lays = laysOf;

export { partnersOf, StemMenu };

/** The chromosomes a hit names: a ribbon's two, a stem's one, a stub's stem's and its partner's. */
export function chromosomesOf(lay: RondoLay, hit: Hit): string[] {
  const i = typeof hit.element?.key === 'string' ? lay.byKey.get(hit.element.key) : undefined;
  if (i === undefined) return [];
  const kind = lay.polyKind[i];
  const ref = lay.polyRef[i];
  const name = (sys: number) => lay.names[sys];
  if (kind === Kind.ribbon)
    return [...new Set([name(lay.stemSystem[lay.ribA[ref]]), name(lay.stemSystem[lay.ribB[ref]])])];
  // a folded stub: its stem's chromosome (its partners are many)
  if (kind === Kind.stub)
    return lay.stubPartner[ref] < 0
      ? [name(lay.stemSystem[lay.stubStem[ref]])]
      : [name(lay.stemSystem[lay.stubStem[ref]]), name(lay.stubPartner[ref])];
  return [name(lay.stemSystem[ref])];
}

/** Where `rondo:partner` stands: the stem it last moved to and the partner it is on (by count), by lay. */
const Walk = new WeakMap<RondoLay, { from: number; next: number }>();

/** A viewer's commands (each viewer its own: they act on the viewer that composed them). */
function viewerParts() {
  let viewer: SeqViewer | undefined;
  let unregister: (() => void) | undefined;
  /** The cursor to the next (or the previous) partner of its stem by count, said. */
  const walk = (args: { back?: boolean }) => {
    const v = viewer;
    const c = v?.selection.cursor.value;
    if (!v || !c) return { message: 'No cursor' };
    for (const lay of lays(v)) {
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
      const partners = partnersOf(lay, from);
      if (!partners.length) return { message: 'No partners' };
      const k0 = w && w.from === from ? w.next : -1;
      const k = (((k0 + (args.back ? -1 : 1)) % partners.length) + partners.length) % partners.length;
      Walk.set(lay, { from, next: k });
      const p = partners[k];
      const position = lay.stemStart[p.stem];
      v.selection.setCursor({ section: c.section, system: lay.systems[lay.stemSystem[p.stem]], position });
      return {
        message: `${lay.names[lay.stemSystem[p.stem]]}, ${position.toLocaleString('en-US')} · ${Math.round(p.count).toLocaleString('en-US')} links${p.strand < 0 ? ', reversed' : ''}`,
      };
    }
    return { message: 'No genome ring' };
  };
  const commands: ViewerCommand[] = [
    {
      id: 'rondo:partner',
      title: 'Go to the next partner',
      category: 'Genome ring',
      // `.` and `,` (the design's; free in the viewer's keys): along a stem's partners by count, back and forth
      shortcut: { key: '.' },
      scope: 'canvas',
      run: () => walk({}),
    },
    {
      id: 'rondo:partner-back',
      title: 'Go to the previous partner',
      category: 'Genome ring',
      shortcut: { key: ',' },
      scope: 'canvas',
      run: () => walk({ back: true }),
    },
    {
      id: 'rondo:select-partners',
      title: 'Select the partners',
      category: 'Genome ring',
      run: (_args: unknown, ctx: { layer?: string; key?: string } = {}) => {
        const v = viewer;
        if (!v) return { message: '' };
        // from a stem's menu: the stem and its partners' stems, one selection
        if (ctx.layer !== undefined && ctx.key !== undefined) return selectStemPartners(v, ctx.layer, ctx.key);
        const shown = v.selection.value.value.shown;
        let added = 0;
        for (const lay of lays(v))
          for (let r = 0; r < lay.nr; r++) {
            if (!lay.ribDrawn[r]) continue;
            for (const [near, far] of [
              [lay.ribA[r], lay.ribB[r]],
              [lay.ribB[r], lay.ribA[r]],
            ]) {
              const sys = lay.systems[lay.stemSystem[near]];
              const hit = shown.parts.some(
                (p) =>
                  p.system === sys &&
                  Array.from({ length: p.ranges.length >> 1 }, (_, k) => k).some(
                    (k) => p.ranges[2 * k] < lay.stemEnd[near] && p.ranges[2 * k + 1] > lay.stemStart[near],
                  ),
              );
              if (!hit) continue;
              v.selection.add(rangesLocation(lay.systems[lay.stemSystem[far]], [lay.stemStart[far], lay.stemEnd[far]]));
              added++;
            }
          }
        return { message: added ? `${added} partners added` : 'No ribbon in the selection' };
      },
    },
    {
      id: 'rondo:open',
      title: 'Show these chromosomes',
      category: 'Genome ring',
      run: (_args: unknown, ctx: { section?: string; hit?: Hit } = {}) => {
        const v = viewer;
        const hit = ctx.hit;
        if (!v || !hit) return { message: '' };
        for (const lay of lays(v)) {
          if (hit.layer !== lay.key || !lay.prep?.opens) continue;
          const names = chromosomesOf(lay, hit);
          if (!names.length) return { message: '' };
          return { message: v.runCommand('section.include', { address: lay.prep.opens, include: names }) };
        }
        return { message: '' };
      },
    },
  ] as unknown as ViewerCommand[];
  return {
    commands,
    composed(v: SeqViewer) {
      viewer = v;
      Viewers.add(v);
      // the stem menu: registered on the viewer's headers at run time (so `viewer.extend` brings it too), once the
      // composition holds this behavior's commands (its entries name them), taken away on dispose
      void Promise.resolve().then(() => {
        if (viewer !== v || unregister) return;
        unregister = v.behavior<{ register(...e: HeaderEntry[]): () => void }>('headers')?.register(...StemMenu);
      });
    },
    dispose() {
      unregister?.();
      unregister = undefined;
      if (viewer) Viewers.delete(viewer);
      viewer = undefined;
    },
  };
}

/** The behavior: its commands and the viewer they act on, one a viewer. */
function behavior(): BehaviorProvider {
  return { name: 'rondo:behavior', create: () => viewerParts() };
}

/** The genome ring extension with the viewer's parts: commands, the `opens` double-click. */
export const rondoViewer: SeqViewerExtension = defineViewerExtension({
  id: 'rondo',
  version: RondoVersion,
  requires: '0.10',
  doc: rondo.doc,
  registry: rondo.registry,
  params: RondoParams as never,
  behaviors: [{ provider: behavior() }],
  bindings: [{ gesture: 'double-click', target: 'item', action: 'rondo:open', layer: 'rondo:ribbons' }],
});
