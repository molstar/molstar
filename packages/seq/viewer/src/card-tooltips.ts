/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Location } from '@molstar/seqx-core/query';
import type { RangesLay, TooltipOut, TooltipPart, TooltipRow, ValuesLay } from '@molstar/seqx-core/render';
import { type TooltipContext, tooltipsAt } from '@molstar/seqx-core/render/tooltips';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer, ResolvedView } from '@molstar/seqx-core/view';
import type { CardAt, CardParts } from './card-provider.js';

// The spec's tooltips on the card (_spec/README.md "Hover card", Tooltips): the body of the built-in `seq:tooltips`
// provider (`card-providers.ts` holds its entry), over the core's `tooltipsAt`, loaded with the tooltips. It answers
// only on drawn marks (an item, a pair, a cell, a letter, a value, a row's cell, or a density bin over its stretch),
// never on an empty band or a gap (a gap letter, a missing value or cell, EC-254); it ignores `whole` (a note on 94 to
// 292 shows over a feature, a cell or a pair the item and region providers said whole). Its rows share the card's
// `hover.cardItems` cap with the place's items, after them; the matches past it are counted (`moreCount`). The live
// region hears the rows shown and the count. A table failed under the mark gives no rows and its issue on the tooltip
// (EC-174, `seq.viewer.card.tooltip`). A warm show allocates what the card without tooltips does, give or take its
// words: the rows are the core's, made once a table row (`tooltips.ts`), the buffers here reused.

/** The layer types that draw a mark at every position they show (their hit has no element). */
const PositionMarks: ReadonlySet<string> = new Set(['letters', 'values', 'grid', 'rows', 'matrix', 'logo']);

/** Letters that are a gap: no mark there (EC-254). */
const GapLetters = new Set(['', '-', '.', ' ']);

/** The lookup's buffers, reused show to show (one card is composed at a time). */
const Out: TooltipOut = { rows: [], more: 0, failed: [] };
const Parts: { system: string; start: number; end: number }[] = [
  { system: '', start: 0, end: 0 },
  { system: '', start: 0, end: 0 },
];
/** The parts' first one alone, both, none: what a place reads (its `parts`), no array made a show. */
const ByCount: readonly (readonly TooltipPart[])[] = [[], [Parts[0]], Parts];
const Item = { key: '' as string | number, start: 0 };
const Place: { layer?: ResolvedLayer; item?: typeof Item; parts: readonly TooltipPart[] } = { parts: [] };
/** The describer the lookup locates addresses by, for this show. */
let locating: CardAt['words'] | undefined;
const Ctx: { -readonly [K in keyof TooltipContext]: TooltipContext[K] } = {
  model: undefined as unknown as StaticModel,
  view: undefined as unknown as ResolvedView,
  data: undefined as unknown as TooltipContext['data'],
  locate: (a: string): Location | undefined => locating?.locate(a),
};

/** Each row's phrase in the live region (`label value`), made once a row. */
const Phrases = new WeakMap<TooltipRow, string>();
function phraseOf(r: TooltipRow): string {
  let p = Phrases.get(r);
  if (p === undefined) Phrases.set(r, (p = r.value ? `${r.label} ${r.value}` : r.label));
  return p;
}

/** Sets part `k` of the reused parts; the count of them. */
function part(k: number, system: string, start: number, end: number): number {
  const p = Parts[k];
  p.system = system;
  p.start = start;
  p.end = end;
  return k + 1;
}

/** The stretches under the hovered mark, into the reused parts: their count (0: no mark there). */
function partsOf(at: CardAt): number {
  const p = at.place;
  const el = p.element;
  const l = at.layer!;
  // a pair: both its ends (a note on either end answers)
  if (el && l.type.name === 'pairs') {
    const owner = l.layer.data?.address;
    const a = owner ? at.model.annotations.get(owner) : undefined;
    const d = a?.pairs;
    const k = d ? d.keys.indexOf(el.key) : -1;
    if (a && d && k >= 0) {
      const as = (a.units && d.aUnit && a.units[d.aUnit[k]]) || a.target;
      const bs = (a.units && d.bUnit && a.units[d.bUnit[k]]) || (a.ySystem ?? a.target);
      part(0, as, d.aStart[k], d.aEnd[k]);
      return part(1, bs, d.bStart[k], d.bEnd[k]);
    }
  }
  if (p.system === undefined || Number.isNaN(at.position)) return 0;
  // a density bin: the positions its pixel column spans
  if (!el && l.type.name === 'ranges' && (l.lay as RangesLay).mode === 'density') {
    const axis = at.section.axis;
    const px = axis ? Math.abs(axis.x(at.position + 1) - axis.x(at.position)) : 1;
    return part(0, p.system, at.position, at.position + Math.max(1, Math.ceil(1 / Math.max(px, 1e-9))));
  }
  let n = part(0, p.system, at.position, at.position + 1);
  if (p.ySystem !== undefined && !Number.isNaN(at.yPosition)) n = part(n, p.ySystem, at.yPosition, at.yPosition + 1);
  return n;
}

/** Whether a position mark is a gap: a gap letter, a missing value, a missing cell (no mark: no tooltip, EC-254). */
function isGap(at: CardAt): boolean {
  const l = at.layer!;
  const p = at.place;
  switch (l.type.name) {
    case 'letters':
      return GapLetters.has(at.words.letterAt(at.section, at.position));
    case 'values': {
      const axis = at.section.axis;
      if (!axis) return true;
      const x = (axis.x(at.position) + axis.x(at.position + 1)) / 2;
      const lay = l.lay as ValuesLay;
      for (let j = 0; j < lay.n; j++) if (lay.x0[j] <= x && x <= lay.x1[j]) return Number.isNaN(lay.v[j]);
      return true;
    }
    case 'grid':
      return (
        p.system === undefined ||
        p.ySystem === undefined ||
        Number.isNaN(at.words.gridValue(l, p.system, at.position, p.ySystem, at.yPosition))
      );
    default:
      return false;
  }
}

export function describeTooltips(at: CardAt, out: CardParts): void {
  const view = at.view;
  const l = at.layer;
  if (!view?.tooltips || !l) return;
  const p = at.place;
  const density = l.type.name === 'ranges' && (l.lay as RangesLay).mode === 'density';
  // drawn marks only: an element under the pointer, a layer that marks every position it shows, a density bin
  if (p.element === undefined && !density && (!PositionMarks.has(l.type.name) || isGap(at))) return;
  const n = partsOf(at);
  const el = p.element;
  const it = el && l.type.items ? at.words.itemsOf(l).find((d) => d.key === el.key) : undefined;
  // the cap shared with the place's items, after them; a single row on the mark always shows
  let items = 0;
  for (const r of out.rows) if (r.item) items++;
  const o = Out;
  o.rows.length = 0;
  o.more = 0;
  o.failed.length = 0;
  Place.layer = l.layer;
  Place.parts = ByCount[n];
  if (it) {
    Item.key = it.key;
    Item.start = it.start;
    Place.item = Item;
  } else Place.item = undefined;
  Ctx.model = at.model;
  Ctx.view = view;
  Ctx.data = at.store.store;
  locating = at.words;
  tooltipsAt(Ctx, Place, o, Math.max(1, at.cardItems - items));
  locating = undefined;
  for (const address of o.failed) at.words.tooltipFailed(address);
  if (!o.rows.length) return;
  let said = '';
  for (const r of o.rows) {
    out.rows.push(r);
    const phrase = phraseOf(r);
    said = said ? `${said}; ${phrase}` : phrase;
  }
  out.moreCount += o.more;
  out.words.push(o.more ? `${said}; and ${o.more.toLocaleString('en-US')} more` : said);
}
