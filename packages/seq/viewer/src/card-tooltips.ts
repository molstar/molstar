/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type TooltipOut, type TooltipPart, tooltipsAt } from '@molstar/seqx-core/render';
import type { CardAt, CardParts, CardProviderEntry } from './card-provider.js';
import type { HoverCardRow } from './hover-card.js';

// The spec's tooltips on the card (_spec/README.md "Hover card", Tooltips): the built-in `seq:tooltips` provider over
// the core's `tooltipsAt`. It answers only on drawn marks (an item, a pair, a cell, a letter, a value, a row's cell, or
// a bin of them), never on an empty band; it ignores `whole` (a note on 94 to 292 shows over a feature, a cell or a
// pair the item and region providers said whole). Its rows share the card's `hover.cardItems` cap with the place's
// items, after them; the matches past it are counted (`moreCount`). The live region hears the rows shown and the count.
// A host replaces it by an entry of its own under its name, over the same `tooltipsAt`.

/** The layer types that draw a mark at every position they show (their hit has no element). */
const PositionMarks: ReadonlySet<string> = new Set(['letters', 'values', 'grid', 'rows', 'matrix', 'logo']);

/** The stretches under the hovered mark: the pointer's position (and a cell's y), a pair's two ends. */
function partsOf(at: CardAt): TooltipPart[] {
  const p = at.place;
  const parts: TooltipPart[] = [];
  const el = p.element;
  // a pair: both its ends (a note on either end answers)
  if (el && at.layer?.type.name === 'pairs') {
    const owner = at.layer.layer.data?.address;
    const a = owner ? at.model.annotations.get(owner) : undefined;
    const d = a?.pairs;
    const k = d ? d.keys.indexOf(el.key) : -1;
    if (a && d && k >= 0) {
      const as = (a.units && d.aUnit && a.units[d.aUnit[k]]) || a.target;
      const bs = (a.units && d.bUnit && a.units[d.bUnit[k]]) || (a.ySystem ?? a.target);
      parts.push(
        { system: as, start: d.aStart[k], end: d.aEnd[k] },
        { system: bs, start: d.bStart[k], end: d.bEnd[k] },
      );
      return parts;
    }
  }
  if (p.system !== undefined && !Number.isNaN(at.position))
    parts.push({ system: p.system, start: at.position, end: at.position + 1 });
  if (p.ySystem !== undefined && !Number.isNaN(at.yPosition))
    parts.push({ system: p.ySystem, start: at.yPosition, end: at.yPosition + 1 });
  return parts;
}

function describeTooltips(at: CardAt, out: CardParts): void {
  const view = at.view;
  const l = at.layer;
  if (!view?.tooltips || !l) return;
  const p = at.place;
  // drawn marks only: an element under the pointer, or a layer that marks every position it shows
  if (p.element === undefined && !PositionMarks.has(l.type.name)) return;
  const parts = partsOf(at);
  const el = p.element;
  const it = el && l.type.items ? at.words.itemsOf(l).find((d) => d.key === el.key) : undefined;
  // the cap shared with the place's items, after them; a single row on the mark always shows
  let items = 0;
  for (const r of out.rows) if (r.item) items++;
  const o: TooltipOut = { rows: [], more: 0, failed: [] };
  tooltipsAt(
    { model: at.model, view, data: at.store.store, locate: (a) => at.words.locate(a) },
    { layer: l.layer, ...(it ? { item: { key: it.key, start: it.start } } : {}), parts },
    o,
    Math.max(1, at.cardItems - items),
  );
  if (!o.rows.length) return;
  const said: string[] = [];
  for (const r of o.rows) {
    const row: HoverCardRow = {
      label: r.label,
      value: r.value,
      ...(r.rich ? { rich: r.rich } : {}),
      ...(r.item ? { item: true } : {}),
    };
    out.rows.push(row);
    said.push(r.value ? `${r.label} ${r.value}` : r.label);
  }
  out.moreCount += o.more;
  out.words.push(`${said.join('; ')}${o.more ? `; and ${o.more.toLocaleString('en-US')} more` : ''}`);
}

/** The spec's tooltips (`seq:tooltips`, order 95: after the place's rows, before the block). */
export const TooltipsCardProvider: CardProviderEntry = {
  kind: 'card-provider',
  name: 'seq:tooltips',
  order: 95,
  describe: describeTooltips,
};
