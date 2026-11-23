/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import {
  continuousLut,
  formatValue,
  type Legend,
  legendTitle,
  LutSize,
  mergeLegends,
  qualifyLegends,
} from '@molstar/seqx-core/theme/index';
import { Anchor, type PassWriter } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import type { Box } from './layer-type.js';
import type { TextMetrics } from './text.js';
import type { SeqTokens } from './tokens.js';
import type { ResolvedLayer, ResolvedView } from '@molstar/seqx-core/view';
import { layerName } from './headers.js';
import { layerKey } from './tracks-columns.js';
import type { LayerLayout, ViewLayout } from './tracks.js';

// Legends (_spec/README.md "Render", Export; "Themes", Legends): derived from each laid-out layer's channels, merged
// when they say the same, never drawn on the canvas over data. An export draws them when asked, laid out here into a
// band outside the data box (a ramp as one row of the `image` primitive, as density strips are).

/**
 * The legends of a layout: each layer's, titled by its track's label, merged; with the notes they say once. A
 * categorical legend lists `maxItems` entries and counts the rest in `more` (the viewer asks for all of them: its part
 * opens a legend in place).
 */
export function viewLegends(layout: ViewLayout, maxItems = 12): { legends: Legend[]; issues: Issue[] } {
  const all: Legend[] = [];
  const add = (l: { type: LayerLayout['type']; lay: unknown; layer: LayerLayout['layer'] }, title: string) => {
    // the figure's own say: a layer's `legend: none` keys nothing (_spec/README.md "Themes", Legends)
    if (l.layer.params.legend === 'none') return;
    const ls = l.type.legend?.(l.lay, { title, maxItems });
    if (ls) all.push(...ls);
  };
  for (const s of layout.sections) {
    for (const t of s.tracks) for (const l of t.layers) add(l, t.track.label);
    // a section's own layers (pairs' arcs and dots), titled by the track that holds them
    for (const l of s.sectionLays) add(l, s.section.tracks.find((t) => t.layers.includes(l.layer))?.label ?? '');
  }
  // a bridge's layers (its ribbons), titled by the bridge's label, after the sections' (_spec/README.md "Render", Bridges)
  for (const b of layout.bridges?.list ?? []) for (const l of b.lays) add(l, b.label);
  // after the merge: legends of one title but of different channels say their channel (`Variants (size)`)
  const legends = qualifyLegends(mergeLegends(all));
  const issues: Issue[] = [];
  for (const l of legends)
    if (l.notes.includes('unlisted'))
      issues.push(
        issueOf(
          RenderIssues,
          'seq.render.legend.unlisted',
          `The legend ${l.title} lists the hashed values painted in view`,
          {
            address: l.layers[0],
            hint: RenderIssues.codes['seq.render.legend.unlisted'].hint,
          },
        ),
      );
  return { legends, issues };
}

/**
 * A view's legends with its hidden layers' kept in place (_spec/README.md "Themes", Legends; the viewer's Legends: a
 * layer's row keeps its place when hidden): a hidden layer draws nothing, so it derives no legend; each legend it last
 * had (`last`, by layer key: what a caller keeps across settles) comes back, its layer listed in `hidden`, joined to a
 * shown legend of the same key (a merged legend, one of its layers hidden) or a row of its own; a hidden layer with none
 * kept (hidden from the start, or its key not derived yet) a row of its name alone, no entries. Each goes before the
 * first legend after it in the view (its tracks' layers, then its bridges'); a layer gone from the view is forgotten,
 * and `names` gets each layer's `layerName`. A layer saying `legend: none` keys nothing, hidden or not.
 */
export function keepHiddenLegends(
  layout: ViewLayout,
  view: ResolvedView | undefined,
  shown: readonly Legend[],
  last: Map<string, readonly Legend[]>,
  names: Map<string, string>,
): readonly Legend[] {
  const place = new Map<string, number>();
  const hidden: string[] = [];
  names.clear();
  const add = (k: string, l: ResolvedLayer, owner: { label: string; layers: readonly ResolvedLayer[] }) => {
    place.set(k, place.size);
    names.set(k, layerName(l, owner));
    if (l.params.hidden !== true) last.set(k, []);
    else if (l.params.legend !== 'none') hidden.push(k);
  };
  for (const s of layout.sections)
    for (const t of s.section.tracks) for (const l of t.layers) add(layerKey(t, l), l, t);
  for (const b of view?.bridges ?? []) for (const l of b.layers) if (l.address) add(l.address, l, b);
  for (const l of shown) for (const k of l.layers) last.set(k, [...(last.get(k) ?? []), l]);
  for (const k of [...last.keys()]) if (!place.has(k)) last.delete(k);
  if (!hidden.length) return shown;
  const out = [...shown];
  const own = new Map<string, { legend: Legend; layers: string[] }>();
  for (const k of hidden) {
    const known = last.get(k);
    const was: readonly Legend[] = known?.length
      ? known
      : [
          {
            key: `hidden|${k}`,
            title: names.get(k)!,
            channel: 'color',
            entries: { kind: 'swatches', items: [] },
            open: false,
            source: 'vocabulary',
            more: 0,
            notes: [],
            layers: [k],
          },
        ];
    for (const w of was) {
      const i = out.findIndex((l) => l.key === w.key);
      if (i >= 0) out[i] = { ...out[i], layers: [...out[i].layers, k], hidden: [...(out[i].hidden ?? []), k] };
      else {
        const g = own.get(w.key);
        if (g) g.layers.push(k);
        else own.set(w.key, { legend: w, layers: [k] });
      }
    }
  }
  for (const { legend, layers } of own.values()) {
    const at = place.get(layers[0])!;
    const i = out.findIndex((l) => (place.get(l.layers[0]) ?? Number.POSITIVE_INFINITY) > at);
    out.splice(i < 0 ? out.length : i, 0, { ...legend, layers, hidden: layers });
  }
  return out;
}

/** The legend band's numbers (CSS px): a row's height, a swatch's side, gaps, the ramp's width. */
export const LegendBand = { row: 18, swatch: 10, gap: 6, ramp: 96, font: '11px sans' } as const;

/**
 * Lays the legends out in `box` (a band below or right of the data) into a pass: a row a legend (its title, then its
 * entries), swatches flowing on, a ramp as an image strip with its ends, steps as adjacent cells. Returns the height
 * used. Text in the table's text inks; the swatches in their colors as resolved.
 */
export function layoutLegends(
  legends: readonly Legend[],
  box: Box,
  tokens: SeqTokens,
  text: TextMetrics,
  out: PassWriter,
): number {
  const B = LegendBand;
  const ink = tokens.colors['msui-text'];
  const muted = tokens.colors['msui-text-muted'];
  let y = box.y;
  for (const l of legends) {
    let x = box.x;
    const label = (s: string, color: number) => {
      const w = text.measure(B.font, s);
      // the baseline of the line the text is on now (a wrapped entry's label follows its swatch down)
      out.text(s, x, y + B.row / 2 + 4, Anchor.start, B.font, color, w);
      x += w + B.gap;
    };
    const wrap = (w: number) => {
      if (x + w > box.x + box.width && x > box.x) {
        y += B.row;
        x = box.x + B.gap * 2;
      }
    };
    label(`${legendTitle(l)}${l.open ? ' (in view)' : ''}`, muted);
    const e = l.entries;
    const top = y + (B.row - B.swatch) / 2;
    if (e.kind === 'swatches') {
      for (const i of e.items) {
        wrap(B.swatch + B.gap + text.measure(B.font, i.label));
        out.fill(x, y + (B.row - B.swatch) / 2, B.swatch, B.swatch, i.color);
        x += B.swatch + 3;
        label(i.label, ink);
      }
      if (l.more) {
        wrap(text.measure(B.font, `+${l.more} more`));
        label(`+${l.more} more`, muted);
      }
    } else if (e.kind === 'ramp') {
      const [lo, hi] = l.domain ?? [0, 1];
      label(formatValue(lo), ink);
      const lut = continuousLut(
        e.stops.map((s) => s.color),
        e.stops.map((s) => s.at),
      );
      const data = new Uint8Array(4 * LutSize);
      for (let i = 0; i < LutSize; i++) {
        const c = lut[i];
        data[4 * i] = c >>> 24;
        data[4 * i + 1] = (c >>> 16) & 0xff;
        data[4 * i + 2] = (c >>> 8) & 0xff;
        data[4 * i + 3] = 0xff;
      }
      out.image(x, top, B.ramp, B.swatch, { width: LutSize, height: 1, data });
      // the ticks where the scale puts them (a log scale's decades, a center at the middle)
      for (const t of e.tickAt) if (t >= 0 && t <= 1) out.fill(x + t * (B.ramp - 1), top + B.swatch - 3, 1, 3, ink);
      x += B.ramp + B.gap;
      label(formatValue(hi), ink);
    } else if (e.kind === 'steps') {
      for (let i = 0; i < e.colors.length; i++) {
        const s = `${formatValue(e.edges[i])}–${formatValue(e.edges[i + 1])}`;
        wrap(B.swatch + B.gap + text.measure(B.font, s));
        out.fill(x, top, B.swatch, B.swatch, e.colors[i]);
        x += B.swatch + 3;
        label(s, ink);
      }
    } else {
      label(e.items.map((i) => `${formatValue(i.value)}: ${formatValue(i.size)}`).join('  '), ink);
    }
    y += B.row;
  }
  return y - box.y;
}

/** The height a band of legends takes at a width (a dry run of `layoutLegends` with no output). */
export function legendsHeight(legends: readonly Legend[], width: number, tokens: SeqTokens, text: TextMetrics): number {
  const sink = { text() {}, fill() {}, image() {} } as unknown as PassWriter;
  return layoutLegends(legends, { x: 0, y: 0, width, height: 0 }, tokens, text, sink);
}
