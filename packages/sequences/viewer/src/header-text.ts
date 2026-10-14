/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A track header's text where the canvas cuts it with an ellipsis (seqspec:interaction#the-legend, the header's
// tooltip): the full text, measured as the composer draws it (draw/composer-headers.ts: a track's label, a rows track's
// row label). Edge-rate: once per hovered header, through the hover's label (labels.ts, `LabelHost.cutText`).

import { Lanes } from "@molstar/seq-core/format";
import type { SeqViewerContext } from "./context.ts";
import { CHEVRON_ROOM, splitChevron } from "./header-toggle.ts";
import { HoverKind } from "./hit.ts";
import type { HoverState } from "./hover.ts";
import { EntryKind } from "./layout.ts";

/** The full text of the hovered header when it is drawn cut; null: none hovered, or drawn whole. */
export function cutHeaderText(ctx: SeqViewerContext, h: HoverState): string | null {
  if (h.kind !== HoverKind.header || h.track < 0 || h.treeNode >= 0) return null;
  const view = ctx.graph.assembled.value?.view;
  const ref = ctx.layout.sectionRef(h.section);
  const l = ref === null ? undefined : ctx.layout.get(ref)?.current;
  if (!view || ref === null || !l) return null;
  const rv = view.resolved;
  const px = ctx.layout.props(ref).labelFontPx;
  const r = l.rows;
  // A rows track's row: its label beside the tree's inset, indented by its depth.
  if (r.layer >= 0 && h.track === r.track && h.lane >= 0) {
    const rend = ctx.registry.get("renderers", view.set.strings.get(rv.layers.type[r.layer]));
    const s = view.prepared[r.layer]?.state;
    if (!rend?.rowModel) return null;
    const inset = rend.header ? Math.min(rend.header.inset(s, rv.layers.params[r.layer]), l.headerWidth * 0.5) : 0;
    const x = 10 + inset + Lanes.depthOf(rv.lanes, h.lane) * 14;
    const font = ctx.text.font("ui", Math.max(8, Math.min(px, Math.floor(r.height * 0.8))));
    return cut(ctx, rend.rowModel.label(s, Lanes.rowOf(rv.lanes, h.lane)), font, l.headerWidth - x - 6);
  }
  // A track's label, indented by its groups' depth (a split row's value is not a label: left out).
  if (h.split >= 0) return null;
  let e = 0;
  while (e < l.entryCount && !(l.entryKind[e] === EntryKind.track && l.entryTrack[e] === h.track)) e++;
  if (e === l.entryCount) return null;
  let depth = 0;
  for (let g = l.entryGroup[e]; g >= 0; g = rv.groups.parent[g]) depth++;
  const t = rv.tracks.label[h.track] >= 0 ? rv.tracks.label[h.track] : rv.tracks.ref[h.track];
  const font = ctx.text.font("ui", Math.max(8, Math.min(px, Math.floor(l.entryHeight[e] * 0.8))), 400);
  // A split track's chevron takes its room from the label (header-toggle.ts).
  const chevron = l.entryRow[e] < 0 && splitChevron(view, l, e) ? CHEVRON_ROOM : 0;
  return cut(ctx, view.set.strings.get(t), font, l.headerWidth - (10 + depth * 14) - 6 - chevron);
}

/**
 * `text` when it is wider than `width` (cut, as `fitIn` reads it), else null. Measured without the label cache (a
 * hover over many rows would fill it) and allocation-free: advances for ASCII, the measurer otherwise.
 */
function cut(ctx: SeqViewerContext, text: string, font: number, width: number): string | null {
  let w = 0;
  const a = ctx.text.advances(font);
  for (let i = 0; i < text.length && w >= 0; i++) {
    const c = text.charCodeAt(i);
    w = c < 128 ? w + a[c] : -1;
  }
  if (w < 0) w = ctx.text.port.measure(ctx.text.fontString(font), text);
  return w > width ? text : null;
}
