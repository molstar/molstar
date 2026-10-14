/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for RN2's feature, pair and logo draws (the Sonnet half; rn2-fixtures.ts holds the hot-path ones): a document
// of every feature layer over a short chain, a logo over letters with known stacks, and the overlay's marks as the
// composer will give them; a recording's items in the shapes the tests read.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { DrawView } from "../src/draw/view.ts";
import { type LayerRenderer, MarkBit, type SeqViewerContext } from "../src/index.ts";
import { type OverlayMarks, SectionMarks } from "../src/marks.ts";
import type { MatrixState } from "../src/rows.ts";
import { LayerFrame, RecordingBackend } from "../src/testing/index.ts";
import { layerOf } from "./draw-fixtures.ts";

export { logoDoc } from "./logo-fixtures.ts";

const PROTEIN = "MSEEKLKELAEKAIRLAKEGNPEEALKLLEEALRLAKESGDPELIARVLELLAEVTKQSVHAYRLDDGKWFLQPGGRVSLEMAPTDLQAAH";

/**
 * Every feature layer over one chain of `n`: "Ranges" (labeled [10, 30) and [50, 56), one lane), "Lanes" (those and an
 * overlapping [20, 44): two lanes), "Brackets" (the first two, as brackets), "Cuts" (boundaries before 20 and 48),
 * "Points" (5, 6, 7 and 30), "Variants" (two at 12, one at 40) and "Secondary" (a helix, a strand, a turn and a
 * one-position bar), with a ruler and letters.
 */
export function featuresDoc(o: { n?: number } = {}) {
  const n = o.n ?? 60;
  const b = createSeqViewSpec({ title: "Features" });
  const a = b.data.sequence({
    ref: "a",
    label: "Chain A",
    value: PROTEIN.repeat(Math.ceil(n / PROTEIN.length)).slice(0, n),
  });
  const ann = (ref: string, value: unknown) =>
    b.data.annotation({ ref, geometry: "ranges", target: a, value } as never);
  const flat = ann("ra", {
    start: [10, 50],
    stop: [30, 56],
    label: ["Alpha domain", "Gamma"],
    category: ["c0", "c1"],
  });
  const lanes = ann("rb", { start: [10, 20, 50], stop: [30, 44, 56], category: ["c0", "c1", "c0"] });
  const cuts = ann("bd", { boundary: [20, 48], category: ["c0", "c1"] });
  const points = ann("pt", { start: [5, 6, 7, 30], category: ["p0", "p1", "p0", "p1"] });
  const variants = ann("vr", { start: [12, 12, 40], alt: ["A", "C", "D"] });
  const secondary = ann("ss", {
    start: [4, 14, 26, 34],
    stop: [12, 24, 32, 35],
    category: ["helix", "strand", "turn", "other"],
  });
  const v = b.view({ ref: "v", label: "Features", height: "box" });
  v.block({ sequence: a });
  v.track({ label: "" }).ruler();
  v.track({ label: "Letters" }).letters();
  const color = { theme: "field_palette", params: { field: "category" } };
  v.track({ label: "Ranges" }).ranges({ data: flat, labels: true, color } as never);
  v.track({ label: "Lanes" }).ranges({ data: lanes, color } as never);
  v.track({ label: "Brackets" }).ranges({ data: flat, display: "brackets", color } as never);
  v.track({ label: "Cuts" }).ranges({ data: cuts, color } as never);
  v.track({ label: "Points" }).points({ data: points, color } as never);
  v.track({ label: "Variants", height: 3 }).variants({ data: variants });
  v.track({ label: "Secondary" }).secondary({ data: secondary });
  return b.getDocument();
}

/** A section's marks with `bit`'s channel lit over block 0's [start, stop) pairs (`runs`): what an overlay reads. */
export function lit(bit: number, runs: readonly (readonly [number, number])[], blocks = 1, length = 1_000_000) {
  const marks = new SectionMarks();
  const set = marks.runs(bit);
  set.reset(blocks);
  for (let b = 0; b < blocks; b++) {
    set.begin(b);
    if (b === 0) for (const [s, t] of runs) set.add(s, t);
    set.close(0, length);
  }
  return marks as OverlayMarks;
}
export { MarkBit };

/** A track's overlay as a frame: draw it into any recording again and again, with `marks` as `v.marks`. */
export function overlayFrame(ctx: SeqViewerContext, label: string, marks: OverlayMarks) {
  const { k, section } = layerOf(ctx, label);
  const entry = ctx.registry.get("renderers", ctx.graph.assembled.value?.view?.prepared[k].renderer ?? "");
  if (!entry?.overlay) throw new Error(`no overlay for "${label}"`);
  const overlay = entry.overlay;
  const wrap = {
    ...entry,
    draw(v: DrawView, s: unknown, p: unknown) {
      v.marks = marks;
      overlay(v, s, p);
    },
  } as LayerRenderer;
  return new LayerFrame(ctx, wrap, section, k);
}
/** A track's overlay drawn once into a recording. */
export function overlayOf(ctx: SeqViewerContext, label: string, marks: OverlayMarks) {
  const rec = new RecordingBackend();
  overlayFrame(ctx, label, marks).draw(rec);
  return rec;
}

/** A recording's paths as [mode, width, color, alpha, ...coordinates]. */
export function paths(r: RecordingBackend): number[][] {
  const b = r.pathItems;
  const out: number[][] = [];
  for (let p = 0; p < b.n; p++) {
    const end = p + 1 < b.n ? b.cstart[p + 1] : b.nc;
    out.push([b.mode[p], b.width[p], b.color[p], b.alpha[p], ...Array.from(b.coords.subarray(b.cstart[p], end))]);
  }
  return out;
}
/** A recording's glyphs as [x, y, radius, glyph, color, alpha]. */
export function glyphs(r: RecordingBackend): number[][] {
  const b = r.glyphItems;
  return Array.from({ length: b.n }, (_, i) => [b.x[i], b.y[i], b.radius[i], b.glyph[i], b.color[i], b.alpha[i]]);
}
export type { MatrixState };
