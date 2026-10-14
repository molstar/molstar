/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for the axis, viewport and focus tests (VC3): views resolved headlessly (no network), a block's piece
// rewritten to a reversed or cyclic one through resolve's `layoutOf` (AX1 and CY1 give the format), and a context on a
// manual clock opened on the managers' document (the old s3 fixtures), sized before the graph settles or after.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import {
  compile,
  layoutOf,
  type Patch,
  type ResolvedView,
  resolveView,
  type SeqViewSpecDoc,
  supportedAll,
} from "@molstar/seq-core/format";
import type { CoordinateMap } from "@molstar/seq-core/model";
import { fixtureLoaders, registries } from "../../core/_test/format/load-fixtures.ts";
import { Axis, type ShownIndex } from "../src/axis.ts";
import { manualScheduler, type SeqViewerContext, type SeqViewerPorts } from "../src/index.ts";
import type { Viewport } from "../src/viewport.ts";
import { viewer } from "./graph-fixtures.ts";

const reg = registries(false);

/** A view of a document, resolved with the kernel's registries and the fixtures' files. */
export async function resolved(doc: SeqViewSpecDoc, view: string) {
  const { set, issues } = await compile(doc, reg, fixtureLoaders());
  const r = resolveView(doc, set, view, reg, supportedAll(reg));
  return { set, view: r.view, issues: [...issues, ...r.issues] };
}
export type Resolved = Awaited<ReturnType<typeof resolved>>;

/** Block b's piece as `[start, stop)` laid out in `dir` over a cyclic domain `[0, length)` when it wraps. */
export function setPiece(view: ResolvedView, b: number, start: number, stop: number, dir: number, length: number) {
  const B = view.blocks;
  const p = layoutOf(B.layout[b].from, start, stop, dir, length);
  [B.start[b], B.stop[b], B.dir[b], B.first[b]] = [start, stop, dir, p.first];
  (B.layout as CoordinateMap[])[b] = p.map;
}

export const axisOf = (r: Resolved, section = 0) => new Axis(r.view, r.set, section, 28, 3);
export const k = (n: number) => n as ShownIndex;

/** Every shown index's block index, and every block index's shown index over [0, domain). */
export function readBack(ax: Axis, b: number, domain: number) {
  return {
    index: Array.from({ length: ax.shown[b] }, (_, s) => ax.indexOf(b, k(s))),
    shown: Array.from({ length: domain }, (_, i) => ax.shownOf(b, i) as number),
    breaks: Array.from({ length: ax.shown[b] }, (_, s) => s).filter((s) => ax.breakAt(b, k(s))),
  };
}
const out = new Int32Array(64);
/** `runToShown` and `shownToRuns` as pairs. */
export function toShown(ax: Axis, b: number, start: number, stop: number) {
  const n = ax.runToShown(b, start, stop, out);
  return Array.from({ length: n }, (_, r) => [out[2 * r], out[2 * r + 1]]);
}
export function toRuns(ax: Axis, b: number, k0: number, k1: number) {
  const n = ax.shownToRuns(b, k0, k1, out);
  return Array.from({ length: Math.min(n, 32) }, (_, r) => [out[2 * r], out[2 * r + 1]]);
}

/** 20 residues; `u` its unit over [5, 15); `site` [7, 10) and [12, 14). View `v` the unit, `f` it filtered by `site`. */
export function rangedDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: "s", value: "MKVLAAGIVGLLLAQRSTVW" });
  const u = b.data.unit({ ref: "u", sequence: s, range: [5, 15] });
  const site = b.data.annotation({
    ref: "site",
    geometry: "ranges",
    target: s,
    value: { start: [7, 12], stop: [10, 14] },
  });
  const v = b.view({ ref: "v", label: "V" });
  v.block({ unit: u }).track({ ref: "v.l" }).letters();
  v.focus({ ref: "v.f", target: u });
  b.view({ ref: "f", label: "F", positions: site }).block({ unit: u }).track({ ref: "f.l" }).letters();
  b.view({ ref: "w", label: "Whole" }).block({ sequence: s }).track({ ref: "w.l" }).letters();
  return b.getDocument();
}

/** `n` sequences of `len` residues side by side in one section (24 chromosomes). */
export function blocksDoc(n: number, len = 10): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const v = b.view({ ref: "v", label: "V" });
  for (let i = 0; i < n; i++) v.block({ sequence: b.data.sequence({ ref: `c${i}`, value: "M".repeat(len) }) });
  v.track({ ref: "v.l" }).letters();
  return b.getDocument();
}

/**
 * The managers' document (the old s3 fixtures): `a` (20 residues) and `b` (10), `b[i]` mapped to `a[10 + i]`; `dom`
 * covers a's [0, 4) and [10, 16). View `v`: one section `v.main` (blocks a and b side by side) opening on a's [2, 12).
 * View `g`: `g.one` and `g.two` over a in scroll group `s`; `g.b` over b (linked to both through the mapping).
 */
export function s3Doc(): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const src = b.data.source({ ref: "pasted", value: ">a\nMKVLATGHEKMKVLATGHEK\n>b\nMKVLATGHEK\n", format: "fasta" });
  const a = b.data.sequence({ ref: "a", source_ref: src, record: "a" });
  const s = b.data.sequence({ ref: "b", source_ref: src, record: "b" });
  b.data.mapping({ ref: "m", from: s, to: a, method: "segments", segments: [[0, 10, 10]] });
  b.data.annotation({ ref: "dom", geometry: "ranges", target: "a", value: { start: [0, 10], stop: [4, 16] } });
  const v = b.view({ ref: "v", label: "V" });
  v.focus({ ref: "v.f", target: "a", start: 2, stop: 12 });
  v.track({ ref: "v.letters" }).letters({ ref: "v.letters.l" });
  v.track({ ref: "v.dom" }).ranges({ ref: "v.dom.l", data: "dom" });
  const g = b.view({ ref: "g", label: "G" });
  g.section({ ref: "g.one", scroll_group: "s" }).block({ sequence: a }).track().letters();
  g.section({ ref: "g.two", scroll_group: "s" }).block({ sequence: a }).track().ranges({ data: "dom" });
  g.section({ ref: "g.b" }).block({ sequence: s }).track().letters();
  return b.getDocument();
}

export const WIDTH = 160;

/**
 * A context on a manual clock with `doc` shown; every viewport sized before the graph first settles, or after it
 * (`late`: the order both UIs use). The layout's header and gutters are 0 px, so the viewport's arithmetic is over the
 * whole width.
 */
export async function open(doc: SeqViewSpecDoc = s3Doc(), o: { late?: boolean; ports?: Partial<SeqViewerPorts> } = {}) {
  const scheduler = manualScheduler();
  const v = viewer({ scheduler, ...o.ports });
  v.ctx.setRenderProps({ headerWidth: 0, gutterLeft: 0, gutterRight: 0 });
  const done = v.ctx.setDocument(doc);
  const size = () => {
    for (const ref of v.ctx.viewports.sections.value) v.ctx.viewports.get(ref)?.setSize(WIDTH, 100, 1);
  };
  if (!o.late) size();
  const r = await done;
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  return { ...v, scheduler, size };
}

export const vp = (ctx: SeqViewerContext, ref: string) => ctx.viewports.get(ref) as Viewport;

/** The block index at the center of a viewport's first block. */
export function center(v: Viewport): number {
  const s = v.live.current;
  const ax = v.axis;
  if (!ax) return -1;
  return ax.indexOf(0, k(Math.floor((s.scrollX + s.width / 2) / (s.cell * ax.span[0]))));
}

export const set = (ref: string, params: Record<string, unknown>): Patch => ({
  format: "seq_view_spec_patch",
  ops: [{ op: "set", ref, params }],
});
