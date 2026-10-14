/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Frames over rows through a mounted surface for the browser heap guards (rows.heap.ts), imported afresh per case: the
// context on a scheduler that never fires (the guard drives `renderNow` itself), a family of 12,000 rows with its tree,
// data columns and sites, or 20,000 lazy rows with their tiles cached, and a step that scrolls and draws or hovers.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { Loaders } from "@molstar/seq-core/format";
import { tilesStream } from "@molstar/seq-core/stream/tiles";
import { tileServer } from "@molstar/seq-core/stream/tiles/server";
import { Inputs, PointerType } from "../../src/index.js";
import { familyWithTree } from "../families.js";
import { close, type Opened, opened } from "./open.js";

/** The procedural tile server (the format's tests) as the viewer's progressive loader. */
function lazy(uri: string, rows: number, columns: number): Loaders {
  const server = tileServer({ rowCount: rows, columnCount: columns });
  return {
    load: (r) => Promise.reject(new Error(`no fixture for ${r.uri ?? r.hostFile}`)),
    stream: tilesStream((u) => (u === uri ? server.transport : null)),
  };
}

const moves = [7, -7, 11, -11];
const pan = (o: Opened, i: number) => {
  const k = i & 3;
  o.ctx.viewports.get(o.section)?.panBy(k < 2 ? moves[k] : 0, k >= 2 ? moves[k] : 0, "pointer");
  o.surface.renderNow();
};
const redraw = (o: Opened) => {
  o.surface.composer.redraw();
  o.surface.invalidate(7);
  o.surface.renderNow();
};

/** A pan and scroll frame over 12,000 rows at row scale `scale`. */
export async function familyFrames(scale: number) {
  const o = await opened(familyWithTree(12_000, 20, { columns: true, sites: true }), {
    width: 1000,
    height: 500,
    dpr: 2,
    still: true,
  });
  const vp = o.ctx.viewports.get(o.section);
  vp?.setView(scale === 1 ? 14 : 3, 0, "api");
  vp?.setRowScale(scale, 0);
  vp?.scrollTo(0, 2000, "api");
  redraw(o);
  return { step: (i: number) => pan(o, i), surface: o.surface, dispose: () => close(o) };
}

/** A pan and scroll frame over 20,000 lazy rows once their tiles are cached. */
export async function lazyFrames() {
  const uri = "https://example.org/lazy";
  const b = createSeqViewSpec({ title: "lazy" });
  const aln = b.data.alignment({
    ref: "lazy",
    row_count: 20_000,
    column_count: 400,
    source: { uri, protocol: "seqstar-tiles/1" },
  } as never);
  b.view({ ref: "v", label: "V", style: "msa", height: "box" })
    .block({ alignment: aln })
    .track({ label: "Rows" })
    .rows();
  const o = await opened(b.getDocument(), {
    width: 1000,
    height: 500,
    dpr: 2,
    still: true,
    loaders: lazy(uri, 20_000, 400),
  });
  const vp = o.ctx.viewports.get(o.section);
  vp?.setView(4, 0, "api");
  vp?.scrollTo(0, 3000, "api");
  for (let i = 0; i < 12; i++) {
    for (let d = 0; d < 4; d++) {
      vp?.panBy(d < 2 ? moves[d] : 0, d >= 2 ? moves[d] : 0, "pointer");
      o.surface.renderNow();
    }
    await new Promise((r) => setTimeout(r, 0));
  }
  redraw(o);
  return { step: (i: number) => pan(o, i), surface: o.surface, dispose: () => close(o) };
}

/** A hover frame over the tree: clades lit, their rows banded, without the label providers (their strings are per element). */
export async function treeHover() {
  const o = await opened(familyWithTree(160, 20, { columns: true }), {
    width: 1000,
    height: 500,
    dpr: 2,
    still: true,
    labels: false,
  });
  const l = o.ctx.layout.get(o.section)?.current;
  if (!l) throw new Error("no layout");
  const p = Inputs.pointer();
  p.section = 0;
  p.type = PointerType.move;
  let nodes = 0;
  const step = (i: number) => {
    p.x = 6 + (i % 13) * 5;
    p.y = l.scrolledTop + 4 + ((i * 7) % 61) * 5;
    o.ctx.input.pointer(p);
    if (o.ctx.hover.live.current.treeNode >= 0) nodes++;
    o.surface.renderNow();
  };
  return { step, surface: o.surface, nodes: () => nodes, dispose: () => close(o) };
}
