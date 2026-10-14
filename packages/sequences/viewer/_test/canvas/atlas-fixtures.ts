/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Frames of a lazy alignment's lettered rows for the atlas's browser heap guard (atlas.heap.ts), imported afresh per
// case: the composer on the Canvas2D backend (a back target, as the surface's main pass), the letters from the glyph
// atlas, the tiles from the in-process tile server, arrived before the first step. A step is a full tracks frame
// (`full`: every letter a cell hit) or a pan through the copy (`pan`: the band copied, its strip's letters drawn), a
// pointer's or (`api`) a program's, committed at once.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { tilesStream } from "@molstar/seq-core/stream/tiles";
import { tileServer } from "@molstar/seq-core/stream/tiles/server";
import { Canvas2DBackend, Canvas2DShared } from "../../src/canvas/index.js";
import { Dirty, SectionComposer } from "../../src/draw/index.js";
import { drawResources } from "../../src/draw/resources.js";
import { createSeqViewer, DefaultSeqViewerSpec, manualScheduler } from "../../src/index.js";
import type { Streams } from "../../src/streams.js";
import { fixtureMeasurer } from "../../src/testing/index.js";

const URI = "https://example.org/lazy";

export async function atlasFrames(kind: "full" | "pan" | "api") {
  const server = tileServer({ rowCount: 5000, columnCount: 4000 });
  const ctx = createSeqViewer(DefaultSeqViewerSpec(), {
    scheduler: manualScheduler(),
    text: fixtureMeasurer(),
    loaders: {
      load: () => Promise.reject(new Error("no fixture")),
      stream: tilesStream((u) => (u === URI ? server.transport : null)),
    },
  });
  const b = createSeqViewSpec();
  const source = { uri: URI, protocol: "seqstar-tiles/1" };
  const aln = b.data.alignment({ ref: "lazy", row_count: 5000, column_count: 4000, source } as never);
  const v = b.view({ ref: "g", label: "G", height: "box" });
  v.block({ alignment: aln });
  v.track({ label: "Rows" }).rows();
  const [W, H, dpr] = [1100, 600, 2];
  const done = ctx.setDocument(b.getDocument());
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(W, H, dpr);
  const r = await done;
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  await ctx.graph.update();
  const section = ctx.viewports.sections.value[0];
  const vp = ctx.viewports.get(section);
  if (!vp) throw new Error("no viewport");
  vp.setView(14, 0, "api");
  vp.setRowScale(1, 0);
  vp.panBy(20_000, 1_000, "api");
  const c = new OffscreenCanvas(W * dpr, H * dpr);
  const shared = new Canvas2DShared(drawResources(ctx));
  const be = new Canvas2DBackend([c.getContext("2d") as OffscreenCanvasRenderingContext2D], shared);
  const composer = new SectionComposer(ctx, section);
  // Every region in view and on either side of the pans arrived.
  for (const dx of [0, 37, -37, 0])
    for (let k = 0; k < 6; k++) {
      vp.panBy(dx, 0, "api");
      (ctx.streams as Streams).pass();
      await new Promise((ok) => setTimeout(ok, 0));
      composer.redraw();
      composer.frame(be, Dirty.tracks);
    }
  const step =
    kind === "full"
      ? () => {
          composer.redraw();
          composer.frame(be, Dirty.tracks);
        }
      : kind === "pan"
        ? (i: number) => {
            // A pointer's pan, as a drag's.
            vp.panBy(i % 2 ? 37 : -37, 0, "pointer");
            composer.frame(be, Dirty.tracks);
          }
        : (i: number) => {
            // A program's (a player, a link): committed each move, its record built only when read.
            vp.panBy(i % 2 ? 37 : -37, 0, "api");
            composer.frame(be, Dirty.tracks);
          };
  return { step, composer, shared, dispose: () => ctx.dispose() };
}
