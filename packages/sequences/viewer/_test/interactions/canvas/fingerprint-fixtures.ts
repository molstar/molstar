/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Frames of the fingerprint through the Canvas2D backend for the browser heap guards (fingerprint.heap.ts), imported
// afresh per case: a scene's context (inline data only: no Node files), the layer's `LayerFrame`, and a step that
// draws it into a canvas, with the kit's composer-in-small filling the view as the real one does.

import { Canvas2DBackend, Canvas2DShared } from "../../../src/canvas/index.ts";
import { Pass, SkinToken } from "../../../src/draw/backend.ts";
import { drawResources } from "../../../src/draw/resources.ts";
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  type LayerRenderer,
  manualScheduler,
} from "../../../src/index.ts";
import { InteractionsViewerRegistry } from "../../../src/interactions/viewer/index.ts";
import { fixtureMeasurer, LayerFrame, Op, RecordingBackend } from "../../../src/testing/index.ts";
import { absentScene, manyScene, mixedScene, type Scene, ShortFingerprint, unitScene, wideScene } from "../scenes.ts";
export const SCENES: Record<string, () => Scene> = {
  unit: () => unitScene(),
  many: () => manyScene(),
  compact: () => manyScene(20, 8, 1),
  mixed: () => mixedScene(),
  absent: () => absentScene(),
  wide: () => wideScene(),
};

/** The scene drawn by `entry` (the fingerprint, or a wrapper of it) into a canvas at 2× per step. */
export async function fingerprintFrames(
  name: keyof typeof SCENES,
  wrap?: (e: LayerRenderer) => LayerRenderer,
  o: { dpr?: number; back?: number } = {},
) {
  const dpr = o.dpr ?? 2;
  const scene = SCENES[name]();
  const spec = { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, InteractionsViewerRegistry, ShortFingerprint] };
  const ctx = createSeqViewer(spec, { scheduler: manualScheduler(), text: fixtureMeasurer() });
  const r = await ctx.setDocument(scene.doc);
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  ctx.views.setCurrent(scene.section.split(".")[0]);
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(900, 500, dpr);
  await ctx.graph.update();
  const m = ctx.graph.assembled.value?.view;
  const k = m?.prepared.findIndex((p) => p.renderer === "interactions:fingerprint") ?? -1;
  if (k < 0) throw new Error("no fingerprint layer");
  const entry = ctx.registry.get("renderers", "interactions:fingerprint") as LayerRenderer;
  const frame = new LayerFrame(ctx, wrap ? wrap(entry) : entry, scene.section, k);
  const canvas = document.createElement("canvas");
  canvas.width = 900 * dpr;
  canvas.height = 500 * dpr;
  const shared = new Canvas2DShared(drawResources(ctx));
  const backend = new Canvas2DBackend(
    [canvas.getContext("2d", { willReadFrequently: true })],
    shared,
    o.back ?? 1 << Pass.main,
  );
  // Mixed: its three blocks of 30 in view (the third the tiles'); wide: cells under a device px.
  const cell = name === "wide" ? 0.05 : name === "mixed" ? 8 : 24;
  if (name !== "unit") ctx.viewports.get(scene.section)?.setView(cell, 0, "api");
  /** What one frame draws, recorded: its rects, the grounds among them, and the tiles' images. */
  const drawn = () => {
    const rec = new RecordingBackend();
    frame.draw(rec);
    const b = rec.rectItems;
    let grounds = 0;
    for (let i = 0; i < b.n; i++) if (b.color[i] === SkinToken.grid) grounds++;
    return { rects: b.n, grounds, images: rec.count(Op.image) };
  };
  return {
    step: () => frame.draw(backend),
    drawn,
    dispose: () => ctx.dispose(),
    frame,
    ctx,
    scene,
    canvas,
    k,
  };
}
