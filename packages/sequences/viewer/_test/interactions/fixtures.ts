/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The scenes (scenes.ts) opened in a viewer context for the Node tests, and their fingerprint layer.

import { KernelRenderers, type SeqViewerContext, type SeqViewerRegistryEntry } from "../../src/index.ts";
import { InteractionsViewerRegistry } from "../../src/interactions/viewer/index.ts";
import { opened } from "../layout-fixtures.ts";
import { type Scene, ShortFingerprint } from "./scenes.ts";

export * from "./scenes.ts";

/**
 * A context showing the scene's section with the kernel renderers and, by default, the extension's entry; `dpr` sizes the
 * section at that pixel ratio.
 */
export async function open(
  scene: Scene,
  o: { registry?: SeqViewerRegistryEntry[]; dpr?: number; size?: readonly [number, number]; render?: object } = {},
) {
  const t = await opened(scene.doc, {
    renderers: [...KernelRenderers],
    registry: o.registry ?? [InteractionsViewerRegistry, ShortFingerprint],
    view: scene.section.split(".")[0],
    render: o.render,
  });
  const [w, h] = o.size ?? [1000, 400];
  t.ctx.viewports.get(scene.section)?.setSize(w, h, o.dpr ?? 1);
  await t.ctx.graph.update();
  return t;
}

/** The fingerprint layer: its index and prepared state. */
export function layerOf(ctx: SeqViewerContext) {
  const m = ctx.graph.assembled.value?.view;
  const k = m?.prepared.findIndex((p) => p.renderer === "interactions:fingerprint") ?? -1;
  return { k, m, state: m?.prepared[k]?.state };
}
