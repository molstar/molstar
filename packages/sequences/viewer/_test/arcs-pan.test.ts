/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A pan's frame draws the `pairs` arcs a full frame draws at that view (seqspec:format#layer-types `pairs`): the arcs'
// scale (the widest span in view, `Arcs.frame`) is the whole view's, never the uncovered strip's. The bug: the strip's
// draw computed the scale over the strip's few pairs (none: every arc flat along the ends, the sashimi's counts as
// thick bars) and kept it, so the layer's whole redraw after the copy reused it. Case 21's junctions (the gene's own
// strand and the + strand) and case 7's contacts, through pans both ways, each frame against a fresh composer's.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { case7, case21 } from "../../core/_test/format/cases.ts";
import { Dirty, type PathBatch, type RenderBackend, SectionComposer } from "../src/draw/index.ts";
import { KernelRenderers, type LayerRenderer } from "../src/index.ts";
import type { PairsState } from "../src/prepare.ts";
import type { RecordingBackend } from "../src/testing/index.ts";
import { COPY_BACK, opNames, recording } from "./composer-fixtures.ts";
import { layerOf } from "./draw-fixtures.ts";
import { opened } from "./layout-fixtures.ts";

/** The pairs renderer, its draws' subpaths noted: per recording, the [from, to) of its path items. */
const spans = new WeakMap<RenderBackend, number[]>();
const kernel = (KernelRenderers as LayerRenderer[]).find((r) => r.name === "pairs") as LayerRenderer;
const noted: LayerRenderer = {
  ...kernel,
  draw(v, s, p) {
    const be = v.backend as RecordingBackend;
    const from = be.pathItems.n;
    kernel.draw?.(v, s, p);
    const list = spans.get(be) ?? [];
    list.push(from, be.pathItems.n);
    spans.set(be, list);
  },
};

/** The pairs' arcs a recording drew, each subpath its verbs and coordinates (to 0.01 px), as a set. */
function arcsOf(r: RecordingBackend): Set<string> {
  const b: PathBatch = r.pathItems;
  const out = new Set<string>();
  const list = spans.get(r) ?? [];
  for (let j = 0; j < list.length; j += 2) {
    for (let i = list[j]; i < list[j + 1]; i++) {
      const v1 = i + 1 < b.n ? b.start[i + 1] : b.nv;
      const c1 = i + 1 < b.n ? b.cstart[i + 1] : b.nc;
      const coords = Array.from(b.coords.subarray(b.cstart[i], c1), (c) => Math.round(c * 100) / 100);
      out.add(`${b.color[i]}:${b.width[i]}:${Array.from(b.verbs.subarray(b.start[i], v1))}:${coords}`);
    }
  }
  return out;
}

/** A pairs state's kept frame forgotten (`Arcs.frame` computes it again at the next draw). */
const forget = (s: unknown) => {
  (s as PairsState).at.axis = null;
};

/** `doc` open at `view`, its sections 1000 × 400; section `at` (a ref, else its track's label) zoomed to [start, stop) of `block`. */
async function panned(doc: SeqViewSpecDoc, view: string, at: string, block: number, start: number, stop: number) {
  const t = await opened(doc, {
    renderers: [...(KernelRenderers as LayerRenderer[]).filter((r) => r !== kernel), noted],
  });
  t.ctx.views.setCurrent(view);
  await t.ctx.graph.update();
  for (const ref of t.ctx.viewports.sections.value) t.ctx.viewports.get(ref)?.setSize(1000, 400, 1);
  await t.ctx.graph.update();
  const section = t.ctx.viewports.sections.value.includes(at) ? at : layerOf(t.ctx, at).section;
  const vp = t.ctx.viewports.get(section);
  vp?.zoomTo(block, start, stop, { durationMs: 0 });
  await t.ctx.graph.update();
  const composer = new SectionComposer(t.ctx, section);
  const be = recording(COPY_BACK);
  /** A frame into the pan's backend, then the same view drawn whole by a fresh composer: their arcs. */
  const step = () => {
    be.reset();
    spans.delete(be);
    composer.frame(be, Dirty.tracks);
    // The arcs' frame is kept on the prepared state, which every composer shares: forgotten, so the full frame's own.
    const m = t.ctx.graph.assembled.value?.view;
    for (let k = 0; m && k < m.prepared.length; k++)
      if (m.prepared[k].renderer === "pairs") forget(m.prepared[k].state);
    const whole = recording(COPY_BACK);
    const fresh = new SectionComposer(t.ctx, section);
    fresh.frame(whole, Dirty.all);
    fresh.dispose();
    return { pan: arcsOf(be), full: arcsOf(whole), drew: be.opCount > 0, copied: opNames(be).includes("copy") };
  };
  step();
  return { vp, step };
}

/** Pans by `dx` CSS px `n` times each way: every frame's arcs are the full frame's, and some frames copied. */
async function sweep(p: Awaited<ReturnType<typeof panned>>, dx: number, n: number) {
  let copies = 0;
  let arcs = 0;
  for (const dir of [1, -1]) {
    for (let i = 0; i < n; i++) {
      p.vp?.panBy(dir * dx, 0, "api");
      const { pan, full, drew, copied } = p.step();
      // At the axis's end the pan moves nothing, and nothing is drawn.
      if (!drew) continue;
      if (copied) copies++;
      arcs += full.size;
      // Every arc drawn in the pan's frame (its strip's and its whole redraw's) is one the full frame draws, and the
      // whole redraw draws them all.
      expect([...pan].filter((a) => !full.has(a))).toEqual([]);
      expect([...full].filter((a) => !pan.has(a))).toEqual([]);
    }
  }
  expect(copies).toBeGreaterThan(n);
  expect(arcs).toBeGreaterThan(0);
}

describe("a pan's frame draws the pairs' arcs of a full frame at its view", () => {
  it("case 21's junctions on the gene's own strand", async () => {
    await sweep(await panned(case21(), "gene2", "gene2.own", 0, 300, 500), 23, 16);
  });
  it("case 21's junctions on the region as written", async () => {
    await sweep(await panned(case21(), "gene2", "gene2.reference", 0, 300, 600), 31, 20);
  });
  it("case 7's contacts, zoomed into barnase", async () => {
    await sweep(await panned(case7(), "complex-one", "Contacts", 0, 20, 70), 19, 24);
  });
});
