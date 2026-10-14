/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1b-b's heap guards (seqspec:sources#the-viewers-side; the memory rule): a pan's pass and composer frame over a
// template's 1,000 instance tracks allocate nothing once laid out (the instances are tracks like any other); a frame
// that allocates an object a step fails the same guard (the positive control); and an instance that leaves (its node
// evicted) is not retained: its resolved layer and prepared state are freed once a later update has run.

import { describe, expect, it, vi } from "vitest";
import type { Streams } from "../src/streams.ts";
import {
  controlCaught,
  describeScaling,
  FRAME_SETTLE_MS,
  FRAME_WARM,
  FRAME_WARM_MS,
  heapScaling,
  isFlat,
  leaking,
  settleFrame,
} from "./heap.ts";

const N = 1_000;
const WARM = FRAME_WARM;
const WARM_MS = FRAME_WARM_MS;
/** The file's warm-up, settle and measurement, as its guards and its control run them. */
const measure = (frame: () => void) => {
  const t0 = performance.now();
  for (let i = 0; i < WARM || performance.now() - t0 < WARM_MS; i++) frame();
  // Settled after the timed warm-up (GS-02): a climb to TurboFan can outlast it.
  settleFrame(frame, FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? "" : describeScaling(r, N);
};
const gc = (globalThis as { gc?: () => void }).gc as () => void;

/** The template document with 1,000 instances (ranges on the protein), composed at 1,000 × 400, DPR 2. */
async function templated() {
  const { composed, recording } = await import("./composer-fixtures.ts");
  const { Dirty } = await import("../src/draw/composer.ts");
  const { range, templateDoc, templateLoaders, P } = await import("./template-fixtures.ts");
  const genes = Array.from({ length: N }, (_, k) => range(`g${k}`, 1 + (k % (P.length - 8))));
  const { loaders } = templateLoaders({ genes, notes: [] });
  const t = await composed(templateDoc(), { width: 1000, height: 400, dpr: 2, loaders });
  const vp = t.vp;
  if (!vp) throw new Error("no viewport");
  vp.zoomTo(0, 20, 120);
  const be = recording();
  for (let k = 0; k < 8; k++) {
    vp.panBy(k & 1 ? -9 : 9, 0, "pointer");
    t.frame(be, Dirty.all);
  }
  let i = 0;
  return {
    t,
    frame: () => {
      vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, "pointer");
      (t.ctx.streams as Streams).pass();
      t.frame(be, Dirty.all);
    },
  };
}

describe("a template's 1,000 instance tracks", () => {
  it("a pan's pass and frame allocate nothing", async () => {
    vi.resetModules();
    const g = await templated();
    expect(g.t.ctx.graph.assembled.value?.view?.resolved.tracks.count).toBeGreaterThan(N);
    expect(guard(g.frame)).toBe("");
  }, 120_000);

  it("the same frame allocating an object a step is caught (positive control)", async () => {
    vi.resetModules();
    const g = await templated();
    expect(controlCaught(() => measure(leaking(g.frame) as () => void))).toBe(true);
  }, 120_000);
});

describe("an instance that leaves", () => {
  it("is not retained once a later update has run", async () => {
    vi.resetModules();
    const { openTemplates, range } = await import("./template-fixtures.ts");
    const t = await openTemplates({ genes: [range("a", 1)] });
    const find = (node: string) => {
      const v = t.ctx.graph.assembled.value?.view;
      for (let k = 0; v && k < v.resolved.layers.count; k++)
        if (v.resolved.layers.node[k] === node)
          return { state: v.prepared[k].state as object, themes: v.prepared[k].themes };
      return null;
    };
    const refs = await (async () => {
      t.manager.nodes.add("genes", [range("b", 9)]);
      await t.frame();
      const b = find("genes/b");
      if (!b) throw new Error("no instance");
      return { state: new WeakRef(b.state), themes: new WeakRef(b.themes) };
    })();
    t.manager.nodes.release(["genes/b"]);
    await t.frame();
    t.manager.nodes.add("genes", [range("c", 12)]);
    await t.frame();
    expect(find("genes/b")).toBeNull();
    for (let k = 0; k < 4 && refs.state.deref(); k++) {
      await new Promise((r) => setTimeout(r, 0));
      gc();
    }
    expect([refs.state.deref(), refs.themes.deref()]).toEqual([undefined, undefined]);
    t.ctx.dispose();
  }, 60_000);
});
