/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ideogram navigator's heap guard (ST5; the memory rule): its strip, redrawn on every frame here (the nav pass is
// drawn by hand, as a layout or skin change would), and the window moved over it allocate nothing once warm. Over hg38's
// committed cytobands on the streamed genome's first three chromosomes (about 200 bands). The positive control: the same
// frames allocating an object each (`leaking`), through the same warm-up and measurement, are caught.

import { describe, expect, it, vi } from "vitest";
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
/** The strip's frame tiers up slowly, as the streamed items' do (streams-items.heap.ts): warmed six times as long. */
const WARM = 6 * FRAME_WARM;
/** The guard's measurement: warmed by count and by time, then `N` frames and `4N` measured. */
const measure = (frame: (i: number) => void) => {
  const t0 = performance.now();
  let i = 0;
  for (; i < WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  // Settled after the timed warm-up (GS-02): a climb to TurboFan can outlast it.
  settleFrame(() => frame(i++), FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? "" : describeScaling(r, N);
};

async function ideogramFrame() {
  const { composed, recording } = await import("./composer-fixtures.ts");
  const { GRCH38, hg38Streams, streamedDoc } = await import("../../core/_test/format/stream-provider.ts");
  const { ideogramNavigator } = await import("../src/index.ts");
  const { readFileSync } = await import("node:fs");
  const { GENOME } = await import("./stream-fixtures.ts");
  const bands = JSON.parse(readFileSync(new URL("./fixtures/hg38-cytobands.json", import.meta.url), "utf8"));
  const s = hg38Streams({ chromosomes: GRCH38.slice(0, 3), letters: true });
  const t = await composed(streamedDoc({ values: false }), {
    width: 1000,
    height: 300,
    dpr: 2,
    loaders: { load: () => Promise.reject(new Error("none")), stream: s.stream },
    registry: [
      {
        navigators: [ideogramNavigator(bands.chromosomes)],
        presets: [{ name: "karyo", geometry: "linear", render: {}, navigator: "ideogram" }],
      },
    ],
  });
  t.ctx.ui.setStyle(GENOME, "karyo");
  await t.ctx.graph.update();
  t.vp?.zoomTo(0, 100_000_000, 100_050_000);
  // The nav pass drawn by hand, as the composer does on a layout or skin change (a pan never redraws the strip).
  const { NavView, navDraw, navDrawWindow } = await import("../src/draw/nav.ts");
  const av = t.ctx.layout.axis(GENOME);
  const entry = t.ctx.layout.navigatorOf(0);
  if (!av || !entry) throw new Error("no navigator");
  const n = new NavView().from(av);
  const be = recording();
  n.backend = be;
  n.text = t.ctx.text;
  n.parts = entry.parts;
  n.dpr = 2;
  n.top = 0;
  n.height = 24;
  n.gripTop = 0;
  n.gripBottom = 24;
  n.names = av.axis?.names ?? [];
  n.from(av);
  let i = 0;
  const at = Float64Array.of(av.navX0, av.navX0 + 9, av.navX0 + 40, av.navX0 + 49);
  // The window moves between two places over the strip; the strip is redrawn each time.
  const frame = () => {
    n.x0 = at[i & 1];
    n.x1 = at[2 + (i++ & 1)];
    be.reset();
    navDraw(entry, n);
    navDrawWindow(entry, n);
  };
  for (let k = 0; k < 8; k++) frame();
  return { frame, backend: be, dispose: () => t.ctx.dispose() };
}

describe("the ideogram navigator allocates nothing per frame", () => {
  it("its strip and the window over it, three chromosomes", async () => {
    vi.resetModules();
    const g = await ideogramFrame();
    g.frame();
    // Not vacuous: the bands, the centromeres' triangles and the window.
    expect(g.backend.rectItems.n).toBeGreaterThan(80);
    expect(g.backend.pathItems.nv).toBeGreaterThanOrEqual(24);
    expect(guard(g.frame)).toBe("");
    g.dispose();
  }, 120_000);

  it("the same frames allocating an object each are caught (positive control)", async () => {
    vi.resetModules();
    const g = await ideogramFrame();
    expect(controlCaught(() => measure(leaking(g.frame)))).toBe(true);
    g.dispose();
  }, 120_000);
});
