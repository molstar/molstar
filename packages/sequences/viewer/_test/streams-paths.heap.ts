/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST7's heap guard (seqspec:sources#requests; the memory rule): a protein's section showing its genome's streamed
// genes, variants and phyloP through its codon map, regions held. A pan's pass (the view's residues taken back through
// the path to its exons' runs, the regions asked) and its frame (items read back on the residues, values over codons,
// a split codon's two parts, placeholders mapped through the path) allocate nothing, over 1,000 residues (three exons)
// and at 100 residues across a split codon; so does the pointer over the genes and variants there (`hitStreamed`
// through the path, every region over a split codon tried), and so do the codon label's reads (a residue's positions,
// `Paths.fromAt`, and their letters held, `lettersAt`) and a pass with letters wanted (the aside layer in its order). A frame allocating an object a step fails the same guard (the
// positive control).

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
/** Six times the shared warm-up, as the other streamed items' guards: colored items tier up slowly. */
const WARM = 6 * FRAME_WARM;
/** The file's warm-up, settle and measurement, as its guards and its control run them. */
const measure = (frame: () => void) => {
  const t0 = performance.now();
  for (let i = 0; i < WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame();
  // Settled after the timed warm-up (GS-02): a climb to TurboFan can outlast it.
  settleFrame(frame, FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? "" : describeScaling(r, N);
};

/** The protein's section (DPR 2, 1,000 × 300) at residues [start, stop) (0, 0: the whole protein), every region held. */
async function proteinFrame(start: number, stop: number) {
  const { composed, recording } = await import("./composer-fixtures.ts");
  const { Dirty } = await import("../src/draw/composer.ts");
  const { GRCH38, codingDoc, hg38Streams } = await import("../../core/_test/format/stream-provider.ts");
  const s = hg38Streams({ chromosomes: GRCH38.slice(0, 3), gated: true });
  const loaders = { load: () => Promise.reject(new Error("none")), stream: s.stream };
  const section = "coding.protein";
  const t = await composed(codingDoc(), { width: 1000, height: 300, dpr: 2, loaders, section });
  const m = t.ctx.streams as import("../src/streams.ts").Streams;
  const vp = t.vp;
  if (!vp) throw new Error("no viewport");
  if (stop > start) vp.zoomTo(0, start, stop);
  else vp.fit();
  const be = recording();
  const asked = () => s.genes.requests.length + s.variants.requests.length + s.phylop.requests.length;
  for (let k = 0, quiet = 0; k < 60 && quiet < 2; k++) {
    const before = asked();
    for (const dx of [9, -9]) {
      vp.panBy(dx, 0, "pointer");
      t.layout(section);
      m.pass();
      for (const p of [s.genome, s.phylop, s.genes, s.variants]) p.release();
      await new Promise((ok) => setTimeout(ok, 0));
    }
    t.frame(be, Dirty.all);
    quiet = asked() === before ? quiet + 1 : 0;
  }
  let i = 0;
  return {
    backend: be,
    asked,
    dispose: () => t.ctx.dispose(),
    frame: () => {
      vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, "pointer");
      m.pass();
      t.frame(be, Dirty.all);
    },
  };
}

/** The protein's section in the input harness at residues 450–550, regions held: a pointer move over its tracks. */
async function proteinHover() {
  const { driven } = await import("./input-fixtures.ts");
  const { KernelRenderers } = await import("../src/index.ts");
  const { GRCH38, codingDoc, hg38Streams } = await import("../../core/_test/format/stream-provider.ts");
  const s = hg38Streams({ chromosomes: GRCH38.slice(0, 3), gated: true, letters: true });
  const section = "coding.protein";
  const t = await driven(codingDoc(), {
    renderers: [...KernelRenderers] as never,
    loaders: { load: () => Promise.reject(new Error("none")), stream: s.stream },
    size: { [section]: [1000, 300] },
    labels: false,
  });
  t.vp(section).zoomTo(0, 450, 550);
  await t.ctx.graph.update();
  for (let k = 0; k < 6; k++) {
    t.layout(section);
    (t.ctx.streams as import("../src/streams.ts").Streams).pass();
    for (const p of [s.genome, s.phylop, s.genes, s.variants]) p.release();
    await new Promise((ok) => setTimeout(ok, 0));
  }
  const l = t.layout(section);
  const release = () => {
    for (const p of [s.genome, s.phylop, s.genes, s.variants]) p.release();
  };
  // One `PointerInput`, as an engine's (the harness's `pointer` makes an event a call).
  const { Inputs, PointerType } = await import("../src/index.ts");
  const e = Inputs.pointer();
  const at = t.ctx.layout.sectionIndex(section);
  const move = PointerType.move;
  return {
    t,
    release,
    move: (i: number) => {
      e.type = move;
      e.section = at;
      e.x = l.dataLeft + ((i * 37) % 900);
      e.y = 40 + ((i * 13) % 120);
      t.ctx.input.pointer(e);
    },
  };
}

describe("a genome's streamed tracks on a protein's axis allocate nothing per frame once held", () => {
  for (const [label, start, stop] of [
    ["1,000 residues", 100, 1100],
    ["residues 450–550", 450, 550],
  ] as const)
    it(`a pass and a frame, ${label}`, async () => {
      vi.resetModules();
      const g = await proteinFrame(start, stop);
      g.frame();
      // Not vacuous: values and items drawn (rects).
      expect(g.backend.rectItems.n).toBeGreaterThan(10);
      const asked = g.asked();
      expect(asked).toBeGreaterThan(0);
      expect(guard(g.frame)).toBe("");
      expect(g.asked()).toBe(asked);
      g.dispose();
    }, 120_000);

  it("the pointer over the genes and variants (a hit through the path)", async () => {
    vi.resetModules();
    const g = await proteinHover();
    let parts = 0;
    for (let i = 0; i < 400; i++) {
      g.move(i);
      if (g.t.ctx.hover.live.current.part >= 0) parts++;
    }
    // Not vacuous: items hit.
    expect(parts).toBeGreaterThan(10);
    // The same 400 moves over and over (a hover that changes item each move, as a sweep does), after 20,000 of them.
    let i = 0;
    for (let w = 0; w < 20_000; w++) g.move(i++ % 400);
    expect(guard(() => g.move(i++ % 400))).toBe("");
    g.t.ctx.dispose();
  }, 120_000);

  it("the codon label's reads: a residue's positions through the path and their letters, held", async () => {
    vi.resetModules();
    const g = await proteinHover();
    const { Paths, CoordinateSystemKind } = await import("@molstar/seq-core/model");
    const { StreamProp } = await import("@molstar/seq-core/format");
    const m = g.t.ctx.streams;
    const v = g.t.ctx.graph.assembled.value?.view;
    let ref = "";
    for (let k = 0; v && k < v.resolved.layers.count; k++)
      if (v.set.strings.get(v.resolved.layers.type[k]) === "values") ref = v.set.strings.get(v.resolved.layers.ref[k]);
    const p = m.get(ref)?.source.paths?.[0];
    const set = v?.set;
    if (!p || !set || p.from.kind !== CoordinateSystemKind.unit) throw new Error("no path");
    const mol = set.molecules[set.units.molecule[p.from.id]];
    const data = mol?.custom.props.get(StreamProp) as never;
    const sys = (data as { systems: Map<string, number> }).systems.get(set.strings.get(mol?.key ?? -1)) ?? -1;
    const at = new Int32Array(3);
    const out = new Uint8Array(3);
    let held = 0;
    const read = (i: number) => {
      const r = 450 + (i % 100);
      for (let k = 0; k < 3; k++) at[k] = Paths.fromAt(p, r, k);
      held += m.lettersAt?.(data, sys, at, 3, out) ?? 0;
    };
    // Each read asks for its own codon's letters (the last read's replace the earlier's); the passes bring them.
    for (const i of [0, 99]) {
      read(i);
      for (let k = 0; k < 3; k++) {
        (m as import("../src/streams.ts").Streams).pass();
        g.release();
        await new Promise((ok) => setTimeout(ok, 0));
      }
    }
    held = 0;
    for (let i = 0; i < 100; i++) read(i);
    // Not vacuous: every base of the 100 residues held.
    expect(held).toBe(300);
    for (let w = 0; w < 20_000; w++) read(w);
    let j = 0;
    expect(guard(() => read(j++))).toBe("");
    // A pass with the aside layer's letters wanted (the third exon's, asked, never answered): pushed into the pass's
    // order, asked, the order cut back.
    const far = (i: number) => {
      const r = 950 + (i % 40);
      for (let k = 0; k < 3; k++) at[k] = Paths.fromAt(p, r, k);
      held += m.lettersAt?.(data, sys, at, 3, out) ?? 0;
      (m as import("../src/streams.ts").Streams).pass();
    };
    for (let w = 0; w < 2_000; w++) far(w);
    expect(guard(() => far(j++))).toBe("");
    g.t.ctx.dispose();
  }, 120_000);

  it("the same frame allocating an object a step is caught (positive control)", async () => {
    vi.resetModules();
    const g = await proteinFrame(450, 550);
    expect(controlCaught(() => measure(leaking(g.frame) as () => void))).toBe(true);
    g.dispose();
  }, 120_000);
});
