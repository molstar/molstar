/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3a: a letters stream's data drawn from the region cache (seqspec:sources#the-viewers-side, `StreamView.letters`;
// seqspec:sources#levels-and-the-zoom-gate), over the streamed genome with a sequence track
// (core/_test/format/stream-provider.ts: synthetic 2bit-like letters, N runs and soft-masked stretches). Its data is asked
// only where a device px holds at most `LETTERS_PER_PX` positions: none at the whole view or at 100 kb, the view's
// regions and then the margin at 10 kb; there `state` is empty past the gate and pending, then loaded, inside it. The
// read gives the upper-case codes held and 0 elsewhere; at letter zoom the letters renderer draws exactly the region's
// letters, thin cells below it, and nothing where no region is held.

import { type StreamedData, StreamProp } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { GRCH38, letterAt } from "../../core/_test/format/stream-provider.ts";
import { RegionState } from "../src/index.ts";
import { LETTERS_PER_PX, type Streams } from "../src/streams.ts";
import { composed, recording, textsOf } from "./composer-fixtures.ts";
import { GENOME, MaxInFlight, openGenome, type StreamedGenome, streamedDoc } from "./stream-fixtures.ts";

const flush = () => new Promise((ok) => setTimeout(ok, 0));
const THREE = GRCH38.slice(0, 3);
const CHR2 = 1;
/** Inside chr2, clear of the fixture's N runs (each 100 kb's first 1 kb). */
const AT = 100_040_100;
const SPAN = 4096;
const R0 = Math.floor(AT / SPAN);
const upper = (s: number, x: number) => letterAt(s, x).toUpperCase().charCodeAt(0);

/** The letters layer's ref in the shown view. */
function lettersRef(t: { ctx: StreamedGenome["ctx"] }): string {
  const m = t.ctx.graph.assembled.value?.view;
  for (let k = 0; m && k < m.resolved.layers.count; k++)
    if (m.set.strings.get(m.resolved.layers.type[k]) === "letters") return m.set.strings.get(m.resolved.layers.ref[k]);
  throw new Error("no letters layer");
}
async function genome(gated = true) {
  const t = await openGenome({ chromosomes: THREE, gated, letters: true, doc: streamedDoc({ sequence: true }) });
  // Every region asked in one pass (the values' share the bound).
  t.ctx.config.set(MaxInFlight, 64);
  return t;
}
function show(t: StreamedGenome, b: number, start: number, stop: number) {
  const vp = t.vp(GENOME);
  if (stop > start) vp.zoomTo(b, start, stop);
  else vp.fit();
  t.layout(GENOME);
  t.manager.pass();
}
async function drain(t: StreamedGenome) {
  for (let k = 0; k < 6; k++) {
    t.streams.genome.release();
    t.streams.phylop.release();
    await flush();
    t.manager.pass();
  }
}

describe("requests: the data only, within the letters gate", () => {
  it("none at the whole view or at 100 kb; at 10 kb the view's regions from the center out, then the margin", async () => {
    const t = await genome();
    const ref = lettersRef(t);
    show(t, 0, 0, 0);
    show(t, CHR2, AT, AT + 100_000);
    expect(t.streams.genome.requests).toEqual([]);
    expect(t.manager.state(ref, CHR2, AT)).toBe(RegionState.empty);
    show(t, CHR2, AT, AT + 10_000);
    const asked = t.streams.genome.requests;
    expect(asked.every((r) => r.level === -1 && r.system === "chr2" && r.stop - r.start === SPAN)).toBe(true);
    const view = asked.filter((r) => r.priority === "view").map((r) => r.start / SPAN);
    const margin = asked.filter((r) => r.priority === "margin").map((r) => r.start / SPAN);
    // The regions over [AT, AT + 10 kb), the one at its middle first; the margin a view width each side.
    const r = (x: number) => Math.floor(x / SPAN);
    expect(view[0]).toBe(r(AT + 5000));
    expect(view.slice().sort()).toEqual(Array.from({ length: r(AT + 9999) - R0 + 1 }, (_, k) => R0 + k));
    expect(Math.min(...margin)).toBe(r(AT - 10_000));
    expect(Math.max(...margin)).toBe(r(AT + 19_999));
    expect(asked.findIndex((r) => r.priority === "margin")).toBe(view.length);
    expect(t.manager.state(ref, CHR2, AT)).toBe(RegionState.pending);
    await drain(t);
    expect(t.manager.state(ref, CHR2, AT)).toBe(RegionState.loaded);
    // Asked once: a pass over the held view asks nothing more.
    const n = t.streams.genome.requests.length;
    show(t, CHR2, AT, AT + 10_000);
    expect(t.streams.genome.requests.length).toBe(n);
  });

  it("`settled` waits for the view's letters", async () => {
    const t = await genome();
    show(t, CHR2, AT, AT + 10_000);
    let done = false;
    void t.manager.settled().then(() => {
      done = true;
    });
    t.streams.phylop.release();
    await flush();
    t.manager.pass();
    await flush();
    expect(done).toBe(false);
    await drain(t);
    await flush();
    expect(done).toBe(true);
  });
});

describe("`StreamView.letters`", () => {
  it("the held codes upper case, 0 where no region is held; the count held", async () => {
    const t = await genome();
    const ref = lettersRef(t);
    show(t, CHR2, AT, AT + 10_000);
    await drain(t);
    const st = (t.manager as Streams).get(ref);
    const out = new Uint8Array(300);
    // Across a region's edge, soft-masked letters included.
    const i0 = (R0 + 1) * SPAN - 150;
    expect(st?.letters?.(CHR2, i0, i0 + 300, out)).toBe(300);
    expect(Array.from(out)).toEqual(Array.from({ length: 300 }, (_, j) => upper(CHR2, i0 + j)));
    expect(new Set(Array.from(out, (c) => String.fromCharCode(c)))).toEqual(new Set("ACGT"));
    // Far from the view: nothing held, zeros written.
    out.fill(7);
    expect(st?.letters?.(CHR2, 50_000_000, 50_000_300, out)).toBe(0);
    expect(out.every((c) => c === 0)).toBe(true);
    // A block that shows another system's letters reads its own.
    expect(st?.letters?.(0, AT, AT + 10, out)).toBe(0);
  });
});

describe("drawn from the cache", () => {
  async function drawn(stop: number, letters = true) {
    const { hg38Streams } = await import("../../core/_test/format/stream-provider.ts");
    const s = hg38Streams({ chromosomes: THREE, letters });
    const loaders = { load: () => Promise.reject(new Error("none")), stream: s.stream };
    const t = await composed(streamedDoc({ sequence: true }), { width: 1000, height: 200, loaders });
    const m = t.ctx.streams as Streams;
    t.vp?.zoomTo(CHR2, AT, AT + stop);
    t.layout(GENOME);
    for (let k = 0; k < 6; k++) {
      m.pass();
      await flush();
    }
    return { t, s, frame: t.frame(recording()) };
  }

  it("at letter zoom, the region's letters exactly; none before the data arrives", async () => {
    const { frame } = await drawn(60);
    // The view's 60 positions and the one the edge cuts.
    const letters = textsOf(frame).filter((x) => /^[ACGTN]$/.test(x));
    expect(letters.length).toBe(61);
    expect(letters.join("")).toBe(Array.from({ length: 61 }, (_, j) => letterAt(CHR2, AT + j).toUpperCase()).join(""));
    const none = await drawn(60, false);
    expect(textsOf(none.frame).filter((x) => /^[ACGTN]$/.test(x))).toEqual([]);
  });

  it("below letter zoom, cells; past the gate nothing is asked or drawn", async () => {
    const thin = await drawn(10_000);
    const failed = await drawn(10_000, false);
    expect(thin.frame.rectItems.n).toBeGreaterThan(failed.frame.rectItems.n + 500);
    const wide = await drawn(1000 * LETTERS_PER_PX * 2);
    expect(wide.s.genome.requests).toEqual([]);
  });
});

describe("the hover label", () => {
  it("names a streamed position's letter once its region is held", async () => {
    const { driven } = await import("./input-fixtures.ts");
    const { hg38Streams } = await import("../../core/_test/format/stream-provider.ts");
    const { KernelRenderers } = await import("../src/index.ts");
    const s = hg38Streams({ chromosomes: THREE, letters: true, gated: true });
    const t = await driven(streamedDoc({ sequence: true }), {
      renderers: [...KernelRenderers] as never,
      loaders: { load: () => Promise.reject(new Error("none")), stream: s.stream },
      size: { [GENOME]: [1000, 200] },
    });
    t.vp(GENOME).zoomTo(CHR2, AT, AT + 60);
    await t.ctx.graph.update();
    const l = t.layout(GENOME);
    const m = t.ctx.graph.assembled.value?.view;
    const h = t.ctx.hover.live.current;
    // Down the data area to the sequence track.
    const label = (k: number) => m?.set.strings.get(m.resolved.tracks.label[k]);
    for (let y = 1; y < 200 && (h.track < 0 || label(h.track) !== "Sequence"); y += 2)
      t.pointer("move", GENOME, l.dataLeft + 5, y);
    const i = h.seq.index as number;
    const title = () => t.ctx.labels.current.value?.title ?? "";
    expect(i).toBeGreaterThanOrEqual(AT);
    // Not held: the position's number, no letter.
    expect(title()).toMatch(/^chr2 \d/);
    for (let k = 0; k < 4; k++) {
      (t.ctx.streams as Streams).pass();
      s.genome.release();
      s.phylop.release();
      await flush();
    }
    // Made again for the arrival: the pointer never moved.
    expect(title()).toBe(`chr2 ${letterAt(CHR2, i).toUpperCase()}${i + 1}`);
    t.ctx.dispose();
  });
});

describe("aside layers (a codon label's letters)", () => {
  it("two that no lazy molecule reads any longer close in one update, and the update completes (VT-01)", async () => {
    const t = await openGenome({ chromosomes: THREE, letters: true, doc: streamedDoc({ sequence: true }) });
    t.vp(GENOME).zoomTo(1, 100_040_100, 100_050_100);
    t.layout(GENOME);
    t.manager.pass();
    const model = t.ctx.graph.assembled.value?.view ?? null;
    let data: StreamedData | undefined;
    for (const m of model?.set.molecules ?? []) if (m?.lazy) data = m.custom.props.get(StreamProp) as StreamedData;
    if (!data) throw new Error("no lazy molecule");
    // Two streams no molecule of the set reads (another document's, as a hover's codon label before a switch).
    const at = Int32Array.of(5);
    const out = new Uint8Array(1);
    t.manager.lettersAt({ ...data }, 0, at, 1, out);
    t.manager.lettersAt({ ...data }, 0, at, 1, out);
    const inner = t.manager as unknown as { asides: unknown[]; open: unknown[] };
    const open = inner.open.length;
    expect(inner.asides.length).toBe(2);
    expect(() => t.manager.update(model)).not.toThrow();
    expect([inner.asides.length, inner.open.length]).toEqual([0, open - 2]);
    t.ctx.dispose();
  });
});
