/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST6: prefetch (seqspec:sources#the-viewers-side, seqspec:viewer#snapshots-and-playback), headless over the streamed
// genome's values (a scripted provider: answers held until released, every request and abort recorded). Past the
// margin, a view width toward the side the view last panned to, asked after the view's and the margin's with `priority:
// "prefetch"`, kept while the view rests and moved with it; a zoom's next level over the view; the focus regions of the
// snapshots next to the one shown; a view's or margin's request takes a prefetch's slot; `Stream.Prefetch` off asks
// none; and a prefetch is never waited for (`settled`).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { GRCH38, hg38Streams, streamedDoc } from "../../core/_test/format/stream-provider.js";
import { SeqViewerConfig } from "../src/index.js";
import { snapshotHints } from "../src/prefetch.js";
import { composed } from "./composer-fixtures.js";
import { GENOME, MaxInFlight, openGenome, type StreamedGenome } from "./stream-fixtures.js";

const flush = () => new Promise((ok) => setTimeout(ok, 0));
const THREE = GRCH38.slice(0, 3);
const CHR2 = 1;
const Prefetch = SeqViewerConfig.Stream.Prefetch;

/** The window [start, stop) of block b in view, then one pass. */
function show(t: StreamedGenome, b: number, start: number, stop: number) {
  t.vp(GENOME).zoomTo(b, start, stop);
  t.layout(GENOME);
  t.manager.pass();
}
const data = (t: StreamedGenome, priority: string) =>
  t.streams.phylop.requests.filter((r) => r.level === -1 && r.priority === priority).map((r) => r.start);

async function opened(prefetch = true) {
  const t = await openGenome({ gated: true, chromosomes: THREE });
  t.ctx.config.set(MaxInFlight, 100);
  t.ctx.config.set(Prefetch, prefetch);
  return t;
}
const AT = 100_000_000;

describe("the lead: past the margin, toward where the view last panned", () => {
  it("nothing before the view has moved; then a view width past the margin on the side it moved toward, after the margin's", async () => {
    const t = await opened();
    show(t, CHR2, AT, AT + 10_000);
    expect(data(t, "prefetch")).toEqual([]);
    // A pan right by half a view: the lead is [hi + w, hi + 2w) at the view's level, last in the order asked.
    show(t, CHR2, AT + 5_000, AT + 15_000);
    const v = t.ctx.layout.axis(GENOME);
    const [lo, hi] = [v?.visLo[CHR2] ?? 0, v?.visHi[CHR2] ?? 0];
    const w = hi - lo;
    const lead = data(t, "prefetch");
    expect(lead.length).toBeGreaterThan(0);
    // Every region of it lies past the margin on the right, and none on the left.
    for (const at of lead) {
      expect(at + 4096).toBeGreaterThan(hi + w);
      expect(at).toBeLessThan(hi + 2 * w + 4096);
    }
    const order = t.streams.phylop.requests.map((r) => r.priority);
    expect(order.lastIndexOf("margin")).toBeLessThan(order.indexOf("prefetch"));
    expect(order.lastIndexOf("view")).toBeLessThan(order.indexOf("prefetch"));
    t.ctx.dispose();
  });

  it("kept while the view rests (not aborted); moved to the other side by a pan back, the old one aborted", async () => {
    const t = await opened();
    show(t, CHR2, AT, AT + 10_000);
    show(t, CHR2, AT + 5_000, AT + 15_000);
    const right = data(t, "prefetch");
    const aborted = t.streams.phylop.aborted.length;
    t.manager.pass();
    t.manager.pass();
    expect(t.streams.phylop.aborted.length).toBe(aborted);
    // Back past where it started: the lead is on the left now, the right's lead (not held) aborted once no pass has
    // needed it for `Stream.CancelGraceMs` (VT-02).
    show(t, CHR2, AT - 5_000, AT + 5_000);
    t.scheduler.advance(t.ctx.config.get(SeqViewerConfig.Stream.CancelGraceMs));
    t.manager.pass();
    const all = data(t, "prefetch");
    expect(all.length).toBeGreaterThan(right.length);
    const left = all.filter((at) => !right.includes(at));
    expect(left.every((at) => at < AT)).toBe(true);
    expect(t.streams.phylop.aborted.filter((r) => r.priority === "prefetch").length).toBeGreaterThan(0);
    t.ctx.dispose();
  });

  it("a prefetch arrives into the cache: panning on finds the region held, nothing more asked for it", async () => {
    const t = await opened();
    show(t, CHR2, AT, AT + 10_000);
    show(t, CHR2, AT + 2_000, AT + 12_000);
    for (let k = 0; k < 4; k++) {
      t.streams.phylop.release();
      await flush();
      t.manager.pass();
    }
    const lead = data(t, "prefetch");
    expect(lead.length).toBeGreaterThan(0);
    const asked = t.streams.phylop.requests.length;
    // A pan the margin's width on lands over the prefetched regions: they are held, so none is asked again.
    const v = t.ctx.layout.axis(GENOME);
    const w = (v?.visHi[CHR2] ?? 0) - (v?.visLo[CHR2] ?? 0);
    show(t, CHR2, AT + 2_000 + w, AT + 12_000 + w);
    const again = t.streams.phylop.requests.slice(asked).filter((r) => r.level === -1 && lead.includes(r.start));
    expect(again).toEqual([]);
    t.ctx.dispose();
  });

  it("`Stream.Prefetch` off asks none", async () => {
    const t = await opened(false);
    show(t, CHR2, AT, AT + 10_000);
    show(t, CHR2, AT + 5_000, AT + 15_000);
    expect(t.streams.phylop.requests.some((r) => r.priority === "prefetch")).toBe(false);
    t.ctx.dispose();
  });

  it("a prefetch is not waited for: `settled` resolves with only prefetches outstanding", async () => {
    const t = await opened();
    show(t, CHR2, AT, AT + 10_000);
    show(t, CHR2, AT + 5_000, AT + 15_000);
    // Answer the view's and the margin's; the prefetches stay held.
    const keep = t.streams.phylop.requests.filter((r) => r.priority === "prefetch").length;
    expect(keep).toBeGreaterThan(0);
    const n = t.streams.phylop.requests.length - keep;
    t.streams.phylop.release(n);
    await flush();
    t.manager.pass();
    let done = false;
    void t.manager.settled().then(() => {
      done = true;
    });
    t.manager.pass();
    await flush();
    expect(done).toBe(true);
    t.ctx.dispose();
  });
});

describe("slots", () => {
  it("prefetches wait while the bound is full of the view's and the margin's requests", async () => {
    const t = await opened();
    show(t, CHR2, AT, AT + 10_000);
    t.ctx.config.set(MaxInFlight, t.manager.stats.current.inFlight);
    show(t, CHR2, AT + 5_000, AT + 15_000);
    expect(t.streams.phylop.requests.some((r) => r.priority === "prefetch")).toBe(false);
    // The slots free as the others land: then the prefetches go.
    t.ctx.config.set(MaxInFlight, 100);
    t.manager.pass();
    expect(t.streams.phylop.requests.some((r) => r.priority === "prefetch")).toBe(true);
    t.ctx.dispose();
  });

  it("a view's request frees a slot by aborting the newest prefetch the view does not need; none to take: it waits", async () => {
    const t = await opened();
    show(t, CHR2, AT, AT + 10_000);
    show(t, CHR2, AT + 5_000, AT + 15_000);
    const prefetches = t.streams.phylop.requests.filter((r) => r.priority === "prefetch").length;
    expect(prefetches).toBeGreaterThan(1);
    t.streams.phylop.release(t.streams.phylop.requests.length - prefetches);
    await flush();
    t.manager.pass();
    const take = () => (t.manager as unknown as { preempt(): boolean }).preempt();
    const [before, flying] = [t.streams.phylop.aborted.length, t.manager.stats.current.inFlight];
    expect(flying).toBe(prefetches);
    for (let k = 1; k <= prefetches; k++) {
      expect(take()).toBe(true);
      expect(t.streams.phylop.aborted.length).toBe(before + k);
    }
    expect(take()).toBe(false);
    expect(t.manager.stats.current.inFlight).toBe(0);
    t.ctx.dispose();
  });
});

describe("the first pan after load", () => {
  it("is followed too, at the default in-flight bound: the lead is asked once the view's and the margin's have landed", async () => {
    const t = await openGenome({ gated: true, chromosomes: THREE });
    t.ctx.config.set(Prefetch, true);
    // Opened as the story is: the view settles on a place and its regions land, a bound of 6 letting a few through a pass.
    const settle = async () => {
      for (let k = 0; k < 60; k++) {
        t.streams.phylop.release();
        await flush();
        t.manager.pass();
      }
    };
    show(t, CHR2, AT, AT + 100_000);
    await settle();
    expect(data(t, "prefetch")).toEqual([]);
    // The first pan: right by a third of the view.
    const w = (t.ctx.layout.axis(GENOME)?.visHi[CHR2] ?? 0) - (t.ctx.layout.axis(GENOME)?.visLo[CHR2] ?? 0);
    show(t, CHR2, AT + w / 3, AT + w / 3 + 100_000);
    await settle();
    const lead = t.streams.phylop.requests.filter((r) => r.priority === "prefetch");
    expect(lead.length).toBeGreaterThan(0);
    t.ctx.dispose();
  });

  it("is followed while the view's regions are still in flight: the lead waits its turn, then goes", async () => {
    const t = await openGenome({ gated: true, chromosomes: THREE });
    show(t, CHR2, AT, AT + 100_000);
    // Nothing has landed when the first pan comes.
    show(t, CHR2, AT + 30_000, AT + 130_000);
    for (let k = 0; k < 60; k++) {
      t.streams.phylop.release();
      await flush();
      t.manager.pass();
    }
    expect(
      data(t, "prefetch").length + t.streams.phylop.requests.filter((r) => r.priority === "prefetch").length,
    ).toBeGreaterThan(0);
    t.ctx.dispose();
  });
});

describe("the first pan with every layer streaming", () => {
  it("genes, variants, letters and a signal at the default bound: the pan comes before the view has landed, the lead follows", async () => {
    const s = hg38Streams({ chromosomes: THREE, letters: true });
    const loaders = { load: () => Promise.reject(new Error("none")), stream: s.stream };
    const t = await composed(streamedDoc({ items: true, sequence: true, values: true }), {
      width: 1000,
      height: 300,
      dpr: 1,
      loaders,
    });
    const frames = async (n: number) => {
      for (let k = 0; k < n; k++) {
        t.scheduler.advance(16);
        await flush();
      }
    };
    t.vp?.zoomTo(1, AT, AT + 100_000);
    await frames(2);
    t.vp?.panBy(300, 0, "api");
    await frames(80);
    const lead = [s.genes, s.variants, s.phylop, s.genome]
      .flatMap((p) => p.requests)
      .filter((r) => r.priority === "prefetch");
    expect(lead.length).toBeGreaterThan(0);
    t.ctx.dispose();
  });
});

describe("a pan before any pass has seen the view", () => {
  it("is followed too: the view's motion is read from its changes, not only from the passes (a hidden page runs no frame between)", async () => {
    const s = hg38Streams({ chromosomes: THREE, letters: true });
    const loaders = { load: () => Promise.reject(new Error("none")), stream: s.stream };
    const t = await composed(streamedDoc({ values: true }), { width: 1000, height: 300, dpr: 1, loaders });
    const frames = async (n: number) => {
      for (let k = 0; k < n; k++) {
        t.scheduler.advance(16);
        await flush();
      }
    };
    const lead = () => s.phylop.requests.filter((r) => r.priority === "prefetch").length;
    // The view and the first pan come before the manager's pass (a frame not yet run: a hidden page's).
    t.vp?.zoomTo(1, AT, AT + 100_000);
    t.vp?.panBy(300, 0, "api");
    await frames(60);
    const first = lead();
    expect(first).toBeGreaterThan(0);
    // And it is toward the side panned: every region of it past the view on the right.
    const hi = t.ctx.layout.axis(GENOME)?.visHi[1] ?? 0;
    expect(s.phylop.requests.filter((r) => r.priority === "prefetch").every((r) => r.stop > AT)).toBe(true);
    expect(hi).toBeGreaterThan(AT);
    t.ctx.dispose();
  });
});

describe("the zoom's next level", () => {
  it("zooming in asks the next finer level over the view; zooming out the next coarser", async () => {
    const t = await opened();
    show(t, CHR2, AT, AT + 200_000);
    show(t, CHR2, AT + 50_000, AT + 100_000);
    const levels = (p: string) =>
      new Set(t.streams.phylop.requests.filter((r) => r.priority === p).map((r) => r.level));
    const mid = t.streams.phylop.requests
      .filter((r) => r.priority === "view" && r.level >= 0 && r.level < 5)
      .map((r) => r.level);
    const finer = [...levels("prefetch")].filter((l) => l < Math.min(...mid));
    expect(finer.length).toBeGreaterThan(0);
    t.ctx.dispose();
  });
});

describe("the snapshots' focus regions", () => {
  /** A genome of three chromosomes with a focus on chr2 and three snapshots moving it (10 kb windows far apart). */
  async function snapshots() {
    const b = createSeqViewSpec({ title: "hg38" });
    const genome = b.data.source({ ref: "genome", stream: "genome", data_kind: "letters" });
    const hg38 = b.data.complex({ ref: "hg38", source_ref: genome, label: "GRCh38" });
    b.data.source({ ref: "phylop_src", stream: "phylop", data_kind: "values", coordinate_system: "hg38" });
    const phylop = b.data.annotation({ ref: "phylop", geometry: "values", source_ref: "phylop_src" });
    const v = b.view({ ref: "browser", label: "hg38" });
    v.block({ complex: hg38 });
    v.focus({ ref: "browser.focus", target: "hg38/chr2", start: AT, stop: AT + 10_000 });
    v.track({ label: "Ruler" }).ruler();
    v.track({ label: "phyloP" }).values({ data: phylop });
    for (const [ref, at] of [
      ["a", AT],
      ["b", 150_000_000],
      ["c", 200_000_000],
    ] as const)
      b.snapshot({
        ref,
        duration_ms: 1000,
        transition_ms: 0,
        ops: [{ op: "set", ref: "browser.focus", params: { start: at, stop: at + 10_000 } }],
      });
    return b.getDocument();
  }

  it("the hints are the focus regions of the snapshots before and after the one shown", async () => {
    const doc = await snapshots();
    const at = (current: string | null) => snapshotHints(doc, current).map((h) => [h.target, h.start] as const);
    expect(at(null)).toEqual([["hg38/chr2", AT]]);
    // Before the first is the base's focus; after the last, none.
    expect(at("a")).toEqual([
      ["hg38/chr2", AT],
      ["hg38/chr2", 150_000_000],
    ]);
    expect(at("b")).toEqual([
      ["hg38/chr2", AT],
      ["hg38/chr2", 200_000_000],
    ]);
    expect(at("c")).toEqual([["hg38/chr2", 150_000_000]]);
  });

  it("the first pan after the document opens (a focus sets the view) is followed by a lead, driven by frames", async () => {
    const s = hg38Streams({ chromosomes: THREE });
    const loaders = { load: () => Promise.reject(new Error("none")), stream: s.stream };
    const t = await composed(await snapshots(), { width: 1000, height: 200, dpr: 1, loaders });
    const frames = async (n: number) => {
      for (let k = 0; k < n; k++) {
        t.scheduler.advance(16);
        await flush();
      }
    };
    await frames(40);
    const before = s.phylop.requests.filter((r) => r.priority === "prefetch").length;
    expect(before).toBe(0);
    t.vp?.panBy(300, 0, "api");
    await frames(60);
    expect(s.phylop.requests.filter((r) => r.priority === "prefetch").length).toBeGreaterThan(0);
    t.ctx.dispose();
  });

  it("entering the next snapshot draws from the cache: its regions were prefetched, none is pending", async () => {
    const s = hg38Streams({ chromosomes: THREE, gated: true });
    const loaders = { load: () => Promise.reject(new Error("none")), stream: s.stream };
    const t = await composed(await snapshots(), { width: 1000, height: 200, dpr: 1, loaders });
    t.ctx.config.set(MaxInFlight, 100);
    const m = t.ctx.streams as import("../src/streams.js").Streams;
    t.layout(GENOME);
    m.pass();
    const answer = async () => {
      for (let k = 0; k < 4; k++) {
        s.phylop.release();
        await flush();
        m.pass();
      }
    };
    await answer();
    await t.ctx.snapshots.apply("a");
    t.layout(GENOME);
    m.pass();
    await answer();
    // The view's own regions asked first; the next snapshot's (b: 150 Mb) as prefetches.
    const ahead = s.phylop.requests.filter(
      (r) => r.priority === "prefetch" && r.level === -1 && r.start >= 150_000_000,
    );
    expect(ahead.length).toBeGreaterThan(0);
    expect(s.phylop.requests.every((r) => r.priority !== "prefetch" || r.start < 160_000_000 || r.level >= 0)).toBe(
      true,
    );
    const asked = s.phylop.requests.length;
    await t.ctx.snapshots.apply("a");
    await t.ctx.snapshots.apply("b");
    t.layout(GENOME);
    m.pass();
    // Entering b asks for no region of its focus: they were held.
    const focus = s.phylop.requests
      .slice(asked)
      .filter((r) => r.level === -1 && r.start >= 150_000_000 && r.start < 150_010_000 && r.priority === "view");
    expect(focus).toEqual([]);
    t.ctx.dispose();
  });
});
