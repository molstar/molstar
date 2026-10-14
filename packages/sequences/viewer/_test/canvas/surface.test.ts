/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The surface's frame scheduling and its resize and DPR path, in Chromium (seqspec:viewer#the-render-backend; CS2):
// events between two frames ask for one frame and their flags are drawn together; a paused context asks for none; a
// section that leaves the view and comes back is drawn on its new viewport and layout; the canvases follow the layout
// at every DPR; a DPR change with no size change is followed. Then, through the real surface: a scripted 30-frame
// resize, as a rail drag makes it, where what the page paints in every frame is a whole frame at the host's size, never
// blank and never a frame behind.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { afterEach, describe, expect, it } from "vitest";
import { Canvas2DBackend, Canvas2DShared, type DprHost, SectionSurface, watchDpr } from "../../src/canvas/index.ts";
import { type DrawResources, SkinToken } from "../../src/draw/backend.ts";
import { Dirty, SectionComposer } from "../../src/draw/index.ts";
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  manualScheduler,
  rafScheduler,
  type Scheduler,
  type SeqViewerContext,
  type Skin,
} from "../../src/index.ts";
import { fixtureMeasurer } from "../../src/testing/index.ts";
import { glyphTableOf } from "../../src/themes.ts";
import { scaleDoc } from "../prepare-fixtures.ts";

const open: { ctx: SeqViewerContext; host: HTMLElement; surface: SectionSurface }[] = [];
afterEach(() => {
  for (const t of open.splice(0)) {
    t.surface.dispose();
    t.ctx.dispose();
    t.host.remove();
  }
});

const pageScheduler = (): Scheduler =>
  rafScheduler({
    requestAnimationFrame: (cb) => window.requestAnimationFrame(cb),
    cancelAnimationFrame: (id) => window.cancelAnimationFrame(id),
    setTimeout: (cb, ms) => window.setTimeout(cb, ms),
    clearTimeout: (id) => window.clearTimeout(id),
    now: () => performance.now(),
  });
const nextFrame = () => new Promise<number>((ok) => requestAnimationFrame(ok));

async function mount(doc: unknown, w: number, h: number, o: { dpr?: number; page?: boolean; section?: number } = {}) {
  const scheduler = o.page ? pageScheduler() : manualScheduler();
  const ctx = createSeqViewer(DefaultSeqViewerSpec(), { scheduler, text: fixtureMeasurer() });
  const r = await ctx.setDocument(doc as never);
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  await ctx.graph.update();
  const host = document.createElement("div");
  host.style.cssText = `position:absolute;left:0;top:0;width:${w}px;height:${h}px`;
  document.body.appendChild(host);
  const section = ctx.viewports.sections.value[o.section ?? 0];
  const surface = new SectionSurface(ctx, section, host, { dpr: o.dpr, readback: true });
  open.push({ ctx, host, surface });
  return { ctx, host, surface, section, scheduler };
}
const canvases = (s: SectionSurface) => Array.from(s.element.querySelectorAll("canvas"));
const pixels = (c: HTMLCanvasElement) =>
  (c.getContext("2d") as CanvasRenderingContext2D).getImageData(0, 0, c.width, c.height).data;

function resources(ctx: SeqViewerContext): DrawResources {
  const glyphs = glyphTableOf(ctx.registry);
  return {
    font: (id) => ctx.text.fontString(id),
    fontPx: (id) => Number.parseFloat(ctx.text.fontString(id).split(" ")[1]),
    token: (t) => {
      const name = Object.entries(SkinToken).find(([, id]) => id === t)?.[0] as keyof Skin;
      return ctx.skin.value[name] as string;
    },
    glyph: (id) => glyphs.glyph(id),
  };
}
/** A fresh composer's whole tracks frame at the context's current view, through a back target, as pixels. */
function reference(ctx: SeqViewerContext, section: string, W: number, H: number) {
  const c = new OffscreenCanvas(W, H);
  const g = c.getContext("2d", { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
  const composer = new SectionComposer(ctx, section);
  composer.frame(new Canvas2DBackend([g], new Canvas2DShared(resources(ctx))), Dirty.tracks);
  composer.dispose();
  return g.getImageData(0, 0, W, H).data;
}
function differing(a: Uint8ClampedArray, b: Uint8ClampedArray): number {
  if (a.length !== b.length) return -1;
  let n = 0;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) n++;
  return n;
}

/** Genes and variants over two blocks: the tracks once seen to blink. */
const genome = () => scaleDoc({ length: 40_000, ranges: 600, variants: 800 });

describe("accessibility of a surface mounted on its own (VB-09)", () => {
  it("is a named group with the role description, its name the section's label, kept up to date", async () => {
    const { ctx, surface, section } = await mount(genome(), 600, 200);
    const el = surface.element;
    expect([el.tabIndex, el.getAttribute("role"), el.getAttribute("aria-roledescription")]).toEqual([
      0,
      "group",
      "sequence view",
    ]);
    const label = el.getAttribute("aria-label");
    expect(label).toBeTruthy();
    await ctx.applyPatch({
      format: "seq_view_spec_patch",
      ops: [{ op: "set", ref: section, params: { label: "Renamed" } }],
    });
    expect(el.getAttribute("aria-label")).toBe("Renamed");
  });
});

describe("frame scheduling", () => {
  it("events between two frames ask for one frame, drawn with all their flags; none while paused", async () => {
    const t = await mount(genome(), 800, 300, { dpr: 1 });
    const sch = t.scheduler as ReturnType<typeof manualScheduler>;
    t.ctx.viewports.get(t.section)?.setView(0.5, 0, "api");
    sch.advance(0);
    const frames = t.surface.stats.frames;
    expect(sch.pendingFrames).toBe(0);
    t.ctx.highlight.set("probe", { systems: ["seqstar:chr0"], runs: Int32Array.of(0, 100, 200) });
    t.ctx.selection.set({ kind: "seq", seq: { chr0: [[300, 310]] } } as never, "api");
    t.ctx.viewports.get(t.section)?.panBy(30, 0, "api");
    t.ctx.ui.setKeyCursor(null);
    expect(sch.pendingFrames).toBe(1);
    const tracks = t.surface.composer.stats.copies + t.surface.composer.stats.fullDraws;
    sch.advance(0);
    expect(t.surface.stats.frames).toBe(frames + 1);
    expect(t.surface.composer.stats.copies + t.surface.composer.stats.fullDraws).toBe(tracks + 1);
    t.ctx.setActive(false);
    t.ctx.viewports.get(t.section)?.panBy(30, 0, "api");
    sch.advance(0);
    expect([sch.pendingFrames, t.surface.stats.frames]).toEqual([0, frames + 1]);
    t.ctx.setActive(true);
    expect(sch.pendingFrames).toBe(1);
    sch.advance(0);
    expect(t.surface.stats.frames).toBe(frames + 2);
  });

  it("a section that leaves the document and comes back is drawn on its new viewport and layout", async () => {
    const doc = (second: boolean) => {
      const b = createSeqViewSpec();
      const a = b.data.sequence({ ref: "a", value: "MKVLATGHEKWRSPQIDNCY".repeat(20) });
      const v = b.view({ ref: "v", label: "V" });
      v.section({ ref: "s1" }).block({ sequence: a }).track({ label: "L" }).letters();
      if (second) v.section({ ref: "s2" }).block({ sequence: a }).track({ label: "L" }).letters();
      return b.getDocument();
    };
    const t = await mount(doc(true), 700, 200, { dpr: 1, section: 1 });
    const sch = t.scheduler as ReturnType<typeof manualScheduler>;
    sch.advance(0);
    const before = [t.ctx.viewports.get("s2"), t.ctx.layout.get("s2")];
    await t.ctx.setDocument(doc(false));
    sch.advance(0);
    await t.ctx.setDocument(doc(true));
    await t.ctx.graph.update();
    sch.advance(0);
    const vp = t.ctx.viewports.get("s2");
    expect([vp, t.ctx.layout.get("s2")].map((x, k) => x !== before[k] && x !== null)).toEqual([true, true]);
    expect(vp?.live.current.width).toBe(700);
    const [main] = canvases(t.surface);
    expect(differing(pixels(main), reference(t.ctx, "s2", main.width, main.height))).toBe(0);
    // Bound to the new viewport: a pan draws.
    vp?.setView(14, 0, "api");
    sch.advance(0);
    const drawn = t.surface.composer.stats.copies + t.surface.composer.stats.fullDraws;
    vp?.panBy(40, 0, "api");
    sch.advance(0);
    expect(t.surface.composer.stats.copies + t.surface.composer.stats.fullDraws).toBe(drawn + 1);
  });
});

describe("size and DPR", () => {
  it("each canvas is the layout's size times the DPR, shown at the layout's CSS size", async () => {
    for (const dpr of [1, 1.5, 2]) {
      const t = await mount(genome(), 701, 263, { dpr });
      const l = t.ctx.layout.get(t.section)?.current;
      if (!l) throw new Error("no layout");
      const [main, over, nav] = canvases(t.surface);
      for (const c of [main, over]) expect([c.width, c.height]).toEqual([Math.ceil(701 * dpr), Math.ceil(263 * dpr)]);
      expect(nav.height).toBe(Math.ceil(l.navHeight * dpr));
      expect(Number.parseFloat(main.style.width) * dpr).toBeCloseTo(main.width, 2);
    }
  });

  it("watchDpr follows changes with no size change, its query re-armed at each", () => {
    const queries: { q: string; on: (() => void) | null }[] = [];
    const win = {
      devicePixelRatio: 1,
      matchMedia(q: string) {
        const m = { q, on: null as (() => void) | null };
        queries.push(m);
        return {
          addEventListener: (_: string, f: () => void) => {
            m.on = f;
          },
          removeEventListener: (_: string, f: () => void) => {
            if (m.on === f) m.on = null;
          },
        } as unknown as MediaQueryList;
      },
    };
    const seen: number[] = [];
    const stop = watchDpr(win as DprHost, (d) => seen.push(d));
    win.devicePixelRatio = 2;
    queries[0].on?.();
    win.devicePixelRatio = 1.25;
    queries[1].on?.();
    expect(seen).toEqual([2, 1.25]);
    expect(queries.map((m) => [m.q, m.on !== null])).toEqual([
      ["(resolution: 1dppx)", false],
      ["(resolution: 2dppx)", false],
      ["(resolution: 1.25dppx)", true],
    ]);
    stop();
    expect(queries[2].on).toBeNull();
  });

  it("a surface redraws at a new DPR with no size change", async () => {
    const real = Object.getOwnPropertyDescriptor(window, "devicePixelRatio");
    const realMatch = window.matchMedia;
    let fire: (() => void) | null = null;
    let dpr = 1;
    Object.defineProperty(window, "devicePixelRatio", { configurable: true, get: () => dpr });
    window.matchMedia = (() =>
      ({
        addEventListener: (_: string, f: () => void) => {
          fire = f;
        },
        removeEventListener: () => {},
      }) as unknown as MediaQueryList) as typeof window.matchMedia;
    try {
      const t = await mount(genome(), 600, 240);
      expect(canvases(t.surface)[0].width).toBe(600);
      dpr = 2;
      (fire as unknown as () => void)();
      const [main] = canvases(t.surface);
      expect([main.width, main.height, t.ctx.viewports.get(t.section)?.live.current.dpr]).toEqual([1200, 480, 2]);
      expect(differing(pixels(main), reference(t.ctx, t.section, 1200, 480))).toBe(0);
    } finally {
      window.matchMedia = realMatch;
      if (real) Object.defineProperty(window, "devicePixelRatio", real);
    }
  });

  for (const dpr of [1, 2])
    it(`a 30-frame resize paints a whole frame at the host's size in every frame (DPR ${dpr})`, async () => {
      const t = await mount(genome(), 900, 300, { dpr, page: true });
      const vp = t.ctx.viewports.get(t.section);
      vp?.setView(0.02, 0, "api");
      await t.ctx.graph.update();
      await nextFrame();
      await nextFrame();
      const widths = Array.from({ length: 30 }, (_, k) => (k < 15 ? 900 - 23 * k : 578 + 31 * (k - 15)));
      let checked = 0;
      for (let k = 0; k <= widths.length; k++) {
        // What the previous frame painted: its canvases sized to that frame's host, the image whole.
        await nextFrame();
        const w = k === 0 ? 900 : widths[k - 1];
        const [main, over] = canvases(t.surface);
        expect([main.width, over.width], `frame ${k}`).toEqual([Math.ceil(w * dpr), Math.ceil(w * dpr)]);
        const d = differing(pixels(main), reference(t.ctx, t.section, main.width, main.height));
        expect(d, `frame ${k}, ${w} px`).toBe(0);
        checked++;
        // As a rail drag does: the host's size changes, the ResizeObserver calls back before this frame's paint.
        if (k < widths.length) t.host.style.width = `${widths[k]}px`;
      }
      expect(checked).toBe(31);
      expect(t.surface.stats.resizes).toBeGreaterThanOrEqual(29);
    });
});
