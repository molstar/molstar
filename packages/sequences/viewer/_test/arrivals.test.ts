/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Per-region redraws on the recording backend's batches, in Node: a frame redrawn because regions arrived equals a full
// redraw, batch for batch. Each device px keeps the sequence of batch items drawn over it since its last clear (an
// opaque item restarts it; a translucent one, a letter or a path is appended), over the frames drawn so far, as a
// canvas keeps their pixels; an image's item at a px is its nearest sample. Fine tiles held back over three columns
// (arrival-fixtures.ts) draw the coarse level in their place alone, and their arrivals one column a frame leave nothing
// of it under the others (reduced tiles, translucent by their share); a jump to rows whose labels and tiles have not
// arrived redraws, pass by pass, what changed. The pixels in Chromium: canvas/band-copy.test.ts.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { tilesStream } from "@molstar/seq-core/stream/tiles";
import { tileServer } from "@molstar/seq-core/stream/tiles/server";
import { describe, expect, it } from "vitest";
import {
  BackendCaps,
  Dirty,
  type GlyphBatch,
  NO_COLOR,
  Pass,
  type PathBatch,
  type RectBatch,
  type RenderBackend,
  SectionComposer,
  type SectorBatch,
  type TextBatch,
} from "../src/draw/index.js";
import { createSeqViewer, DefaultSeqViewerSpec, manualScheduler } from "../src/index.js";
import type { Streams } from "../src/streams.js";
import { fixtureMeasurer, RecordedImages } from "../src/testing/index.js";
import { holdBack } from "./arrival-fixtures.js";

const URI = "https://example.org/lazy";
const [W, H] = [700, 300];

const mix = (h: number, k: number) => (Math.imul(h ^ k, 0x9e3779b1) + 0x7f4a7c15) | 0;
function hashOf(s: string, h: number): number {
  for (let i = 0; i < s.length; i++) h = mix(h, s.charCodeAt(i));
  return h;
}

/** The main pass's batch items per device px, as a canvas composites them (see above); copies and clips applied. */
class Stacks implements RenderBackend {
  readonly caps = BackendCaps.Copy | BackendCaps.Back;
  readonly images = new RecordedImages();
  w = 0;
  h = 0;
  px = new Int32Array(0);
  private on = false;
  private snap: Int32Array | null = null;
  private clips: number[][] = [];
  private get clip() {
    return this.clips[this.clips.length - 1];
  }
  beginPass(pass: number, width: number, height: number): void {
    this.on = pass === Pass.main;
    if (!this.on) return;
    if (width !== this.w || height !== this.h)
      [this.w, this.h, this.px] = [width, height, new Int32Array(width * height)];
    this.clips = [[0, 0, width, height]];
    this.snap = null;
  }
  endPass(): void {}
  /** Each px of the rect within the clip. */
  private each(x: number, y: number, w: number, h: number, f: (p: number, x: number, y: number) => void) {
    if (!this.on) return;
    const [cx0, cy0, cx1, cy1] = this.clip;
    for (let j = Math.max(y, cy0); j < Math.min(y + h, cy1); j++)
      for (let i = Math.max(x, cx0); i < Math.min(x + w, cx1); i++) f(j * this.w + i, i, j);
  }
  private put(p: number, key: number, opaque: boolean) {
    this.px[p] = opaque ? key : mix(this.px[p], key);
  }
  clear(x: number, y: number, w: number, h: number): void {
    this.each(x, y, w, h, (p) => {
      this.px[p] = 0;
    });
  }
  pushClip(x: number, y: number, w: number, h: number): void {
    const [cx0, cy0, cx1, cy1] = this.clip;
    this.clips.push([Math.max(cx0, x), Math.max(cy0, y), Math.min(cx1, x + w), Math.min(cy1, y + h)]);
  }
  popClip(): void {
    this.clips.pop();
  }
  copy(dx: number, dy: number, x: number, y: number, w: number, h: number): void {
    if (!this.on) return;
    this.snap ??= this.px.slice();
    const s = this.snap;
    this.each(x, y, w, h, (p, i, j) => {
      const [a, b] = [i - dx, j - dy];
      this.px[p] = a >= 0 && a < this.w && b >= 0 && b < this.h ? s[b * this.w + a] : 0;
    });
  }
  present(): void {}
  stale(): void {}
  rects(b: RectBatch): void {
    for (let k = 0; k < b.n; k++) {
      if (b.color[k] === NO_COLOR) continue;
      const key = mix(b.color[k], b.alpha[k]);
      const op = b.alpha[k] === 255;
      const [x, y, w, h, s] = [b.x[k], b.y[k], b.w[k], b.h[k], b.stroke[k]];
      if (s === 0) this.each(x, y, w, h, (p) => this.put(p, key, op));
      else
        for (const [ex, ey, ew, eh] of [
          [x, y, w, s],
          [x, y + h - s, w, s],
          [x, y + s, s, h - 2 * s],
          [x + w - s, y + s, s, h - 2 * s],
        ])
          this.each(ex, ey, ew, eh, (p) => this.put(p, key, op));
    }
  }
  /** An item drawn at its anchor px (letters, glyphs, sectors, paths: their identity, not their shape). */
  private at(x: number, y: number, key: number) {
    this.each(Math.round(x), Math.round(y), 1, 1, (p) => this.put(p, key, false));
  }
  sectors(b: SectorBatch): void {
    for (let k = 0; k < b.n; k++) this.at(b.cx[k], b.cy[k], mix(mix(b.color[k], b.r1[k] * 64), b.a0[k] * 64));
  }
  paths(b: PathBatch): void {
    for (let k = 0; k < b.n; k++)
      this.at(b.coords[b.cstart[k]], b.coords[b.cstart[k] + 1], mix(mix(b.color[k], b.alpha[k]), b.mode[k]));
  }
  glyphs(b: GlyphBatch): void {
    for (let k = 0; k < b.n; k++) this.at(b.x[k], b.y[k], mix(mix(b.glyph[k], b.color[k]), b.radius[k] * 64));
  }
  text(b: TextBatch): void {
    for (let k = 0; k < b.n; k++)
      this.at(b.x[k], b.y[k], hashOf(b.text[k], mix(mix(b.color[k], b.font), b.align * 8 + b.baseline)));
  }
  image(img: number, r: Float64Array): void {
    const pix = this.images.content[img] as Uint32Array;
    const iw = this.images.w[img];
    const [sx, sy, sw, sh, dx, dy, dw, dh] = r;
    this.each(dx, dy, dw, dh, (p, i, j) => {
      const u = sx + Math.floor(((i + 0.5 - dx) * sw) / dw);
      const v = sy + Math.floor(((j + 0.5 - dy) * sh) / dh);
      const c = pix[v * iw + u];
      const a = c >>> 24;
      if (a > 0) this.put(p, c | 0, a === 255);
    });
  }
}

async function lazy(dpr: number, cell: number, rowScale: number) {
  const server = tileServer({ rowCount: 3000, columnCount: 4000 });
  const ctx = createSeqViewer(DefaultSeqViewerSpec(), {
    scheduler: manualScheduler(),
    text: fixtureMeasurer(),
    loaders: {
      load: () => Promise.reject(new Error("no fixture")),
      stream: tilesStream((u) => (u === URI ? server.transport : null)),
    },
  });
  const b = createSeqViewSpec();
  const aln = b.data.alignment({
    ref: "lazy",
    row_count: 3000,
    column_count: 4000,
    source: { uri: URI, protocol: "seqstar-tiles/1" },
  } as never);
  const v = b.view({ ref: "g", label: "G", height: "box" });
  v.block({ alignment: aln });
  v.track({ label: "Rows" }).rows();
  const done = ctx.setDocument(b.getDocument());
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(W, H, dpr);
  const r = await done;
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  await ctx.graph.update();
  const section = ctx.viewports.sections.value[0];
  const vp = ctx.viewports.get(section);
  if (!vp) throw new Error("no viewport");
  vp.setView(cell, 0, "api");
  vp.setRowScale(rowScale, 0);
  vp.panBy(1500 * cell, 600, "api");
  await ctx.graph.update();
  const arrive = async (each?: (k: number) => void) => {
    const st = ctx.streams as Streams;
    for (let k = 0; k < 8; k++) {
      st.pass();
      await new Promise((ok) => setTimeout(ok, 0));
      each?.(k);
      if (st.stats.current.inFlight === 0 && k > 1) break;
    }
  };
  await arrive();
  const surface = () => {
    const be = new Stacks();
    const composer = new SectionComposer(ctx, section);
    return { be, composer, frame: () => composer.frame(be, Dirty.tracks) };
  };
  return { ctx, section, vp, arrive, surface };
}

/** The drawn frame against a full redraw: device px whose item sequences differ, and the first. */
function compare(a: Stacks, b: Stacks, where: string) {
  let n = 0;
  let first = -1;
  for (let p = 0; p < a.px.length; p++)
    if (a.px[p] !== b.px[p]) {
      n++;
      if (first < 0) first = p;
    }
  expect(n, `${where}, first at (${first % a.w}, ${Math.floor(first / a.w)})`).toBe(0);
}

describe("a frame redrawn for arrivals alone equals a full redraw, batch for batch", () => {
  it.each([
    [1, 0.5, 0.5 / 18],
    [2, 0.3, 0.6 / 18],
  ])(
    "fine tiles held back, then arriving one column a frame, at DPR %i, cell %f, row scale %f",
    async (dpr, cell, rs) => {
      const t = await lazy(dpr, cell, rs);
      const h = holdBack(t.ctx, t.section);
      const lv = Array.from(h.levels);
      expect(lv[2] !== lv[0] || lv[3] !== lv[1], `levels ${lv}`).toBe(true);
      const [drawn, full] = [t.surface(), t.surface()];
      const c = h.center();
      for (const tx of [c - 1, c, c + 1]) h.hold(tx);
      const check = (where: string) => {
        drawn.frame();
        full.composer.redraw();
        full.frame();
        compare(drawn.be, full.be, where);
      };
      check("held");
      for (const tx of [c, c + 1, c - 1]) {
        h.release(tx);
        check(`column ${tx} arrived`);
      }
      // The arrivals were redrawn, not drawn whole.
      expect(drawn.composer.stats.fullDraws).toBe(1);
      t.ctx.dispose();
    },
  );

  it("a jump to rows not arrived: each pass's arrivals redraw what changed, the row labels with their chunk", async () => {
    const t = await lazy(1, 14, 1);
    const [drawn, full] = [t.surface(), t.surface()];
    const check = (where: string) => {
      drawn.frame();
      full.composer.redraw();
      full.frame();
      compare(drawn.be, full.be, where);
    };
    check("arrived");
    for (const row of [1300, 2600]) {
      t.vp.scrollTo(t.vp.live.current.scrollX, row * 18, "api");
      check(`rows from ${row}, pending`);
      await t.arrive((p) => check(`rows from ${row}, pass ${p}`));
    }
    t.ctx.dispose();
  });
});
