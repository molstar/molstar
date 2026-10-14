/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The recording backend (seqspec:viewer#the-render-backend): every call in an op log, `present` and `stale` included;
// every batch's items in its stores, so tests read colors, counts and positions; images and layers kept; a replay
// reproduces the frame; its hash sees a change to any one field, and a composer's frames hash alike on fresh contexts.

import { describe, expect, it } from "vitest";
import {
  BackendCaps,
  GlyphBatch,
  Pass,
  PathBatch,
  PathMode,
  PathVerb,
  RectBatch,
  SectorBatch,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
} from "../../src/draw/api.js";
import { Dirty } from "../../src/draw/index.js";
import { hashRecording, Op, RecordingBackend } from "../../src/testing/index.js";
import { COPY_BACK, composed, recording } from "../composer-fixtures.js";
import { dotPlotDoc } from "../dot-plot-fixtures.js";
import { family } from "../overlay-docs.js";
import { allTypes, linesDoc } from "../prepare-fixtures.js";
import { drawFrame, Frame, PIXELS } from "../support/frame.js";

const codes = (r: RecordingBackend) => Array.from({ length: r.opCount }, (_, k) => r.ops[k * 4]);
const args = (r: RecordingBackend, k: number, n: number) => [
  ...r.args.subarray(r.ops[k * 4 + 3], r.ops[k * 4 + 3] + n),
];

describe("recording backend", () => {
  it("draws whole frames: no copy, no back target", () => {
    expect(new RecordingBackend().caps & (BackendCaps.Copy | BackendCaps.Back)).toBe(0);
  });

  it("logs the composer's calls with their arguments, presentation included", () => {
    const r = new RecordingBackend();
    r.beginPass(Pass.overlay, 300, 100, 2);
    r.clear(0, 0, 300, 100);
    r.pushClip(10, 0, 50, 100);
    r.popClip();
    r.copy(-4, 0, 0, 0, 300, 100);
    r.stale(0, 0, 100, 20, 0, 0, 200, 20);
    r.present(0, 20, 300, 30);
    r.endPass();
    expect(codes(r)).toEqual([
      Op.beginPass,
      Op.clear,
      Op.pushClip,
      Op.popClip,
      Op.copy,
      Op.stale,
      Op.present,
      Op.endPass,
    ]);
    expect(args(r, 0, 4)).toEqual([Pass.overlay, 300, 100, 2]);
    expect(args(r, 4, 6)).toEqual([-4, 0, 0, 0, 300, 100]);
    expect(args(r, 5, 8)).toEqual([0, 0, 100, 20, 0, 0, 200, 20]);
    expect(args(r, 6, 4)).toEqual([0, 20, 300, 30]);
    expect(r.count(Op.present)).toBe(1);
  });

  it("copies each batch's items, so the renderer may reuse it at once", () => {
    const r = new RecordingBackend();
    const b = new RectBatch();
    b.begin(0xff0000, 128);
    b.push(1, 2, 3, 4);
    r.rects(b);
    b.reset();
    b.begin(SkinToken.selectionBorder, 255, 1);
    b.push(5, 6, 7, 8);
    b.push(9, 10, 11, 12);
    r.rects(b);
    const t = new TextBatch();
    t.style(4, TextAlign.center, TextBaseline.middle, 0.5);
    t.begin(0x00ff00);
    t.push(1.5, 2.5, "A");
    r.text(t);
    expect(r.count(Op.rects)).toBe(2);
    expect([r.ops[4 + 1], r.ops[4 + 2]]).toEqual([1, 2]);
    const s = r.rectItems;
    expect(s.n).toBe(3);
    expect([...s.x.subarray(0, 3)]).toEqual([1, 5, 9]);
    expect([...s.color.subarray(0, 3)]).toEqual([0xff0000, SkinToken.selectionBorder, SkinToken.selectionBorder]);
    expect([s.alpha[0], s.stroke[0], s.stroke[2]]).toEqual([128, 0, 1]);
    expect(args(r, 2, 4)).toEqual([4, TextAlign.center, TextBaseline.middle, 0.5]);
    expect([r.textItems.text[0], r.textItems.x[0]]).toEqual(["A", 1.5]);
    r.reset();
    expect([r.opCount, r.rectItems.n, r.textItems.n]).toEqual([0, 0, 0]);
  });

  it("keeps images' pixels and layers' recordings", () => {
    const r = new RecordingBackend();
    const px = r.images.fromPixels(2, 2, PIXELS);
    const l = r.images.beginLayer(10, 4);
    const b = new RectBatch();
    b.begin(0x0000ff);
    b.push(0, 0, 10, 4);
    l.rects(b);
    const layer = r.images.endLayer();
    expect([px, layer]).toEqual([0, 1]);
    expect(r.images.content[0]).toEqual(PIXELS);
    expect(r.images.content[0]).not.toBe(PIXELS);
    expect((r.images.content[1] as RecordingBackend).rectItems.n).toBe(1);
    expect([r.images.w[1], r.images.h[1]]).toEqual([10, 4]);
    r.images.dispose(0);
    expect(r.images.content[0]).toBeNull();
  });

  it("replays a frame into another backend as it was recorded", () => {
    const a = new RecordingBackend();
    const f = new Frame(1000, 300);
    const img = a.images.fromPixels(2, 2, PIXELS);
    drawFrame(a, f, 3, img);
    const b = new RecordingBackend();
    a.replay(b);
    expect(codes(b)).toEqual(codes(a));
    expect([...b.args.subarray(0, b.argCount)]).toEqual([...a.args.subarray(0, a.argCount)]);
    for (const k of ["rectItems", "sectorItems", "pathItems", "glyphItems", "textItems"] as const) {
      expect(b[k].n).toBe(a[k].n);
      expect(b[k].color.subarray(0, b[k].n)).toEqual(a[k].color.subarray(0, a[k].n));
    }
    expect(b.rectItems.x.subarray(0, 1000)).toEqual(a.rectItems.x.subarray(0, 1000));
    expect(b.textItems.text.slice(0, 300)).toEqual(a.textItems.text.slice(0, 300));
    expect(b.pathItems.coords.subarray(0, b.pathItems.nc)).toEqual(a.pathItems.coords.subarray(0, a.pathItems.nc));
    // The image was made in the target for the replay, and kept there: a recording keeps what it is given.
    expect([b.images.w[0], b.images.h[0], b.images.content[0]]).toEqual([2, 2, PIXELS]);
    expect(hashRecording(b)).toBe(hashRecording(a));
  });

  it("hashes everything it holds: a change to any one field of any item, argument or image changes the hash", () => {
    const base = hashRecording(synthetic(SYNTH));
    // Deterministic: the same frame drawn again, and a replayed copy, hash alike.
    expect(hashRecording(synthetic(SYNTH))).toBe(base);
    const copy = new RecordingBackend();
    synthetic(SYNTH).replay(copy);
    expect(hashRecording(copy)).toBe(base);
    const same = CHANGES.filter(([, change]) => hashRecording(synthetic({ ...SYNTH, ...change })) === base);
    expect(same.map(([name]) => name)).toEqual([]);
  });

  it("a composer's frames hash alike on fresh contexts, and differ between documents", async () => {
    const docs = [allTypes(), linesDoc(400, { style: "circular" }), dotPlotDoc(), family()];
    const run = async (k: number) => {
      const t = await composed(docs[k], { width: 900 });
      const be = recording(COPY_BACK);
      const first = hashRecording(t.frame(be, Dirty.all));
      t.vp?.setView(12, 40, "api");
      t.frame(be, Dirty.all);
      t.vp?.panBy(30, 0, "api");
      return [first, hashRecording(t.frame(be, Dirty.all))];
    };
    const a = [];
    for (let k = 0; k < docs.length; k++) a.push(await run(k));
    for (let k = 0; k < docs.length; k++) expect(await run(k)).toEqual(a[k]);
    expect(new Set(a.map((x) => x[0])).size).toBe(docs.length);
  });
});

/** One item of every kind, each field from `s`; an image of pixels and a layer of one rect, disposed or not. */
const SYNTH = {
  pass: Pass.main,
  passW: 300,
  clipX: 10,
  copyDx: -4,
  staleSx: 3,
  presentY: 20,
  rectX: 1,
  rectY: 2,
  rectW: 3,
  rectH: 4,
  rectColor: 0xff0000,
  rectAlpha: 128,
  rectStroke: 1,
  sectorCx: 5.5,
  sectorCy: 6.5,
  sectorR0: 1.25,
  sectorR1: 9.75,
  sectorA0: 0.1,
  sectorA1: 1.2,
  sectorColor: 0x00ff00,
  sectorAlpha: 200,
  pathVerb: PathVerb.line as number,
  pathX: 3.5,
  pathY: 4.5,
  pathColor: 0x0000ff,
  pathAlpha: 77,
  pathMode: PathMode.stroke as number,
  pathWidth: 3,
  glyphX: 7.5,
  glyphY: 8.5,
  glyphRadius: 2.5,
  glyphId: 2,
  glyphColor: 0x222222,
  glyphAlpha: 90,
  glyphFont: 3,
  textX: 11.5,
  textY: 12.5,
  textSx: 1,
  textSy: 1,
  text: "Ab",
  textColor: 0x123456,
  textAlpha: 250,
  textFont: 4,
  textAlign: TextAlign.center as number,
  textBaseline: TextBaseline.middle as number,
  textAngle: 0.5,
  pixel: 7,
  imageW: 2,
  imageX: 0.5,
  smooth: 0,
  layerColor: 0x0000ff,
  disposed: 0,
};
/** Each change of one field: integers by one, doubles by 1e-9, the string by a letter. */
const CHANGES: [string, Partial<typeof SYNTH>][] = Object.entries(SYNTH).map(([k, v]) => [
  k,
  {
    [k]:
      typeof v === "string"
        ? `${v}~`
        : k === "pathVerb"
          ? PathVerb.move
          : k === "pathMode"
            ? PathMode.fill
            : k === "disposed" || k === "smooth" || k === "rectStroke"
              ? 1 - v
              : Number.isInteger(v)
                ? v + 1
                : v + 1e-9,
  },
]);

function synthetic(s: typeof SYNTH): RecordingBackend {
  const r = new RecordingBackend();
  r.beginPass(s.pass, s.passW, 100, 2);
  r.clear(0, 0, s.passW, 100);
  r.pushClip(s.clipX, 0, 50, 100);
  r.copy(s.copyDx, 0, 0, 0, 300, 100);
  r.stale(s.staleSx, 0, 100, 20, 0, 0, 200, 20);
  const rb = new RectBatch();
  rb.begin(s.rectColor, s.rectAlpha, s.rectStroke);
  rb.push(s.rectX, s.rectY, s.rectW, s.rectH);
  r.rects(rb);
  const sb = new SectorBatch();
  sb.begin(s.sectorColor, s.sectorAlpha);
  sb.push(s.sectorCx, s.sectorCy, s.sectorR0, s.sectorR1, s.sectorA0, s.sectorA1);
  r.sectors(sb);
  const pb = new PathBatch();
  pb.begin(s.pathColor, s.pathAlpha, s.pathMode, s.pathWidth);
  pb.moveTo(1, 1);
  const o = pb.verb(s.pathVerb);
  pb.coords[o] = s.pathX;
  pb.coords[o + 1] = s.pathY;
  r.paths(pb);
  const gb = new GlyphBatch();
  gb.font = s.glyphFont;
  gb.begin(s.glyphColor, s.glyphAlpha);
  gb.push(s.glyphX, s.glyphY, s.glyphRadius, s.glyphId);
  r.glyphs(gb);
  const tb = new TextBatch();
  tb.style(s.textFont, s.textAlign, s.textBaseline, s.textAngle);
  tb.begin(s.textColor, s.textAlpha);
  const i = tb.add(s.text);
  tb.x[i] = s.textX;
  tb.y[i] = s.textY;
  tb.sx[i] = s.textSx;
  tb.sy[i] = s.textSy;
  r.text(tb);
  const px = PIXELS.slice(0, 4);
  px[1] = s.pixel;
  const a = r.images.fromPixels(s.imageW, 2, px);
  const layer = r.images.beginLayer(10, 4);
  const lb = new RectBatch();
  lb.begin(s.layerColor);
  lb.push(0, 0, 10, 4);
  layer.rects(lb);
  const b = r.images.endLayer();
  const at = new Float64Array([s.imageX, 0, 2, 2, 0, 0, 2, 2]);
  r.image(a, at, s.smooth === 1);
  r.image(b, at, false);
  if (s.disposed === 1) r.images.dispose(b);
  r.popClip();
  r.present(0, s.presentY, 300, 30);
  r.endPass();
  return r;
}
