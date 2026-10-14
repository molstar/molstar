/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The recording backend (seqspec:viewer#the-render-backend), for tests in Node: every batch's items are copied into
// grow-only stores of the same batch classes, and every call into an op log (code, start, count, args), `present` and
// `stale` included, so a test asserts on colors, counts and positions. It draws whole frames (`caps` none). `replay`
// draws a recorded frame into another backend. `reset` starts a frame; a warm recording allocates nothing.

import {
  type BackendImages,
  GlyphBatch,
  PathBatch,
  RectBatch,
  type RenderBackend,
  SectorBatch,
  TextBatch,
} from "../draw/backend.js";

export const Op = {
  beginPass: 0,
  endPass: 1,
  clear: 2,
  pushClip: 3,
  popClip: 4,
  copy: 5,
  present: 6,
  stale: 7,
  rects: 8,
  sectors: 9,
  paths: 10,
  glyphs: 11,
  text: 12,
  image: 13,
} as const;

export class RecordedImages implements BackendImages {
  /** Per image: its size, and its pixels or the recording of its layer. */
  readonly w: number[] = [];
  readonly h: number[] = [];
  readonly content: (Uint32Array | RecordingBackend | null)[] = [];
  private layer: RecordingBackend | null = null;
  private lw = 0;
  private lh = 0;
  fromPixels(w: number, h: number, rgba: Uint32Array): number {
    return this.keep(w, h, rgba.slice(0, w * h));
  }
  beginLayer(w: number, h: number): RenderBackend {
    this.layer = new RecordingBackend();
    this.lw = w;
    this.lh = h;
    return this.layer;
  }
  endLayer(): number {
    return this.keep(this.lw, this.lh, this.layer);
  }
  dispose(img: number): void {
    this.content[img] = null;
  }
  private keep(w: number, h: number, c: Uint32Array | RecordingBackend | null): number {
    this.w.push(w);
    this.h.push(h);
    return this.content.push(c) - 1;
  }
}

export class RecordingBackend implements RenderBackend {
  readonly caps = 0;
  readonly images = new RecordedImages();
  /** Every item drawn this frame, in order, per kind. */
  readonly rectItems = new RectBatch();
  readonly sectorItems = new SectorBatch();
  readonly pathItems = new PathBatch();
  readonly glyphItems = new GlyphBatch();
  readonly textItems = new TextBatch();
  /** Per op: its code, the first item in its store, its item count, and the offset of its arguments in `args`. */
  ops = new Int32Array(256);
  opCount = 0;
  args = new Float64Array(256);
  argCount = 0;
  private readonly r = new RectBatch();
  private readonly s = new SectorBatch();
  private readonly p = new PathBatch();
  private readonly g = new GlyphBatch();
  private readonly t = new TextBatch();
  private readonly at = new Float64Array(8);

  /** A new frame: the log and the stores emptied, their room kept. */
  reset(): void {
    this.opCount = this.argCount = 0;
    this.rectItems.reset();
    this.sectorItems.reset();
    this.pathItems.reset();
    this.glyphItems.reset();
    this.textItems.reset();
  }
  /** The ops of a kind this frame. */
  count(code: number): number {
    let n = 0;
    for (let k = 0; k < this.opCount; k++) if (this.ops[k * 4] === code) n++;
    return n;
  }

  beginPass(pass: number, width: number, height: number, dpr: number): void {
    const a = this.log(Op.beginPass, 0, 0, 4);
    this.args[a] = pass;
    this.args[a + 1] = width;
    this.args[a + 2] = height;
    this.args[a + 3] = dpr;
  }
  endPass(): void {
    this.log(Op.endPass, 0, 0, 0);
  }
  clear(x: number, y: number, w: number, h: number): void {
    this.rect4(Op.clear, x, y, w, h);
  }
  pushClip(x: number, y: number, w: number, h: number): void {
    this.rect4(Op.pushClip, x, y, w, h);
  }
  popClip(): void {
    this.log(Op.popClip, 0, 0, 0);
  }
  copy(dx: number, dy: number, x: number, y: number, w: number, h: number): void {
    const a = this.log(Op.copy, 0, 0, 6);
    this.args[a] = dx;
    this.args[a + 1] = dy;
    this.args[a + 2] = x;
    this.args[a + 3] = y;
    this.args[a + 4] = w;
    this.args[a + 5] = h;
  }
  present(x: number, y: number, w: number, h: number): void {
    this.rect4(Op.present, x, y, w, h);
  }
  stale(sx: number, sy: number, sw: number, sh: number, dx: number, dy: number, dw: number, dh: number): void {
    const a = this.log(Op.stale, 0, 0, 8);
    this.args[a] = sx;
    this.args[a + 1] = sy;
    this.args[a + 2] = sw;
    this.args[a + 3] = sh;
    this.args[a + 4] = dx;
    this.args[a + 5] = dy;
    this.args[a + 6] = dw;
    this.args[a + 7] = dh;
  }
  rects(b: RectBatch): void {
    this.log(Op.rects, this.rectItems.n, b.n, 0);
    this.rectItems.append(b, 0, b.n);
  }
  sectors(b: SectorBatch): void {
    this.log(Op.sectors, this.sectorItems.n, b.n, 0);
    this.sectorItems.append(b, 0, b.n);
  }
  paths(b: PathBatch): void {
    this.log(Op.paths, this.pathItems.n, b.n, 0);
    this.pathItems.append(b, 0, b.n);
  }
  glyphs(b: GlyphBatch): void {
    this.args[this.log(Op.glyphs, this.glyphItems.n, b.n, 1)] = b.font;
    this.glyphItems.append(b, 0, b.n);
  }
  text(b: TextBatch): void {
    const a = this.log(Op.text, this.textItems.n, b.n, 4);
    this.args[a] = b.font;
    this.args[a + 1] = b.align;
    this.args[a + 2] = b.baseline;
    this.args[a + 3] = b.angle;
    this.textItems.append(b, 0, b.n);
  }
  image(img: number, r: Float64Array, smooth: boolean): void {
    const a = this.log(Op.image, img, 0, 9);
    for (let k = 0; k < 8; k++) this.args[a + k] = r[k];
    this.args[a + 8] = +smooth;
  }

  /**
   * The recorded frame drawn into `into`, its images made there and, unless `into` is a recording, disposed after. Not
   * for a hot path.
   */
  replay(into: RenderBackend): void {
    const made = new Map<number, number>();
    const A = this.args;
    for (let k = 0; k < this.opCount; k++) {
      const code = this.ops[k * 4];
      const start = this.ops[k * 4 + 1];
      const count = this.ops[k * 4 + 2];
      const a = this.ops[k * 4 + 3];
      switch (code) {
        case Op.beginPass:
          into.beginPass(A[a], A[a + 1], A[a + 2], A[a + 3]);
          break;
        case Op.endPass:
          into.endPass();
          break;
        case Op.clear:
          into.clear(A[a], A[a + 1], A[a + 2], A[a + 3]);
          break;
        case Op.pushClip:
          into.pushClip(A[a], A[a + 1], A[a + 2], A[a + 3]);
          break;
        case Op.popClip:
          into.popClip();
          break;
        case Op.copy:
          into.copy(A[a], A[a + 1], A[a + 2], A[a + 3], A[a + 4], A[a + 5]);
          break;
        case Op.present:
          into.present(A[a], A[a + 1], A[a + 2], A[a + 3]);
          break;
        case Op.stale:
          into.stale(A[a], A[a + 1], A[a + 2], A[a + 3], A[a + 4], A[a + 5], A[a + 6], A[a + 7]);
          break;
        case Op.rects:
          this.r.reset();
          this.r.append(this.rectItems, start, count);
          into.rects(this.r);
          break;
        case Op.sectors:
          this.s.reset();
          this.s.append(this.sectorItems, start, count);
          into.sectors(this.s);
          break;
        case Op.paths:
          this.p.reset();
          this.p.append(this.pathItems, start, count);
          into.paths(this.p);
          break;
        case Op.glyphs:
          this.g.reset();
          this.g.append(this.glyphItems, start, count);
          this.g.font = A[a];
          into.glyphs(this.g);
          break;
        case Op.text:
          this.t.reset();
          this.t.style(A[a], A[a + 1], A[a + 2], A[a + 3]);
          this.t.append(this.textItems, start, count);
          into.text(this.t);
          break;
        default: {
          let id = made.get(start);
          if (id === undefined) {
            id = this.made(into, start);
            made.set(start, id);
          }
          for (let j = 0; j < 8; j++) this.at[j] = A[a + j];
          into.image(id, this.at, A[a + 8] === 1);
        }
      }
    }
    // A recording keeps what it is given (a replayed copy hashes as the original); another backend frees them.
    if (!(into instanceof RecordingBackend)) for (const id of made.values()) into.images.dispose(id);
  }

  private made(into: RenderBackend, img: number): number {
    const c = this.images.content[img];
    const w = this.images.w[img];
    const h = this.images.h[img];
    if (c instanceof Uint32Array) return into.images.fromPixels(w, h, c);
    const layer = into.images.beginLayer(w, h);
    c?.replay(layer);
    return into.images.endLayer();
  }
  private rect4(code: number, x: number, y: number, w: number, h: number): void {
    const a = this.log(code, 0, 0, 4);
    this.args[a] = x;
    this.args[a + 1] = y;
    this.args[a + 2] = w;
    this.args[a + 3] = h;
  }
  /** An op logged with room for `n` arguments: their offset. */
  private log(code: number, start: number, count: number, n: number): number {
    if ((this.opCount + 1) * 4 > this.ops.length) {
      const o = new Int32Array(this.ops.length * 2);
      o.set(this.ops);
      this.ops = o;
    }
    if (this.argCount + n > this.args.length) {
      const g = new Float64Array(Math.max(this.args.length * 2, this.argCount + n));
      g.set(this.args);
      this.args = g;
    }
    const k = this.opCount++ * 4;
    this.ops[k] = code;
    this.ops[k + 1] = start;
    this.ops[k + 2] = count;
    this.ops[k + 3] = this.argCount;
    this.argCount += n;
    return this.ops[k + 3];
  }
}

/**
 * A hash of everything a recording holds (16 hex digits): the op log with each image op's image by its content, the
 * arguments, every item of every store (each field, path verbs and coordinates, text strings), and the image store
 * (each image's size and its pixels or its layer's hash, or a disposed one's mark), layers hashed alike. Two frames
 * that draw anything differently, down to a glyph's x, hash differently; for proving a change draws as before. Not for
 * a hot path.
 */
export function hashRecording(r: RecordingBackend): string {
  const h = new Hash();
  hashInto(r, h);
  return h.hex();
}

const contentHashes = new WeakMap<object, string>();
/** An image's content hashed once: pixels are copied in and a layer is closed when kept. */
function contentHash(c: Uint32Array | RecordingBackend): string {
  let s = contentHashes.get(c);
  if (s === undefined) {
    s = c instanceof Uint32Array ? new Hash().words(c, c.length).hex() : hashRecording(c);
    contentHashes.set(c, s);
  }
  return s;
}

function hashInto(r: RecordingBackend, h: Hash): void {
  const im = r.images;
  h.word(r.opCount);
  for (let k = 0; k < r.opCount; k++) {
    const code = r.ops[4 * k];
    h.word(code);
    const img = r.ops[4 * k + 1];
    const c = code === Op.image ? im.content[img] : undefined;
    if (code !== Op.image) h.word(img);
    else if (c) h.str(contentHash(c));
    else h.word(-1);
    h.word(r.ops[4 * k + 2]);
  }
  h.bytes(r.args, r.argCount * 8);
  const R = r.rectItems;
  for (const a of [R.x, R.y, R.w, R.h, R.color]) h.words(a, R.n);
  for (const a of [R.alpha, R.stroke]) h.bytes(a, R.n);
  const S = r.sectorItems;
  for (const a of [S.cx, S.cy, S.r0, S.r1, S.a0, S.a1]) h.bytes(a, S.n * 8);
  h.words(S.color, S.n);
  h.bytes(S.alpha, S.n);
  const P = r.pathItems;
  for (const a of [P.start, P.cstart, P.color]) h.words(a, P.n);
  for (const a of [P.alpha, P.mode, P.width]) h.bytes(a, P.n);
  h.bytes(P.verbs, P.nv);
  h.bytes(P.coords, P.nc * 8);
  const G = r.glyphItems;
  for (const a of [G.x, G.y, G.radius]) h.bytes(a, G.n * 8);
  for (const a of [G.glyph, G.color]) h.words(a, G.n);
  h.bytes(G.alpha, G.n);
  const T = r.textItems;
  for (const a of [T.x, T.y, T.sx, T.sy]) h.bytes(a, T.n * 8);
  h.words(T.color, T.n);
  h.bytes(T.alpha, T.n);
  for (let i = 0; i < T.n; i++) h.str(T.text[i]);
  // The images made, in order, and which were disposed.
  h.word(im.content.length);
  for (let i = 0; i < im.content.length; i++) {
    const c = im.content[i];
    h.word(im.w[i]);
    h.word(im.h[i]);
    if (c) h.str(contentHash(c));
    else h.word(-1);
  }
}

/** MurmurHash3's 32-bit mix over words, in two lanes of different seeds (64 bits). */
class Hash {
  private a = 0x9747b28c;
  private b = 0x2a5f1c3d;
  word(x: number): void {
    let k = Math.imul(x | 0, 0xcc9e2d51);
    k = Math.imul((k << 15) | (k >>> 17), 0x1b873593);
    let a = this.a ^ k;
    a = (a << 13) | (a >>> 19);
    this.a = (Math.imul(a, 5) + 0xe6546b64) | 0;
    let b = this.b ^ k;
    b = (b << 17) | (b >>> 15);
    this.b = (Math.imul(b, 9) + 0x85ebca6b) | 0;
  }
  /** The first `n` elements of a 32-bit array, then `n`. */
  words(v: Int32Array | Uint32Array, n: number): this {
    for (let i = 0; i < n; i++) this.word(v[i]);
    this.word(n);
    return this;
  }
  /** The first `n` bytes of a typed array, then `n`. */
  bytes(v: ArrayBufferView, n: number): void {
    const u = new Uint8Array(v.buffer, v.byteOffset, n);
    const m = n >> 2;
    if ((v.byteOffset & 3) === 0) {
      const w = new Uint32Array(v.buffer, v.byteOffset, m);
      for (let i = 0; i < m; i++) this.word(w[i]);
    } else
      for (let i = 0; i < m; i++)
        this.word(u[4 * i] | (u[4 * i + 1] << 8) | (u[4 * i + 2] << 16) | (u[4 * i + 3] << 24));
    for (let i = 4 * m; i < n; i++) this.word(u[i]);
    this.word(n);
  }
  str(s: string): void {
    for (let i = 0; i < s.length; i++) this.word(s.charCodeAt(i));
    this.word(s.length);
  }
  hex(): string {
    return (this.a >>> 0).toString(16).padStart(8, "0") + (this.b >>> 0).toString(16).padStart(8, "0");
  }
}
