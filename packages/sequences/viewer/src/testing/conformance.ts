/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The renderer conformance kit (seqspec:viewer#extension-points): `rendererConformance(entry, fixtures)`, run by every
// renderer, built-in or a host's, as it lands. Per declared geometry it opens the fixture's document, draws the entry's
// layers at three zooms (wrapped: its lines) through a small composer of its own (CS1's is the real one), and checks:
// the same batches on the recording backend as on a tee's two sides, a recording and the SVG backend (the SVG side saw
// the same calls, wrote every item, and wrote only finite numbers); no
// allocation per frame after warm-up (the fixture's heap guard); every drawn rect and glyph hits back to its layer, an
// item one of its data's; the legend lists the fixture's labels; a context without the entry (or a geometry it does
// not declare) draws the fallback with one issue per layer. An entry with no draw yet (LY3b's headless halves) has its
// hits sampled over its entries instead. With `reversed` (AX1), a reversed piece draws the forward draw mirrored within
// each block (`mirror`). On rings (CG1b; seqspec:viewer#the-render-backend) an entry that declares `circular` is drawn
// through the engine's ring map and its sectors and glyphs hit back as above; one that does not is drawn as bars by the
// composer, which says so once per layer (`svs.render.fallback`).

import { Lanes } from "@molstar/seq-core/format";
import { Annotations } from "@molstar/seq-core/model";
import type { SeqViewerContext } from "../context.ts";
import {
  type Axis,
  AxisGeometry,
  type BackendImages,
  DrawView,
  type GlyphBatch,
  Hit,
  HoverKind,
  type LayerRenderer,
  type LegendSink,
  NO_COLOR,
  Pass,
  type PathBatch,
  type RectBatch,
  type RenderBackend,
  type SectorBatch,
  type ShownIndex,
  SkinToken,
  TextAlign,
  type TextBatch,
} from "../draw/api.ts";
import { Dirty, SectionComposer } from "../draw/composer.ts";
import { Miniature } from "../draw/miniature.ts";
import { drawResources } from "../draw/resources.ts";
import { RingBackend } from "../draw/ring.ts";
import type { Skin } from "../graph.ts";
import { EntryKind, entryBorder, facingOf, SectionLayouts } from "../layout.ts";
import { Ring } from "../ring.ts";
import { lanesInView } from "../streams.ts";
import { SvgBackend } from "../svg/backend.ts";
import { type EdgeKind, EdgeKinds } from "./edges.ts";
import { Op, RecordingBackend } from "./recording.ts";

export interface ConformanceFixtures {
  /**
   * A context showing a document with layers of the entry's type, sized; with the entry registered, or without it;
   * `direction` -1: every section's pieces reversed (asked only with `reversed`).
   */
  open(o: {
    /** `circular`: every section on rings (`style: "circular"`). */
    readonly geometry: "linear" | "wrapped" | "circular";
    readonly entry: boolean;
    readonly direction?: 1 | -1;
    /** The document in an edge geometry (`edgeDocument`; asked only with `edges`, linear, with the entry). */
    readonly edge?: EdgeKind;
  }): Promise<SeqViewerContext>;
  /**
   * Whether `open` takes `direction`: then, linear at the fit, every rect, glyph and path of the reversed draw is the
   * forward one mirrored within its block, to 2 device px (centered text where labels were not thinned; images by
   * the pixels each device column shows). Items
   * across blocks (a pair's arc between two) and cell gaps (on each cell's right either way) are left out.
   */
  readonly reversed?: boolean;
  /**
   * Whether `open` takes `edge`: then, linear, every layer of the type is drawn and hit back in each edge geometry
   * (`EdgeKinds`: reversed, a cyclic origin, stacked, filtered, copies), at the fit and at 16 px a position from
   * either end.
   */
  readonly edges?: boolean;
  /**
   * The edge geometries the fixture's document cannot take, each with why (an alignment's columns have no cyclic
   * origin): reported as skipped. Every other one must apply (a reversed piece, an origin past 0, a stacked section, a
   * filter, more blocks for copies), else its check fails: a rewrite that changed nothing checks nothing.
   */
  readonly edgesSkipped?: Readonly<Partial<Record<EdgeKind, string>>>;
  /**
   * The zooms that may draw nothing (`fit`: the fit, the edge geometries' too; `middle`: √(16·fit); `close`: 16 px a
   * position): at any other `batches` fails a draw with no item (a rect, sector, path, glyph, text or image), and a
   * layer of the type must draw an item at one zoom at least, so a renderer that regressed to drawing nothing does
   * not pass (VB-04). A fixture that fits at less than a px a position exercises the density and LOD paths too.
   */
  readonly empty?: readonly ConformanceZoom[];
  /** The labels the legend lists for the data, in order; omitted: listed, not compared. */
  readonly legend?: readonly string[];
  /** A guard that `frame` allocates nothing after warm-up (Node: the tests' `heapScaling`); omitted: not checked. */
  heap?(frame: () => void): void | Promise<void>;
  /**
   * Whether a drawn element hits an item of the layer's data, for every geometry or per geometry (default: when the
   * entry has its own `hit`).
   */
  readonly items?: boolean | Readonly<Partial<Record<"linear" | "wrapped" | "circular", boolean>>>;
}
/** The kit's zooms, linear: the fit, √(16·fit) and 16 px a position. */
export type ConformanceZoom = "fit" | "middle" | "close";
const ZOOMS: readonly ConformanceZoom[] = ["fit", "middle", "close"];
export interface ConformanceCheck {
  readonly check: "batches" | "heap" | "hits" | "legend" | "fallback" | "mirror" | "edges";
  readonly geometry: string;
  readonly ok: boolean;
  readonly detail: string;
}

/** Every call to both backends, in order; `a`'s images ids are returned (both are fresh, so they agree). */
export class TeeBackend implements RenderBackend {
  readonly a: RenderBackend;
  readonly b: RenderBackend;
  readonly caps: number;
  readonly images: BackendImages;
  constructor(a: RenderBackend, b: RenderBackend) {
    this.a = a;
    this.b = b;
    this.caps = a.caps & b.caps;
    this.images = {
      fromPixels: (w, h, rgba) => {
        b.images.fromPixels(w, h, rgba);
        return a.images.fromPixels(w, h, rgba);
      },
      beginLayer: (w, h) => new TeeBackend(a.images.beginLayer(w, h), b.images.beginLayer(w, h)),
      endLayer: () => {
        b.images.endLayer();
        return a.images.endLayer();
      },
      dispose: (img) => {
        a.images.dispose(img);
        b.images.dispose(img);
      },
    };
  }
  beginPass(pass: number, w: number, h: number, dpr: number) {
    this.a.beginPass(pass, w, h, dpr);
    this.b.beginPass(pass, w, h, dpr);
  }
  endPass() {
    this.a.endPass();
    this.b.endPass();
  }
  clear(x: number, y: number, w: number, h: number) {
    this.a.clear(x, y, w, h);
    this.b.clear(x, y, w, h);
  }
  pushClip(x: number, y: number, w: number, h: number) {
    this.a.pushClip(x, y, w, h);
    this.b.pushClip(x, y, w, h);
  }
  popClip() {
    this.a.popClip();
    this.b.popClip();
  }
  copy(dx: number, dy: number, x: number, y: number, w: number, h: number) {
    this.a.copy(dx, dy, x, y, w, h);
    this.b.copy(dx, dy, x, y, w, h);
  }
  present(x: number, y: number, w: number, h: number) {
    this.a.present(x, y, w, h);
    this.b.present(x, y, w, h);
  }
  stale(sx: number, sy: number, sw: number, sh: number, dx: number, dy: number, dw: number, dh: number) {
    this.a.stale(sx, sy, sw, sh, dx, dy, dw, dh);
    this.b.stale(sx, sy, sw, sh, dx, dy, dw, dh);
  }
  rects(r: RectBatch) {
    this.a.rects(r);
    this.b.rects(r);
  }
  sectors(s: SectorBatch) {
    this.a.sectors(s);
    this.b.sectors(s);
  }
  paths(p: PathBatch) {
    this.a.paths(p);
    this.b.paths(p);
  }
  glyphs(g: GlyphBatch) {
    this.a.glyphs(g);
    this.b.glyphs(g);
  }
  text(t: TextBatch) {
    this.a.text(t);
    this.b.text(t);
  }
  image(img: number, r: Float64Array, smooth: boolean) {
    this.a.image(img, r, smooth);
    this.b.image(img, r, smooth);
  }
}

/** The first difference between two recorded frames (ops, arguments, every item store); null: the same. */
export function frameDiff(a: RecordingBackend, b: RecordingBackend): string | null {
  const same = (x: unknown, y: unknown): boolean =>
    ArrayBuffer.isView(x) || Array.isArray(x)
      ? (x as ArrayLike<unknown>).length === (y as ArrayLike<unknown>).length &&
        Array.prototype.every.call(x, (v, i) => Object.is(v, (y as ArrayLike<unknown>)[i]))
      : Object.is(x, y);
  if (a.opCount !== b.opCount || !same(a.ops, b.ops) || !same(a.args, b.args)) return "ops";
  for (const store of ["rectItems", "sectorItems", "pathItems", "glyphItems", "textItems"] as const)
    for (const [k, v] of Object.entries(a[store])) if (!same(v, (b[store] as never)[k])) return `${store}.${k}`;
  return null;
}

const OPS = Object.fromEntries(Object.entries(Op).map(([name, code]) => [code, name]));
/** The arguments each call logs. */
const CALL_ARGS: Record<number, number> = {
  [Op.beginPass]: 4,
  [Op.clear]: 4,
  [Op.pushClip]: 4,
  [Op.copy]: 6,
  [Op.present]: 4,
  [Op.stale]: 8,
};
/** A recording's calls as the SVG backend's trace lines (`SvgOptions.trace`): the op, and its arguments or item count. */
export function traceOf(r: RecordingBackend): string[] {
  const out: string[] = [];
  for (let k = 0; k < r.opCount; k++) {
    const [code, start, count, at] = r.ops.subarray(4 * k, 4 * k + 4);
    const name = OPS[code];
    const a = (n: number) => Array.from(r.args.subarray(at, at + n));
    out.push(
      code === Op.rects || code === Op.sectors || code === Op.paths
        ? `${name} ${count}`
        : code === Op.glyphs
          ? `glyphs ${count} ${r.args[at]}`
          : code === Op.text
            ? `text ${count} ${a(4).join(" ")}`
            : code === Op.image
              ? `image ${start} ${a(9).join(" ")}`
              : [name, ...a(CALL_ARGS[code] ?? 0)].join(" "),
    );
  }
  return out;
}
/**
 * What the SVG side did with the frame the recording saw: the same calls in the same order, an element for every item
 * (NO_COLOR ones are left out by both), and no number that was not finite. Null: all of it.
 */
export function svgDiff(r: RecordingBackend, svg: SvgBackend): string | null {
  const want = traceOf(r);
  const got = svg.trace ?? [];
  for (let k = 0; k < Math.max(want.length, got.length); k++)
    if (want[k] !== got[k]) return `svg call ${k}: ${got[k] ?? "(none)"}, recorded ${want[k] ?? "(none)"}`;
  const kept = (c: Int32Array, n: number) => c.subarray(0, n).filter((v) => v !== NO_COLOR).length;
  const items = {
    rects: kept(r.rectItems.color, r.rectItems.n),
    sectors: kept(r.sectorItems.color, r.sectorItems.n),
    paths: kept(r.pathItems.color, r.pathItems.n),
    glyphs: kept(r.glyphItems.color, r.glyphItems.n),
    texts: kept(r.textItems.color, r.textItems.n),
  };
  // Muted passes (the overlay and navigator by default) write nothing, so only an SVG that draws every pass compares.
  for (const [kind, n] of Object.entries(items))
    if (svg.stats[kind as keyof typeof items] !== n)
      return `svg wrote ${svg.stats[kind as keyof typeof items]} ${kind}, recorded ${n}`;
  return svg.stats.nonFinite > 0 ? `svg wrote ${svg.stats.nonFinite} number(s) that were not finite` : null;
}

const reset = (b: RenderBackend) => {
  if (b instanceof RecordingBackend || b instanceof SvgBackend) b.reset();
};

/** Caches per prepared state, held across frames and kits as the composer holds them (dropped with the state). */
const caches = new WeakMap<object, { skin: unknown; c: unknown }>();
function cacheOf(entry: LayerRenderer, state: unknown, skin: Skin): unknown {
  const held = state !== null && typeof state === "object" ? caches.get(state) : undefined;
  if (held && held.skin === skin) return held.c;
  const c = entry.cache?.create(state, skin);
  if (state !== null && typeof state === "object") caches.set(state, { skin, c });
  return c;
}

/**
 * The kit's composer: the layer drawn on each of its layout entries, as CS1 will (`DrawView` filled per entry, the
 * section's numberings from the published view); `tiny` draws its `tiny` member at a height below `minHeight`. Tests
 * draw a layer into any backend through it.
 */
export class LayerFrame {
  readonly v = new DrawView();
  readonly direct = new RecordingBackend();
  /** A recording and the SVG backend (every pass, a call trace), given the same calls. */
  readonly tee: TeeBackend;
  readonly ctx: SeqViewerContext;
  readonly entry: LayerRenderer;
  readonly section: string;
  readonly k: number;
  constructor(ctx: SeqViewerContext, entry: LayerRenderer, section: string, k: number) {
    this.ctx = ctx;
    this.entry = entry;
    this.section = section;
    this.k = k;
    this.tee = new TeeBackend(
      new RecordingBackend(),
      new SvgBackend(drawResources(ctx), { overlay: true, nav: true, trace: true }),
    );
  }
  tiny = false;
  /** The layer's cache (made once per prepared state) and the tracks' miniature host, as the composer holds them. */
  private cache: unknown = null;
  private mini: Miniature | null = null;
  /** On rings, the map and the backend through it, as the composer holds them (mapped at a new layout only). */
  private readonly ring = new Ring();
  private readonly ringBe = new RingBackend();
  private ringAt = -1;
  draw(to: RenderBackend): number {
    if (to instanceof TeeBackend) {
      reset(to.a);
      reset(to.b);
    } else reset(to);
    const { ctx, v, k, section } = this;
    const m = ctx.graph.assembled.value?.view;
    const live = ctx.layout.get(section);
    const l = live?.current;
    const ax = ctx.layout.axis(section);
    const draw = this.tiny ? this.entry.tiny : this.entry.draw;
    if (!m || !l || !live || !ax?.axis || !draw) return 0;
    const rings = l.geometry === AxisGeometry.circular;
    if (rings && this.ringAt !== live.version.value) {
      this.ringAt = live.version.value;
      this.ringBe.map(this.ring.from(l, ax, ax.dpr));
    }
    this.ringBe.to = to;
    const backend = rings ? this.ringBe : to;
    const t = m.resolved.layers.track[k];
    const p = m.prepared[k];
    backend.beginPass(Pass.main, Math.ceil(l.width * ax.dpr), Math.ceil(l.height * ax.dpr), ax.dpr);
    this.cache = this.entry.cache ? cacheOf(this.entry, p.state, ctx.skin.value) : null;
    const rows = l.rows;
    const banded = rows.layer === k && (rows.lazy || rows.height < ctx.layout.props(section).denseRowHeight);
    let n = 0;
    for (let e = 0; e < l.entryCount; e++) {
      if (l.entryKind[e] !== EntryKind.track || l.entryTrack[e] !== t) continue;
      v.from(ax);
      v.ring = rings ? this.ring : null;
      if (rings) {
        v.geometry = AxisGeometry.linear;
        this.ring.drawn = v;
      }
      v.backend = backend;
      v.skin = ctx.skin.value;
      v.props = ctx.layout.props(section);
      v.text = ctx.text;
      v.layer = k;
      v.setThemes(p.themes);
      v.stream = ctx.streams.get(ctx.graph.assembled.value?.set.strings.get(m.resolved.layers.ref[k]) ?? "");
      v.top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - ax.scrollY;
      // The entry's content stands above its border, as the composer draws it (`SectionLayouts.entryBorder`).
      v.height = this.tiny
        ? Math.max(2, (this.entry.minHeight ?? 6) - 1)
        : l.entryHeight[e] -
          (entryBorder(l, e, v.props.rowBorders) ? (Math.max(1, Math.round(ax.dpr)) | 0) / ax.dpr : 0);
      v.numberingFrom(m.resolved, ctx.layout.sectionIndex(section));
      v.lanes[0] = v.lanes[1] = 0;
      lanesInView(l, v, v.lanes);
      // Only what is in view, as the composer draws.
      if ((e >= l.pinnedCount && v.top + v.height <= l.scrolledTop) || v.top >= l.height) continue;
      v.facing = facingOf(l, ax.axis, k, t, e);
      v.cache = this.cache;
      v.tracks = this.mini ??= new Miniature(ctx, section);
      v.dataRow = v.row = -1;
      if (rows.layer >= 0 && t === rows.track && e >= rows.entry && !rows.lazy) {
        const lane = SectionLayouts.laneOf(rows, l.entryRow[e]);
        v.row = lane >= 0 ? lane - rows.laneStart : -1;
        v.dataRow = lane >= 0 ? Lanes.rowOf(m.resolved.lanes, lane) : -1;
      }
      if (l.geometry === AxisGeometry.wrapped) {
        const line = l.entryRow[e];
        v.line = line;
        v.onlyBlock = l.lines.block[line];
        v.lineStart = l.lines.start[line];
        v.lineStop = l.lines.stop[line];
      } else v.onlyBlock = ax.axis.stacked ? l.entryBlock[e] : -1;
      if (banded && !this.tiny && this.entry.rows) {
        // Dense or lazy rows: one `rows` call over the display rows in view (the composer's band).
        const top = l.scrolledTop + l.entryTop[rows.entry] - ax.scrollY;
        const d0 = Math.max(0, Math.floor((l.scrolledTop - top) / rows.height));
        const d1 = Math.min(rows.count - 1, Math.ceil((l.height - top) / rows.height) - 1);
        v.top = top + d0 * rows.height;
        v.height = rows.height;
        v.row = v.dataRow = v.onlyBlock = -1;
        if (d1 >= d0) this.entry.rows(v, p.state, m.resolved.layers.params[k], d0, d1);
        n++;
        break;
      }
      draw(v, p.state, m.resolved.layers.params[k]);
      n++;
    }
    backend.endPass();
    return n;
  }
}

/** The layers of the entry's type (declared and drawn as it) in the current view: [section ref, layer index]. */
function layersOf(ctx: SeqViewerContext, name: string, drawn: boolean): [string, number, string][] {
  const m = ctx.graph.assembled.value?.view;
  const out: [string, number, string][] = [];
  if (!m) return out;
  const { layers, tracks, sections } = m.resolved;
  const s = (id: number) => m.set.strings.get(id);
  for (let k = 0; k < layers.count; k++)
    if (s((drawn ? layers.type : layers.declared)[k]) === name)
      out.push([s(sections.ref[tracks.section[layers.track[k]]]), k, s(layers.ref[k])]);
  return out;
}

/** The items a recording drew: rects, sectors, paths, glyphs, texts and images. */
function drawnCount(r: RecordingBackend): number {
  let n = r.rectItems.n + r.sectorItems.n + r.pathItems.n + r.glyphItems.n + r.textItems.n;
  for (let o = 0; o < r.opCount; o++) if (r.ops[4 * o] === Op.image) n++;
  return n;
}

/** Item hits must name the layer's data. */
function inData(ctx: SeqViewerContext, k: number, annotation: number): boolean {
  const rv = ctx.graph.assembled.value?.view?.resolved;
  if (!rv) return false;
  for (let d = rv.layers.dataStart[k]; d < rv.layers.dataStart[k] + rv.layers.dataCount[k]; d++)
    if (rv.data[d] === annotation) return true;
  return false;
}

/**
 * Whether the item a hit names covers the position it names (LC-17: a hit naming another item of its data passed): a
 * `ranges` row of the item, on the block's own unit, over the index, within the device px's positions (below a device px
 * a cell, `w`, a px holds several); a boundary, a position either side. Other geometries and placements through a path
 * pass.
 */
function covers(ctx: SeqViewerContext, hit: Hit, ax: Axis | null, w: number): boolean {
  const set = ctx.graph.assembled.value?.view?.set;
  const ann = set?.annotations[hit.annotation];
  if (!set || !ann || !ax || ann.geometry !== "ranges" || hit.block < 0 || hit.index < 0) return true;
  const per = w < 1 ? Math.ceil(1 / w) : 0;
  const { start, rows } = Annotations.itemRows(ann);
  const unit = ax.unit[hit.block];
  let mine = false;
  for (let j = start[hit.item]; j < start[hit.item + 1]; j++) {
    const r = rows[j];
    const u = ann.slots.unit[ann.slot[r]];
    // On the block's own unit (a slot over a molecule: each of its units).
    if (unit < 0 || (u !== unit && (u >= 0 || ann.slots.molecule[ann.slot[r]] !== set.units.molecule[unit]))) continue;
    mine = true;
    const a = ann.start[r];
    const e = Math.max(ann.stop[r], a + 1);
    if (hit.index >= a - 1 - per && hit.index <= e + per) return true;
  }
  return !mine;
}

type Say = (check: ConformanceCheck["check"], ok: boolean, detail: string) => void;

/**
 * The format's fallback (linear, wrapped): without the entry, and in a geometry it does not declare, one issue per
 * layer; a declared geometry its format entry's layouts exclude, a mismatch. False: nothing more to check in it.
 */
async function formatFallback(
  entry: LayerRenderer,
  fx: ConformanceFixtures,
  geometry: "linear" | "wrapped",
  say: Say,
): Promise<boolean> {
  const declared = entry.geometries.includes(geometry);
  const ctx0 = await fx.open({ geometry, entry: !declared });
  const fell = layersOf(ctx0, entry.name, false);
  const issues = ctx0.graph.issues.value.filter(
    (i) => (i.code === "svs.render.fallback" || i.code === "svs.render.omitted") && fell.some((f) => f[2] === i.ref),
  );
  const kept = fell.filter(([, k]) => layersOf(ctx0, entry.name, true).some((d) => d[1] === k));
  say(
    "fallback",
    fell.length > 0 && issues.length === fell.length && kept.length === 0,
    `${fell.length} layer(s), ${issues.length} issue(s) [${issues.map((i) => i.code).join(", ")}]`,
  );
  const layouts = ctx0.registry.svs.get("layerTypes", entry.name)?.layouts ?? [];
  ctx0.dispose();
  if (declared && !layouts.includes(geometry === "wrapped" ? "text" : "tracks")) {
    say("fallback", false, `declares ${geometry}, but its type's layouts are ${layouts.join(", ")}`);
    return false;
  }
  return declared;
}

/**
 * On rings without `circular`: the composer draws each layer of the type as bars, its sectors in the layer's ring, and
 * reports `svs.render.fallback` once per layer, naming it.
 */
async function ringless(entry: LayerRenderer, fx: ConformanceFixtures): Promise<[boolean, string]> {
  const ctx = await fx.open({ geometry: "circular", entry: true });
  const said: string[] = [];
  const sub = ctx.events.issue.subscribe((i) => {
    if (i.code === "svs.render.fallback") said.push(i.message);
  });
  // A layer in a collapsed group or a hidden track is not laid out, so not drawn.
  const laid = (section: string, k: number) => {
    const l = ctx.layout.get(section)?.current;
    const t = ctx.graph.assembled.value?.view?.resolved.layers.track[k] ?? -2;
    return !!l?.entryTrack.subarray(0, l.entryCount).includes(t);
  };
  const layers = layersOf(ctx, entry.name, true).filter(([section, k]) => laid(section, k));
  let bars = 0;
  try {
    for (const section of new Set(layers.map((x) => x[0]))) {
      if (ctx.layout.get(section)?.current.geometry !== AxisGeometry.circular)
        return [false, `${section} is not on rings`];
      const c = new SectionComposer(ctx, section);
      const r = new RecordingBackend();
      c.frame(r, Dirty.tracks);
      c.frame(r, Dirty.tracks);
      c.dispose();
      bars += r.sectorItems.n;
    }
  } finally {
    sub.unsubscribe();
    ctx.dispose();
  }
  // Another type's layers in the fixture fall back too: only this type's are counted.
  const mine = said.filter((m) => layers.some(([, , ref]) => m.includes(`"${ref}"`)));
  const named = layers.filter(([, , ref]) => mine.some((m) => m.includes(`"${ref}"`)));
  const ok = layers.length > 0 && mine.length === layers.length && named.length === layers.length && bars > 0;
  return [ok, `${layers.length} layer(s), ${mine.length} issue(s), ${bars} sectors`];
}

/** Runs the kit; one result per check and geometry (linear, wrapped and circular). */
export async function rendererConformance(entry: LayerRenderer, fx: ConformanceFixtures): Promise<ConformanceCheck[]> {
  const out: ConformanceCheck[] = [];
  for (const geometry of ["linear", "wrapped", "circular"] as const) {
    const say = (check: ConformanceCheck["check"], ok: boolean, detail: string) =>
      out.push({ check, geometry, ok, detail });
    const declared = entry.geometries.includes(geometry);
    if (geometry === "circular" && !declared) {
      say("fallback", ...(await ringless(entry, fx)));
      continue;
    }
    if (geometry !== "circular" && !(await formatFallback(entry, fx, geometry, say))) continue;

    const ctx = await fx.open({ geometry, entry: true });
    const layers = layersOf(ctx, entry.name, true);
    // A layer of another type drawn as this one (its fallback) may have nothing this type draws.
    const own = new Set(layersOf(ctx, entry.name, false).map((x) => x[1]));
    if (ctx.registry.get("renderers", entry.name) !== entry || layers.length === 0) {
      say("batches", false, "the entry is not registered, or no layer of its type is drawn");
      ctx.dispose();
      continue;
    }
    const hit = Hit.create();
    const labels: string[] = [];
    const sink: LegendSink = {
      item: (label) => labels.push(label),
      ramp: () => labels.push("(ramp)"),
      size: (label) => labels.push(label),
    };
    let shown = 0;
    for (const [section, k] of layers) {
      // A layer in a collapsed group or a hidden track has no entries to draw or hit.
      const l = ctx.layout.get(section)?.current;
      const t = ctx.graph.assembled.value?.view?.resolved.layers.track[k];
      if (!l?.entryTrack.subarray(0, l.entryCount).includes(t ?? -2)) continue;
      shown++;
      if ((geometry === "circular") !== (l.geometry === AxisGeometry.circular)) {
        say("batches", false, `${section}: the fixture's section is ${l.geometry === 2 ? "" : "not "}on rings`);
        continue;
      }
      const f = new LayerFrame(ctx, entry, section, k);
      const vp = ctx.viewports.get(section);
      vp?.fit();
      const fit = vp?.live.current.cell ?? 1;
      const cells = geometry === "linear" ? [fit, Math.sqrt(fit * 16), 16] : [fit];
      let drawn = 0;
      for (let z = 0; z < cells.length; z++) {
        const cell = cells[z];
        if (geometry === "linear") vp?.setView(cell, 0, "api");
        const at = `${section} layer ${k} at ${cell.toFixed(2)} px`;
        // The draw, then below the row's minimum (`tiny`); the hits on the draw's elements.
        for (const tiny of [true, false]) {
          f.tiny = tiny;
          if (!(tiny ? entry.tiny : entry.draw)) continue;
          const where = tiny ? `${at} (tiny)` : at;
          const n = f.draw(f.direct);
          f.draw(f.tee);
          const diff =
            frameDiff(f.direct, f.tee.a as RecordingBackend) ??
            svgDiff(f.tee.a as RecordingBackend, f.tee.b as SvgBackend);
          const items = drawnCount(f.direct);
          if (!tiny) drawn += items;
          say(
            "batches",
            n > 0 && diff === null && (items > 0 || fx.empty?.includes(ZOOMS[z]) === true || !own.has(k)),
            `${where}: ${n} entries, ${f.direct.opCount} ops, ${items} items${diff ? `; ${diff}` : ""}`,
          );
          if (fx.heap) {
            try {
              await fx.heap(() => f.draw(f.direct));
              say("heap", true, where);
            } catch (e) {
              say("heap", false, `${where}: ${e instanceof Error ? e.message : String(e)}`);
            }
          }
        }
        const items = typeof fx.items === "object" ? fx.items[geometry] : fx.items;
        out.push({ check: "hits", geometry, ...hits(ctx, f, hit, at, items ?? entry.hit !== undefined) });
      }
      if (entry.draw && own.has(k) && drawn === 0)
        say("batches", false, `${section} layer ${k}: no item drawn at any zoom`);
      const p = ctx.graph.assembled.value?.view?.prepared[k];
      if (p) entry.legend?.(ctx.graph.assembled.value?.view?.resolved.layers.params[k], p.themes, sink, p.state);
    }
    if (shown === 0) say("hits", false, "no layer of its type is laid out");
    const ok = fx.legend ? JSON.stringify(labels) === JSON.stringify(fx.legend) : true;
    say("legend", ok, labels.join(", "));
    ctx.dispose();
    if (fx.reversed && geometry === "linear") {
      const rev = await fx.open({ geometry, entry: true, direction: -1 });
      const fwd = await fx.open({ geometry, entry: true, direction: 1 });
      for (const [section, k] of layersOf(fwd, entry.name, true)) {
        for (const c of [fwd, rev]) {
          const vp = c.viewports.get(section);
          vp?.fit();
        }
        const a = drawnItems(fwd, entry, section, k, true);
        const b = drawnItems(rev, entry, section, k, false);
        const miss = unmatched(a, b);
        const extra = unmatched(b, a);
        const at = `${section} layer ${k}: ${a.length} items`;
        const show = (x: Item[]) => x.slice(0, 2).map((i) => `${i.kind} ${i.x}+${i.w} at ${i.y}`);
        say(
          "mirror",
          miss.length === 0 && extra.length === 0,
          `${at}, ${miss.length} unmatched ${show(miss)}, ${extra.length} extra ${show(extra)}`,
        );
      }
      rev.dispose();
      fwd.dispose();
    }
    if (fx.edges && geometry === "linear") {
      // Each section's blocks without a rewrite: copies must add some.
      const base = await fx.open({ geometry, entry: true });
      const blocks = new Map(base.viewports.sections.value.map((s) => [s, base.layout.axis(s)?.axis?.count ?? 0]));
      base.dispose();
      for (const edge of EdgeKinds) {
        const why = fx.edgesSkipped?.[edge];
        if (why) say("edges", true, `${edge}: skipped, ${why}`);
        else await edges(entry, fx, edge, blocks, say);
      }
    }
  }
  return out;
}

/**
 * One edge geometry: each layer of the type drawn (the same batches on the tee) and hit back, at the fit and at 16 px
 * from either end.
 */
async function edges(
  entry: LayerRenderer,
  fx: ConformanceFixtures,
  edge: EdgeKind,
  blocks: ReadonlyMap<string, number>,
  say: Say,
): Promise<void> {
  const ctx = await fx.open({ geometry: "linear", entry: true, edge });
  try {
    const layers = layersOf(ctx, entry.name, true);
    if (layers.length === 0) say("edges", false, `${edge}: no layer of its type is drawn`);
    // A layer of another type drawn as this one (its fallback) may have nothing this type draws.
    const own = new Set(layersOf(ctx, entry.name, false).map((x) => x[1]));
    const hit = Hit.create();
    for (const section of new Set(layers.map((x) => x[0]))) {
      const ax = ctx.layout.axis(section)?.axis;
      const applied =
        !!ax &&
        (edge === "reversed"
          ? ax.dir.subarray(0, ax.count).includes(-1)
          : edge === "origin"
            ? ax.first.subarray(0, ax.count).some((f, b) => f > 0 && ax.dir[b] === 1)
            : edge === "stacked"
              ? ax.stacked
              : edge === "filtered"
                ? ax.positions.some((p) => p !== null)
                : ax.count > (blocks.get(section) ?? 0));
      if (!applied) say("edges", false, `${edge}: ${section}: the geometry did not apply (its document kept it out)`);
    }
    for (const [section, k] of layers) {
      const l = ctx.layout.get(section)?.current;
      const t = ctx.graph.assembled.value?.view?.resolved.layers.track[k];
      if (!l?.entryTrack.subarray(0, l.entryCount).includes(t ?? -2)) continue;
      const f = new LayerFrame(ctx, entry, section, k);
      const vp = ctx.viewports.get(section);
      vp?.fit();
      const fit = vp?.live.current.cell ?? 1;
      // At the fit, then 16 px a position from the start and at the end (a seam, a copy, a stacked block's end).
      for (const [cell, x] of [
        [fit, 0],
        [16, 0],
        [16, 1e9],
      ]) {
        vp?.setView(cell, x, "api");
        const at = `${edge}: ${section} layer ${k} at ${cell.toFixed(2)} px${x > 0 ? ", the end" : ""}`;
        if (entry.draw) {
          f.tiny = false;
          const n = f.draw(f.direct);
          f.draw(f.tee);
          const diff =
            frameDiff(f.direct, f.tee.a as RecordingBackend) ??
            svgDiff(f.tee.a as RecordingBackend, f.tee.b as SvgBackend);
          // At the fit something is drawn (VB-04), unless the layer is another type's fallback or the fixture names
          // the fit empty; nearer, the view may hold no item.
          const items = drawnCount(f.direct);
          const need = cell === fit && x === 0 && own.has(k) && !fx.empty?.includes("fit");
          say(
            "edges",
            n > 0 && diff === null && (items > 0 || !need),
            `${at}: ${n} entries, ${items} items${diff ? `; ${diff}` : ""}`,
          );
        }
        const items = typeof fx.items === "object" ? fx.items.linear : fx.items;
        const h = hits(ctx, f, hit, at, items ?? entry.hit !== undefined);
        say("edges", h.ok, h.detail);
      }
    }
    // Reversed, the draw at the fit is the forward one mirrored within each block (`mirror`, as with `reversed`
    // fixtures): an entry that ignores the direction draws and hits back alike (LC-17).
    if (edge === "reversed" && entry.draw) {
      const fwd = await fx.open({ geometry: "linear", entry: true });
      try {
        for (const [section, k] of layersOf(fwd, entry.name, true)) {
          if (!layers.some((x) => x[0] === section && x[1] === k)) continue;
          fwd.viewports.get(section)?.fit();
          ctx.viewports.get(section)?.fit();
          const a = drawnItems(fwd, entry, section, k, true);
          const b = drawnItems(ctx, entry, section, k, false);
          const miss = unmatched(a, b);
          const extra = unmatched(b, a);
          say(
            "mirror",
            miss.length === 0 && extra.length === 0,
            `reversed: ${section} layer ${k}: ${a.length} items, ${miss.length} unmatched, ${extra.length} extra`,
          );
        }
      } finally {
        fwd.dispose();
      }
    }
  } catch (e) {
    say("edges", false, `${edge}: ${e instanceof Error ? e.message : String(e)}`);
  } finally {
    ctx.dispose();
  }
}

/** A drawn rect, glyph, path (its x hull), centered letter or image column, device px; `key` its style. */
export interface MirrorItem {
  kind: number;
  x: number;
  w: number;
  y: number;
  h: number;
  key: string;
}
type Item = MirrorItem;
/** A layer's items drawn at the current view; `mirror`: x mirrored within its block. */
function drawnItems(ctx: SeqViewerContext, entry: LayerRenderer, section: string, k: number, mirror: boolean): Item[] {
  const f = new LayerFrame(ctx, entry, section, k);
  f.draw(f.direct);
  const v = ctx.layout.axis(section);
  const ax = v?.axis;
  if (!v || !ax) return [];
  const ends = Array.from({ length: ax.count }, (_, b) => [
    Math.round(v.x(b, 0 as ShownIndex) * v.dpr),
    Math.round(v.x(b, ax.shown[b] as ShownIndex) * v.dpr),
  ]);
  return mirrorItems(f.direct, ends, mirror, (font) => ctx.text.advances(font), v.dpr);
}
/**
 * What a recording drew (in `pass`; -1: every pass) as items, each in the block extent of `ends` ([x0, x1] device px
 * per block) that holds it (items across blocks left out), its x mirrored there when `mirror`: rects (less cell gaps,
 * a device px in the surface color on each cell's right in either direction), glyphs, paths by their x hull, centered
 * text, single letters stretched from their left by their center (`advances`: the font's), images by device column.
 */
export function mirrorItems(
  rec: RecordingBackend,
  ends: readonly number[][],
  mirror: boolean,
  advances: (font: number) => ArrayLike<number>,
  dpr: number,
  pass = -1,
): MirrorItem[] {
  const out: Item[] = [];
  const add = (kind: number, x0: number, x1: number, y: number, h: number, key: string) => {
    const slack = kind === 4 ? 0 : 2;
    const e = ends.find(([a, z]) => x0 >= a - slack && x1 <= z + slack);
    if (!e) return;
    const [m0, m1] = mirror ? [e[0] + e[1] - x1, e[0] + e[1] - x0] : [x0, x1];
    out.push({ kind, x: m0, w: m1 - m0, y, h, key });
  };
  const { rectItems: r, glyphItems: g, pathItems: p, textItems: t } = rec;
  let inPass = pass < 0;
  for (let o = 0; o < rec.opCount; o++) {
    const [code, start, count, at] = rec.ops.subarray(4 * o, 4 * o + 4);
    if (code === Op.beginPass) inPass = pass < 0 || rec.args[at] === pass;
    if (!inPass) continue;
    const end = start + count;
    if (code === Op.rects)
      for (let i = start; i < end; i++)
        if (r.color[i] !== SkinToken.surface || r.w[i] > Math.max(1, Math.round(dpr)))
          add(0, r.x[i], r.x[i] + r.w[i], r.y[i], r.h[i], `${r.color[i]} ${r.alpha[i]} ${r.stroke[i]}`);
    if (code === Op.glyphs)
      for (let i = start; i < end; i++) add(1, g.x[i], g.x[i], g.y[i], g.radius[i], `${g.glyph[i]} ${g.color[i]}`);
    if (code === Op.paths)
      for (let i = start; i < end; i++) {
        const c1 = i + 1 < p.n ? p.cstart[i + 1] : p.nc;
        let [lo, hi] = [Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY];
        for (let c = p.cstart[i]; c < c1; c += 2) [lo, hi] = [Math.min(lo, p.coords[c]), Math.max(hi, p.coords[c])];
        if (c1 > p.cstart[i]) add(2, lo, hi, c1 - p.cstart[i], p.mode[i], `${p.color[i]}`);
      }
    // Centered text (letters, numbers); a label left of its anchor reads from the left in either direction.
    if (code === Op.text && rec.args[at + 1] === TextAlign.center)
      for (let i = start; i < end; i++) add(3, t.x[i], t.x[i], t.y[i], 0, `${t.text[i]} ${t.color[i]}`);
    // A letter stretched from its left (a logo's): by its center.
    if (code === Op.text && rec.args[at + 1] === TextAlign.left) {
      const adv = advances(rec.args[at]);
      for (let i = start; i < end; i++) {
        if (t.text[i].length !== 1) continue;
        const c = t.x[i] + ((adv[t.text[i].charCodeAt(0) & 127] || 12) * t.sx[i] * dpr) / 2;
        add(3, c, c, t.y[i], t.sy[i], `${t.text[i]} ${t.color[i]}`);
      }
    }
    // An image of pixels: each device column it covers, keyed by the source column's pixels (nearest sample).
    const px = code === Op.image ? rec.images.content[start] : null;
    if (!(px instanceof Uint32Array)) continue;
    const [sx, , sw, sh, dx, dy, dw, dh] = rec.args.subarray(at, at + 8);
    const iw = rec.images.w[start];
    for (let x = Math.ceil(dx); x < dx + dw; x++) {
      const u = Math.min(iw - 1, Math.floor(sx + ((x + 0.5 - dx) / dw) * sw));
      let col = "";
      for (let y = 0; y < sh; y++) col += `${px[y * iw + u].toString(36)} `;
      add(4, x, x + 1, dy, dh, col);
    }
  }
  return out;
}
/**
 * The items of `a` unmatched in `b`: within each group of one kind, key, y and height, taken in x order, pairs more
 * than 2 device px apart in x or width, and those past the other group's count. Text whose counts differ was thinned
 * (labels where they fit, from the left): not compared. An image column matches the same pixels within 2 columns.
 */
export function unmatched(a: readonly MirrorItem[], b: readonly MirrorItem[]): MirrorItem[] {
  const groups = (xs: Item[]) => {
    const m = new Map<string, Item[]>();
    for (const i of xs) {
      const k = `${i.kind}|${i.key}|${i.y}|${i.h}`;
      m.set(k, [...(m.get(k) ?? []), i]);
    }
    for (const g of m.values()) g.sort((p, q) => p.x - q.x || p.w - q.w);
    return m;
  };
  // Image columns by place: a device column across two positions samples either.
  const cols = new Set(b.filter((i) => i.kind === 4).map((i) => `${i.y}|${i.h}|${i.x}|${i.key}`));
  const near = (p: Item) => [0, -1, 1, -2, 2].some((d) => cols.has(`${p.y}|${p.h}|${p.x + d}|${p.key}`));
  const out: Item[] = a.filter((p) => p.kind === 4 && !near(p));
  const gb = groups(b.filter((i) => i.kind !== 4));
  for (const [k, g] of groups(a.filter((i) => i.kind !== 4))) {
    const h = gb.get(k) ?? [];
    if (g[0].kind === 3 && h.length !== g.length) continue;
    g.forEach((p, j) => {
      if (j >= h.length || Math.abs(h[j].x - p.x) > 2 || Math.abs(h[j].w - p.w) > 2) out.push(p);
    });
  }
  return out;
}

/**
 * Drawn rects, sectors, glyphs and (off rings) centered text hit back to their layer, an item that covers the position
 * hit where the layer hits items; without a draw, the layer's entries sampled across the data area.
 */
function hits(ctx: SeqViewerContext, f: LayerFrame, hit: Hit, at: string, itemHits: boolean) {
  const l = ctx.layout.get(f.section)?.current;
  const ax = ctx.layout.axis(f.section);
  if (!l || !ax) return { ok: false, detail: `${at}: no layout` };
  const points: number[] = [];
  const t = ctx.graph.assembled.value?.view?.resolved.layers.track[f.k];
  const r = f.direct.rectItems;
  const g = f.direct.glyphItems;
  const sc = f.direct.sectorItems;
  if (f.entry.draw) {
    for (let i = 0; i < r.n; i++) points.push((r.x[i] + r.w[i] / 2) / ax.dpr, (r.y[i] + r.h[i] / 2) / ax.dpr);
    // On rings, each sector's middle (rects drawn through the ring map).
    for (let i = 0; i < sc.n; i++) {
      const a = (sc.a0[i] + sc.a1[i]) / 2;
      const rr = (sc.r0[i] + sc.r1[i]) / 2;
      points.push((sc.cx[i] + rr * Math.cos(a)) / ax.dpr, (sc.cy[i] + rr * Math.sin(a)) / ax.dpr);
    }
    for (let i = 0; i < g.n; i++) points.push(g.x[i] / ax.dpr, g.y[i] / ax.dpr);
    // Off rings, each centered text's anchor (a letter's, a number's): on the track, an item's where the layer hits
    // items (LC-17). A path's vertices are no hit target (an arc's ends lie on cell corners).
    const r0 = f.direct;
    for (let o = 0; l.geometry !== AxisGeometry.circular && o < r0.opCount; o++) {
      const [code, start, count, at] = r0.ops.subarray(4 * o, 4 * o + 4);
      if (code !== Op.text || r0.args[at + 1] !== TextAlign.center) continue;
      for (let i = start; i < start + count; i++) points.push(r0.textItems.x[i] / ax.dpr, r0.textItems.y[i] / ax.dpr);
    }
  } else {
    for (let e = 0; e < l.entryCount; e++) {
      if (l.entryKind[e] !== EntryKind.track || l.entryTrack[e] !== t) continue;
      const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - ax.scrollY;
      if (top + l.entryHeight[e] / 2 >= l.height || (e >= l.pinnedCount && top < l.scrolledTop)) continue;
      for (let s = 0; s < 32; s++)
        points.push(l.dataLeft + ((s + 0.5) * (l.dataRight - l.dataLeft)) / 32, top + l.entryHeight[e] / 2);
    }
  }
  let bad = 0;
  let items = 0;
  for (let i = 0; i < points.length; i += 2) {
    // A glyph's center past the data area's edge is clipped there: nothing to hover.
    if (points[i] < l.dataLeft || points[i] >= l.dataRight) continue;
    ctx.layout.hitTest(f.section, points[i], points[i + 1], hit);
    const item = hit.kind === HoverKind.item;
    if (item) items++;
    // Drawn: its own layer (an item where it hits items; another layer's item on its track is over it, as a rows
    // layer's lanes); sampled: its track, an item of this layer's in its data.
    if (f.entry.draw ? (hit.layer !== f.k && !(item && hit.track === t)) || (itemHits && !item) : hit.track !== t)
      bad++;
    else if (item && hit.layer === f.k && !inData(ctx, f.k, hit.annotation)) bad++;
    else if (item && hit.layer === f.k && !covers(ctx, hit, ax.axis, ax.cell * ax.span[hit.block] * ax.dpr)) bad++;
  }
  // Drawn, nothing in view is no failure; sampled, the layer's entries must be hit.
  const ok = (points.length > 0 || f.entry.draw !== undefined) && bad === 0;
  return { ok, detail: `${at}: ${points.length / 2} points, ${items} items, ${bad} bad` };
}
