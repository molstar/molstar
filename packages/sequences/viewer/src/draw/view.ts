/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The engine's view of one layer in one pass: the section's headless axis view, extended with the
// pass's backend, the skin, the text manager, the layer's channel themes and what the overlay reads. One per surface
// and pass, preallocated; the composer (CS1) fills it per entry, the conformance kit (testing) in small.

import type { ResolvedView } from "@molstar/seq-core/format";
import type { ComposedPath, Numbering } from "@molstar/seq-core/model";
import { AxisGeometry, SectionAxisView, type ShownIndex } from "../axis.js";
import type { Skin } from "../graph.js";
import { DefaultRenderProps, type RenderProps } from "../layout.js";
import { MarkSet, type OverlayMarks } from "../marks.js";
import { runXInto } from "../runs.js";
import type { PendingSelection } from "../selection.js";
import type { TextManagerApi } from "../text.js";
import { type Theme, ThemeLocation } from "../themes.js";
import type { RenderBackend } from "./backend.js";

/**
 * The circular geometry's view of an entry (seqspec:viewer#the-render-backend): its ring and angles; `Ring` (ring.ts)
 * is it.
 */
export interface RingView {
  /** The center and the angle of layout offset 0 (−π/2: top), device px and radians. */
  cx: number;
  cy: number;
  theta0: number;
  /** The entry's ring, device px. */
  readonly rInner: number;
  readonly rOuter: number;
  /** Its angular columns: its circumference in device px. */
  columns: number;
  /**
   * The unrolled frame, CSS px (ring.ts): the outer radius, the radius per px of depth, the frame's top, its x period
   * (the circumference): a draw reads them to place a point at a radius (a chord's control point).
   */
  readonly r: number;
  readonly scale: number;
  readonly top: number;
  readonly w: number;
  /** Radians, clockwise from `theta0`; block gaps as angular gaps. */
  angleOf(block: number, shown: ShownIndex): number;
  /** [block, shown] under an angle; false: a gap. */
  shownAt(angle: number, out: Int32Array): boolean;
}

/** No runs in any block: what an overlay drawn without marks reads (VL-22: `null` there threw in `litPairs`). */
const NO_RUNS = new MarkSet();
const NO_MARKS: OverlayMarks = {
  bitsAt: () => 0,
  any: () => false,
  runs: () => NO_RUNS,
  primary: () => -1,
  item: () => -1,
  itemAnnotation: () => -1,
  element: () => false,
};

/** `seqv:map`'s host, the composer's (draw/miniature.ts): the section's tracks drawn in miniature into an image's layer. */
export interface TracksHost {
  /** Into `out[at..at + 3]` what the miniature depends on besides the map's geometry: layout, regions, fonts, model and skin. */
  stamp(out: Float64Array, at: number): void;
  /** Tracks `refs` (none: all but `miniature: false` ones) over `v`'s map, `w` × `h` device px, into `into`; false: none. */
  draw(into: RenderBackend, refs: readonly string[], v: DrawView, w: number, h: number): boolean;
}

let frames = 0;
/** A new composer frame's id (`DrawView.frameId`): never one another frame had (wraps after 2^30). */
export function nextFrameId(): number {
  frames = (frames + 1) & 0x3fffffff;
  return frames;
}

export class DrawView extends SectionAxisView {
  /** The pass's target. */
  backend: RenderBackend = null as unknown as RenderBackend;
  skin: Skin = null as unknown as Skin;
  props: RenderProps = DefaultRenderProps;
  text: TextManagerApi = null as unknown as TextManagerApi;
  /** The layer drawn: its index in the resolved view (`ThemeLocation.at` reads it). */
  layer = -1;
  /**
   * Per block, the whole view's shown range [viewLo, viewHi): `visLo`/`visHi` narrowed to a pan's uncovered strip
   * leaves these as they are, so a decision over the view (items or density, marks or none) is the same in a strip as
   * in a full frame, as a `stable` renderer's must be (VL-07).
   */
  viewLo: Int32Array = new Int32Array(0);
  viewHi: Int32Array = new Int32Array(0);
  /** The layer's channel themes in its type's channel order; the first color and glyph channels. */
  themes: readonly (Theme | null)[] = [];
  color: Theme<"color"> | null = null;
  glyph: Theme<"glyph"> | null = null;
  /** Reused by every draw. */
  readonly loc: ThemeLocation = ThemeLocation.create();
  /** On rings, the map of this unrolled draw onto the circle (`geometry` is then linear); else null. */
  ring: RingView | null = null;
  /** The renderer's cache for this prepared state and skin. */
  cache: unknown = null;
  miniature = false;
  /** The tracks in miniature (`seqv:map`); null: none. */
  tracks: TracksHost | null = null;
  /** A rows layer's data row in the entry drawn; -1: not a rows entry. */
  dataRow = -1;
  columnX = 0;
  columnWidth = 0;
  /**
   * The section's rows layer's lanes in view (`lanesInView`, the stream manager's): lane offsets [k0, k1), display rows
   * [d0, d1], display row 0's surface y and a display row's height (CSS px); k1 <= k0: none. `blockLevels` reads it.
   */
  readonly lanes = new Float64Array(6);
  /** What the channels light: the overlay and navigator passes. */
  marks: OverlayMarks = NO_MARKS;
  /** The selection being dragged (the overlay's view): a dot plot draws its rectangle; null: none. */
  pending: PendingSelection | null = null;
  /**
   * Set by a draw that left work for another frame (a dot plot's tiles past its budget): the engine draws that layer's
   * entries again next frame. Reset by the engine.
   */
  more = false;
  /** The composer frame drawing (unique across sections and passes; `nextFrameId`): a draw's per-frame work keys on it. */
  frameId = 0;
  /**
   * Per layer and stacked block, the density mode last drawn (`densityMode`; 1: density): the composer's overlay view
   * shares its main view's, so an overlay reads the main pass's decision (`densityShown`; seqspec:interaction#hover).
   */
  densityModes: Map<number, number> = new Map();
  /**
   * Per block of the section, its numbering and the path from its coordinate system to the block's, read live from
   * the published view (`numberingFrom`): a numbering switch reruns no prepare.
   */
  numbering: readonly (Numbering | null)[] = [];
  numberingPath: readonly (ComposedPath | null)[] = [];
  private numberedBy: ResolvedView | null = null;
  private numberedAt = -1;
  private readonly xs2 = new Float64Array(1);

  /** The section's axis view under the current layout and viewport, copied (arrays shared). Allocation-free. */
  from(a: SectionAxisView): this {
    this.axis = a.axis;
    this.geometry = a.geometry;
    this.first = a.first;
    this.dir = a.dir;
    this.blockX = a.blockX;
    this.blockX0 = a.blockX0;
    this.scrollDev = a.scrollDev;
    this.span = a.span;
    this.positions = a.positions;
    this.visLo = this.viewLo = a.visLo;
    this.visHi = this.viewHi = a.visHi;
    this.cell = a.cell;
    this.lod = a.lod;
    this.reductions = a.reductions;
    this.dpr = a.dpr;
    this.top = a.top;
    this.height = a.height;
    this.left = a.left;
    this.right = a.right;
    this.gutterLeft = a.gutterLeft;
    this.gutterRight = a.gutterRight;
    this.regionTop = a.regionTop;
    this.regionBottom = a.regionBottom;
    this.scrollY = a.scrollY;
    this.navX0 = a.navX0;
    this.navX1 = a.navX1;
    this.navWidth = a.navWidth;
    this.navX0s = a.navX0s;
    this.navX1s = a.navX1s;
    this.navCell = a.navCell;
    this.navGap = a.navGap;
    this.vnavWidth = a.vnavWidth;
    this.runs = a.runs;
    this.gap = a.gap;
    this.fontPx = a.fontPx;
    this.laneHeight = a.laneHeight;
    this.trackGap = a.trackGap;
    this.onlyBlock = this.row = this.lane = this.line = -1;
    return this;
  }

  /** Section s's numberings from a published view; refilled only when the view or section changed (edge-rate). */
  numberingFrom(view: ResolvedView, s: number): this {
    if (view === this.numberedBy && s === this.numberedAt) return this;
    return this.renumber(view, s);
  }
  /** Apart from `numberingFrom`: a closure there would make its context on every call (V8 makes it on entry). */
  private renumber(view: ResolvedView, s: number): this {
    const b0 = view.sections.blockStart[s];
    this.numberedBy = view;
    this.numberedAt = s;
    this.numbering = view.blocks.numbering.slice(b0, b0 + view.sections.blockCount[s]);
    this.numberingPath = this.numbering.map((_, b) => view.paths[view.blocks.numberingPath[b0 + b]] ?? null);
    return this;
  }

  /**
   * The left edge of shown k of block b in device px, an int32 (no double crosses a call, in any tier: `| 0`): linear,
   * `x` times the DPR exactly; wrapped, on the view's line. For k within a cell of the view (`visLo` − 1 to `visHi` +
   * 1): far off it the px pass 2^31 and wrap (VL-12).
   */
  xd(b: number, k: number): number {
    if (this.geometry !== AxisGeometry.wrapped)
      return (Math.round((this.blockX0[b] + k * this.span[b] * this.cell) * this.dpr) - this.scrollDev) | 0;
    runXInto(this.runs[b], this.lineStart, k, this.gap, this.xs2, 0);
    return Math.round((this.left + this.xs2[0]) * this.dpr) | 0;
  }
  /** The right edge of shown k (device px): the next edge, or wrapped its advance (a group gap may follow). */
  xe(b: number, k: number): number {
    if (this.geometry !== AxisGeometry.wrapped) return this.xd(b, k + 1);
    runXInto(this.runs[b], this.lineStart, k, this.gap, this.xs2, 0);
    return Math.round((this.left + this.xs2[0] + this.runs[b].advance) * this.dpr) | 0;
  }

  /** The layer's themes: the first color and glyph channels from its type's channel kinds. */
  setThemes(themes: readonly (Theme | null)[]): void {
    this.themes = themes;
    this.color = this.glyph = null;
    for (let i = 0; i < themes.length; i++) {
      const t = themes[i];
      if (t?.kind === "color" && !this.color) this.color = t as Theme<"color">;
      if (t?.kind === "glyph" && !this.glyph) this.glyph = t as Theme<"glyph">;
    }
  }
}
