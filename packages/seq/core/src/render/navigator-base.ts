/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { DefaultStoreSettings, type StoreSettings, type Want } from '@molstar/seqx-core/data';
import {
  inversePath,
  type Path,
  pathDir,
  projectEnd,
  projectPoint,
  projectRanges,
  type StaticModel,
} from '@molstar/seqx-core/resolve';
import {
  blockPositionAt,
  blockSpan,
  cellOffset,
  type LayerKinds,
  nearestOffsetOf,
  offsetOf,
  placeOfUnit,
  positionAt,
  type ResolvedNavigator,
  type ResolvedSection,
  type ResolvedView,
  sameBlocks,
  unitOfPlace,
} from '@molstar/seqx-core/view';
import type { AxisTransform, Viewport } from './axis.js';
import type { ThemeContext } from './channels.js';
import { DefaultRenderConfig, type RenderConfig } from './config.js';
import { type DrawLists } from './draw-lists.js';
import type { ImageCache } from './images.js';
import { RenderIssues } from './issues.js';
import { layerSource } from './layer-data.js';
import type { Emphasis, Overlay, StoreReader } from './layer-type.js';
import { BlockPlacement } from './layout/blocks.js';
import type { WindowPart } from './layout/types.js';
import { type NavigatorPicture, navigatorPicture } from './navigator-picture.js';
import { Renderer, type RendererOptions } from './renderer.js';
import type { RowsLay } from './rows.js';
import type { TextMetrics } from './text.js';
import type { SeqTokens } from './tokens.js';
import { type LayoutState, type SectionLayout, type ViewLayout } from './tracks.js';
import { DefaultZoomStrip, type ZoomOverBlocks, floorExtent, NavWindowFlags } from './window.js';

import {
  arcText,
  blockIndexAt,
  foldsOf,
  NavigatorWindow,
  type Picture,
  PictureId,
  pictureConfig,
  pictureOf,
  rowsLayOf,
  sectionIn,
  spanText,
  trackHasRows,
} from './navigator-parts.js';

export { blockPlaceText, type NavigatorHit, NavigatorWindow, navigatorBox, windowText } from './navigator-parts.js';

// A navigator (_spec/README.md "Render", Navigators): one part over a section's axis, its rows or both, drawing its
// system whole (the fit window, in the bound section's display order) through the layout interface, as a view of its
// own over the spec's data (a section of the bound section's tracks small, a chromosome's bands, the rows summarized),
// and over it the bound sections' windows: each section's viewport carried along the path from its system to the
// navigator's (identity, a region into its chromosome, a protein into its gene). Its draw lists are its own: the base
// (`background`, `data`) is keyed by its axis, size, tokens and the store versions it read, so a pan redraws only the
// emphasis passes, the `window` pass first (the windows, their wash, their grips), then the marks (the selection's and
// the highlights' ticks, the keyboard cursor and the pointer's hairline). Never in a section's lists.

/** The gap between the context and the zoom strip, where the connectors run (CSS px, fixed: a look). */
const ZoomGapPx = 6;

/** The zoom strip's section in the picture's view. */
export interface NavigatorRendererOptions {
  readonly model: StaticModel;
  /** The view the navigator belongs to, and the navigator (a node of it, or a host's definition resolved alike). */
  readonly view: ResolvedView;
  readonly navigator: ResolvedNavigator;
  readonly data: StoreReader;
  readonly tokens: SeqTokens;
  readonly text: TextMetrics;
  readonly config?: RenderConfig;
  /** The layer types and section layouts its picture draws with: the view's renderer's (the composition's registry). */
  readonly types: RendererOptions['types'];
  readonly layouts: RendererOptions['layouts'];
  /** The pictures `content` names (the composition's `navigator-picture` entries): one the table lacks draws the axis. */
  readonly pictures: Readonly<Record<string, NavigatorPicture>>;
  /**
   * The layer kinds its picture is resolved with (the composition's, `SeqContext.kinds`): a picture's track of an
   * extension's type is drawn as that type; the built-ins' when absent.
   */
  readonly kinds?: LayerKinds;
  /** The track sizings and the host's port its picture sizes tracks with: the view's renderer's. */
  readonly sizings?: RendererOptions['sizings'];
  readonly trackSizing?: RendererOptions['trackSizing'];
  readonly theme?: Partial<Pick<ThemeContext, 'schemes' | 'measures'>>;
  /** The shared image cache (the view's renderer's): a navigator's image tiles take the view's budget. */
  readonly images?: ImageCache;
  /**
   * The view's renderer's layout state (`Renderer.layoutState`): the picture colors by the view's held domains, never
   * by a domain of its own over the whole axis, so a value has one color in both.
   */
  readonly domains?: LayoutState;
  /** The store's settings (its tile bound: what a layer read at the fit would cost); the defaults when absent. */
  readonly storeSettings?: StoreSettings;
  /**
   * Which way it runs (`navigator.direction`): the bound section's display order (`display`, the default: a drag right
   * moves the view right), or `genomic`, positions up left to right whatever the section's direction.
   */
  readonly direction?: 'display' | 'genomic';
  /**
   * The zoom strip (`navigator.zoomStrip`, `DefaultZoomStrip` when absent): under a reachable span narrower than
   * `minPx`, `heightPx` more.
   */
  readonly zoomStrip?: { readonly minPx: number; readonly heightPx: number; readonly overBlocks?: ZoomOverBlocks };
}

export interface NavigatorFrameInput {
  /** The navigator's box; a zoom strip, when it has one, draws `zoomHeight` more below it. */
  readonly size: { readonly width: number; readonly height: number };
  readonly dpr?: number;
  readonly settled?: boolean;
  readonly overlay?: Overlay;
  /** The view as last laid out: the bound sections' windows are read from it (undefined: none drawn yet). */
  readonly view?: ViewLayout;
  /**
   * The axes each bound section's window moves along, by address (`NavWindowFlags.x | .y`; 0: idle); absent, the
   * windows are an export's outlines (no grips, no idle state).
   */
  readonly moves?: Readonly<Record<string, number>>;
  /** The lit part (hovered or dragged) and its section. */
  readonly lit?: { readonly section: string; readonly part: WindowPart };
  /** The section the pointer last acted on: its window drawn on top and hit first. */
  readonly top?: string;
  /** The pointer's hairline (x, y in the navigator; NaN: none). */
  readonly pointerX?: number;
  readonly pointerY?: number;
  readonly emphasis?: Emphasis;
  readonly cvdSafe?: boolean;
  /** A sweep in progress: the region swept, x on the navigator (drawn as a window outline). */
  readonly sweep?: { readonly x0: number; readonly x1: number };
}

export interface NavigatorFrame {
  readonly lists: DrawLists;
  /** Whether the picture was laid out again (else its base was kept and only the emphasis passes drawn). */
  readonly rebuilt: boolean;
  readonly base: boolean;
  /** What the picture reads, `view` priority only (a navigator's tiles are few and never thrash as margins would). */
  readonly wants: readonly Want[];
  readonly windows: readonly NavigatorWindow[];
  readonly issues: readonly Issue[];
}

/** The passes the navigator draws itself over its picture's. */

/** A window's x over blocks, on the work section and the context (scratch, module-private). */
const BlockXs = new Float64Array(4);
const Ext = new Float64Array(8);

/**
 * Draws a navigator: its picture through a renderer of its own (the layout interface at the fit window), its windows
 * from the view's last layout. Allocates nothing once warm on a window move (the picture's key holds).
 */
export class NavigatorBase {
  readonly navigator: ResolvedNavigator;
  /** The picture it draws (the entry `content` names, else the axis). */
  readonly entry: NavigatorPicture;
  readonly config: RenderConfig;
  protected readonly pictureConfig: RenderConfig;
  protected zoomConfig: RenderConfig | undefined;
  protected readonly zoomStrip: NonNullable<NavigatorRendererOptions['zoomStrip']>;
  /** The zoom strip's track and the gap above it (their sum `heightPx`, never more). */
  protected readonly zoomTrack: number;
  protected readonly zoomGap: number;
  /** The gap the picture's sections are laid out with: the strip at the box's bottom whatever the context's content. */
  protected sectionGap = 0;
  /** The box's height the frame was given (the strip below it). */
  protected boxHeight = 0;
  /** The picture has a zoom strip (made with one); its viewport, the reachable span in the picture's units. */
  protected zoomed = false;
  protected zoomU0 = 0;
  protected zoomU1 = 0;
  /**
   * The reachable span (the last known: kept while a bound section is pending): its units on the picture, `[0]` to
   * `[1]` (past the axis's total over a cyclic origin: the shorter arc), and its ends as positions on the navigator's
   * system (`[2]`, `[3]`; over the origin the arc's start and end).
   */
  protected readonly reach = new Float64Array([0, 0, 0, 0]);
  /** The span's pieces this frame (units on the picture, `[a, b)` pairs; `pieceN` of them) and their positions' hull. */
  protected pieces = new Float64Array(16);
  protected pieceN = 0;
  protected readonly piecePos = new Float64Array(2);
  /** The mark's drawn pieces on the context (a second over a cyclic origin, else NaN): a hover over them names it. */
  protected readonly markX = new Float64Array([Number.NaN, Number.NaN, Number.NaN, Number.NaN]);
  /** The token table last set (`setTokens`): a picture made again draws with it. */
  protected tokens: SeqTokens;
  protected picture: Picture | undefined;
  protected renderer: Renderer | undefined;
  protected left = new Set<string>();
  protected readonly windows: NavigatorWindow[];
  protected readonly issues: Issue[] = [];
  protected layout: ViewLayout | undefined;
  protected readonly viewports: Record<string, Viewport> = {};
  protected viewportsFor: Readonly<Record<string, Viewport>> = {};
  protected rowScale = 0;
  protected wants: readonly Want[] = [];
  protected readonly frameRecord: {
    lists: DrawLists;
    rebuilt: boolean;
    base: boolean;
    wants: readonly Want[];
    windows: readonly NavigatorWindow[];
    issues: readonly Issue[];
  };
  protected readonly frameInput: {
    size: { width: number; height: number };
    viewports?: Readonly<Record<string, Viewport>>;
    dpr: number;
    settled: boolean;
    overlay?: Overlay;
    cvdSafe: boolean;
  } = { size: { width: 0, height: 0 }, dpr: 1, settled: false, cvdSafe: false };
  /** Each bound section's path to the navigator's system and back (undefined: identity). */
  protected readonly paths = new Map<string, { to?: Path; back?: Path }>();
  protected readonly scratch: number[] = [];
  /** A navigator over blocks (its picture a section of blocks: windows by places, never along paths). */
  protected readonly overBlocks: boolean;
  /** Over blocks, per bound section (made again with the picture): equal block lists (identity in units), its reach. */
  protected readonly identical = new Map<string, boolean>();
  protected readonly reaches = new Map<string, readonly [number, number] | undefined>();
  /** Over blocks: the window's width on the context this frame (what turns the zoom strip on), the span's blocks. */
  protected windowPx = Number.POSITIVE_INFINITY;
  protected spanFirst = -1;
  protected spanLast = -1;

  constructor(readonly options: NavigatorRendererOptions) {
    this.navigator = options.navigator;
    this.overBlocks = !!options.navigator.blocks;
    this.entry = navigatorPicture(options.navigator, options.pictures);
    this.config = options.config ?? DefaultRenderConfig;
    this.pictureConfig = pictureConfig(this.config);
    this.zoomStrip = options.zoomStrip ?? DefaultZoomStrip;
    // a small `heightPx` takes its gap first: the strip never lays out taller than it
    const h = Math.max(0, this.zoomStrip.heightPx);
    this.zoomGap = Math.min(ZoomGapPx, Math.floor(h / 2));
    // over blocks the strip is a genome's ruler, its labels in bases above its ticks: at least a ruler's height (the
    // view's, `heights.ruler`), so its labels and body are inside it, never clipped (`heightPx: 0` still off)
    this.zoomTrack =
      h > 0 && this.overBlocks ? Math.max(h - this.zoomGap, this.config.heights.ruler) : h - this.zoomGap;
    this.sectionGap = this.zoomGap;
    this.tokens = options.tokens;
    this.windows = this.navigator.of.map(() => new NavigatorWindow());
    const nav = this.navigator;
    for (const a of nav.missing)
      this.say('seq.render.navigator.of', 'The navigator names a section the view does not have', a);
    if (!nav.of.length) this.say('seq.render.navigator.of', 'The navigator binds no section of the view', nav.address);
    for (const a of nav.of) {
      const p = nav.paths[a];
      if (!p) this.say('seq.render.navigator.path', 'No path joins the section’s system to the navigator’s', a);
      else this.paths.set(a, p.edges.length ? { to: p, back: inversePath(p) } : {});
    }
    const bound = options.view.sections.find((s) => s.address === nav.of[0]);
    if ((nav.axes !== 'positions' || nav.content === 'rows') && bound && !bound.tracks.some(trackHasRows))
      this.say('seq.render.navigator.rows', 'The navigator draws rows and its section has none', nav.address);
    this.frameRecord = {
      lists: undefined as never,
      rebuilt: false,
      base: false,
      wants: [],
      windows: this.windows,
      issues: this.issues,
    };
    this.build();
  }

  protected say(code: keyof typeof RenderIssues.codes, message: string, address: string): void {
    if (this.issues.some((i) => i.code === code && i.address === address)) return;
    this.issues.push(issueOf(RenderIssues, code, message, { address, hint: RenderIssues.codes[code].hint }));
  }

  /** A key of the picture (a layer's open domain) as the view's: its tracks are the bound section's, copied. */
  protected readonly viewKey = (key: string): string => {
    const at = key.indexOf(this.picturePrefix);
    return at < 0 ? key : key.slice(0, at) + this.boundPrefix + key.slice(at + this.picturePrefix.length);
  };
  protected get picturePrefix(): string {
    return `${this.navigator.address}/${PictureId}/`;
  }
  protected get boundPrefix(): string {
    return `${this.navigator.of[0]}/`;
  }

  /** Whether the picture fills the navigator's height (its entry's `fills`: a `bands` track, the rows). */
  protected get fills(): boolean {
    return !!this.entry.fills || this.navigator.axes !== 'positions';
  }

  /** The height the picture was made for (a `bands` track or the rows fill it). */
  protected builtHeight = 0;

  protected build(height = this.builtHeight): void {
    const o = this.options;
    this.builtHeight = height;
    // the zoom strip's track: its height less the connectors' gap
    const zoom = this.zoomed ? this.zoomTrack : 0;
    this.picture = pictureOf(
      o.model,
      o.view,
      this.navigator,
      this.entry,
      o.kinds,
      this.config,
      this.left,
      height,
      o.direction ?? 'display',
      zoom,
    );
    if (this.zoomed && this.zoomConfig?.sectionGap !== this.sectionGap)
      this.zoomConfig = pictureConfig(this.config, this.sectionGap);
    this.renderer = this.picture
      ? new Renderer({
          model: o.model,
          view: this.picture.view,
          data: o.data,
          // the tokens last set: a theme change survives the picture made again
          tokens: this.tokens,
          text: o.text,
          config: this.zoomed ? this.zoomConfig! : this.pictureConfig,
          types: o.types,
          layouts: o.layouts,
          ...(o.sizings ? { sizings: o.sizings } : {}),
          ...(o.trackSizing ? { trackSizing: o.trackSizing } : {}),
          picture: true,
          ...(o.theme ? { theme: o.theme } : {}),
          ...(o.images ? { images: o.images } : {}),
          ...(o.domains ? { domains: { state: o.domains, map: this.viewKey } } : {}),
        })
      : undefined;
    this.layout = undefined;
    this.rowScale = 0;
    this.identical.clear();
    this.reaches.clear();
  }

  /** x of a unit on a section of the picture: a row of blocks by its row's arithmetic (gaps, chips), else its axis. */
  protected xOn(s: SectionLayout, u: number): number {
    const p = s.placement;
    return p instanceof BlockPlacement && !p.stacked ? p.row.xOf(u) : s.axis!.xOfUnit(u);
  }

  /** `xOn` into `out[k]` (a window move's: nothing boxed). */
  /** `xInto` of the unit `src[i]` (read here: a unit passed on would be boxed). */
  protected xIntoFrom(s: SectionLayout, src: Float64Array, i: number, out: Float64Array, k: number): void {
    const p = s.placement;
    if (p instanceof BlockPlacement && !p.stacked) p.row.xInto(src[i], out, k);
    else out[k] = s.axis!.xOfUnit(src[i]);
  }

  protected xInto(s: SectionLayout, u: number, out: Float64Array, k: number): void {
    const p = s.placement;
    if (p instanceof BlockPlacement && !p.stacked) p.row.xInto(u, out, k);
    else out[k] = s.axis!.xOfUnit(u);
  }

  /** The unit at x on a section of the picture (a row of blocks: in a gap the next block's first). */
  protected unitOn(s: SectionLayout, x: number): number {
    const p = s.placement;
    if (p instanceof BlockPlacement && !p.stacked) return p.row.unitAt(x);
    const a = s.axis!;
    return a.u0 + (x - a.x0) / a.scale;
  }

  /** Over blocks: whether a bound section shows the picture's block list (identity in units). */
  protected same(bound: ResolvedSection): boolean {
    let v = this.identical.get(bound.address);
    if (v === undefined) {
      const pic = this.picture?.section;
      v = !!pic && sameBlocks(bound, pic);
      this.identical.set(bound.address, v);
    }
    return v;
  }

  /** Over blocks: a bound section's unit on the picture, through its place (`{ system, position }`). */
  protected toPicture(bound: ResolvedSection, u: number, side: 'start' | 'end'): number | undefined {
    if (this.same(bound)) return u;
    const pic = this.picture?.section;
    const at = pic && placeOfUnit(bound, u, side);
    return pic && at ? unitOfPlace(pic, at) : undefined;
  }

  /** Over blocks: a picture's unit in a bound section's units; undefined where the section shows no block of it. */
  protected fromPicture(bound: ResolvedSection, u: number, side: 'start' | 'end'): number | undefined {
    if (this.same(bound)) return u;
    const pic = this.picture?.section;
    const at = pic && placeOfUnit(pic, u, side);
    if (!at) return undefined;
    if (bound.blocks) return unitOfPlace(bound, at);
    if (at.system !== bound.axis.system) return undefined;
    return offsetOf(bound.axis, at.position) ?? nearestOffsetOf(bound.axis, at.position);
  }

  /** Over blocks: what a bound section reaches on the picture (its blocks' or its axis's places, their hull). */
  protected reachOf(bound: ResolvedSection): readonly [number, number] | undefined {
    if (this.reaches.has(bound.address)) return this.reaches.get(bound.address);
    const pic = this.picture?.section;
    let lo = Number.POSITIVE_INFINITY;
    let hi = Number.NEGATIVE_INFINITY;
    const take = (system: string, r: ArrayLike<number> | undefined) => {
      for (let k = 0; pic && r && k < r.length; k++) {
        const u = unitOfPlace(pic, { system, position: r[k] });
        if (u === undefined) continue;
        lo = Math.min(lo, u);
        hi = Math.max(hi, u);
      }
    };
    if (bound.blocks) for (const b of bound.blocks) take(b.axis.system, b.axis.ranges);
    else take(bound.axis.system, bound.axis.ranges);
    const r = lo < hi ? ([lo, hi] as const) : undefined;
    this.reaches.set(bound.address, r);
    return r;
  }

  /** Sets the token table (a theme change): the picture draws again with it. */
  setTokens(tokens: SeqTokens): void {
    this.tokens = tokens;
    if (this.renderer) this.renderer.tokens = tokens;
  }

  /** The picture's last layout (what hits read). */
  get pictureLayout(): ViewLayout | undefined {
    return this.layout;
  }

  /** The picture's section as last laid out (with a zoom strip, the context). */
  get pictureSection(): SectionLayout | undefined {
    return this.layout?.sections[0];
  }

  /** What the zoom strip adds below the navigator's box (CSS px): 0 without one. */
  get zoomHeight(): number {
    return this.zoomed ? this.zoomTrack + this.zoomGap : 0;
  }

  /** The section the windows are on: the zoom strip, else the picture. */
  protected get work(): SectionLayout | undefined {
    return this.layout?.sections[this.zoomed ? 1 : 0];
  }

  /** What the mark says (a hover over it): the reachable span on the navigator's system; undefined without a zoom strip. */
  markText(): string | undefined {
    if (!this.zoomed || !(this.reach[1] > this.reach[0])) return undefined;
    // over blocks the span is blocks: their names
    if (this.overBlocks) {
      const blocks = this.picture?.section.blocks;
      const name = (b: number) => (blocks?.[b] ? blocks[b].label : '');
      if (this.spanFirst < 0) return undefined;
      return this.spanFirst === this.spanLast
        ? name(this.spanFirst)
        : `${name(this.spanFirst)} – ${name(this.spanLast)}`;
    }
    // over a cyclic origin the arc, from its start on past the origin to its end
    if (this.wraps) return arcText(this.options.model, this.navigator, this.reach[2], this.reach[3]);
    return spanText(this.options.model, this.navigator, this.reach[2], this.reach[3]);
  }

  /**
   * x on the picture (the context, with a zoom strip) of positions `[start, end)` of the navigator's system, in its
   * display order (a reversed axis's, collapsed introns' nearest): `[x0, x1]`, x0 ≤ x1; undefined before the picture is
   * laid out, or where none is shown (`end` past `start`: across a cyclic origin, two calls). For a host's marks.
   */
  xRangeOf(start: number, end: number): readonly [number, number] | undefined {
    const axis = this.pictureSection?.axis;
    if (!axis || !(end > start)) return undefined;
    // each end's cell, its left edge in display order (a collapsed intron's: the nearest shown)
    const a = cellOffset(axis, start) ?? nearestOffsetOf(axis, start);
    const b = cellOffset(axis, end - 1) ?? nearestOffsetOf(axis, end - 1);
    if (a === undefined || b === undefined) return undefined;
    const x0 = axis.xOfUnit(Math.min(a, b));
    const x1 = axis.xOfUnit(Math.max(a, b) + 1);
    return [Math.min(x0, x1), Math.max(x0, x1)];
  }

  /** The position of the navigator's system at x on the picture (the context), clamped to it; undefined before a layout. */
  positionAtX(x: number): number | undefined {
    const axis = this.pictureSection?.axis;
    if (!axis) return undefined;
    const p = axis.positionAt(Math.max(axis.x0, Math.min(axis.x0 + axis.width - 1e-9, x)));
    return p >= 0 ? p : undefined;
  }

  /** What the picture says at a position of the navigator's system (its entry's `readout`); undefined: the label's. */
  readout(position: number): string | undefined {
    return this.entry.readout?.(position, { model: this.options.model, navigator: this.navigator });
  }

  /** The context's axis's total (the picture's units), and whether the reachable span runs over its cyclic origin. */
  protected get turn(): number {
    return this.layout?.sections[0]?.axis?.total ?? 0;
  }
  protected get wraps(): boolean {
    const t = this.turn;
    return t > 0 && this.reach[1] > t + 1e-9;
  }

  /**
   * A unit on the zoom strip's axis: one before its start goes on a turn (its viewport past a cyclic origin, its pieces
   * twice), so the strip runs on over the origin.
   */
  protected onStrip(u: number): number {
    return this.zoomed && this.zoomU1 > this.turn + 1e-9 && u < this.zoomU0 - 1e-9 ? u + this.turn : u;
  }

  /**
   * The reachable span from this frame's pieces (none: a bound section pending, the last kept): their hull, or on a
   * cyclic axis the shorter arc that covers them (its end past the total over the origin); a bound section on the
   * navigator's own axis is the whole axis.
   */
  protected readReach(): void {
    const n = this.pieceN;
    const ctx = this.layout?.sections[0]?.axis;
    if (!n || !ctx) return;
    const p = this.pieces;
    const t = ctx.total;
    const reach = this.reach;
    if (!ctx.cyclic) {
      let lo = Number.POSITIVE_INFINITY;
      let hi = Number.NEGATIVE_INFINITY;
      for (let i = 0; i < n; i++) {
        lo = Math.min(lo, p[2 * i]);
        hi = Math.max(hi, p[2 * i + 1]);
      }
      reach[0] = lo;
      reach[1] = hi;
      reach[2] = this.piecePos[0];
      reach[3] = this.piecePos[1];
      return;
    }
    // the shorter arc: from each piece's start, the farthest end met going on (a piece across that start: no arc)
    let best = t;
    let start = 0;
    for (let i = 0; i < n; i++) {
      const s = p[2 * i];
      let len = 0;
      for (let j = 0; j < n && len <= t; j++) {
        let a = p[2 * j] - s;
        if (a < -1e-9) a += t;
        const b = a + (p[2 * j + 1] - p[2 * j]);
        len = b > t + 1e-9 ? Number.POSITIVE_INFINITY : Math.max(len, b);
      }
      if (len < best - 1e-9) {
        best = len;
        start = s;
      }
    }
    reach[0] = start;
    reach[1] = start + best;
    if (reach[1] > t + 1e-9) {
      // over the origin: the arc's ends as positions, its start then its end
      reach[2] = positionAt(ctx, start, 'start') ?? this.piecePos[0];
      reach[3] = positionAt(ctx, reach[1] - t, 'end') ?? this.piecePos[1];
    } else {
      reach[2] = this.piecePos[0];
      reach[3] = this.piecePos[1];
    }
  }

  /** Adds a piece of the reachable span (units on the picture) and its ends as positions. */
  protected addPiece(ua: number, ub: number, pa: number, pb: number): void {
    if (2 * this.pieceN + 2 > this.pieces.length) {
      const grown = new Float64Array(2 * this.pieces.length);
      grown.set(this.pieces);
      this.pieces = grown;
    }
    this.pieces[2 * this.pieceN] = Math.min(ua, ub);
    this.pieces[2 * this.pieceN + 1] = Math.max(ua, ub);
    this.pieceN++;
    this.piecePos[0] = Math.min(this.piecePos[0], pa, pb);
    this.piecePos[1] = Math.max(this.piecePos[1], pa, pb);
  }

  /** Fits the rows to their band; true when the scale moved (laid out again). */
  protected fitRows(layout: ViewLayout, view: ViewLayout | undefined): boolean {
    const lay = rowsLayOf(layout.sections[0]);
    if (!lay || lay.total === undefined) return false;
    const count = Math.max(1, lay.display.count);
    const scale = lay.box.height / count;
    if (!(scale > 0) || Math.abs(scale - this.rowScale) < 1e-9) return false;
    this.rowScale = scale;
    this.setViewports(view);
    return true;
  }

  /**
   * What no summary serves is left out, said once (`seq.render.navigator.rows`, `.layer`): a lazy alignment's rows (its
   * every row would be read), an annotation's ranges on a source with no coarse level whose fit read would pass the
   * store's bound on tiles. True when the picture was made again without them.
   */
  protected leaveOut(layout: ViewLayout): boolean {
    const s = layout.sections[0];
    if (!s?.axis) return false;
    const left = new Set(this.left);
    const settings = this.options.storeSettings ?? DefaultStoreSettings;
    for (const t of s.tracks) {
      const viewTrack = `${this.navigator.of[0]}${t.track.address.slice(t.track.address.lastIndexOf('/'))}`;
      for (const l of t.layers) {
        // the rows layer by its registered name, never its class (a composition with no rows bundles none)
        const lay = l.type.name === 'rows' ? (l.lay as RowsLay) : undefined;
        if (lay) {
          if (lay.total !== undefined && !lay.held && lay.total * lay.columns > settings.maxBlockValues) {
            left.add(viewTrack);
            this.say('seq.render.navigator.rows', 'The alignment is not summarized: its rows are left out', l.key);
          }
          continue;
        }
        if (l.layer.shape !== 'ranges') continue;
        const src = layerSource(this.options.model, l.layer);
        if (!src || !this.options.model.sources.has(src.source)) continue;
        const levels = this.options.data.coverage(src.source)?.levels ?? [1];
        const coarse = levels[levels.length - 1] ?? 1;
        if (s.axis.total / (settings.tileBins * coarse) > settings.maxTilesPerWant) {
          left.add(viewTrack);
          this.say('seq.render.navigator.layer', 'A layer with no coarse level is left out of the navigator', l.key);
        }
      }
    }
    if (left.size === this.left.size) return false;
    this.left = left;
    this.build();
    return true;
  }

  /** Each bound section's window from the view's last layout, carried to the navigator's frame. */
  protected placeWindows(view: ViewLayout | undefined): void {
    const nav = this.navigator;
    const pic = this.work;
    const axis = pic?.axis;
    this.pieceN = 0;
    this.piecePos[0] = Number.POSITIVE_INFINITY;
    this.piecePos[1] = Number.NEGATIVE_INFINITY;
    for (let i = 0; i < nav.of.length; i++) {
      const w = this.windows[i];
      w.section = nav.of[i];
      w.shown = false;
      w.axes = 0;
      if (!pic || !view) continue;
      const s = sectionIn(view, nav.of[i]);
      const box = pic.placement.box;
      // across the picture's data box (the zoom strip's, below the context, with one); rows: their band
      w.lo = box.x;
      w.hi = box.x + box.width;
      w.ly = this.zoomed ? pic.y : 0;
      w.hy = this.zoomed ? pic.y + pic.height : this.frameInput.size.height;
      w.x0 = w.lo;
      w.x1 = w.hi;
      w.y0 = w.ly;
      w.y1 = w.hy;
      w.wrap = false;
      if (nav.axes !== 'rows') {
        if (!axis || !s?.axis || !this.placeX(w, s, axis)) continue;
        w.axes |= NavWindowFlags.x;
      }
      if (nav.axes !== 'positions') {
        const vr = s && rowsLayOf(s);
        const nr = rowsLayOf(pic);
        if (!vr || !nr || vr.total === undefined || nr.total === undefined) {
          if (nav.axes === 'rows') continue;
        } else {
          const rows = vr.display.count - vr.display.pinned;
          const vis = Math.max(0, vr.box.height - vr.pinnedHeight) / vr.rowScale;
          const top = nr.box.y + nr.pinnedHeight;
          w.ly = top;
          w.hy = top + Math.max(0, nr.display.count - nr.display.pinned) * nr.rowScale;
          w.y0 = top + vr.top * nr.rowScale;
          w.y1 = top + Math.min(rows, vr.top + vis) * nr.rowScale;
          w.r0 = vr.top;
          w.r1 = Math.min(rows, vr.top + vis);
          w.rows = rows;
          w.axes |= NavWindowFlags.y;
        }
      }
      const c = this.config.navigator.window;
      // floored within the picture, not the reachable span: a region's window on its chromosome is seen, the hairline
      // at its truth
      Ext[0] = w.x0;
      Ext[1] = w.x1;
      Ext[2] = box.x;
      Ext[3] = box.x + box.width;
      Ext[4] = w.y0;
      Ext[5] = w.y1;
      Ext[6] = w.ly;
      Ext[7] = w.hy;
      floorExtent(Ext, 0, c);
      floorExtent(Ext, 4, c);
      w.dx0 = Ext[0];
      w.dx1 = Ext[1];
      w.dy0 = Ext[4];
      w.dy1 = Ext[5];
      w.shown = true;
    }
  }

  /** The window along x: the section's view as positions, carried along the path, as x on the picture. */
  protected placeX(w: NavigatorWindow, s: SectionLayout, axis: AxisTransform): boolean {
    if (this.overBlocks) return this.placeBlocksX(w, s);
    const va = s.axis!;
    const path = this.paths.get(w.section);
    if (!path) return false;
    let p0 = positionAt(va, va.u0, 'start');
    let p1 = positionAt(va, va.u1, 'end');
    if (p0 === undefined || p1 === undefined) return false;
    if (path.to) {
      p0 = projectPoint(path.to, p0);
      p1 = projectEnd(path.to, p1);
      if (p0 === undefined || p1 === undefined) return false;
    }
    w.p0 = p0;
    w.p1 = p1;
    w.rev = s.section.axis.dir * (path.to ? pathDir(path.to) : 1) < 0;
    const u0 = offsetOf(axis, p0) ?? nearestOffsetOf(axis, p0);
    const u1 = offsetOf(axis, p1) ?? nearestOffsetOf(axis, p1);
    if (u0 === undefined || u1 === undefined) return false;
    let x0 = axis.xOfUnit(this.onStrip(u0));
    let x1 = axis.xOfUnit(this.onStrip(u1));
    // the section's display order carried along the path against the picture's (a reversed section on a genomic
    // picture, a protein on a minus-strand gene): its start is the picture's right end
    if (s.section.axis.dir * (path.to ? pathDir(path.to) : 1) !== this.picture!.section.axis.dir) [x0, x1] = [x1, x0];
    // a window over a cyclic origin: two pieces, from its start to the end and from the start to its end
    if (x1 < x0 && axis.cyclic) {
      w.wrap = true;
      w.w0 = axis.x0;
      w.w1 = x1;
      x1 = axis.x0 + axis.width;
    } else if (x1 < x0) [x0, x1] = [x1, x0];
    w.x0 = x0;
    w.x1 = x1;
    // what the section can reach: its own extent here (a region within its chromosome), the rest washed
    if (path.to || !this.navigator.sectionAxis) {
      // its own pieces (the layout's, once a source gave the length; a follower's own, never the image it shows)
      const own = s.section.axis.ranges ?? (s.section.follow ? [] : va.pieces);
      const r = path.to ? projectRanges(path.to, own, this.scratch) : own;
      let lo = Number.POSITIVE_INFINITY;
      let hi = Number.NEGATIVE_INFINITY;
      for (let k = 0; k + 1 < r.length; k += 2) {
        const ua = nearestOffsetOf(axis, r[k]) ?? 0;
        const ub = nearestOffsetOf(axis, r[k + 1]) ?? 0;
        // on a strip over a cyclic origin, a piece past the origin goes on a turn (its start the turn's when the
        // piece starts at the origin itself)
        const sa = this.onStrip(ua);
        const a = axis.xOfUnit(sa);
        const b = axis.xOfUnit(sa + (ub - ua));
        lo = Math.min(lo, a, b);
        hi = Math.max(hi, a, b);
        // the span's pieces (the zoom strip's viewport, the mark's text)
        this.addPiece(ua, ub, r[k], r[k + 1]);
      }
      if (lo < hi) {
        w.lo = Math.max(w.lo, lo);
        w.hi = Math.min(w.hi, hi);
      }
    } else {
      // on the navigator's own axis: the section reaches the whole of it
      const t = this.turn;
      if (t > 0) this.addPiece(0, t, 0, t);
    }
    return true;
  }

  /**
   * Over blocks (_spec/README.md "Render", Navigators): the bound section's window in the picture's units (identity
   * where the block lists are equal, else by its ends' places), its ends' systems and positions (what its text says),
   * the wash outside what the section reaches, and the zoom strip's span: the blocks the window touches, their hull.
   */
  protected placeBlocksX(w: NavigatorWindow, s: SectionLayout): boolean {
    const work = this.work;
    const ctx = this.layout?.sections[0];
    const pic = this.picture?.section;
    const blocks = pic?.blocks;
    const va = s.axis;
    if (!work || !ctx || !pic || !blocks?.length || !va) return false;
    const bound = s.section;
    const a = this.toPicture(bound, va.u0, 'start');
    const b = this.toPicture(bound, va.u1, 'end');
    if (a === undefined || b === undefined) return false;
    const i = blockIndexAt(blocks, a, 'start');
    const j = blockIndexAt(blocks, b, 'end');
    w.s0 = blocks[i].axis.system;
    w.p0 = blockPositionAt(blocks[i], a, 'start') ?? 0;
    w.s1 = blocks[j].axis.system;
    w.p1 = blockPositionAt(blocks[j], b, 'end') ?? 0;
    // through a typed array: an x returned from the row's arithmetic would be boxed on each window move
    const X = BlockXs;
    this.xInto(work, a, X, 0);
    this.xInto(work, b, X, 1);
    this.xInto(ctx, a, X, 2);
    this.xInto(ctx, b, X, 3);
    w.x0 = X[0] < X[1] ? X[0] : X[1];
    w.x1 = X[0] < X[1] ? X[1] : X[0];
    this.windowPx = X[3] > X[2] ? X[3] - X[2] : X[2] - X[3];
    // what the section reaches: a section of other blocks (one chromosome) only its own, the rest washed
    if (!this.same(bound)) {
      const r = this.reachOf(bound);
      if (r) {
        this.xInto(work, r[0], X, 0);
        this.xInto(work, r[1], X, 1);
        w.lo = Math.max(w.lo, Math.min(X[0], X[1]));
        w.hi = Math.min(w.hi, Math.max(X[0], X[1]));
      }
    }
    const first = Math.min(i, j);
    const last = Math.max(i, j);
    this.spanFirst = first;
    this.spanLast = last;
    this.addPiece(blocks[first].offset ?? 0, (blocks[last].offset ?? 0) + blockSpan(blocks[last]), 0, 0);
    return true;
  }

  /**
   * The picture's viewports (one record while they hold): the whole axis, its rows band fitted, the view's folds; a
   * zoom strip's the reachable span (`zoom`: it moved).
   */
  protected setViewports(view: ViewLayout | undefined, zoom = false): void {
    const p = this.picture;
    const zoomAt = this.zoomed ? this.picture?.view.sections[1]?.address : undefined;
    if (zoomAt && (zoom || !this.viewports[zoomAt])) {
      this.viewports[zoomAt] = { start: this.zoomU0, end: this.zoomU1 };
      this.viewportsFor = { ...this.viewports };
      this.frameInput.viewports = this.viewportsFor;
    }
    if (!p?.rows || !this.rowScale) {
      if (!zoomAt && this.frameInput.viewports) delete this.frameInput.viewports;
      return;
    }
    const folds = foldsOf(view, this.navigator.of[0]);
    const v = this.viewports[p.section.address];
    const band = v?.bands?.[p.rows.picture];
    if (band && band.rowScale === this.rowScale && band.folds === folds) return;
    this.viewports[p.section.address] = {
      start: 0,
      end: Number.MAX_SAFE_INTEGER,
      bands: { [p.rows.picture]: { bandTop: 0, rowScale: this.rowScale, ...(folds ? { folds } : {}) } },
    };
    this.viewportsFor = { ...this.viewports };
    this.frameInput.viewports = this.viewportsFor;
  }
}
