/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf, type Node, type Value } from '@molstar/core/repr/tree';
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
  LayerTypes,
  nearestOffsetOf,
  offsetOf,
  positionAt,
  type ResolvedNavigator,
  type ResolvedSection,
  type ResolvedView,
  resolveViewNode,
} from '@molstar/seqx-core/view';
import type { AxisTransform, Viewport } from './axis.js';
import type { ThemeContext } from './channels.js';
import { DefaultRenderConfig, type RenderConfig, renderConfig } from './config.js';
import type { DrawLists, PassName, PassWriter } from './draw-lists.js';
import type { ImageCache } from './images.js';
import { RenderIssues } from './issues.js';
import { layerSource } from './layer-data.js';
import type { Emphasis, Overlay, StoreReader } from './layer-type.js';
import type { WindowPart } from './layout/types.js';
import { Renderer } from './renderer.js';
import { RowsLay } from './rows.js';
import type { TextMetrics } from './text.js';
import type { SeqTokens } from './tokens.js';
import { hitView, type LayoutState, type SectionLayout, sectionLayoutName, type ViewLayout } from './tracks.js';
import { drawNavigatorWindow, floorExtent, NavWindow, NavWindowFlags, navigatorPart } from './window.js';

// A navigator (_spec/README.md "Render", Navigators): one part over a section's axis, its rows or both, drawing its
// system whole (the fit window, in the bound section's display order) through the layout interface, as a view of its
// own over the spec's data (a section of the bound section's tracks small, a chromosome's bands, the rows summarized),
// and over it the bound sections' windows: each section's viewport carried along the path from its system to the
// navigator's (identity, a region into its chromosome, a protein into its gene). Its draw lists are its own: the base
// (`background`, `data`) is keyed by its axis, size, tokens and the store versions it read, so a pan redraws only the
// emphasis passes, the `window` pass first (the windows, their wash, their grips), then the marks (the selection's and
// the highlights' ticks, the keyboard cursor and the pointer's hairline). Never in a section's lists.

/** A navigator's default thickness by its form (CSS px), and a miniature's box. */
export const NavigatorSizes = { axis: 20, tracks: 48, bands: 24, rows: 96, both: [160, 96] } as const;

/** A navigator's box for a view of `width` × `height` (its `size`, else its form's, along or across the view). */
export function navigatorBox(
  nav: ResolvedNavigator,
  width: number,
  height: number,
): { readonly width: number; readonly height: number } {
  const size = nav.size;
  if (Array.isArray(size)) return { width: size[0], height: size[1] };
  if (size === undefined && nav.axes === 'both')
    return { width: NavigatorSizes.both[0], height: NavigatorSizes.both[1] };
  const t =
    (size as number | undefined) ??
    (nav.axes === 'rows'
      ? NavigatorSizes.rows
      : nav.content === 'axis'
        ? NavigatorSizes.axis
        : NavigatorSizes[nav.content]);
  return nav.place === 'left' || nav.place === 'right' ? { width: t, height } : { width, height: t };
}

/** What a navigator's picture is: its view, the section that draws it, the bound section's rows track there. */
interface Picture {
  readonly view: ResolvedView;
  readonly section: ResolvedSection;
  /** The bound section's first rows track and its copy here (their addresses), when the picture draws rows. */
  readonly rows?: { readonly view: string; readonly picture: string };
}

const PictureId = 'picture';

function nodeAt(root: Node, address: string): Node | undefined {
  let n: Node | undefined = root;
  for (const s of address.split('/').filter(Boolean)) n = n?.children?.find((c) => c.id === s);
  return n;
}

const isRuler = (t: Node) => (t.children ?? []).every((l) => l.kind !== 'layer' || l.params?.type === 'ruler');
const hasRows = (t: Node) => (t.children ?? []).some((l) => l.kind === 'layer' && l.params?.type === 'rows');

/**
 * A track copied small: its height scaled (`navigator.trackScale`, at least `trackMinPx`), no label, a rows layer
 * with no header (its tree and numbers are the view's) filling what is left.
 */
function smallTrack(t: Node, config: RenderConfig): Node {
  const h = t.params?.height;
  const params: Record<string, Value> = { ...t.params, label: '' };
  if (typeof h === 'number') params.height = Math.max(config.navigator.trackMinPx, h * config.navigator.trackScale);
  const children = (t.children ?? []).map((l) =>
    l.kind === 'layer' && l.params?.type === 'rows'
      ? { ...l, params: { ...l.params, tree_width: 0, numbers: 'none', height: 'fill' } }
      : l,
  );
  return { ...t, params, children };
}

/** A small track's height: its own, else its first layer type's default, both at the picture's scale. */
function smallHeight(t: Node, config: RenderConfig): number {
  const h = t.params?.height;
  if (typeof h === 'number') return h;
  const type = t.children?.find((l) => l.kind === 'layer')?.params?.type;
  const heights = config.heights as Readonly<Record<string, number>>;
  const base = (typeof type === 'string' ? heights[type] : undefined) ?? config.heights.ranges;
  return Math.max(config.navigator.trackMinPx, base * config.navigator.trackScale);
}

/**
 * The navigator's picture (_spec/README.md "Render", Navigators): a view of one section over its system whole, in the
 * bound section's display order, holding what its `content` draws; `left` the tracks it leaves out (what no summary
 * serves). Resolved like any view, so its layers are the spec's layers, read and colored as the view's.
 */
function pictureOf(
  model: StaticModel,
  view: ResolvedView,
  nav: ResolvedNavigator,
  config: RenderConfig,
  left: ReadonlySet<string>,
  height: number,
  direction: 'display' | 'genomic',
): Picture | undefined {
  const bound = view.sections.find((s) => s.address === nav.of[0]);
  if (!bound || nav.system === undefined) return undefined;
  const node = nodeAt(model.spec, bound.address);
  const same = nav.sectionAxis && nav.system === bound.axis.system;
  const dir = direction === 'genomic' ? 1 : bound.axis.dir;
  const params: Record<string, Value> = { of: nav.system, direction: dir === -1 ? 'reversed' : 'genomic' };
  const blocks = same ? (node?.children ?? []).filter((c) => c.kind === 'block') : [];
  const tracks: Node[] = [];
  let rows: Picture['rows'];
  const content = nav.axes === 'rows' ? 'rows' : nav.content;
  if (content === 'bands' && nav.bands)
    tracks.push({
      kind: 'track',
      id: 'bands',
      // the bands fill the navigator
      params: { label: '', height: Math.max(1, height) },
      children: [{ kind: 'layer', id: 'bands', params: { type: 'bands', data: nav.bands } }],
    });
  else if (content === 'tracks' || content === 'rows')
    for (const t of node?.children ?? []) {
      if (t.kind !== 'track' || t.id === undefined || t.params?.hidden === true) continue;
      const address = `${bound.address}/${t.id}`;
      if (left.has(address)) continue;
      const wanted = content === 'rows' ? hasRows(t) && !rows : nav.tracks ? nav.tracks.includes(address) : !isRuler(t);
      if (!wanted) continue;
      if (hasRows(t) && !rows) rows = { view: address, picture: `${nav.address}/${PictureId}/${t.id}` };
      tracks.push(smallTrack(t, config));
    }
  // the rows fill what the other tracks leave of the navigator's height (an overview's whole column)
  const rowsAt = tracks.findIndex(hasRows);
  if (rowsAt >= 0 && height > 0) {
    let rest = height;
    for (let i = 0; i < tracks.length; i++) if (i !== rowsAt) rest -= smallHeight(tracks[i], config) + 1;
    const t = tracks[rowsAt];
    tracks[rowsAt] = { ...t, params: { ...t.params, height: Math.max(config.navigator.trackMinPx, rest) } };
  }
  const viewNode: Node = {
    kind: 'view',
    id: nav.id,
    params: { layout: 'tracks' },
    children: [{ kind: 'section', id: PictureId, params, children: [...blocks, ...tracks] }],
  };
  const resolved = resolveViewNode(model, viewNode, nav.address);
  let section = resolved.sections[0];
  if (!section) return undefined;
  // the same system: the bound section's own axis (its pieces, collapsed introns and direction), never a guess
  if (same && !bound.blocks)
    section = { ...section, axis: dir === bound.axis.dir ? bound.axis : genomicAxis(bound.axis) };
  return { view: { ...resolved, sections: [section], navigators: [] }, section, ...(rows ? { rows } : {}) };
}

/** An axis read up left to right (`navigator.direction: genomic`): a reversed axis's pieces in genomic order. */
function genomicAxis(axis: ResolvedSection['axis']): ResolvedSection['axis'] {
  const r = axis.ranges;
  if (!r) return { ...axis, dir: 1 };
  const out: number[] = [];
  for (let k = r.length - 2; k >= 0; k -= 2) out.push(r[k], r[k + 1]);
  return { ...axis, dir: 1, ranges: out };
}

/** The config a navigator's picture lays out with: no inset, label column or direction mark; its tracks and lanes small. */
function pictureConfig(base: RenderConfig): RenderConfig {
  const s = base.navigator.trackScale;
  const min = base.navigator.trackMinPx;
  const heights = {} as Record<keyof RenderConfig['heights'], number>;
  for (const k of Object.keys(base.heights) as (keyof RenderConfig['heights'])[])
    heights[k] = Math.max(min, base.heights[k] * s);
  return renderConfig({
    ...base,
    inset: 0,
    labelWidth: 0,
    labelGap: 0,
    trackGap: 1,
    sectionGap: 0,
    axis: { directionMark: 'none' },
    heights,
    // a lane at the tracks' scale: a small track keeps its lanes (its items drawn, not only their density)
    laneHeight: Math.max(min, base.laneHeight * s),
    laneGap: base.laneGap * s,
  });
}

/** A bound section's window in the navigator's frame (reused frame to frame). */
export class NavigatorWindow {
  section = '';
  /** Drawn: the section is placed and the navigator reaches it. */
  shown = false;
  /** The true ends along x (a second piece over a cyclic origin: `wrap`), along y; the drawn ends after the floor. */
  x0 = 0;
  x1 = 0;
  wrap = false;
  w0 = 0;
  w1 = 0;
  y0 = 0;
  y1 = 0;
  dx0 = 0;
  dx1 = 0;
  dy0 = 0;
  dy1 = 0;
  /** What the window can reach along x and y (the bound section's extent here; the rows' band). */
  lo = 0;
  hi = 0;
  ly = 0;
  hy = 0;
  /** The axes it is drawn over (`NavWindowFlags.x`, `.y`). */
  axes = 0;
  /** Its ends as positions of the navigator's system, and as display rows (what its title and slider say). */
  p0 = 0;
  p1 = 0;
  r0 = 0;
  r1 = 0;
  rows = 0;
}

/**
 * The sections that draw a ruler rail by default (_spec/README.md "Render", Navigators): every linear section (tracks)
 * that no shown navigator binds along its positions: the view's declared ones and a host's (`hosted`) replace it; a
 * hidden one brings it back.
 */
export function railSections(view: ResolvedView, hosted: readonly ResolvedNavigator[] = []): string[] {
  const bound = new Set<string>();
  for (const list of [view.navigators, hosted])
    for (const n of list) if (!n.hidden && n.axes !== 'rows') for (const a of n.of) bound.add(a);
  return view.sections
    .filter((s) => sectionLayoutName(s, view) === 'tracks' && !bound.has(s.address))
    .map((s) => s.address);
}

/** What a window says (an export's title, a slider's value): its positions on the navigator's system, its rows. */
export function windowText(model: StaticModel, nav: ResolvedNavigator, w: NavigatorWindow): string {
  const name = (nav.system && model.systems.get(nav.system)?.name) ?? '';
  const f = (n: number) => Math.round(n).toLocaleString('en-US');
  const parts: string[] = [];
  if (w.axes & NavWindowFlags.x)
    parts.push(`${name}:${f(Math.min(w.p0, w.p1) + 1)}–${f(Math.max(w.p0, w.p1))} of ${name}`);
  if (w.axes & NavWindowFlags.y) parts.push(`rows ${f(w.r0 + 1)}–${f(w.r1)} of ${f(w.rows)}`);
  return parts.join(', ');
}

/** A point on a navigator: the window part under it, its section and axis, the position, row and band there. */
export interface NavigatorHit {
  /** The bound section whose window part is under the point; the first bound section outside every window. */
  readonly section: string;
  readonly part: WindowPart;
  /** The axis the part moves (`positions` or `rows`; the body both, where the navigator covers both). */
  readonly axis: 'positions' | 'rows';
  /** The position at the navigator's scale on its system, the display row of the bound rows, a band's name. */
  readonly position?: number;
  readonly row?: number;
  readonly band?: string;
}

export interface NavigatorRendererOptions {
  readonly model: StaticModel;
  /** The view the navigator belongs to, and the navigator (a node of it, or a host's definition resolved alike). */
  readonly view: ResolvedView;
  readonly navigator: ResolvedNavigator;
  readonly data: StoreReader;
  readonly tokens: SeqTokens;
  readonly text: TextMetrics;
  readonly config?: RenderConfig;
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
}

export interface NavigatorFrameInput {
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
const OwnPasses: readonly PassName[] = ['window', 'highlight', 'selection', 'hover', 'cursor'];

const Ext = new Float64Array(8);
const FlagNoStart: number = NavWindowFlags.noStart;
const FlagNoEnd: number = NavWindowFlags.noEnd;
const Seen: number[] = [];

/**
 * Draws a navigator: its picture through a renderer of its own (the layout interface at the fit window), its windows
 * from the view's last layout. Allocates nothing once warm on a window move (the picture's key holds).
 */
export class NavigatorRenderer {
  readonly navigator: ResolvedNavigator;
  readonly config: RenderConfig;
  private readonly pictureConfig: RenderConfig;
  private picture: Picture | undefined;
  private renderer: Renderer | undefined;
  private left = new Set<string>();
  private readonly windows: NavigatorWindow[];
  private readonly issues: Issue[] = [];
  private layout: ViewLayout | undefined;
  private readonly viewports: Record<string, Viewport> = {};
  private viewportsFor: Readonly<Record<string, Viewport>> = {};
  private rowScale = 0;
  private wants: readonly Want[] = [];
  private readonly frameRecord: {
    lists: DrawLists;
    rebuilt: boolean;
    base: boolean;
    wants: readonly Want[];
    windows: readonly NavigatorWindow[];
    issues: readonly Issue[];
  };
  private readonly frameInput: {
    size: { width: number; height: number };
    viewports?: Readonly<Record<string, Viewport>>;
    dpr: number;
    settled: boolean;
    overlay?: Overlay;
    cvdSafe: boolean;
  } = { size: { width: 0, height: 0 }, dpr: 1, settled: false, cvdSafe: false };
  /** Each bound section's path to the navigator's system and back (undefined: identity). */
  private readonly paths = new Map<string, { to?: Path; back?: Path }>();
  private readonly scratch: number[] = [];

  constructor(readonly options: NavigatorRendererOptions) {
    this.navigator = options.navigator;
    this.config = options.config ?? DefaultRenderConfig;
    this.pictureConfig = pictureConfig(this.config);
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

  private say(code: keyof typeof RenderIssues.codes, message: string, address: string): void {
    if (this.issues.some((i) => i.code === code && i.address === address)) return;
    this.issues.push(issueOf(RenderIssues, code, message, { address, hint: RenderIssues.codes[code].hint }));
  }

  /** A key of the picture (a layer's open domain) as the view's: its tracks are the bound section's, copied. */
  private readonly viewKey = (key: string): string => {
    const at = key.indexOf(this.picturePrefix);
    return at < 0 ? key : key.slice(0, at) + this.boundPrefix + key.slice(at + this.picturePrefix.length);
  };
  private get picturePrefix(): string {
    return `${this.navigator.address}/${PictureId}/`;
  }
  private get boundPrefix(): string {
    return `${this.navigator.of[0]}/`;
  }

  /** Whether the picture fills the navigator's height (a `bands` track, the rows). */
  private get fills(): boolean {
    return (
      this.navigator.content === 'bands' || this.navigator.content === 'rows' || this.navigator.axes !== 'positions'
    );
  }

  /** The height the picture was made for (a `bands` track or the rows fill it). */
  private builtHeight = 0;

  private build(height = this.builtHeight): void {
    const o = this.options;
    this.builtHeight = height;
    this.picture = pictureOf(o.model, o.view, this.navigator, this.config, this.left, height, o.direction ?? 'display');
    this.renderer = this.picture
      ? new Renderer({
          model: o.model,
          view: this.picture.view,
          data: o.data,
          tokens: o.tokens,
          text: o.text,
          config: this.pictureConfig,
          ...(o.theme ? { theme: o.theme } : {}),
          ...(o.images ? { images: o.images } : {}),
          ...(o.domains ? { domains: { state: o.domains, map: this.viewKey } } : {}),
        })
      : undefined;
    this.layout = undefined;
    this.rowScale = 0;
  }

  /** Sets the token table (a theme change): the picture draws again with it. */
  setTokens(tokens: SeqTokens): void {
    if (this.renderer) this.renderer.tokens = tokens;
  }

  /** The picture's last layout (what hits read). */
  get pictureLayout(): ViewLayout | undefined {
    return this.layout;
  }

  /** The picture's section as last laid out. */
  get pictureSection(): SectionLayout | undefined {
    return this.layout?.sections[0];
  }

  /** Lays the picture out (or keeps it), then draws the windows and marks over it. Valid until the next call. */
  frame(input: NavigatorFrameInput): NavigatorFrame {
    const f = this.frameRecord;
    // the bands and the rows fill the navigator's height: the picture is made again for another
    if (this.fills && input.size.height !== this.builtHeight) this.build(input.size.height);
    const r = this.renderer;
    if (!r) {
      f.rebuilt = false;
      f.base = false;
      f.wants = [];
      return f;
    }
    const fi = this.frameInput;
    fi.size.width = input.size.width;
    fi.size.height = input.size.height;
    fi.dpr = input.dpr ?? 1;
    fi.settled = input.settled ?? false;
    fi.cvdSafe = input.cvdSafe ?? false;
    if (input.overlay) fi.overlay = input.overlay;
    else delete fi.overlay;
    this.setViewports(input.view);
    let frame = r.frame(fi);
    // the rows fit the band: their px a row follows the count the picture placed (a second layout, on a change only)
    if (frame.rebuilt && this.fitRows(frame.layout, input.view)) frame = r.frame(fi);
    if (frame.rebuilt && this.leaveOut(frame.layout)) return this.frame(input);
    this.layout = frame.layout;
    if (frame.rebuilt) {
      this.wants = frame.wants.filter((w) => w.priority !== 'margin');
      for (const i of frame.layout.issues) if (!this.issues.includes(i)) this.issues.push(i);
    }
    this.placeWindows(input.view);
    this.drawOver(frame.lists, input);
    f.lists = frame.lists;
    f.rebuilt = frame.rebuilt;
    f.base = frame.base;
    f.wants = this.wants;
    return f;
  }

  /** The picture's viewport: the whole axis, its rows band fitted, the view's folds (one record while they hold). */
  private setViewports(view: ViewLayout | undefined): void {
    const p = this.picture;
    if (!p?.rows || !this.rowScale) {
      if (this.frameInput.viewports) delete this.frameInput.viewports;
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

  /** Fits the rows to their band; true when the scale moved (laid out again). */
  private fitRows(layout: ViewLayout, view: ViewLayout | undefined): boolean {
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
  private leaveOut(layout: ViewLayout): boolean {
    const s = layout.sections[0];
    if (!s?.axis) return false;
    const left = new Set(this.left);
    const settings = this.options.storeSettings ?? DefaultStoreSettings;
    for (const t of s.tracks) {
      const viewTrack = `${this.navigator.of[0]}${t.track.address.slice(t.track.address.lastIndexOf('/'))}`;
      for (const l of t.layers) {
        if (l.lay instanceof RowsLay) {
          if (l.lay.total !== undefined && !l.lay.held && l.lay.total * l.lay.columns > settings.maxBlockValues) {
            left.add(viewTrack);
            this.say('seq.render.navigator.rows', 'The alignment is not summarized: its rows are left out', l.key);
          }
          continue;
        }
        if (LayerTypes[l.layer.type]?.shape !== 'ranges') continue;
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
  private placeWindows(view: ViewLayout | undefined): void {
    const nav = this.navigator;
    const pic = this.layout?.sections[0];
    const axis = pic?.axis;
    for (let i = 0; i < nav.of.length; i++) {
      const w = this.windows[i];
      w.section = nav.of[i];
      w.shown = false;
      w.axes = 0;
      if (!pic || !view) continue;
      const s = sectionIn(view, nav.of[i]);
      const box = pic.placement.box;
      // across the picture's data box; rows: their band
      w.lo = box.x;
      w.hi = box.x + box.width;
      w.ly = 0;
      w.hy = this.frameInput.size.height;
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
      Ext[0] = w.x0;
      Ext[1] = w.x1;
      Ext[2] = w.lo;
      Ext[3] = w.hi;
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
  private placeX(w: NavigatorWindow, s: SectionLayout, axis: AxisTransform): boolean {
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
    const u0 = offsetOf(axis, p0) ?? nearestOffsetOf(axis, p0);
    const u1 = offsetOf(axis, p1) ?? nearestOffsetOf(axis, p1);
    if (u0 === undefined || u1 === undefined) return false;
    let x0 = axis.xOfUnit(u0);
    let x1 = axis.xOfUnit(u1);
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
        const a = axis.xOfUnit(nearestOffsetOf(axis, r[k]) ?? 0);
        const b = axis.xOfUnit(nearestOffsetOf(axis, r[k + 1]) ?? 0);
        lo = Math.min(lo, a, b);
        hi = Math.max(hi, a, b);
      }
      if (lo < hi) {
        w.lo = Math.max(w.lo, lo);
        w.hi = Math.min(w.hi, hi);
      }
    }
    return true;
  }

  /** The windows (the top one last), then the marks, in the navigator's own passes. */
  private drawOver(lists: DrawLists, input: NavigatorFrameInput): void {
    lists.reset(OwnPasses);
    const out = lists.writer('window');
    const tokens = this.renderer!.tokens;
    const c = this.config.navigator.window;
    const dpr = input.dpr ?? 1;
    const moves = input.moves;
    // the reachable span's wash: fixed, where a window cannot go
    for (let i = 0; i < this.windows.length; i++) {
      const w = this.windows[i];
      if (!w.shown || !(w.axes & NavWindowFlags.x)) continue;
      const box = this.layout!.sections[0].placement.box;
      const wash = tokens.colors['seq-wash'];
      if (w.lo > box.x) {
        const g = out.args;
        g[0] = box.x;
        g[1] = w.ly;
        g[2] = w.lo - box.x;
        g[3] = w.hy - w.ly;
        g[4] = wash;
        g[5] = 0;
        g[6] = 0;
        g[7] = 0;
        g[8] = 0;
        out.rectArgs();
      }
      if (w.hi < box.x + box.width) {
        const g = out.args;
        g[0] = w.hi;
        g[1] = w.ly;
        g[2] = box.x + box.width - w.hi;
        g[3] = w.hy - w.ly;
        g[4] = wash;
        g[5] = 0;
        g[6] = 0;
        g[7] = 0;
        g[8] = 0;
        out.rectArgs();
      }
    }
    for (let pass = 0; pass < 2; pass++)
      for (let i = 0; i < this.windows.length; i++) {
        const w = this.windows[i];
        const isTop = w.section === input.top;
        if (!w.shown || (pass === 0) === isTop) continue;
        const m = moves ? (moves[w.section] ?? 0) : -1;
        const flags = m < 0 ? NavWindowFlags.outline : m === 0 ? NavWindowFlags.idle : m & w.axes;
        const lit = input.lit?.section === w.section ? input.lit.part : undefined;
        this.drawWindow(out, w, flags, m > 0 ? lit : undefined, dpr, c, tokens);
      }
    // a sweep: the region it would zoom to, an outline over the picture
    const sw = input.sweep;
    if (sw) {
      const W = NavWindow;
      W[0] = W[4] = Math.min(sw.x0, sw.x1);
      W[1] = W[5] = Math.max(sw.x0, sw.x1);
      W[2] = W[6] = W[10] = 0;
      W[3] = W[7] = W[11] = this.frameInput.size.height;
      W[8] = W[0];
      W[9] = W[1];
      drawNavigatorWindow(out, NavWindowFlags.outline, 'body', dpr, c, tokens);
    }
    this.drawMarks(lists, input, tokens, dpr);
  }

  private drawWindow(
    out: PassWriter,
    w: NavigatorWindow,
    flags: number,
    lit: WindowPart | undefined,
    dpr: number,
    c: RenderConfig['navigator']['window'],
    tokens: SeqTokens,
  ): void {
    const W = NavWindow;
    W[2] = w.dy0;
    W[3] = w.dy1;
    W[6] = w.y0;
    W[7] = w.y1;
    W[8] = w.lo;
    W[9] = w.hi;
    W[10] = w.ly;
    W[11] = w.hy;
    if (w.wrap) {
      // over the origin: its two pieces, the wash between them only, a grip at each true end (the first piece's start,
      // the second's end): where a piece meets the axis's end it is body
      W[0] = w.w0;
      W[1] = w.w1;
      W[4] = w.w0;
      W[5] = w.w1;
      W[8] = w.w0;
      W[9] = w.w1;
      drawNavigatorWindow(out, flags | FlagNoStart, lit, dpr, c, tokens);
      W[8] = w.w1;
      W[9] = w.x0;
    }
    W[0] = w.wrap ? w.x0 : w.dx0;
    W[1] = w.wrap ? w.x1 : w.dx1;
    W[4] = w.x0;
    W[5] = w.x1;
    drawNavigatorWindow(out, w.wrap ? flags | FlagNoEnd : flags, lit, dpr, c, tokens);
  }

  /** The selection's and the highlights' ticks at their places, the cursor's hairline, the pointer's hairline. */
  private drawMarks(lists: DrawLists, input: NavigatorFrameInput, tokens: SeqTokens, dpr: number): void {
    const pic = this.layout?.sections[0];
    const axis = pic?.axis;
    if (!pic || !axis) return;
    const marks = this.config.navigator.marks;
    const e = input.emphasis;
    const system = axis.system;
    const box = pic.placement.box;
    const h = Math.min(4, box.height);
    const y = box.y + box.height - h;
    const px = 1 / dpr;
    if (e && marks.includes('selection')) {
      const part = partOn(e.selection, system) ?? partOn(e.projected?.selection, system);
      if (part) ticks(lists.writer('selection'), axis, part.ranges, y, h, px, tokens.colors['msui-accent']);
    }
    if (e?.highlights && marks.includes('highlights')) {
      const out = lists.writer('highlight');
      for (let k = 0; k < e.highlights.length; k++) {
        const part = partOn(e.highlights[k].location, system);
        if (part) ticks(out, axis, part.ranges, y, h, px, e.highlights[k].color);
      }
    }
    if (e?.cursor && e.cursor.system === system && marks.includes('cursor')) {
      const x = axis.x(e.cursor.position);
      if (Number.isFinite(x)) {
        const g = lists.writer('cursor').args;
        g[0] = x - px / 2;
        g[1] = box.y;
        g[2] = px;
        g[3] = box.height;
        g[4] = tokens.colors['msui-text'];
        g[5] = 0;
        g[6] = 0;
        g[7] = 0;
        g[8] = 0;
        lists.writer('cursor').rectArgs();
      }
    }
    const px0 = input.pointerX ?? Number.NaN;
    if (Number.isFinite(px0) && px0 >= box.x && px0 < box.x + box.width)
      lists
        .writer('hover')
        .fill(px0 - px / 2, box.y, px, this.frameInput.size.height - box.y, tokens.colors['seq-hover']);
    const py0 = input.pointerY ?? Number.NaN;
    if (Number.isFinite(py0) && this.navigator.axes !== 'positions') {
      const g = lists.writer('hover').args;
      g[0] = box.x;
      g[1] = py0 - px / 2;
      g[2] = box.width;
      g[3] = px;
      g[4] = tokens.colors['seq-hover'];
      g[5] = 0;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      lists.writer('hover').rectArgs();
    }
  }

  /** The part under a point of the navigator, the position, row and band there; undefined off the picture. */
  hit(x: number, y: number, touch = false, top?: string): NavigatorHit | undefined {
    const pic = this.layout?.sections[0];
    if (!pic) return undefined;
    const c = this.config.navigator.window;
    const nav = this.navigator;
    let best: { w: NavigatorWindow; part: WindowPart; axis: 'positions' | 'rows' } | undefined;
    Seen.length = 0;
    for (let k = -1; k < this.windows.length; k++) {
      // the window the pointer last acted on first (overlapping windows: it takes the press)
      const i = k < 0 ? this.windows.findIndex((w) => w.section === top) : k;
      if (i < 0 || Seen.includes(i)) continue;
      Seen.push(i);
      const w = this.windows[i];
      if (!w.shown) continue;
      const along = w.axes & NavWindowFlags.x ? partOf(x, w, c, touch) : 'body';
      const across = w.axes & NavWindowFlags.y ? navigatorPart(y, w.dy0, w.dy1, w.ly, w.hy, c, touch) : 'body';
      if (along === undefined && across === undefined) continue;
      const a = along ?? 'outside';
      const b = across ?? 'outside';
      let part: WindowPart;
      let axis: 'positions' | 'rows' = nav.axes === 'rows' ? 'rows' : 'positions';
      if ((b === 'start' || b === 'end') && a !== 'outside') {
        part = b;
        axis = 'rows';
      } else if ((a === 'start' || a === 'end') && b !== 'outside') part = a;
      else if (a === 'body' && b === 'body') part = 'body';
      else part = 'outside';
      if (part !== 'outside') {
        best = { w, part, axis };
        break;
      }
      best ??= { w, part, axis };
    }
    const section = best?.w.section ?? nav.of[0];
    if (section === undefined) return undefined;
    const h = this.layout ? hitView(this.layout, x, y) : undefined;
    const axis = pic.axis;
    const p = axis ? axis.positionAt(Math.max(axis.x0, Math.min(axis.x0 + axis.width - 1e-9, x))) : -1;
    const position = p >= 0 ? p : undefined;
    const row = this.rowAt(y);
    return {
      section,
      part: best?.part ?? 'outside',
      axis: best?.axis ?? (nav.axes === 'rows' ? 'rows' : 'positions'),
      ...(position !== undefined && nav.axes !== 'rows' ? { position } : {}),
      ...(row !== undefined ? { row } : {}),
      ...(h?.state !== undefined ? { band: h.state } : {}),
    };
  }

  /** The display row (below the pinned band, fractional) at y on the picture's rows; undefined with none. */
  rowAt(y: number): number | undefined {
    const nr = this.layout ? rowsLayOf(this.layout.sections[0]) : undefined;
    if (!nr || nr.total === undefined || !(nr.rowScale > 0)) return undefined;
    const top = nr.box.y + nr.pinnedHeight;
    const rows = nr.display.count - nr.display.pinned;
    return Math.max(0, Math.min(rows, (y - top) / nr.rowScale));
  }

  /**
   * The bound section's display unit at x (its own axis's offset), through the inverse of the path from its system:
   * what a drag on its window moves it to. Past the picture's ends the drag keeps going at the navigator's scale, as the
   * path carries the edge (a codon a protein position). Undefined while the picture is not placed or the place maps
   * nowhere.
   */
  sectionUnitAt(section: string, x: number): number | undefined {
    const axis = this.layout?.sections[0]?.axis;
    if (!axis) return undefined;
    const u = axis.u0 + (x - axis.x0) / axis.scale;
    const inside = Math.max(0, Math.min(axis.total, u));
    const o = this.unitOfPicture(section, axis, inside);
    if (o === undefined || inside === u) return o;
    const step = u < 0 ? 1 : -1;
    const near = this.unitOfPicture(section, axis, inside + step);
    return near === undefined ? o : o + (u - inside) * ((near - o) / step);
  }

  /** The bound section's display unit at the picture's unit `u` (within the picture). */
  private unitOfPicture(section: string, axis: AxisTransform, u: number): number | undefined {
    const s = this.options.view.sections.find((v) => v.address === section);
    const path = this.paths.get(section);
    if (!s || !path) return undefined;
    // the picture's far end is a window's end edge: carried as one
    const end = u >= axis.total;
    let p = positionAt(axis, u, end ? 'end' : 'start');
    if (p === undefined) return undefined;
    if (path.back) {
      p = end ? projectEnd(path.back, p) : projectPoint(path.back, p);
      if (p === undefined) return undefined;
    }
    return offsetOf(s.axis, p) ?? nearestOffsetOf(s.axis, p);
  }

  /** x on the navigator of a bound section's display unit (through the path from its system): a sweep's ends. */
  xOfSectionUnit(section: string, u: number): number {
    const axis = this.layout?.sections[0]?.axis;
    const s = this.options.view.sections.find((v) => v.address === section);
    const path = this.paths.get(section);
    if (!axis || !s || !path) return Number.NaN;
    let p = positionAt(s.axis, u);
    if (p === undefined) return Number.NaN;
    if (path.to) p = projectPoint(path.to, p);
    if (p === undefined) return Number.NaN;
    const o = offsetOf(axis, p) ?? nearestOffsetOf(axis, p);
    return o === undefined ? Number.NaN : axis.xOfUnit(o);
  }

  /** The issues the navigator met: its bindings, its layers left out, its picture's. */
  get allIssues(): readonly Issue[] {
    return this.issues;
  }
}

/**
 * The part along x of a window that may wrap a cyclic origin (two pieces: from its start to the axis's end, from the
 * axis's start to its end): where a piece meets the axis's end it is body, so its grips are at its true ends only.
 */
function partOf(
  x: number,
  w: NavigatorWindow,
  c: RenderConfig['navigator']['window'],
  touch: boolean,
): WindowPart | undefined {
  if (!w.wrap) return navigatorPart(x, w.dx0, w.dx1, w.lo, w.hi, c, touch);
  const a = navigatorPart(x, w.x0, w.x1, w.lo, w.hi, c, touch);
  if (a === 'end') return 'body';
  if (a !== 'outside') return a;
  const b = navigatorPart(x, w.w0, w.w1, w.lo, w.hi, c, touch);
  return b === 'start' ? 'body' : b;
}

function trackHasRows(t: ResolvedSection['tracks'][number]): boolean {
  return t.layers.some((l) => l.type === 'rows');
}

function sectionIn(view: ViewLayout, address: string): SectionLayout | undefined {
  const s = view.sections;
  for (let i = 0; i < s.length; i++) if (s[i].section.address === address) return s[i];
  return undefined;
}

/** The first rows layer's lay of a section's first strip. */
function rowsLayOf(s: SectionLayout | undefined): RowsLay | undefined {
  if (!s) return undefined;
  const lays = s.lays;
  for (let i = 0; i < lays.length; i++) {
    const l = lays[i];
    for (let j = 0; j < l.length; j++) if (l[j].lay instanceof RowsLay) return l[j].lay as RowsLay;
  }
  return undefined;
}

/** The folds the view's rows band shows (its session's, or its layer's `collapsed`), by identity. */
function foldsOf(view: ViewLayout | undefined, section: string | undefined) {
  const s = view && section ? sectionIn(view, section) : undefined;
  const lay = rowsLayOf(s);
  return lay?.foldsApply ? lay.folds : undefined;
}

function partOn(loc: Emphasis['selection'], system: string) {
  if (!loc) return undefined;
  for (let i = 0; i < loc.parts.length; i++) if (loc.parts[i].system === system) return loc.parts[i];
  return undefined;
}

/** Ticks over position ranges at the picture's foot, each at least a device pixel. */
function ticks(
  out: PassWriter,
  axis: AxisTransform,
  ranges: Int32Array,
  y: number,
  h: number,
  px: number,
  color: number,
): void {
  for (let k = 0; k + 1 < ranges.length; k += 2) {
    const a = axis.x(ranges[k]);
    const b = axis.x(ranges[k + 1]);
    if (!Number.isFinite(a) || !Number.isFinite(b)) continue;
    const x0 = Math.min(a, b);
    const w = Math.max(px, Math.abs(b - a));
    {
      const g = out.args;
      g[0] = x0;
      g[1] = y;
      g[2] = w;
      g[3] = h;
      g[4] = color;
      g[5] = 0;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
    }
  }
}
