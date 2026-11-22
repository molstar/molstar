/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashString } from '@molstar/core/data/util/hash-functions';
import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { RenderIssues } from './issues.js';
import type { Coverage, DataEntry, Want } from '@molstar/seqx-core/data';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedView } from '@molstar/seqx-core/view';
import type { Viewport } from './axis.js';
import type { ThemeContext } from './channels.js';
import { type LayerColors, layerColors, viewTheme } from './layer-color.js';
import { DefaultRenderConfig, type RenderConfig } from './config.js';
import { DrawLists, EmphasisPasses, type PassName } from './draw-lists.js';
import { ImageCache } from './images.js';
import type { RowsZoom } from './rail.js';
import { EmptyOverlay, type Emphasis, type LayerType, type Overlay, type StoreReader } from './layer-type.js';
import type { SectionLayoutType } from './layout/types.js';
import type { TextMetrics } from './text.js';
import type { SeqTokens } from './tokens.js';
import { drawView, type HeldDomains, type LayoutInput, LayoutState, layoutView, type ViewLayout } from './tracks.js';

// Frames keyed by what they read (_spec/README.md "Render", The key): the model's hash, the view, the size, the
// viewport, the device pixel ratio, the token table's version, the text metrics' version, the overlay's version, whether
// the view settled, the config, and the versions of the store entries the last layout read. An equal key reuses the
// layout and the lists and redraws only the emphasis passes (a hover); anything else lays out again into the same
// buffers. The check compares the parts as they are, building nothing: a hover allocates nothing once warm
// ("Nothing allocated once warm"); the key's string is made when the lists are rebuilt.

/** The passes an emphasis redraw draws when the selection's tint is under the data and the selection changed. */
const BaseAndEmphasisPasses: readonly PassName[] = ['background', 'data', ...EmphasisPasses];

export interface RendererOptions {
  readonly model: StaticModel;
  readonly view: ResolvedView;
  readonly data: StoreReader;
  readonly tokens: SeqTokens;
  readonly text: TextMetrics;
  readonly config?: RenderConfig;
  /**
   * The layer types it draws, by name: the composition's registry (_spec/README.md "Context"), never a built-in table
   * (`BuiltInLayerTypes` of `render/built-ins` is the default compositions' list). A layer of no type here is said.
   */
  readonly types: Readonly<Record<string, LayerType<any, any>>>;
  /**
   * The section layouts by name, the composition's registry (`BuiltInLayouts` the default list): the viewer's navigation
   * reads the same (`LayoutSource`). A layout not here is said; its section is laid out as tracks.
   */
  readonly layouts: Readonly<Record<string, SectionLayoutType<any>>>;
  /** The track sizings a section's `heights` may name past the built-ins (the composition's registry), and the host's port. */
  readonly sizings?: LayoutInput['sizings'];
  readonly trackSizing?: LayoutInput['trackSizing'];
  /**
   * A navigator's picture (no label lane).
   * @internal
   */
  readonly picture?: boolean;
  /** The schemes and measures the themes read (the built-ins by default): the composition's registry entries. */
  readonly theme?: Partial<Pick<ThemeContext, 'schemes' | 'measures'>>;
  /** The shared image cache (_spec/README.md "Render", The shared reduction): another renderer's (a navigator's view). */
  readonly images?: ImageCache;
  /**
   * Another renderer's held domains (its `layoutState`, keys mapped): a navigator's picture colors by its view's
   * (_spec/README.md "Render", Navigators). Their version is part of the key.
   */
  readonly domains?: HeldDomains;
}

export interface FrameInput {
  readonly size: { readonly width: number; readonly height: number };
  readonly viewport?: Viewport;
  /** A section's own window by its address, before `viewport` (an export placing sections one by one). */
  readonly viewports?: Readonly<Record<string, Viewport>>;
  readonly dpr?: number;
  readonly settled?: boolean;
  readonly overlay?: Overlay;
  readonly emphasis?: Emphasis;
  /** The reader's color-blind substitution (`theme.cvdSafe`): part of the key, never the spec. */
  readonly cvdSafe?: boolean;
  /** Lays out a text section's every line, its box its content (an export: the figure). */
  readonly figure?: boolean;
  /** Draws the windows a figure leaves out by default (a ring's, in an export: `export.window`). */
  readonly windows?: boolean;
  /**
   * The clades each rows layer folds, by its track's address, where its band gives none (an export's `folds`): part of
   * the key by identity. A renderer keeps them for its life (an export's are one set).
   */
  readonly folds?: Readonly<Record<string, readonly string[]>>;
  /**
   * The sections that draw a ruler rail, by address (_spec/README.md "Render", Navigators): part of the key by
   * identity (a caller keeps the record while it holds). Which of them can move (`railsActive`) is not: it is set on
   * the rails each frame, an emphasis redraw included, so idle wakes with no relayout.
   */
  readonly rails?: Readonly<Record<string, unknown>>;
  readonly railsActive?: Readonly<Record<string, boolean>>;
  /**
   * The host's rows zoom (the viewer's `rows.zoom`): the bounds and the on/off of the rows rails' grips; part of the
   * key by identity (a caller keeps the record while it holds). Absent: `DefaultRowsZoom`.
   */
  readonly rowsZoom?: RowsZoom;
  /**
   * Each header's parts (_spec/README.md "Render", Track headers: `LayoutInput.headers`): part of the key by identity (a
   * caller makes a new one when its entries or their state change).
   */
  readonly headers?: LayoutInput['headers'];
}

/**
 * A frame's lists, layout and key. Valid until the next `frame` call: the lists are the renderer's, reused, and an
 * emphasis redraw returns one kept record, its fields set again (a hover makes no object). Copy what you keep.
 */
export interface Frame {
  readonly lists: DrawLists;
  readonly layout: ViewLayout;
  /** The draw lists' key. */
  readonly key: string;
  /** Whether this frame laid out again (else only the emphasis passes were drawn). */
  readonly rebuilt: boolean;
  /**
   * Whether the base passes (`background`, `data`) were drawn again: on a rebuild, and on an emphasis redraw whose
   * selection changed while its tint is drawn under the data (`emphasis.selectionPlace: 'under'`).
   */
  readonly base: boolean;
  /** What the layout needs from the store. */
  readonly wants: readonly Want[];
}

/** The parts of a frame's key that come from its input and the renderer (the store's reads apart). */
class KeyParts {
  model = '';
  view = '';
  width = Number.NaN;
  height = Number.NaN;
  dpr = Number.NaN;
  viewport = false;
  viewports: Readonly<Record<string, Viewport>> | undefined = undefined;
  start = Number.NaN;
  end = Number.NaN;
  tokens = Number.NaN;
  text = Number.NaN;
  overlay = Number.NaN;
  settled = false;
  cvdSafe = false;
  figure = false;
  windows = false;
  folds: FrameInput['folds'] = undefined;
  rails: FrameInput['rails'] = undefined;
  rowsZoom: FrameInput['rowsZoom'] = undefined;
  headers: FrameInput['headers'] = undefined;
  /** Bumped when `headers` changes: the key's text. */
  headersVersion = 0;
  domains = Number.NaN;

  // the parts are read from their objects here, never passed as numbers (a call boxes the numbers it is passed)

  /** Whether they equal a frame's. */
  equals(r: Renderer, input: FrameInput) {
    const vp = input.viewport;
    return (
      this.model === r.options.model.hash &&
      this.view === r.options.view.address &&
      this.width === input.size.width &&
      this.height === input.size.height &&
      this.dpr === (input.dpr ?? 1) &&
      this.viewport === (vp !== undefined) &&
      (vp === undefined || (this.start === vp.start && this.end === vp.end)) &&
      this.viewports === input.viewports &&
      this.tokens === r.tokens.version &&
      this.text === r.options.text.version &&
      this.overlay === (input.overlay ?? EmptyOverlay).version &&
      this.settled === (input.settled ?? false) &&
      this.cvdSafe === (input.cvdSafe ?? false) &&
      this.figure === (input.figure ?? false) &&
      this.windows === (input.windows ?? false) &&
      this.folds === input.folds &&
      this.rails === input.rails &&
      this.rowsZoom === input.rowsZoom &&
      this.headers === input.headers &&
      this.domains === (r.options.domains?.state.domainsVersion ?? 0)
    );
  }

  set(r: Renderer, input: FrameInput) {
    this.model = r.options.model.hash;
    this.view = r.options.view.address;
    this.width = input.size.width;
    this.height = input.size.height;
    this.dpr = input.dpr ?? 1;
    this.viewport = input.viewport !== undefined;
    this.viewports = input.viewports;
    this.start = input.viewport?.start ?? Number.NaN;
    this.end = input.viewport?.end ?? Number.NaN;
    this.tokens = r.tokens.version;
    this.text = r.options.text.version;
    this.overlay = (input.overlay ?? EmptyOverlay).version;
    this.settled = input.settled ?? false;
    this.cvdSafe = input.cvdSafe ?? false;
    this.folds = input.folds;
    this.figure = input.figure ?? false;
    this.windows = input.windows ?? false;
    this.rails = input.rails;
    this.rowsZoom = input.rowsZoom;
    if (this.headers !== input.headers) {
      this.headers = input.headers;
      this.headersVersion++;
    }
    this.domains = r.options.domains?.state.domainsVersion ?? 0;
  }

  /** The key's text (the store's reads appended by the caller). */
  toString(configKey: string): string {
    let vp = this.viewport ? `${this.start}:${this.end}` : 'all';
    if (this.viewports)
      for (const [a, v] of Object.entries(this.viewports).sort(([x], [y]) => (x < y ? -1 : 1)))
        vp += `|${a}=${v.start}:${v.end}${v.scale !== undefined ? `@${v.scale}` : ''}${v.top ? `^${v.top}` : ''}`;
    return [
      this.model,
      this.view,
      `${this.width}x${this.height}@${this.dpr}`,
      vp,
      this.tokens,
      this.text,
      this.overlay,
      this.settled ? 's' : 'm',
      this.cvdSafe ? 'cvd' : '',
      this.figure ? 'figure' : '',
      this.windows ? 'windows' : '',
      this.rails ? `rails:${Object.keys(this.rails).sort().join(',')}` : '',
      this.rowsZoom ? `rows:${this.rowsZoom.minPx}:${this.rowsZoom.maxPx}:${this.rowsZoom.enabled}` : '',
      this.headers ? `headers:${this.headersVersion}` : '',
      configKey,
    ].join('|');
  }
}

/** Sets which rails can move (outside the key: an idle flip redraws the window pass alone). */
function railsActive(layout: ViewLayout, active: FrameInput['railsActive']): void {
  const s = layout.sections;
  for (let i = 0; i < s.length; i++) {
    const rs = s[i].rails;
    for (let j = 0; rs && j < rs.length; j++) rs[j].active = active?.[rs[j].key] === true;
    // a rows rail is idle while every row is shown, whatever the view can do
    const rows = s[i].rowsRails;
    for (let j = 0; rows && j < rows.length; j++) rows[j].active = active?.[rows[j].key] === true && !rows[j].whole;
  }
}

/**
 * What a polar warp reduced in the frame's data pass (an image's rows on a ring: `seq.render.ring.image`), said
 * once a section among the layout's issues (_spec/README.md "Render", Rings layout).
 */
function refusedImages(layout: ViewLayout): void {
  for (const s of layout.sections) {
    const w = s.placement.warpOf?.(0);
    if (!w?.reduced || layout.issues.some((i) => i.code === 'seq.render.ring.image' && i.address === s.section.address))
      continue;
    (layout.issues as Issue[]).push(
      issueOf(RenderIssues, 'seq.render.ring.image', 'Rows reduced on the ring', {
        address: s.section.address,
        hint: RenderIssues.codes['seq.render.ring.image'].hint,
      }),
    );
  }
}

/** A frame the renderer returns again on each emphasis redraw (a hover makes no object). */
interface ReusedFrame {
  lists: DrawLists;
  layout: ViewLayout;
  key: string;
  rebuilt: boolean;
  base: boolean;
  wants: readonly Want[];
}

export class Renderer {
  tokens: SeqTokens;
  readonly lists: DrawLists;
  private readonly state = new LayoutState();

  /** What the layouts keep between frames (its held domains: a navigator over this view reads them). */
  get layoutState(): LayoutState {
    return this.state;
  }

  /**
   * Its view's layers' color channels as its frames compile them (_spec/README.md "Render", Layer colors): its schemes
   * and measures, the reader's substitution of its last frame, its held domains (another renderer's when it reads
   * those). What a host colors a structure by a track with (`@molstar/seqx-core/structure`'s `colorRuns`).
   */
  layerColors(): LayerColors {
    const d = this.options.domains;
    return layerColors(
      this.options.view,
      viewTheme(
        this.options.view,
        { ...this.options.theme, cvdSafe: this.parts.cvdSafe },
        d ? (key) => d.state.domains.get(d.map(key)) : (key) => this.state.domains.get(key),
      ),
    );
  }
  private readonly config: RenderConfig;
  /** The shared image cache of its frames (every layer's image tiles, one budget). */
  readonly images: ImageCache;
  private readonly configKey: string;
  private layout: ViewLayout | undefined;
  private readonly parts = new KeyParts();
  private readonly reads = new Map<string, number>();
  private readonly coverages = new Map<string, Coverage | undefined>();
  /** The store's refusals when laid out (`StoreReader.refusals`): a change lays the view out again. */
  private refusals: number | undefined;
  /** The last layout's reads as arrays: the check walks them with no iterator. */
  private readKeys: string[] = [];
  private readVersions: number[] = [];
  private coverageKeys: string[] = [];
  private coverageValues: (Coverage | undefined)[] = [];
  private emphasisFrame: ReusedFrame | undefined;
  /** The selection the base passes were last drawn with (its tint under the data). */
  private baseSelection: unknown;
  private baseProjected: unknown;
  private baseHighlights: unknown;
  private readonly reader: StoreReader;

  constructor(readonly options: RendererOptions) {
    this.tokens = options.tokens;
    this.config = options.config ?? DefaultRenderConfig;
    this.images = options.images ?? new ImageCache(this.config.imageBytes);
    this.lists = new DrawLists(this.config.internedStrings);
    this.configKey = (hashString(JSON.stringify(this.config)) >>> 0).toString(16);
    const data = options.data;
    this.reader = {
      ...(data.settings ? { settings: data.settings } : {}),
      coverage: (s) => {
        const c = data.coverage(s);
        this.coverages.set(s, c);
        return c;
      },
      keysOf: (w) => data.keysOf(w),
      ...(data.tooMany ? { tooMany: (w: Want) => data.tooMany!(w) } : {}),
      ...(data.tilesOf ? { tilesOf: (w: Want) => data.tilesOf!(w) } : {}),
      ...(data.layerRefused ? { layerRefused: (l: string) => data.layerRefused!(l) } : {}),
      ...(data.layerLoad ? { layerLoad: (w: readonly Want[]) => data.layerLoad!(w) } : {}),
      get refusals() {
        return data.refusals;
      },
      get: (k): DataEntry | undefined => {
        const e = data.get(k);
        this.reads.set(k, e?.version ?? -1);
        return e;
      },
    };
  }

  /** Whether the entries and coverages the last layout read are as they were. */
  private readsHold(): boolean {
    const data = this.options.data;
    const keys = this.readKeys;
    const versions = this.readVersions;
    for (let i = 0; i < keys.length; i++) if ((data.get(keys[i])?.version ?? -1) !== versions[i]) return false;
    const sources = this.coverageKeys;
    const coverages = this.coverageValues;
    for (let i = 0; i < sources.length; i++) if (data.coverage(sources[i]) !== coverages[i]) return false;
    if (data.refusals !== this.refusals) return false;
    return true;
  }

  /**
   * Lays the view out for the input (or, with an equal key, redraws the emphasis passes alone) into the renderer's
   * lists. The returned frame is valid until the next call: copy what you keep.
   */
  frame(input: FrameInput): Frame {
    const { model, view } = this.options;
    const dpr = input.dpr ?? 1;
    const overlay = input.overlay ?? EmptyOverlay;
    const parts = this.parts;
    if (this.layout && parts.equals(this, input) && this.readsHold()) {
      // the tints under the data are in the background pass: a changed selection or highlight draws the base again
      const e = input.emphasis;
      const base =
        this.config.emphasis.selectionPlace === 'under' &&
        (e?.selection !== this.baseSelection ||
          e?.projected?.selection !== this.baseProjected ||
          e?.highlights !== this.baseHighlights);
      if (base) {
        this.baseSelection = e?.selection;
        this.baseProjected = e?.projected?.selection;
        this.baseHighlights = e?.highlights;
      }
      railsActive(this.layout, input.railsActive);
      drawView(this.layout, this.lists, e, base ? BaseAndEmphasisPasses : EmphasisPasses);
      let f = this.emphasisFrame;
      if (!f)
        f = this.emphasisFrame = {
          lists: this.lists,
          layout: this.layout,
          key: '',
          rebuilt: false,
          base: false,
          wants: [],
        };
      f.base = base;
      f.layout = this.layout;
      f.key = this.lists.key;
      f.wants = this.layout.wants;
      return f;
    }
    this.reads.clear();
    this.coverages.clear();
    const layout = layoutView(
      {
        model,
        view,
        size: input.size,
        ...(input.viewport ? { viewport: input.viewport } : {}),
        ...(input.viewports ? { viewports: input.viewports } : {}),
        dpr,
        tokens: this.tokens,
        text: this.options.text,
        config: this.config,
        data: this.reader,
        overlay,
        settled: input.settled ?? false,
        types: this.options.types,
        layouts: this.options.layouts,
        ...(this.options.sizings ? { sizings: this.options.sizings } : {}),
        ...(this.options.trackSizing ? { trackSizing: this.options.trackSizing } : {}),
        ...(this.options.picture ? { picture: true } : {}),
        theme: { ...this.options.theme, cvdSafe: input.cvdSafe ?? false },
        ...(input.figure ? { figure: true } : {}),
        ...(input.windows ? { windows: true } : {}),
        ...(input.folds ? { folds: input.folds } : {}),
        ...(input.rails ? { rails: input.rails } : {}),
        ...(input.rowsZoom ? { rowsZoom: input.rowsZoom } : {}),
        ...(input.headers ? { headers: input.headers } : {}),
        ...(this.options.domains ? { domains: this.options.domains } : {}),
        images: this.images,
      },
      this.state,
    );
    this.layout = layout;
    parts.set(this, input);
    railsActive(layout, input.railsActive);
    drawView(layout, this.lists, input.emphasis);
    refusedImages(layout);
    this.baseSelection = input.emphasis?.selection;
    this.baseProjected = input.emphasis?.projected?.selection;
    this.baseHighlights = input.emphasis?.highlights;
    const reads = [...this.reads].sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0));
    let sig = '';
    for (const [k, v] of reads) sig += `${k}=${v};`;
    sig += `refusals=${this.options.data.refusals ?? 0};`;
    this.readKeys = reads.map((r) => r[0]);
    this.readVersions = reads.map((r) => r[1]);
    this.coverageKeys = [...this.coverages.keys()];
    this.coverageValues = [...this.coverages.values()];
    this.refusals = this.options.data.refusals;
    this.lists.key = `${parts.toString(this.configKey)}|${(hashString(sig) >>> 0).toString(16)}`;
    return { lists: this.lists, layout, key: this.lists.key, rebuilt: true, base: true, wants: layout.wants };
  }
}
