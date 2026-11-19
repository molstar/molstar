/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashString } from '@molstar/core/data/util/hash-functions';
import { hasErrors, type Issue, issueOf, type Node } from '@molstar/core/repr/tree';
import {
  type AssetReader,
  DataStore,
  type FormatRuntime,
  type OpenedSource,
  RealClock,
  type StoreClock,
  SourceOpenError,
  type Want,
} from '@molstar/seqx-core/data';
import { checkSeqSpec, type SpecBounds } from '@molstar/seqx-core/lang';
import {
  isEmptyLocation,
  type Location,
  locationOfResult,
  resolveAddressIn,
  unionLocations,
} from '@molstar/seqx-core/query';
import {
  OverlayBuilder,
  overlayQueryModel,
  type ResolveConfig,
  resolveSpec,
  type StaticModel,
  sourceLengths,
} from '@molstar/seqx-core/resolve';
import {
  nearestOffsetOf,
  projectOntoView,
  type ResolvedView,
  resolveView,
  resolveViews,
} from '@molstar/seqx-core/view';
import type { Viewport } from './axis.js';
import { type RenderConfigOverrides, renderConfig } from './config.js';
import { accessibleOf, locationText } from './accessibility.js';
import { DrawLists } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { layerSource } from './layer-data.js';
import type { ThemeContext } from './channels.js';
import type { Emphasis } from './layer-type.js';
import { LegendBand, layoutLegends, legendsHeight, viewLegends } from './legends.js';
import { OverlaySources } from './overlay-sources.js';
import {
  type NavigatorFrame,
  NavigatorRenderer,
  type NavigatorRendererOptions,
  navigatorBox,
  railSections,
  windowText,
} from './navigator.js';
import { Renderer } from './renderer.js';
import type { ViewLayout } from './tracks.js';
import { drawListsToSvg, escapeXml, type SvgOptions } from './svg.js';
import { bundledTextMetrics, type TextMetrics } from './text.js';
import { builtInTokens, type SeqTokens } from './tokens.js';

// Headless rendering (_spec/README.md "Render", Export): a spec to an SVG string, asynchronously. It resolves the spec,
// opens the sources the view reads, wants what the export's size and viewport need through a want set of its own
// (never the live view's, EC-125), waits for those entries to settle or a timeout (then the export says it is
// unsettled), and draws once more with the view settled (open domains computed, lanes packed). Legends are drawn in a
// band outside the data box only when asked (`legend`, else the view's `export.legend`, else none), and each shared
// scale's domain as drawn is recorded (`data-seq-scale-<id>`): an open end unioned over what the export draws.

/**
 * The store's sources a view reads: its placed layers' data, its sections' systems' sources (for a length), and its
 * navigators' systems' and bands' (_spec/README.md "Render", Navigators: a chromosome under a region).
 */
export function viewSources(model: StaticModel, view: ResolvedView): Set<string> {
  const out = new Set<string>();
  for (const n of view.navigators) {
    const source = n.system ? model.systems.get(n.system)?.source : undefined;
    if (source) out.add(source);
    const a = n.bands ? model.annotations.get(n.bands)?.access : undefined;
    if (a) out.add(a.kind === 'inline' ? a.source.address : a.source);
  }
  for (const s of view.sections) {
    // the systems it lays out: its axis's, its blocks', a grid's y
    for (const system of [s.axis.system, ...(s.blocks ?? []).map((b) => b.axis.system), s.y?.axis.system]) {
      const source = system === undefined ? undefined : model.systems.get(system)?.source;
      if (source) out.add(source);
    }
    for (const t of s.tracks) {
      if (t.template) out.add(t.template.source);
      for (const l of t.layers) {
        const src = layerSource(model, l);
        if (src) out.add(src.source);
        // a rows layer's joined fields: the sources whose own system is its rows system (Render, Row columns)
        const rows = l.type === 'rows' && l.data ? model.annotations.get(l.data.address)?.rows : undefined;
        if (rows) for (const b of model.sources.values()) if (b.entry.system === rows) out.add(b.address);
      }
    }
  }
  return out;
}

/** An asset reader that reads nothing: a headless render reads no file or URL unless given a reader (security). */
export const NoAssets: AssetReader = {
  readText: (uri) => Promise.reject(new Error(`This render reads no assets: ${uri}`)),
  readBytes: (uri) => Promise.reject(new Error(`This render reads no assets: ${uri}`)),
};

export interface AttachOptions {
  readonly assets?: AssetReader;
  readonly signal?: AbortSignal;
  /** The store's sources to open (all the model names when absent). */
  readonly only?: ReadonlySet<string>;
  /** Collects the addresses this call attached. */
  readonly attached?: Set<string>;
}

/**
 * Opens the model's sources into the store: inline data through its `seq` runtime, a source node through the runtime
 * bound to it. One already attached is kept. What cannot be opened is an issue, and its layers stay pending.
 */
export async function attachSources(
  model: StaticModel,
  store: DataStore,
  options: AttachOptions = {},
): Promise<Issue[]> {
  const ctx = {
    assets: options.assets ?? NoAssets,
    signal: options.signal ?? new AbortController().signal,
    limits: { maxRows: store.settings.maxRows, maxBytes: store.settings.maxBytes },
  };
  const opens = new Map<string, () => Promise<OpenedSource>>();
  for (const plan of [
    ...[...model.systems.values()].map((s) => s.access),
    ...[...model.annotations.values()].map((a) => a.access),
  ])
    if (plan?.kind === 'inline') opens.set(plan.source.address, () => plan.runtime.open(plan.source, ctx));
  for (const b of model.sources.values()) {
    const runtime = b.runtime;
    if (runtime) opens.set(b.address, () => runtime.open(b.entry, ctx));
  }
  const issues: Issue[] = [];
  await Promise.all(
    [...opens].map(async ([address, open]) => {
      if ((options.only && !options.only.has(address)) || store.coverage(address)) return;
      try {
        const opened = await open();
        // a caller that gave up (a timeout) attaches nothing late
        if (ctx.signal.aborted) {
          opened.dispose();
          return;
        }
        // what opening said without refusing the source (rows left out) is reported with the refusals
        if (opened.issues) issues.push(...opened.issues.map((i) => ({ ...i, address })));
        if (store.attach(address, opened)) options.attached?.add(address);
      } catch (e) {
        if (ctx.signal.aborted) return;
        if (e instanceof SourceOpenError) issues.push(...e.issues.map((i) => ({ ...i, address })));
        else
          issues.push(
            issueOf(
              RenderIssues,
              'seq.render.source',
              `The source could not be opened: ${String((e as Error)?.message ?? e)}`,
              {
                address,
                hint: RenderIssues.codes['seq.render.source'].hint,
              },
            ),
          );
      }
    }),
  );
  return issues;
}

export interface OpenModelOptions extends AttachOptions {
  /** The format runtimes the application composed (inline data always reads through `seq`). */
  readonly runtimes?: readonly FormatRuntime[];
  readonly config?: Partial<ResolveConfig>;
  /** The view's address (the first view when absent): what is opened is what it reads. */
  readonly view?: string;
  /** Opens what the view reads (`view`, an export's), or every source of the model (`all`, the default). */
  readonly scope?: 'view' | 'all';
  /** The spec's first resolve, when the caller made it already. */
  readonly model?: StaticModel;
}

export interface OpenedModel {
  /** The model resolved with the lengths the opened sources gave (the first resolve when none did). */
  readonly model: StaticModel;
  /** Its view (undefined when the spec has none of the address). */
  readonly view?: ResolvedView;
  /** The model's, the view's, and what opening the sources said. */
  readonly issues: readonly Issue[];
}

/**
 * A spec's model with its sources open (_spec/README.md "Resolve", A length from coverage): resolved; the sources
 * that hold a system opened first; resolved again with the lengths their coverages gave; then the rest opened over
 * that model. Each source is opened once (one that failed is not tried again). The one implementation of the
 * protocol: an export, a test's harness and a host's one-shot render use it; the viewer keeps its store in step itself.
 */
export async function openModel(spec: Node, store: DataStore, options: OpenModelOptions = {}): Promise<OpenedModel> {
  const resolveOptions = {
    ...(options.runtimes ? { runtimes: options.runtimes } : {}),
    ...(options.config ? { config: options.config } : {}),
  };
  const viewOf = (m: StaticModel) => (options.view ? resolveView(m, options.view) : resolveViews(m)[0]);
  const scoped = (m: StaticModel, v: ResolvedView | undefined) =>
    options.scope === 'view' && v ? viewSources(m, v) : undefined;
  let model = options.model ?? resolveSpec(spec, resolveOptions);
  let view = viewOf(model);
  const shown = scoped(model, view);
  const first = new Set(
    [...model.systems.values()].flatMap((s) => (s.source && (!shown || shown.has(s.source)) ? [s.source] : [])),
  );
  const attach = (m: StaticModel, only: Set<string> | undefined) => {
    const { only: _o, ...rest } = options;
    return attachSources(m, store, { ...rest, ...(only ? { only } : {}) });
  };
  const said = first.size ? await attach(model, first) : [];
  const lengths = sourceLengths(model, (source) => store.coverage(source));
  if (lengths.size) {
    model = resolveSpec(spec, { ...resolveOptions, lengths });
    view = viewOf(model);
  }
  // the rest, over the model with the lengths: each opened once
  let rest = scoped(model, view);
  if (!rest) {
    rest = new Set<string>();
    for (const b of model.sources.values()) rest.add(b.address);
    for (const plan of [
      ...[...model.systems.values()].map((s) => s.access),
      ...[...model.annotations.values()].map((a) => a.access),
    ])
      if (plan?.kind === 'inline') rest.add(plan.source.address);
  }
  for (const s of first) rest.delete(s);
  const opened = await attach(model, rest);
  return { model, ...(view ? { view } : {}), issues: [...model.issues, ...(view?.issues ?? []), ...said, ...opened] };
}

/**
 * A viewport in positions, by section address (`[start, end)` on each section's system): read through each section's
 * axis (its direction and pieces), so a caller needs no knowledge of them.
 */
export interface PositionsViewport {
  readonly positions: Readonly<Record<string, readonly [number, number]>>;
}

export interface RenderSvgOptions {
  /** CSS px; a height `auto` is the sections' content's (a text section every line). */
  readonly size: { readonly width: number; readonly height: number | 'auto' };
  /** Draws a ring's window (the view's `export.window` when absent; false: an exported map is a figure). */
  readonly window?: boolean;
  /** The window of display units of every leading section, or positions per section. */
  readonly viewport?: Viewport | PositionsViewport;
  /** A leading section's own window of display units by its address, before `viewport` (a viewer's sections). */
  readonly viewports?: Readonly<Record<string, Viewport>>;
  /** The format runtimes the application composed (inline data always reads through `seq`). */
  readonly runtimes?: readonly FormatRuntime[];
  /** The bounds the spec's check holds (`DefaultSpecBounds` when absent). */
  readonly bounds?: SpecBounds;
  /** The view's address; the spec's first view when absent. */
  readonly view?: string;
  /** The built-in token table of a theme (`light` by default), or a table. */
  readonly theme?: 'light' | 'dark';
  readonly tokens?: SeqTokens;
  /** The text metrics port: the bundled table by default. */
  readonly text?: TextMetrics;
  readonly config?: RenderConfigOverrides;
  readonly resolve?: Partial<ResolveConfig>;
  /** What sources may read: nothing by default. */
  readonly assets?: AssetReader;
  /** How long the export waits for its sources to open and its data to settle, ms: 10 s by default. */
  readonly timeoutMs?: number;
  /** The time the timeout runs on; default the given store's clock, else the real one (a fake one in tests). */
  readonly clock?: StoreClock;
  /**
   * Prefixes the SVG's ids (clips, masks, the title and description an `aria-labelledby` names): by default `s` and a
   * hash of the spec, the title, the size, the theme, the view, the highlights and the selection, so exports that
   * differ in them can be inlined in one page. Two equal exports share their ids: a host that inlines both passes each
   * its own `idPrefix`.
   */
  readonly idPrefix?: string;
  readonly dpr?: number;
  readonly emphasis?: Emphasis;
  /**
   * A selection as an address, resolved by the export itself once its data settled (its own store's nodes, so a ref
   * to a dynamic node resolves when the export has read it) and drawn with `emphasis.selection`.
   */
  readonly selection?: string;
  /** A store to read through (a viewer's): the export wants through a want set of its own. A new one by default. */
  readonly store?: DataStore;
  /** Legends in a band `below` the data or to its `right`, or `none`: the view's `export.legend` when absent. */
  readonly legend?: 'none' | 'below' | 'right';
  /** A reader's color-blind substitution (off: an export draws the figure's colors unless asked). */
  readonly cvdSafe?: boolean;
  /**
   * The clades each rows layer folds, by its track's address (a viewer passes its own, as it passes the selection);
   * the layer's `collapsed` when absent.
   */
  readonly folds?: Readonly<Record<string, readonly string[]>>;
  /** The schemes and measures the themes read (the built-ins by default). */
  readonly registry?: Partial<Pick<ThemeContext, 'schemes' | 'measures'>>;
  /**
   * The view's navigators (_spec/README.md "Render", Export): `omit`ted (a window is session, not a figure), or
   * `include`d, stacked by `place` around the view, their windows as outlines and wash, and the ruler rails as plain
   * windows (no grips, no idle state). The view's `export.navigators` when absent, else `omit`.
   */
  readonly navigators?: 'omit' | 'include';
  /**
   * How the included navigators draw (a viewer passes its own `navigator.zoomStrip` and `navigator.direction`, so its
   * export is what its parts show): `DefaultZoomStrip` and the display order when absent.
   */
  readonly navigator?: Pick<NavigatorRendererOptions, 'zoomStrip' | 'direction'>;
}

/** The provisional height a layout for a height `auto` is measured at. */
const AutoHeight = 100_000;

/** The sections' content's height: the last one's bottom and the inset. */
function contentHeight(layout: ViewLayout): number {
  const last = layout.sections[layout.sections.length - 1];
  return Math.ceil(last ? last.y + last.height + layout.config.inset : 2 * layout.config.inset);
}

/** A positions viewport as each section's window of display units (an axis waiting for its length: from 0). */
function offsetsOf(view: ResolvedView, v: PositionsViewport): Record<string, Viewport> {
  const out: Record<string, Viewport> = {};
  for (const s of view.sections) {
    const r = v.positions[s.address];
    if (!r) continue;
    const a = s.axis.ranges ? (nearestOffsetOf(s.axis, r[0]) ?? r[0]) : r[0];
    const b = s.axis.ranges ? (nearestOffsetOf(s.axis, r[1]) ?? r[1]) : r[1];
    out[s.address] = { start: Math.min(a, b), end: Math.max(a, b) };
  }
  return out;
}

export interface RenderSvgResult {
  readonly svg: string;
  /** Whether everything the export needed arrived before the timeout. */
  readonly settled: boolean;
  /** The spec's, the view's, the sources' and the render's issues. */
  readonly issues: readonly Issue[];
  readonly lists: DrawLists;
  /** What `selection` named, as drawn; absent when none was given or it named nothing (then its issues are kept). */
  readonly selection?: Location;
}

function wantsKey(wants: readonly Want[]): string {
  return wants.map((w) => `${w.source}|${w.system}|${w.start}|${w.end}|${w.level}|${w.priority}`).join(';');
}

/** A spec drawn to an SVG string once the data its view needs has settled (or the timeout passed). */
export async function renderSvg(spec: Node, options: RenderSvgOptions): Promise<RenderSvgResult> {
  const issues: Issue[] = [...checkSeqSpec(spec, options.bounds)];
  const blank = (settled: boolean): RenderSvgResult => {
    const lists = new DrawLists();
    lists.width = options.size.width;
    lists.height = options.size.height === 'auto' ? 0 : options.size.height;
    return { svg: drawListsToSvg(lists), settled, issues, lists };
  };
  if (hasErrors(issues)) return blank(true);
  const resolveOptions = {
    ...(options.runtimes ? { runtimes: options.runtimes } : {}),
    config: options.resolve ?? {},
  };
  let model = resolveSpec(spec, resolveOptions);
  const firstView = options.view ? resolveView(model, options.view) : resolveViews(model)[0];
  if (!firstView) {
    issues.push(...model.issues);
    issues.push(
      issueOf(RenderIssues, 'seq.render.view', 'The spec has no such view', {
        ...(options.view ? { address: options.view } : {}),
        hint: RenderIssues.codes['seq.render.view'].hint,
      }),
    );
    return blank(true);
  }
  let view: ResolvedView = firstView;
  const own = !options.store;
  const clock = options.clock ?? options.store?.clock ?? RealClock;
  const store = options.store ?? new DataStore({ clock });
  const controller = new AbortController();
  const wanter = store.wanter();
  const timeout = options.timeoutMs ?? 10_000;
  // the deadline covers opening the sources too
  const until = clock.now() + timeout;
  const attached = new Set<string>();
  // tier 2: the store's dynamic nodes, and the sources they need (their inline rows, the nested sources templates
  // reach), attached as they arrive and let go at the end
  const overlays = new OverlayBuilder({
    ...(options.runtimes ? { runtimes: options.runtimes } : {}),
    ...(options.resolve ? { config: options.resolve } : {}),
  });
  const nested = new OverlaySources(store, options.assets ?? NoAssets);
  let settled = true;
  try {
    let timer: unknown;
    const expired = new Promise<undefined>((resolve) => {
      timer = clock.setTimeout(() => resolve(undefined), timeout);
    });
    // the sources opened and the model resolved again with their lengths; taken only when that wins the race (a late
    // one never swaps the model under the render)
    const opened = await Promise.race([
      openModel(spec, store, {
        ...resolveOptions,
        ...(options.assets ? { assets: options.assets } : {}),
        ...(options.view ? { view: options.view } : {}),
        signal: controller.signal,
        attached,
        scope: 'view',
        model,
      }),
      expired,
    ]);
    clock.clearTimeout(timer);
    if (opened?.view) {
      model = opened.model;
      view = opened.view;
      issues.push(...opened.issues);
    } else {
      issues.push(...model.issues, ...view.issues);
      settled = false;
      controller.abort();
    }
    const tokens = options.tokens ?? builtInTokens(options.theme ?? 'light');
    const config = renderConfig(options.config);
    const renderer = new Renderer({
      model,
      view,
      data: store,
      tokens,
      text: options.text ?? bundledTextMetrics(),
      config,
      ...(options.registry ? { theme: options.registry } : {}),
    });
    const given = options.viewport;
    const auto = options.size.height === 'auto';
    const windows = options.window ?? view.export?.window === true;
    // the navigators an export includes (shown ones), drawn as outlines; the rails plain (a figure's)
    const included = (options.navigators ?? view.export?.navigators ?? 'omit') === 'include';
    const navs = included
      ? view.navigators
          .filter((n) => !n.hidden)
          .map(
            (navigator) =>
              new NavigatorRenderer({
                model,
                view,
                navigator,
                data: store,
                tokens,
                text: options.text ?? bundledTextMetrics(),
                config,
                ...(options.registry ? { theme: options.registry } : {}),
                images: renderer.images,
                domains: renderer.layoutState,
                storeSettings: store.settings,
                ...(options.navigator?.zoomStrip ? { zoomStrip: options.navigator.zoomStrip } : {}),
                ...(options.navigator?.direction ? { direction: options.navigator.direction } : {}),
              }),
          )
      : [];
    const rails = included ? Object.fromEntries(railSections(view).map((a) => [a, true])) : undefined;
    const frameInput = {
      // a height `auto` is measured from a first layout (`fitHeight`)
      size: { width: options.size.width, height: auto ? AutoHeight : (options.size.height as number) },
      // a text section exports every line (the figure), never only those in view
      figure: true,
      ...(windows ? { windows: true } : {}),
      ...(given && 'positions' in given
        ? { viewports: { ...offsetsOf(view, given), ...options.viewports } }
        : {
            ...(given ? { viewport: given } : {}),
            ...(options.viewports ? { viewports: options.viewports } : {}),
          }),
      dpr: options.dpr ?? 1,
      ...(options.emphasis ? { emphasis: options.emphasis } : {}),
      ...(options.cvdSafe ? { cvdSafe: true } : {}),
      ...(options.folds ? { folds: options.folds } : {}),
      ...(rails ? { rails } : {}),
    };
    /** Each navigator's frame over the view's layout, at its box around the view. */
    const navFrames = (layout: ViewLayout, settled: boolean) =>
      navs.map((n) =>
        n.frame({
          size: navigatorBox(n.navigator, layout.width, layout.height),
          dpr: options.dpr ?? 1,
          settled,
          overlay: overlays.of(model, store),
          view: layout,
          ...(options.cvdSafe ? { cvdSafe: true } : {}),
        }),
      );
    let last = '';
    // what a layout wants can follow from what arrived (a length, a follower's range): ask until it holds; the
    // export waits for its view alone, never its margins
    for (let round = 0; settled && round < 8; round++) {
      const overlay = overlays.of(model, store);
      await nested.sync(overlay, view);
      let f = renderer.frame({ ...frameInput, overlay: overlays.of(model, store), settled: false });
      if (auto) {
        frameInput.size = { width: options.size.width, height: contentHeight(f.layout) };
        f = renderer.frame({ ...frameInput, overlay: overlays.of(model, store), settled: false });
      }
      // what the layout shows stays while the export waits
      nested.hold(f.layout.shown);
      const wants = [
        ...f.wants.filter((w) => w.priority !== 'margin'),
        ...navFrames(f.layout, false).flatMap((n) => n.wants),
      ];
      const key = `${wantsKey(wants)}#${store.nodeVersion.value}`;
      if (key === last) break;
      last = key;
      wanter.want(wants);
      if (!(await wanter.whenSettled(Math.max(0, until - clock.now())))) {
        settled = false;
        break;
      }
    }
    // the selection as an address: resolved over what the export holds now, as a viewer resolves it over its own
    let selection: Location | undefined;
    if (options.selection !== undefined) {
      const q = overlayQueryModel(model, overlays.of(model, store), { data: store });
      const r = resolveAddressIn(spec, q, options.selection);
      if (r.result) selection = locationOfResult(spec, q, r.result).location;
      else issues.push(...r.issues);
    }
    const withSelection =
      selection && !isEmptyLocation(selection)
        ? {
            ...options.emphasis,
            selection: options.emphasis?.selection ? unionLocations(options.emphasis.selection, selection) : selection,
          }
        : options.emphasis;
    // the selection carried onto the linked sections, as a viewer draws it (unless the caller gave its own)
    const drawn = withSelection?.selection;
    const emphasis =
      withSelection && drawn && !withSelection.projected?.selection
        ? {
            ...withSelection,
            projected: { ...withSelection.projected, selection: projectOntoView(model, view, drawn) },
          }
        : withSelection;
    let frame = renderer.frame({
      ...frameInput,
      ...(emphasis ? { emphasis } : {}),
      overlay: overlays.of(model, store),
      settled: true,
    });
    if (auto && frame.layout.height !== contentHeight(frame.layout)) {
      frameInput.size = { width: options.size.width, height: contentHeight(frame.layout) };
      frame = renderer.frame({
        ...frameInput,
        ...(emphasis ? { emphasis } : {}),
        overlay: overlays.of(model, store),
        settled: true,
      });
    }
    issues.push(...frame.layout.issues);
    // what the overlay's sources could not open, while their nodes are held
    issues.push(...nested.issues);
    if (!settled)
      issues.push(
        issueOf(RenderIssues, 'seq.render.unsettled', 'The export did not settle before its timeout', {
          hint: RenderIssues.codes['seq.render.unsettled'].hint,
        }),
      );
    // a shared store's issues are its other consumers' too: only those of the sources this export attached
    issues.push(...store.issues.value.filter((i) => own || (i.address !== undefined && attached.has(i.address))));
    const domains: string[] = [];
    for (const s of frame.layout.sections)
      for (const t of s.tracks)
        for (const l of t.layers) {
          const d = l.type.domain?.(l.lay);
          if (d) domains.push(`${l.key}=${d[0]},${d[1]}`);
        }
    const { legends, issues: legendIssues } = viewLegends(frame.layout);
    issues.push(...legendIssues);
    // each shared scale's domain as drawn: its legend's (a color scale) or its layers' y
    const scales: Record<string, string> = {};
    for (const l of legends)
      if (l.scale && l.domain) scales[`seq-scale-${l.scale.slice(l.scale.lastIndexOf('/') + 1)}`] = l.domain.join(',');
    for (const s of frame.layout.sections)
      for (const t of s.tracks)
        for (const l of t.layers) {
          const sc = l.layer.params.scale;
          const d = l.type.domain?.(l.lay);
          if (typeof sc === 'string' && d) scales[`seq-scale-${sc.slice(sc.lastIndexOf('/') + 1)}`] = d.join(',');
        }
    const place = options.legend ?? (view.export?.legend as RenderSvgOptions['legend']) ?? 'none';
    if (place !== 'none' && legends.length) {
      const lists = frame.lists;
      const text = options.text ?? bundledTextMetrics();
      const inset = config.inset;
      const band =
        place === 'below'
          ? { x: inset, y: lists.height + LegendBand.gap, width: lists.width - 2 * inset, height: 0 }
          : { x: lists.width + LegendBand.gap, y: inset, width: 240, height: 0 };
      const h = legendsHeight(legends, band.width, tokens, text);
      const w0 = lists.width;
      const h0 = lists.height;
      if (place === 'below') lists.height = h0 + h + 2 * LegendBand.gap;
      else {
        lists.width = w0 + band.width + 2 * LegendBand.gap;
        lists.height = Math.max(h0, h + 2 * inset);
      }
      const bg = lists.writer('background');
      if (place === 'below') bg.fill(0, h0, lists.width, lists.height - h0, tokens.colors['msui-bg']);
      else bg.fill(w0, 0, lists.width - w0, lists.height, tokens.colors['msui-bg']);
      if (lists.height > h0 && place === 'right') bg.fill(0, h0, w0, lists.height - h0, tokens.colors['msui-bg']);
      layoutLegends(legends, band, tokens, text, lists.writer('labels'));
    }
    // the ids (clips, masks, the title and description) tell two exports of one model apart: by what differs in them
    // (the title, size, tokens, view, highlights, selection); two equal exports share their ids, so a host that inlines
    // both gives each its own `idPrefix`
    const titleText = String(spec.params?.title ?? '') || 'Sequence view';
    const idPrefix =
      options.idPrefix ??
      `s${(
        hashString(
          [
            model.hash,
            titleText,
            options.size.width,
            options.size.height,
            tokens.version,
            view.address,
            ...(emphasis?.highlights ?? []).map((h) => locationText(model, h.location)),
            emphasis?.selection ? locationText(model, emphasis.selection) : '',
          ].join('|'),
        ) >>> 0
      ).toString(36)}-`;
    let svg = drawListsToSvg(frame.lists, {
      fonts: tokens.fonts,
      fade: config.fade,
      idPrefix,
      accessible: accessibleOf(model, frame.layout, titleText, emphasis?.highlights),
      metadata: {
        'seq-settled': String(settled),
        ...(domains.length ? { 'seq-domains': domains.join(';') } : {}),
        ...scales,
      },
      sections: frame.layout.sections.map((s) => ({
        address: s.section.address,
        direction: s.section.axis.dir === -1 ? 'reversed' : 'genomic',
        pieces: s.section.axis.pieces,
        y: s.y,
        height: s.height,
      })),
    });
    if (navs.length) {
      const drawn = navFrames(frame.layout, true);
      for (const n of navs) issues.push(...n.allIssues);
      svg = withNavigators(svg, frame.lists, model, navs, drawn, {
        fonts: tokens.fonts,
        fade: config.fade,
        idPrefix,
      });
    }
    return { svg, settled, issues, lists: frame.lists, ...(selection ? { selection } : {}) };
  } finally {
    nested.dispose();
    wanter.dispose();
    controller.abort();
    if (own) store.retire();
  }
}

/**
 * The view's SVG with its navigators stacked by `place` around it (top and bottom across its width, left and right
 * along its height, a miniature at its own size): each a group with its windows' `<title>` and `data-navigator` its
 * id; the view's metadata on the outer root.
 */
function withNavigators(
  view: string,
  lists: DrawLists,
  model: StaticModel,
  navs: readonly NavigatorRenderer[],
  frames: readonly NavigatorFrame[],
  svg: Pick<SvgOptions, 'fonts' | 'fade' | 'idPrefix'>,
): string {
  const vw = lists.width;
  const vh = lists.height;
  // each navigator's box, a zoom strip's height below it
  const boxes = navs.map((n) => {
    const b = navigatorBox(n.navigator, vw, vh);
    return { width: b.width, height: b.height + n.zoomHeight };
  });
  const sum = (p: string, k: 'width' | 'height') =>
    boxes.reduce((t, b, i) => t + (navs[i].navigator.place === p ? b[k] : 0), 0);
  const top = sum('top', 'height');
  const left = sum('left', 'width');
  const width = left + vw + sum('right', 'width');
  const height = top + vh + sum('bottom', 'height');
  const at = { top: 0, bottom: top + vh, left: 0, right: left + vw };
  const nest = (inner: string, x: number, y: number) => inner.replace(/^<svg /, `<svg x="${x}" y="${y}" `);
  const meta = /^<svg [^>]*?viewBox="[^"]*"([^>]*)>/.exec(view)?.[1] ?? '';
  const parts = [nest(view, left, top)];
  navs.forEach((n, i) => {
    const b = boxes[i];
    const place = n.navigator.place;
    const x = place === 'left' ? at.left : place === 'right' ? at.right : left;
    const y = place === 'top' ? at.top : place === 'bottom' ? at.bottom : top;
    if (place === 'top') at.top += b.height;
    else if (place === 'bottom') at.bottom += b.height;
    else if (place === 'left') at.left += b.width;
    else at.right += b.width;
    const titles = frames[i].windows
      .filter((w) => w.shown)
      .map((w) => `Window: ${windowText(model, n.navigator, w)}`)
      .join('; ');
    const inner = drawListsToSvg(frames[i].lists, { ...svg, idPrefix: `${svg.idPrefix ?? ''}n${i}-` });
    parts.push(
      `<g data-navigator="${escapeXml(n.navigator.id)}"><title>${escapeXml(titles || n.navigator.id)}</title>${nest(inner, x, y)}</g>`,
    );
  });
  return `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}"${meta}>${parts.join('')}</svg>`;
}
