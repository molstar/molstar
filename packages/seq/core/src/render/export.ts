/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashString } from '@molstar/core/data/util/hash-functions';
import { type Issue, issueOf, type Node } from '@molstar/core/repr/tree';
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
import { channelSource } from '@molstar/seqx-core/data/channels';
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
  axisTotal,
  BuiltInLayerKinds,
  type LayerKinds,
  nearestOffsetOf,
  projectOntoView,
  type ResolvedView,
  resolveView,
  resolveViews,
  type SystemPosition,
  unitOfPlace,
} from '@molstar/seqx-core/view';
import type { Viewport } from './axis.js';
import { type RenderConfigOverrides, renderConfig } from './config.js';
import { accessibleOf, locationText } from './accessibility.js';
import { itemTitles, type TitleOptions } from './export-titles.js';
import { DrawLists } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { layerSource } from './layer-data.js';
import type { ThemeContext } from './channels.js';
import type { Emphasis } from './layer-type.js';
import { LegendBand, layoutLegends, legendsHeight, viewLegends } from './legends.js';
import { OverlaySources } from './overlay-sources.js';
import type { NavigatorRendererOptions } from './navigator.js';
import type { NavigatorExport } from './navigator-export.js';
import { type PositionsRailMode, railSections } from './rail.js';
import { Renderer, type RendererOptions } from './renderer.js';
import type { ViewLayout } from './tracks.js';
import type { HeaderPartsOf } from './headers.js';
import { drawListsToSvg, type SvgAccessible } from './svg.js';
import { bundledTextMetrics, type TextMetrics } from './text.js';
import { builtInTokens, type SeqTokens } from './tokens.js';
import { tooltipSources } from './tooltip-wants.js';

// Headless rendering (_spec/README.md "Render", Export): a spec to an SVG string, asynchronously. It resolves the spec,
// opens the sources the view reads, wants what the export's size and viewport need through a want set of its own
// (never the live view's, EC-125), waits for those entries to settle or a timeout (then the export says it is
// unsettled), and draws once more with the view settled (open domains computed, lanes packed). Legends are drawn in a
// band outside the data box only when asked (`legend`, else the view's `export.legend`, else none), and each shared
// scale's domain as drawn is recorded (`data-seq-scale-<id>`): an open end unioned over what the export draws.

/**
 * The store's sources a view reads: its placed layers' data, its sections' systems' sources (for a length), the sources
 * of the numberings its layers name (a ruler's: no layer's data names them), and its navigators' systems' and bands'
 * (_spec/README.md "Render", Navigators: a chromosome under a region; Export).
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
    const systems = [s.axis.system, ...(s.blocks ?? []).map((b) => b.axis.system)];
    for (const t of s.tracks) {
      if (t.template) out.add(t.template.source);
      for (const l of t.layers) {
        const src = layerSource(model, l);
        if (src) out.add(src.source);
        // a numbering a layer names (a ruler's, a text's numbers): a source's rows when the system declares it so, read
        // at layout, which wants them only once the source is open (Render, Numberings from a source)
        const id = l.params.numbering;
        if (typeof id === 'string' && id)
          for (const system of [l.data?.system, ...systems]) {
            const n = system === undefined ? undefined : model.numberingOf(system, id)?.source;
            if (n) out.add(n);
          }
        // a rows layer's joined fields: the sources whose own system is its rows system (Render, Row columns)
        const rows = l.type === 'rows' && l.data ? model.annotations.get(l.data.address)?.rows : undefined;
        if (rows) for (const b of model.sources.values()) if (b.entry.system === rows) out.add(b.address);
      }
    }
  }
  tooltipSources(model, view, out);
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
    ...[...model.annotations.values()].flatMap((a) => (a.ends ? [...a.ends.built()] : [])),
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
  /** The layer kinds its view is resolved with (the built-ins' when absent): a composition's (`SeqContext.kinds`). */
  readonly kinds?: LayerKinds;
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
  const viewOf = (m: StaticModel) =>
    options.view ? resolveView(m, options.view, options.kinds) : resolveViews(m, options.kinds)[0];
  const scoped = (m: StaticModel, v: ResolvedView | undefined) =>
    options.scope === 'view' && v ? viewSources(m, v) : undefined;
  let model = options.model ?? resolveSpec(spec, resolveOptions);
  // a lazy part the spec needs (the genome chunk): resolved again once it lands, never drawn without it
  if (model.loading) {
    // never rejects (a failed load is the next model's `seq.resolve.genome.load`); caught all the same
    await model.loading.catch(() => undefined);
    model = resolveSpec(spec, resolveOptions);
  }
  let view = viewOf(model);
  const shown = scoped(model, view);
  const first = new Set(
    [...model.systems.values()].flatMap((s) => (s.source && (!shown || shown.has(s.source)) ? [s.source] : [])),
  );
  // an assembly's sources list its sequences: opened first, so its systems derive at the next resolve
  for (const a of model.assemblies.values()) for (const s of [a.source, a.aliasSource]) if (s) first.add(s);
  // data over an assembly lists the names it reads each sequence by: opened first too, its plans named at the next
  // resolve (_spec/README.md "Resolve", Assemblies)
  for (const a of model.annotations.values())
    if (a.bySystem && a.access?.kind === 'source') {
      const s = channelSource(a.access.source);
      if (!shown || shown.has(s)) first.add(s);
    }
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
      ...[...model.annotations.values()].flatMap((a) => (a.ends ? [...a.ends.built()] : [])),
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
  /**
   * Each section's window: `[start, end)` on its system, or for a section of blocks each end a system and a position
   * (`{ system, position }`, _spec/README.md "Resolved view", Places): never units.
   */
  readonly positions: Readonly<Record<string, readonly [number, number] | readonly [SystemPosition, SystemPosition]>>;
}

export interface RenderSvgOptions extends TitleOptions {
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
  /** The layer kinds the view is resolved with (the built-ins' when absent): a composition's (`SeqContext.kinds`). */
  readonly kinds?: LayerKinds;
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
   * Each header's parts (_spec/README.md "Render", Track headers): a viewer's entries, so the export draws the states
   * it shows (the slots never are); the built-in states when absent.
   */
  readonly headers?: HeaderPartsOf;
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
  readonly navigator?: Pick<NavigatorRendererOptions, 'zoomStrip' | 'direction'> & {
    /** Whether the included rails keep beside a navigator (`always`, the default), `unlessNavigator` or `none`. */
    readonly positionsRail?: PositionsRailMode;
  };
  /**
   * The layer types drawn (_spec/README.md "Context", Registry entries): the composition's registry. `exportSvg` draws
   * none when absent (each layer said); `renderSvg` draws the built-ins.
   */
  readonly types?: RendererOptions['types'];
  /** The section layouts (the composition's registry): tracks for every section in `exportSvg` when absent. */
  readonly layouts?: RendererOptions['layouts'];
  /** The track sizings past the built-ins and the host's `trackSizing` port (the composition's). */
  readonly sizings?: RendererOptions['sizings'];
  readonly trackSizing?: RendererOptions['trackSizing'];
  /**
   * The navigators' pictures (the composition's `navigator-picture` entries): what a navigator's `content` may name,
   * the one source of the view's pictures (`kinds.pictures`) and of their drawing. Given, or with `navigatorExport`
   * (none: each navigator the axis, said), the view is resolved with these; else with `kinds` as given.
   */
  readonly pictures?: NavigatorRendererOptions['pictures'];
  /**
   * How included navigators are drawn (`NavigatorSvgExport` of `render/navigator-export`): `renderSvg`'s default; none
   * in `exportSvg` when absent (an export asked to include navigators then says it left them out).
   */
  readonly navigatorExport?: NavigatorExport;
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
    if (typeof r[0] === 'object' || typeof r[1] === 'object') {
      // ends by system and position: each through the blocks' places (a position no block shows: the nearest edge)
      const a = typeof r[0] === 'object' ? unitOfPlace(s, r[0]) : undefined;
      const b = typeof r[1] === 'object' ? unitOfPlace(s, r[1], a) : undefined;
      if (a === undefined || b === undefined) continue;
      out[s.address] = { start: Math.min(a, b), end: Math.max(a, b) };
      continue;
    }
    const a = s.axis.ranges ? (nearestOffsetOf(s.axis, r[0]) ?? r[0]) : r[0];
    const b = s.axis.ranges ? (nearestOffsetOf(s.axis, r[1]) ?? r[1]) : r[1];
    // over a cyclic axis's end the second position is before the first: the window runs on, its end past the total
    const total = s.axis.cyclic ? axisTotal(s.axis) : undefined;
    out[s.address] =
      total !== undefined && b < a ? { start: a, end: b + total } : { start: Math.min(a, b), end: Math.max(a, b) };
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

/** The most errors a description names; the rest are counted. */
const MaxDescribedErrors = 10;

/** A figure's description with the spec's errors said, the first ten and a count of the rest (none: unchanged). */
export function withIssues(a: SvgAccessible, errors: readonly Issue[]): SvgAccessible {
  if (!errors.length) return a;
  const said = errors
    .slice(0, MaxDescribedErrors)
    .map((i) => `${i.code}${i.address ? ` at ${i.address}` : ''}: ${i.message}`)
    .join('; ');
  const more = errors.length > MaxDescribedErrors ? ` and ${errors.length - MaxDescribedErrors} more` : '';
  return {
    ...a,
    description: `${a.description}. The spec has ${errors.length === 1 ? 'an error' : 'errors'}: ${said}${more}`,
  };
}

function wantsKey(wants: readonly Want[]): string {
  return wants.map((w) => `${w.source}|${w.system}|${w.start}|${w.end}|${w.level}|${w.priority}`).join(';');
}

/**
 * A spec drawn to an SVG string once the data its view needs has settled (or the timeout passed), with what the options
 * give and nothing else: the layer types, layouts and navigators of its composition (`renderSvg` gives the built-ins).
 */
export async function exportSvg(spec: Node, given: RenderSvgOptions): Promise<RenderSvgResult> {
  // one source of the pictures: the view resolved with those the navigators are drawn from
  const pictures = given.pictures ?? (given.navigatorExport ? {} : undefined);
  const options = pictures
    ? { ...given, pictures, kinds: { ...(given.kinds ?? BuiltInLayerKinds), pictures: Object.keys(pictures) } }
    : given;
  const issues: Issue[] = [...checkSeqSpec(spec, options.bounds)];
  const blank = (settled: boolean): RenderSvgResult => {
    const lists = new DrawLists();
    lists.width = options.size.width;
    lists.height = options.size.height === 'auto' ? 0 : options.size.height;
    return { svg: drawListsToSvg(lists), settled, issues, lists };
  };
  // a spec with an error is drawn as the canvas draws it (what resolves); its errors go in the picture's description
  const specErrors = issues.filter((i) => i.severity === 'error');
  // except a spec past the bounds (`seq.lang.bound.*`): they guard what resolving costs, so it is not resolved
  if (specErrors.some((i) => i.code.startsWith('seq.lang.bound.'))) return blank(true);
  const resolveOptions = {
    ...(options.runtimes ? { runtimes: options.runtimes } : {}),
    config: options.resolve ?? {},
  };
  let model = resolveSpec(spec, resolveOptions);
  const firstView = options.view
    ? resolveView(model, options.view, options.kinds)
    : resolveViews(model, options.kinds)[0];
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
        ...(options.kinds ? { kinds: options.kinds } : {}),
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
    // a figure shows data: every lane of a capped track at its content height, to `export.maxTrackPx`
    const asked = renderConfig(options.config);
    const config =
      asked.export.lanes === 'all' ? { ...asked, track: { ...asked.track, maxPx: asked.export.maxTrackPx } } : asked;
    const renderer = new Renderer({
      model,
      view,
      data: store,
      tokens,
      text: options.text ?? bundledTextMetrics(),
      config,
      // none of a kind not given (each layer said, each section laid out as tracks)
      types: options.types ?? {},
      layouts: options.layouts ?? {},
      ...(options.sizings ? { sizings: options.sizings } : {}),
      ...(options.trackSizing ? { trackSizing: options.trackSizing } : {}),
      ...(options.registry ? { theme: options.registry } : {}),
    });
    const given = options.viewport;
    const auto = options.size.height === 'auto';
    const windows = options.window ?? view.export?.window === true;
    // the navigators an export includes (shown ones), drawn as outlines; the rails plain (a figure's)
    const included = (options.navigators ?? view.export?.navigators ?? 'omit') === 'include';
    const navs =
      included && options.navigatorExport
        ? options.navigatorExport({
            // the view's renderer's: its model, view, store, tokens, text, config, types, layouts, sizings, theme
            ...renderer.options,
            pictures: options.pictures!,
            ...(options.kinds ? { kinds: options.kinds } : {}),
            images: renderer.images,
            domains: renderer.layoutState,
            storeSettings: store.settings,
            ...options.navigator,
            dpr: options.dpr ?? 1,
            cvdSafe: options.cvdSafe ?? false,
            overlay: () => overlays.of(model, store),
          })
        : undefined;
    if (included && !navs && view.navigators.some((n) => !n.hidden))
      issues.push(
        issueOf(RenderIssues, 'seq.render.navigator.export', 'This export draws no navigators: they are left out', {
          hint: RenderIssues.codes['seq.render.navigator.export'].hint,
        }),
      );
    // the rails plain; a navigator this export cannot draw replaces none (its section keeps its rail)
    const railed = navs ? view : { ...view, navigators: [] };
    const rails = included
      ? Object.fromEntries(railSections(railed, [], options.navigator?.positionsRail).map((a) => [a, true]))
      : undefined;
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
      ...(options.headers ? { headers: options.headers } : {}),
    };
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
      const wants = [...f.wants.filter((w) => w.priority !== 'margin'), ...(navs?.frames(f.layout, false).wants ?? [])];
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
      accessible: withIssues(
        accessibleOf(model, frame.layout, titleText, emphasis?.highlights, itemTitles(options, view, model, store)),
        specErrors,
      ),
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
        ...(s.section.blocks
          ? {
              blocks: s.section.blocks.map((b) => {
                const g = b.group !== undefined ? s.section.groups?.[b.group] : s.section.groups?.[0];
                return {
                  address: b.address,
                  system: b.axis.system,
                  direction: b.axis.dir === -1 ? 'reversed' : 'genomic',
                  scale: b.scale ?? 1,
                  ...(g ? { assembly: g.assembly } : {}),
                };
              }),
            }
          : {}),
      })),
    });
    // a bridge a group with its sides (its ribbons are the passes' cubic paths)
    if (frame.layout.bridges) svg = frame.layout.bridges.svg(svg);
    if (navs?.count) {
      const drawn = navs.frames(frame.layout, true).frames;
      issues.push(...navs.issues());
      svg = navs.wrap(svg, frame.lists, model, drawn, { fonts: tokens.fonts, fade: config.fade, idPrefix });
    }
    return { svg, settled, issues, lists: frame.lists, ...(selection ? { selection } : {}) };
  } finally {
    nested.dispose();
    wanter.dispose();
    controller.abort();
    if (own) store.retire();
  }
}
