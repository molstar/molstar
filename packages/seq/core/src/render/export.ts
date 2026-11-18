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
import { DrawLists } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { layerSource } from './layer-data.js';
import type { Emphasis } from './layer-type.js';
import { OverlaySources } from './overlay-sources.js';
import { Renderer } from './renderer.js';
import { drawListsToSvg } from './svg.js';
import { bundledTextMetrics, type TextMetrics } from './text.js';
import { builtInTokens, type SeqTokens } from './tokens.js';

// Headless rendering (_spec/README.md "Render", Export): a spec to an SVG string, asynchronously. It resolves the spec,
// opens the sources the view reads, wants what the export's size and viewport need through a want set of its own
// (never the live view's, EC-125), waits for those entries to settle or a timeout (then the export says it is
// unsettled), and draws once more with the view settled (open domains computed, lanes packed).

/** The store's sources a view reads: its placed layers' data, and its sections' systems' sources (for a length). */
export function viewSources(model: StaticModel, view: ResolvedView): Set<string> {
  const out = new Set<string>();
  for (const s of view.sections) {
    const source = model.systems.get(s.axis.system)?.source;
    if (source) out.add(source);
    for (const t of s.tracks) {
      if (t.template) out.add(t.template.source);
      for (const l of t.layers) {
        const src = layerSource(model, l);
        if (src) out.add(src.source);
      }
    }
  }
  return out;
}

/** An asset reader that reads nothing: a headless render reads no file or URL unless given a reader (security). */
export const NoAssets: AssetReader = {
  readText: (uri) => Promise.reject(new Error(`This render reads no assets: ${uri}`)),
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
  const ctx = { assets: options.assets ?? NoAssets, signal: options.signal ?? new AbortController().signal };
  const opens = new Map<string, () => Promise<OpenedSource>>();
  for (const plan of [
    ...[...model.systems.values()].map((s) => s.access),
    ...[...model.annotations.values()].map((a) => a.access),
  ])
    if (plan.kind === 'inline') opens.set(plan.source.address, () => plan.runtime.open(plan.source, ctx));
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

/**
 * A viewport in positions, by section address (`[start, end)` on each section's system): read through each section's
 * axis (its direction and pieces), so a caller needs no knowledge of them.
 */
export interface PositionsViewport {
  readonly positions: Readonly<Record<string, readonly [number, number]>>;
}

export interface RenderSvgOptions {
  /** CSS px. */
  readonly size: { readonly width: number; readonly height: number };
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
  /** Prefixes the SVG's ids: by default `s` and a hash of the spec, so two exports inlined in one page differ. */
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
    lists.height = options.size.height;
    return { svg: drawListsToSvg(lists), settled, issues, lists };
  };
  if (hasErrors(issues)) return blank(true);
  const model = resolveSpec(spec, {
    ...(options.runtimes ? { runtimes: options.runtimes } : {}),
    config: options.resolve ?? {},
  });
  issues.push(...model.issues);
  const view = options.view ? resolveView(model, options.view) : resolveViews(model)[0];
  if (!view) {
    issues.push(
      issueOf(RenderIssues, 'seq.render.view', 'The spec has no such view', {
        ...(options.view ? { address: options.view } : {}),
        hint: RenderIssues.codes['seq.render.view'].hint,
      }),
    );
    return blank(true);
  }
  issues.push(...view.issues);
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
    const opened = await Promise.race([
      attachSources(model, store, {
        ...(options.assets ? { assets: options.assets } : {}),
        signal: controller.signal,
        only: viewSources(model, view),
        attached,
      }),
      expired,
    ]);
    clock.clearTimeout(timer);
    if (opened) issues.push(...opened);
    else {
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
    });
    const given = options.viewport;
    const frameInput = {
      size: options.size,
      ...(given && 'positions' in given
        ? { viewports: { ...offsetsOf(view, given), ...options.viewports } }
        : {
            ...(given ? { viewport: given } : {}),
            ...(options.viewports ? { viewports: options.viewports } : {}),
          }),
      dpr: options.dpr ?? 1,
      ...(options.emphasis ? { emphasis: options.emphasis } : {}),
    };
    let last = '';
    // what a layout wants can follow from what arrived (a length, a follower's range): ask until it holds; the
    // export waits for its view alone, never its margins
    for (let round = 0; settled && round < 8; round++) {
      const overlay = overlays.of(model, store);
      await nested.sync(overlay, view);
      const f = renderer.frame({ ...frameInput, overlay: overlays.of(model, store), settled: false });
      // what the layout shows stays while the export waits
      nested.hold(f.layout.shown);
      const wants = f.wants.filter((w) => w.priority !== 'margin');
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
      const q = overlayQueryModel(model, overlays.of(model, store));
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
    const frame = renderer.frame({
      ...frameInput,
      ...(emphasis ? { emphasis } : {}),
      overlay: overlays.of(model, store),
      settled: true,
    });
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
          const legend = l.type.legend?.(l.lay);
          if (legend?.domain) domains.push(`${l.key}=${legend.domain[0]},${legend.domain[1]}`);
        }
    const svg = drawListsToSvg(frame.lists, {
      fonts: tokens.fonts,
      fade: config.fade,
      idPrefix: options.idPrefix ?? `s${(hashString(model.hash) >>> 0).toString(36)}-`,
      metadata: { 'seq-settled': String(settled), ...(domains.length ? { 'seq-domains': domains.join(';') } : {}) },
      sections: frame.layout.sections.map((s) => ({
        address: s.section.address,
        direction: s.section.axis.dir === -1 ? 'reversed' : 'genomic',
        pieces: s.section.axis.pieces,
        y: s.y,
        height: s.height,
      })),
    });
    return { svg, settled, issues, lists: frame.lists, ...(selection ? { selection } : {}) };
  } finally {
    nested.dispose();
    wanter.dispose();
    controller.abort();
    if (own) store.retire();
  }
}
