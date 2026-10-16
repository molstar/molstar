/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `renderSvg` (seqspec:viewer#the-render-backend): a whole view as one SVG string, headless. It makes a context from a
// document (the default registry, a manual clock, the port's text measurer) or reuses one, sizes each section from
// `width` (its content's height), has a `SectionComposer` (and, with `nav`, the `VerticalNav` at its box) draw each
// into one `SvgBackend`, and stacks the sections. It draws the view as it stands (zoom, scroll, selection): a context
// to export is positioned first, and is left as it was.
// A document's streamed sources are awaited (`streams.settled`, ST1): its sections sized as they are drawn (their
// content's height), the manual clock run until every region in view has loaded or failed, again while an arrival
// changes a height; `timeout` bounds the wait. A live context draws its regions as they stand (a host awaits
// `ctx.streams.settled()` first to have them all).
//
// A live context is borrowed, not copied (a fresh one would load and prepare its data again): every section is resized,
// drawn and put back (size, DPR, zoom, scroll, row scale) in one synchronous run with no await, so no frame and no
// ResizeObserver callback sees the export's size, and a clamp at that size (a scroll past the narrower end) does not
// stay. What precedes it, a measurer, a skin or another view and the graph update they need, can span tasks: while the
// export changes any of them the viewer is paused (`setActive(false)`), so its surfaces draw nothing of it. The run
// commits no viewport (`Viewport.quiet`): a host syncing a URL from `committed` sees no export.

import type { Loaders, SeqViewSpecDoc } from '@molstar/seq-core/format';
import { ChannelMarks } from '../channels.js';
import { createSeqViewer, type SeqIssue, type SeqViewerContext } from '../context.js';
import { Dirty, SectionComposer } from '../draw/composer.js';
import { drawResources } from '../draw/resources.js';
import { VerticalNav } from '../draw/vnav.js';
import type { Skin } from '../graph.js';
import { Hover } from '../hover.js';
import type { RenderProps } from '../layout.js';
import type { SeqViewerRegistryEntry } from '../registry.js';
import { type ManualScheduler, manualScheduler } from '../scheduler.js';
import { DefaultRegistry } from '../default-spec.js';
import { Viewport } from '../viewport.js';
import { SvgBackend, type SvgOptions, type SvgPorts } from './backend.js';
import { hasMetrics, type MetricsMeasurer, metricsMeasurer } from './metrics.js';
import { DefaultFontMetrics } from './metrics-tables.js';

export interface RenderSvgOptions extends SvgOptions {
  /** CSS px, every section's width. */
  width: number;
  /** The view to show (`ctx.views.setCurrent`); default the context's. */
  view?: string;
  /** The sections to draw, stacked in this order; default all of the view's. */
  sections?: readonly string[];
  /** CSS px between stacked sections. */
  gap?: number;
  ports?: SvgPorts;
  /** For a document: registry entries beyond `DefaultRegistry`, loaders for its sources, render props. */
  registry?: readonly SeqViewerRegistryEntry[];
  loaders?: Loaders;
  render?: Partial<RenderProps>;
  /** The skin to draw with; a context passed in gets its own back. */
  skin?: Skin;
  /**
   * For a document: at most this many ms (wall clock) waiting for its sources (a static load, the streamed regions),
   * then what stands is drawn and `unsettled` is called. Default: no bound.
   */
  timeout?: number;
  /**
   * For a document: called once when the export draws before what it waits for arrived (`timeout`, or the heights
   * still moving after the rounds), with the stream requests still in flight.
   */
  unsettled?: (inFlight: number) => void;
  /** For a document: its context's issues as it was drawn (a source not read, a provider not there, a section left out). */
  issues?: (issues: readonly SeqIssue[]) => void;
}

let fallback: MetricsMeasurer | null = null;
/**
 * The measurer of a context `renderSvg` makes when the host gives none: the committed tables
 * (seqspec:viewer#text-measuring).
 */
const defaultMeasurer = () => (fallback ??= metricsMeasurer(DefaultFontMetrics));

/** A section as the export found it: size, DPR, zoom, scroll and row scale. */
type Kept = [ref: string, width: number, height: number, dpr: number, cell: number, x: number, y: number, rows: number];

/** The SVG of `source` (a context, positioned as it should print; or a document, drawn fitted) at `opts.width`. */
export async function renderSvg(source: SeqViewerContext | SeqViewSpecDoc, opts: RenderSvgOptions): Promise<string> {
  const ports = opts.ports ?? {};
  if (!('setDocument' in source)) {
    const scheduler = manualScheduler();
    const ctx = createSeqViewer(
      { registry: [DefaultRegistry, ...(opts.registry ?? [])], behaviors: [], render: opts.render },
      { scheduler, loaders: opts.loaders, text: ports.measurer ?? defaultMeasurer() },
    );
    let timer: ReturnType<typeof setTimeout> | undefined;
    try {
      if (opts.skin) ctx.setSkin(opts.skin);
      const end = opts.timeout === undefined ? Number.POSITIVE_INFINITY : Date.now() + opts.timeout;
      // The load within the bound too: a source that never answers is drawn as it stands (nothing of it).
      const late = new Promise<null>((ok) => {
        if (Number.isFinite(end)) timer = setTimeout(() => ok(null), Math.max(0, end - Date.now()));
      });
      const r = await Promise.race([ctx.setDocument(source), late]);
      if (r === null) opts.unsettled?.(ctx.streams.stats.current.inFlight);
      else {
        if (!r.ok) throw new Error(`renderSvg: the source does not load: ${JSON.stringify(r.issues)}`);
        if (opts.view) ctx.views.setCurrent(opts.view);
        await ctx.graph.update();
        await settle(ctx, opts, scheduler, end);
      }
      const svg = await draw(ctx, opts, null).toString(ports);
      opts.issues?.(ctx.issues.current.value);
      return svg;
    } finally {
      clearTimeout(timer);
      ctx.dispose();
    }
  }
  const ctx = source;
  const was = {
    measurer: ctx.text.port,
    mono: ctx.text.families.mono,
    ui: ctx.text.families.ui,
    skin: ctx.skin.value,
    view: ctx.views.current.value,
  };
  const pause = ctx.ui.active.value && (!!ports.measurer || !!opts.skin || (!!opts.view && opts.view !== was.view));
  const kept: Kept[] = [];
  const quiet: (() => void)[] = [];
  let be: SvgBackend;
  try {
    if (pause) ctx.setActive(false);
    if (ports.measurer) ctx.setTextMeasurer(ports.measurer);
    if (opts.skin) ctx.setSkin(opts.skin);
    if (opts.view) ctx.views.setCurrent(opts.view);
    await ctx.graph.update();
    // From here to the restore, synchronous (the header). No viewport commits: each is put back as it was, so a host
    // that syncs from `committed` sees no export.
    for (const ref of ctx.viewports.sections.value) {
      const vp = ctx.viewports.get(ref);
      if (vp instanceof Viewport) quiet.push(vp.quiet());
    }
    be = draw(ctx, opts, kept);
  } finally {
    if (!ctx.disposed) {
      // The text first: the layouts rebuilt below are the live ones, whose rows' height a scroll is clamped against.
      if (ctx.text.port !== was.measurer) ctx.setTextMeasurer(was.measurer);
      ctx.text.setFamilies(was.mono, was.ui);
      for (let k = kept.length - 1; k >= 0; k--) {
        const [ref, w, h, dpr, cell, x, y, rows] = kept[k];
        const vp = ctx.viewports.get(ref);
        if (!vp) continue;
        vp.setSize(w, h, dpr);
        // A section with no width keeps its view for the next size; a fitted one refits to the same.
        if (w <= 0) continue;
        ctx.layout.get(ref);
        const s = vp.live.current;
        if (s.cell !== cell || s.scrollX !== x) vp.setView(cell, x, 'api');
        if (s.rowScale !== rows || s.scrollY !== y) vp.setRowView(rows, y, 'api');
      }
      for (const q of quiet) q();
      if (was.view !== null && ctx.views.current.value !== was.view) ctx.views.setCurrent(was.view);
      if (opts.skin) {
        // `setSkin(null)` shows the reader's skin; a host's override held before is set again.
        ctx.setSkin(null);
        if (ctx.skin.value !== was.skin) ctx.setSkin(was.skin);
      }
      if (pause) ctx.setActive(true);
    }
  }
  return be.toString(ports);
}

/** Section `ref` sized as the export draws it: once to build the layout, then to its content's height; that height. */
function sized(ctx: SeqViewerContext, ref: string, width: number, scale: number): number {
  const vp = ctx.viewports.get(ref);
  if (!vp) return -1;
  vp.setSize(width, 1000, scale);
  const h = Math.max(1, Math.ceil(ctx.layout.contentHeight(ref)?.value ?? 1000));
  vp.setSize(width, h, scale);
  return h;
}

/**
 * A document's streams settled at the size the export draws (`sized`), the clock run: again while an arrival changed a
 * section's height (a streamed alignment's rows), at most `ROUNDS` times; within `opts.timeout`, else `unsettled`.
 */
async function settle(ctx: SeqViewerContext, opts: RenderSvgOptions, scheduler: ManualScheduler, end: number) {
  const refs = opts.sections ?? ctx.viewports.sections.value;
  const heights = refs.map(() => -1);
  /** Each section sized at the export's width: whether a height moved since the last round. */
  const resized = () => {
    let moved = false;
    for (let k = 0; k < refs.length; k++) {
      const h = sized(ctx, refs[k], opts.width, opts.scale ?? 1);
      if (h !== heights[k]) moved = true;
      heights[k] = h;
    }
    return moved;
  };
  for (let round = 0; resized(); round++) {
    // Heights still moving after the rounds: drawn as they stand, and said so.
    if (round === ROUNDS) return opts.unsettled?.(ctx.streams.stats.current.inFlight);
    let done = false;
    const settled = ctx.streams.settled().then(() => {
      done = true;
    });
    // The clock is read after a first wait, so what has nothing to wait for is settled whatever the bound.
    for (let first = true; !done; first = false) {
      if (!first && Date.now() >= end) return opts.unsettled?.(ctx.streams.stats.current.inFlight);
      scheduler.advance(16);
      await Promise.race([settled, new Promise((ok) => setTimeout(ok, 1))]);
    }
    await ctx.graph.update();
  }
}
const ROUNDS = 4;

/** Each section sized and drawn into one backend, synchronously; `kept` collects the sections as they were. */
function draw(ctx: SeqViewerContext, opts: RenderSvgOptions, kept: Kept[] | null): SvgBackend {
  ctx.text.setFamilies(ctx.skin.value.fontMono, ctx.skin.value.fontUi);
  const scale = opts.scale ?? 1;
  const measurer = ctx.text.port;
  const be = new SvgBackend(drawResources(ctx), {
    ...opts,
    scale,
    metrics: opts.metrics ?? (hasMetrics(measurer) ? (f) => measurer.metrics(f) : undefined),
  });
  const marks =
    opts.overlay || opts.nav
      ? new ChannelMarks({
          graph: ctx.graph,
          hover: opts.hover ? ctx.hover : new Hover(),
          selection: ctx.selection,
          highlight: ctx.highlight,
          systems: ctx.systems,
          layout: ctx.layout,
        })
      : null;
  const flags = Dirty.tracks | (opts.overlay ? Dirty.overlay : 0) | (opts.nav ? Dirty.nav : 0);
  let y = 0;
  try {
    for (const ref of opts.sections ?? ctx.viewports.sections.value) {
      const vp = ctx.viewports.get(ref);
      if (!vp) throw new Error(`renderSvg: no section "${ref}" in the view.`);
      const s = vp.live.current;
      kept?.push([ref, s.width, s.height, s.dpr, s.cell, s.scrollX, s.scrollY, s.rowScale]);
      // The height is the content's: sized once to build the layout, then to what it asks.
      sized(ctx, ref, opts.width, scale);
      const composer = new SectionComposer(ctx, ref, marks ?? undefined);
      // The vertical overview and the scrollbar, navigator passes: drawn with `nav` (its passes are muted otherwise).
      const vnav = marks && opts.nav ? new VerticalNav(ctx, ref, marks) : null;
      try {
        be.beginSection(0, y);
        composer.frame(be, flags);
        const l = ctx.layout.get(ref)?.current;
        if (vnav && l) {
          // In device px from its box's origin, as the surface's two canvases are placed.
          const box = vnav.box;
          vnav.place(l, ctx.layout.props(ref), scale);
          if (box[2] > 0 && box[3] > 0) {
            be.beginOffset(box[0] / scale, box[1] / scale);
            vnav.frame(be, flags);
            be.endOffset();
          }
        }
        be.endSection();
      } finally {
        composer.dispose();
        vnav?.dispose();
      }
      y += (ctx.layout.get(ref)?.current.height ?? 0) + (opts.gap ?? 8);
    }
  } finally {
    marks?.dispose();
  }
  return be;
}
