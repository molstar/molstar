/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import type { Want } from '@molstar/seqx-core/data';
import {
  expandTemplate,
  type ResolvedLayer,
  type ResolvedSection,
  type ResolvedTrack,
  sameBlocks,
  templateNodes,
} from '@molstar/seqx-core/view';
import { AxisTransform, type Viewport } from './axis.js';
import { DefaultRenderConfig, fontSize, type RenderConfig } from './config.js';
import { RenderIssues } from './issues.js';
import { type Box, EmptyOverlay, type LabelSink, type LayerType, type LayoutContext, sectionId } from './layer-type.js';
import type { LettersLay } from './letters.js';
import type { ThemeContext } from './channels.js';
import { colorKey, defaultColors, viewTheme } from './layer-color.js';
import { LabelLane } from './sizing.js';
import { blockLength, chipText } from './layout/blocks.js';
import { TracksLayout } from './layout/tracks.js';
import type { GridPlacement } from './layout/grid.js';
import type { PlaceContext, Placement, SectionLayoutType } from './layout/types.js';
import { narrowColumns, treeKept } from './narrow.js';
import { placeRail, placeRowsRail, Rail, SlimRailPx } from './rail.js';
import type { RowsLay } from './rows.js';
import { rowsCount, rowsFillHeight, rowsShown } from './rows-fill.js';
import {
  type LayerLayout,
  type LayoutInput,
  LayoutState,
  type SectionLayout,
  type TrackLayout,
  type ViewLayout,
} from './tracks-types.js';
import { followAxis, followedWindow, multiSystem, ownPieces, sectionLayoutName, turned } from './tracks-axes.js';
import { type FrameSinks, reusable, sameGridY, SectionCache, viewStamp } from './tracks-cache.js';
import { acrossSpan, pendingPlacement, stripColumn, stripSpan } from './tracks-strips.js';
import {
  columnsFit,
  elidedLabel,
  labelColumnWidth,
  layerKey,
  measured,
  numbersGutter,
  rowsHeader,
  titleRowOf,
  viewColumns,
  widened,
} from './tracks-columns.js';
import { instanceInView, wantTemplate } from './tracks-templates.js';
import { bandHeights, firstType, rowsFill, rulerTrack } from './tracks-heights.js';
import { sectionHeaderOf, trackHeaderOf } from './tracks-headers.js';
import { collapsedSection } from './headers.js';

export { multiSystem, sectionLayoutName } from './tracks-axes.js';
export { labelColumnWidth, layerKey } from './tracks-columns.js';
export { BreakSlop, drawView } from './tracks-draw.js';
export { hitView } from './tracks-hit.js';
export {
  type HeldDomains,
  type LayerLayout,
  type LayoutInput,
  LayoutState,
  type SectionLayout,
  type TrackLayout,
  type ViewLayout,
} from './tracks-types.js';

// The tracks layout (_spec/README.md "Render", Tracks layout): sections stacked, each with its axis transform from its
// system and range (a follower's from its leader's view), tracks with heights in px under one label column, each
// track's layers in order over its box. Tracks below the page are not laid out: the draw lists are bounded by the
// screen. The layout keeps each layer's last layout, so buffers and lanes carry over frame to frame. This module is the
// layout pass (`layoutView`) and the tracks' public names; its records, axes, reuse, strips, label column, templates
// and heights, and the view's draw, emphasis and hit-test, are the `tracks-*.ts` modules beside it.

/** No layer types, no layouts: what a layout is given when its input names none (it draws nothing and says so). */
const NoTypes: Readonly<Record<string, LayerType<any, any>>> = {};
const NoLayouts: Readonly<Record<string, SectionLayoutType<any>>> = {};

/** A section's direction mark, or none: a reversed axis says so (`reversed`, or `along` what it reads along). */
function directionMark(section: ResolvedSection, config: RenderConfig): string | undefined {
  if (config.axis.directionMark === 'none' || section.axis.dir !== -1) return undefined;
  const along = section.axis.along;
  return along === undefined ? '← reversed' : `← along ${along.slice(along.lastIndexOf('/') + 1)}`;
}

const NoLays: readonly LayerLayout[] = [];
const NoSectionLays: LayerLayout[] = [];

export function layoutView(input: LayoutInput, state = new LayoutState(), again = false): ViewLayout {
  const header = rowsHeader(input.view, input.rails);
  const given = input.config ?? DefaultRenderConfig;
  const base = measured(given, labelColumnWidth(input.view, given, input.text));
  // the first numbers take their own room in the label column, the labels keep theirs
  const first = header.numbers ? Math.ceil(input.text.measure(base.itemFont, '00000')) + 4 : 0;
  // the last numbers' gutter, and beside it a rows rail's column (fixed: a rail idle or live keeps the data's width)
  const gutter = (header.numbers ? numbersGutter(input.text, base) : 0) + (header.rail ? SlimRailPx : 0);
  // the row columns: the left ones widen the label column, the right ones the right gutter (and a gap before a rows
  // rail), yielding to a narrow view; kept while the view, the width and the room they read are
  const fixed = 2 * base.inset + base.labelWidth + base.labelGap + header.tree + first + gutter;
  const fit = columnsFit(state, viewColumns(input.view), input.size.width, fixed, base, header.rail);
  // a narrow frame: the label column to its floor, then the tree collapsed, the plot kept (narrow.ts)
  const narrow = header.tracks
    ? narrowColumns(
        input.size.width,
        fixed + (fit ? fit.left + fit.right : 0),
        base.labelWidth,
        header.tree,
        base,
        fixed + (fit?.own ?? 0),
        treeKept(input.viewports),
      )
    : undefined;
  const config = widened(
    narrow ? measured(given, narrow.labelWidth) : base,
    (narrow?.tree ?? header.tree) + first + (fit?.left ?? 0),
  );
  const right = gutter + (fit?.right ?? 0);
  const grants = fit?.grants;
  // a text section's lines leave the right columns their room (the text layout known by its registered name, as the
  // rows and grid ones below: never by the built-in's identity, which would bundle it into every composition)
  const textWidth = input.size.width - (fit?.right ?? 0);
  const extents = new Map<string, [number, number]>();
  const theme: ThemeContext = viewTheme(
    input.view,
    input.theme,
    input.domains
      ? (key) => input.domains!.state.domains.get(input.domains!.map(key))
      : (key) => state.domains.get(key),
    (key, lo, hi) => {
      const e = extents.get(key);
      if (!e) extents.set(key, [lo, hi]);
      else {
        if (lo < e[0]) e[0] = lo;
        if (hi > e[1]) e[1] = hi;
      }
    },
  );
  const dpr = input.dpr ?? 1;
  const types = input.types ?? NoTypes;
  const overlay = input.overlay ?? EmptyOverlay;
  const { width, height } = input.size;
  const wants: Want[] = [];
  const issues: Issue[] = [];
  const seen = new Set<string>();
  const report = (i: Issue) => {
    const k = `${i.code}|${i.address ?? ''}`;
    if (seen.has(k)) return;
    seen.add(k);
    issues.push(i);
  };
  const dataX = config.inset + config.labelWidth + config.labelGap;
  const dataW = Math.max(1, width - dataX - config.inset - right);
  // axes first, leaders before their followers (a follower may come first in the view)
  const bySection = new Map(input.view.sections.map((s) => [s.address, s]));
  const layouts = input.layouts ?? NoLayouts;
  /** A section's layout: its name's in the layouts given, else tracks (a layout not given is said, once). */
  const layoutOf = (section: ResolvedSection) => {
    const name = sectionLayoutName(section, input.view);
    const layout = layouts[name];
    if (layout) return layout;
    if (name !== TracksLayout.name)
      report(
        issueOf(
          RenderIssues,
          'seq.render.layout',
          `The layout ${name} is not drawn here: the section is laid out as tracks`,
          {
            address: section.address,
            hint: RenderIssues.codes['seq.render.layout'].hint,
          },
        ),
      );
    return TracksLayout;
  };
  /**
   * What a follower reads of its leader (_spec/README.md "Render", Section layouts): the leader's viewport as its layout
   * reads it (`followWindow`: a ring's window, a text's lines in view), else its transform's window (tracks).
   */
  const leaderWindow = (section: ResolvedSection, axis: AxisTransform): Viewport => {
    const layout = layoutOf(section);
    if (!layout.followWindow) return { start: axis.u0, end: axis.u1 };
    const viewport = input.viewports?.[section.address] ?? input.viewport;
    const heights = bandHeights(
      input,
      state,
      config,
      types,
      dataW,
      section,
      section.tracks,
      layout,
      true,
      report,
    ).heights;
    const pads = section.tracks.flatMap((t) => {
      const p = t.padding ?? config.trackPadding;
      return [p.top, p.bottom];
    });
    return layout.followWindow(section, axis, viewport, {
      width: layout.name === 'text' ? textWidth : width,
      config,
      text: input.text,
      heights,
      pads,
    });
  };
  const axes = new Map<string, AxisTransform | undefined>();
  const visiting = new Set<string>();
  /** A following ring's window: its leader's range on its own axis (made when a ring follows). */
  let followed: Map<string, Viewport> | undefined;
  const axisOf = (section: ResolvedSection): AxisTransform | undefined => {
    if (axes.has(section.address)) return axes.get(section.address);
    const at = section.address;
    let made: { pieces: readonly number[]; viewport?: Viewport } | undefined;
    if (section.follow) {
      const leaderSection = bySection.get(section.follow.leader);
      if (visiting.has(at) || !leaderSection) {
        const cycle = visiting.has(at);
        const code = cycle ? 'seq.render.follow.cycle' : 'seq.render.follow.leader';
        report(
          issueOf(
            RenderIssues,
            code,
            cycle ? 'The section follows itself through its leaders' : 'The section’s leader is not in the view',
            {
              address: at,
              hint: RenderIssues.codes[code].hint,
            },
          ),
        );
      } else if (sameBlocks(section, leaderSection) && section.axis.system !== leaderSection.axis.system) {
        // equal block lists follow by identity in units (_spec/README.md "Render", Blocks): the leader's window is the
        // follower's (a ring of the genome and its linear detail, two sections of one assembly)
        visiting.add(at);
        const leader = axisOf(leaderSection);
        visiting.delete(at);
        const own = ownPieces(section, input.model, input.data);
        if (leader && own) {
          const window = leaderWindow(leaderSection, leader);
          made = { pieces: own, viewport: window };
          if (layoutOf(section).window) {
            followed ??= new Map();
            followed.set(at, window);
            made = { pieces: own };
          }
        }
      } else if (
        (multiSystem(section) || multiSystem(leaderSection)) &&
        section.axis.system !== leaderSection.axis.system
      ) {
        // which block's map a follow would carry is a figure's choice no rule makes (_spec/README.md "Render",
        // Blocks); a grid over a section of blocks is on the blocks' own units, its follow the identity ("The grid
        // layout", A grid over blocks)
        report(
          issueOf(RenderIssues, 'seq.render.follow.blocks', 'Sections of unequal block lists cannot follow', {
            address: at,
            hint: RenderIssues.codes['seq.render.follow.blocks'].hint,
          }),
        );
      } else {
        visiting.add(at);
        const leader = axisOf(leaderSection);
        visiting.delete(at);
        const own = ownPieces(section, input.model, input.data);
        if (leader) made = followAxis(section, own, leader, leaderWindow(leaderSection, leader), config);
        // a follower that draws its whole axis with a window (a ring): the leader's range is its window
        if (made && own && layoutOf(section).window) {
          followed ??= new Map();
          followed.set(at, followedWindow(section, own, made));
          made = { pieces: own };
        }
      }
    }
    if (!made) {
      const own = ownPieces(section, input.model, input.data);
      // a layout that reads its viewport its own way (text: an anchor; a ring: its window) places over the whole axis
      const viewport = layoutOf(section).followWindow ? undefined : (input.viewports?.[at] ?? input.viewport);
      if (own) made = { pieces: turned(section, own, viewport), ...(viewport ? { viewport } : {}) };
    }
    const axis = made?.pieces.length
      ? AxisTransform.on(section.axis, made.pieces, dataX, dataW, dpr, made.viewport)
      : undefined;
    if (!axes.has(at)) axes.set(at, axis);
    return axes.get(at);
  };
  for (const section of input.view.sections) axisOf(section);
  /**
   * A grid's y (_spec/README.md "Render", The grid layout): x's pieces and window when it follows x on x's system; a
   * leader's image through its path (x's own, or a section's); else its own window (`viewport.y`), all of it by default.
   */
  const yOf = (section: ResolvedSection, axis: AxisTransform | undefined): PlaceContext['y'] => {
    const y = section.y;
    const own = y?.axis.ranges;
    if (!y || !own || !axis) return undefined;
    if (y.followsX && !y.follow)
      return { axis: { ...y.axis, dir: axis.dir }, pieces: axis.pieces, viewport: { start: axis.u0, end: axis.u1 } };
    if (y.follow) {
      const leaderSection = y.followsX ? section : bySection.get(y.follow.leader);
      const leader = y.followsX ? axis : leaderSection && axes.get(leaderSection.address);
      if (leaderSection && leader) {
        const window = y.followsX ? { start: axis.u0, end: axis.u1 } : leaderWindow(leaderSection, leader);
        const made = followAxis({ ...section, axis: y.axis, follow: y.follow }, own, leader, window, config);
        if (made) return { axis: y.axis, pieces: made.pieces, ...(made.viewport ? { viewport: made.viewport } : {}) };
      }
      return { axis: y.axis, pieces: own };
    }
    const v = (input.viewports?.[section.address] ?? input.viewport)?.y;
    return { axis: y.axis, pieces: own, ...(v ? { viewport: v } : {}) };
  };
  const labelSize = fontSize(config.labelFont);
  const sections: SectionLayout[] = [];
  const used = new Set<string>();
  const showing = new Set<string>();
  // what a section's reuse compares besides its own viewport and reads (_spec/README.md "Render", The key)
  const stamp = viewStamp(input, config, dpr, overlay, state);
  const sinks: FrameSinks = { wants, report, theme, used, showing };
  const colors = defaultColors(input.view);
  let y = config.inset;
  // a section whose every track is hidden collapses fully (the application lists them): no gap, no line, no rail
  let anyShown = false;
  for (const section of input.view.sections) {
    const collapsed = collapsedSection(section);
    if (anyShown && !collapsed) y += config.sectionGap;
    if (!collapsed) anyShown = true;
    const axis = axes.get(section.address);
    if (!axis)
      report(
        issueOf(RenderIssues, 'seq.render.axis.pending', 'The section’s length is not known yet', {
          address: section.address,
          hint: RenderIssues.codes['seq.render.axis.pending'].hint,
        }),
      );
    const layout = layoutOf(section);
    const viewport = input.viewports?.[section.address] ?? input.viewport;
    // an unchanged section is not laid out again: its placement and lays are drawn again (_spec/README.md "Render",
    // The key: per section)
    const kept = state.sections.get(section.address);
    // a grid's y: part of its key (a leader it follows may move it alone)
    const gridY = section.y ? yOf(section, axis) : undefined;
    if (kept && reusable(kept, section, layout, stamp, y, axis, viewport, input.data) && sameGridY(kept.gridY, gridY)) {
      kept.replay(sinks);
      if (layout.window) setWindow(layout, kept.layout.placement, section.address, viewport, followed, input, config);
      sections.push(kept.layout);
      y += kept.layout.height;
      continue;
    }
    // the kept layout's placement, read before the cache starts again
    const prevPlacement = kept?.layout?.layout === layout ? kept.layout.placement : undefined;
    const cache = kept ?? new SectionCache();
    cache.begin(input.data, sinks);
    const mine = cache.reads;
    const sectionReport = cache.report;
    const sectionTheme = cache.theme;
    const want = cache.want;
    const sectionTop = y;
    // a ruler rail on a linear section (tracks, not a stack of blocks); a slim row of its own when no ruler is shown
    const railable = layout === TracksLayout && !(section.blocks && section.params.arrange === 'stack');
    const railOn = input.rails?.[section.address] !== undefined && railable && !collapsed;
    // its header line above its first track, in a layout that draws track headers (_spec/README.md "Render", Track
    // headers)
    const sectionHeader =
      layout === TracksLayout && !collapsed
        ? sectionHeaderOf(input, section, sectionId(section), y, config)
        : undefined;
    const top = sectionHeader ? y + config.header.sectionPx : y;
    // a template is expanded here, over the overlay, into a track per node it selects; its source's region in view is
    // wanted, so the nodes it brings arrive
    const shown: ResolvedTrack[] = [];
    // each shown track's index in the section's (a template's instance: −1), what a block's layers are kept by, and
    // the section's track it came from (an instance's template)
    const trackOf: number[] = [];
    const sourceOf: ResolvedTrack[] = [];
    for (let ti = 0; ti < section.tracks.length; ti++) {
      const track = section.tracks[ti];
      if (!track.template) {
        trackOf[shown.length] = ti;
        sourceOf[shown.length] = track;
        shown.push(track);
        continue;
      }
      if (axis) wantTemplate(input, mine, section, track.template, axis, overlay, want);
      const instances = expandTemplate(input.model, overlay, section, track, config.maxTemplateTracks);
      if (instances.length >= config.maxTemplateTracks && templateNodes(track, overlay).length > instances.length)
        sectionReport(
          issueOf(
            RenderIssues,
            'seq.render.template.limit',
            `The template selects more than ${instances.length} nodes`,
            {
              address: track.address,
              hint: RenderIssues.codes['seq.render.template.limit'].hint,
            },
          ),
        );
      for (const t of instances) {
        // an instance whose node's region is out of view is not laid out: its track goes, its hold with it
        if (axis && !instanceInView(input, mine, section, track.template, axis, overlay, t, state)) continue;
        trackOf[shown.length] = -1;
        sourceOf[shown.length] = track;
        shown.push(t);
      }
    }
    const byType = layout.heights?.(config) ?? config.heights;
    const sizedBands = bandHeights(input, state, config, types, dataW, section, shown, layout, false, sectionReport);
    // the kept heights, copied only where this layout changes them (a rows band's)
    let heights = sizedBands.heights as number[];
    const trackContents = sizedBands.contents;
    // each track's padding inside its well, `[top, bottom]` a track: its own, else the config's (_spec/README.md
    // "Render", Tracks layout)
    const pads: number[] = [];
    let padded = false;
    for (const t of shown) {
      const p = t.padding ?? config.trackPadding;
      pads.push(p.top, p.bottom);
      if (p.top > 0 || p.bottom > 0) padded = true;
    }
    const ruler = railOn ? rulerTrack(shown, heights) : -1;
    // the rail draws just below the ruler's bottom border: the ruler's well grows by the rail's row, so it never covers
    // the ticks or labels; only where a rail is drawn live (the canvas with rails on): an export has no row, a figure's
    // rails (navigators included) lie on the border (_spec/README.md "Render", Navigators)
    if (ruler >= 0 && input.figure !== true) {
      pads[2 * ruler + 1] += SlimRailPx;
      padded = true;
    }
    // a rows band that fills (_spec/README.md "Render", The rows layer): the section's room left on the page, at least
    // its default height (a figure, laid out for its content, takes the default)
    // a band of every row in a layout whose bands do not scroll (text: an alignment figure's lines repeat its rows)
    if (layout !== TracksLayout)
      for (let j = 0; j < shown.length; j++) {
        if (shown[j].height !== undefined || rowsFill(shown[j]) === undefined) continue;
        const rows = rowsCount(shown[j], input.model);
        if (rows === undefined) continue;
        if (heights === sizedBands.heights) heights = heights.slice();
        heights[j] = rows * byType.rows + titleRowOf(shown[j], config);
      }
    const fill = shown.findIndex((t) => t.height === undefined && rowsFill(t) === 'fill');
    if (fill >= 0 && !input.figure && layout === TracksLayout) {
      let used = 0;
      for (let j = 0; j < heights.length; j++)
        used += (j !== fill ? heights[j] + config.trackGap : 0) + pads[2 * j] + pads[2 * j + 1];
      if (heights === sizedBands.heights) heights = heights.slice();
      const n = rowsShown(shown[fill], input.model, config.rows.pinInPlace);
      heights[fill] = rowsFillHeight(
        heights[fill],
        height - top - used - config.inset,
        n,
        titleRowOf(shown[fill], config),
      );
    }
    const nb = shown.length;
    const lays: LayerLayout[][] = [];
    const stripBreaks: number[][] = [];
    // one layer's layout into its kept lay (its cache key: the layer's, by strip where a layout has several)
    const layOne = (
      track: ResolvedTrack,
      layer: ResolvedLayer,
      type: LayerType<any, any>,
      cacheKey: string,
      ctx: LayoutContext,
    ): LayerLayout => {
      let prep = state.preps.get(layer);
      if (prep === undefined) state.preps.set(layer, (prep = type.prepare(layer, input.model)));
      const keptLay = state.lays.get(cacheKey);
      const lay = type.layout(prep, ctx, mine, overlay, keptLay?.type === type ? keptLay.lay : undefined);
      state.lays.set(cacheKey, { type, lay });
      cache.addUsed(cacheKey);
      // what a shown layer draws and goes through (a node, an item system)
      if (layer.data) cache.addShowing(layer.data.address, layer.data.system);
      const edges = layer.path?.edges;
      if (edges) for (let e = 0; e < edges.length; e++) cache.addShowing(edges[e].from, edges[e].to);
      return { key: ctx.key, layer, type, lay };
    };
    /** A layer's type: the first of its fallbacks drawn here; none is said. */
    const typeOf = (layer: ResolvedLayer, key: string): LayerType<any, any> | undefined => {
      const type = firstType(layer, types);
      if (!type)
        sectionReport(
          issueOf(RenderIssues, 'seq.render.layer.type', 'No type of the layer’s fallbacks is drawn here', {
            address: key,
            hint: RenderIssues.codes['seq.render.layer.type'].hint,
          }),
        );
      return type;
    };
    /** Band `j`'s label lane, its layers' labels beside sharing it (kept in the state, one a band). */
    const labelLaneOf = (j: number): LabelLane => (state.labelLanes[j] ??= new LabelLane());
    /** A band's own scroll and row zoom (the viewport's `bands`, by track address). */
    const bandOf = (track: ResolvedTrack) => viewport?.bands?.[track.address];
    // one rows layer a section (the cursor has one row order): the first shown, a second said and not drawn
    let firstRows: string | undefined;
    for (const t of shown)
      for (const l of t.layers)
        if (firstRows === undefined && l.type === 'rows' && l.placed && l.params.hidden !== true)
          firstRows = layerKey(t, l);
    // the lays of a strip's bands (_spec/README.md "Render", Section layouts: the lay cache bounded by the strips in view);
    // a `section` layer is laid out once a section, after
    const layStrip = (placement: Placement, i: number, labels?: LabelSink): boolean => {
      let pill = false;
      const strips = placement.strips;
      const a = strips.axes[i];
      const keyed = layout !== TracksLayout || strips.n > 1 || section.blocks !== undefined;
      // a block's strip draws its layers resolved onto the block's system (a path a block)
      const block = section.blocks?.[strips.block[i]];
      const span = stripSpan(strips, i);
      for (let j = 0; j < nb; j++) {
        const k = i * nb + j;
        const out: LayerLayout[] = [];
        lays[k] = out;
        strips.used[k] = 0;
        const h = strips.bandH[k];
        // a band below the page is not laid out: nothing of it is drawn or asked for
        if (!(strips.bandY[k] + strips.dy[i] < height) || (h <= 0 && layout !== TracksLayout)) continue;
        const track = shown[j];
        const band = bandOf(track);
        if (trackContents[j]) labelLaneOf(j).reset();
        const box: Box = { x: a.x0, y: strips.bandY[k], width: a.width, height: h };
        const layers = block && trackOf[j] >= 0 ? block.layers[trackOf[j]] : track.layers;
        for (const layer of layers) {
          if (!layer.placed || layer.params.hidden === true) continue; // a hidden layer is not drawn
          const key = layerKey(track, layer);
          const type = typeOf(layer, key);
          // a layout of two axes lays every layer out a strip at a time, `section` layers too (pairs as dots)
          if (!type || (type.scope === 'section' && !layout.axesY)) continue;
          if (type.name === 'rows' && key !== firstRows) {
            sectionReport(
              issueOf(RenderIssues, 'seq.render.rows.second', 'A section draws one rows layer', {
                address: key,
                hint: RenderIssues.codes['seq.render.rows.second'].hint,
              }),
            );
            continue;
          }
          const ck = colorKey(sourceOf[j], layer);
          const colorIndex = ck === undefined ? undefined : colors.get(ck);
          const ctx: LayoutContext = {
            box,
            axis: a,
            section,
            key,
            dpr,
            config,
            tokens: input.tokens,
            text: input.text,
            settled: input.settled ?? false,
            theme: sectionTheme,
            want,
            report: sectionReport,
            viewWidth: width,
            ...(colorIndex !== undefined ? { colorIndex } : {}),
            ...(trackContents[j] ? { content: trackContents[j], labelLane: labelLaneOf(j) } : {}),
            strip: i,
            span: Array.isArray(layer.params.across) ? acrossSpan(strips, i, shown, layer.params.across) : span,
            place: placement,
            labels: labels ?? 'inside',
            ...(band ? { bandTop: band.bandTop } : {}),
            ...(band?.rowScale !== undefined ? { rowScale: band.rowScale } : {}),
            ...(band?.folds !== undefined
              ? { folds: band.folds }
              : input.folds?.[track.address]
                ? { folds: input.folds[track.address] }
                : {}),
            ...(input.images ? { images: input.images } : {}),
            ...(strips.ys[i] ? { axisY: strips.ys[i] } : {}),
            // a chip's pill: its label elided to its slot (its positions' px on its strip: the whole chip's)
            ...(block && strips.chip[i]
              ? { chip: chipText(input.text, config, block.label, a.scale * blockLength(block)) }
              : {}),
            ...(grants && type.name === 'rows' && grants.has(key) ? { rowColumns: grants.get(key)! } : {}),
            ...(narrow && narrow.tree !== header.tree && type.name === 'rows' ? { treePx: narrow.tree } : {}),
          };
          // kept by the strip's first unit, and its block's (a stack's blocks all start at 0)
          let at = block ? `${key}#${strips.block[i]}@${strips.first[i]}` : `${key}@${strips.first[i]}`;
          // a grid over blocks: a strip a pair of blocks, kept by its y block too
          if (block && strips.ys[i]) at += `|${strips.yBlock[i]}@${strips.yFirst[i]}`;
          const l = layOne(track, layer, type, keyed ? at : key, ctx);
          const u = type.used ? type.used(l.lay) : h;
          if (u > strips.used[k]) strips.used[k] = u;
          out.push(l);
          if (type.name === 'letters' && strips.chip[i]) {
            const ll = l.lay as LettersLay;
            if (ll.isChip && ll.chip !== '' && ll.n > 0) pill = true;
          }
        }
      }
      return pill;
    };
    /** The `section` layers: once a section, over its box on the page, the placement their reader. */
    const laySection = (placement: Placement, a: AxisTransform): LayerLayout[] => {
      let out = NoSectionLays;
      if (layout.axesY) return out;
      for (let band = 0; band < shown.length; band++)
        for (const layer of shown[band].layers) {
          const track = shown[band];
          if (!layer.placed || layer.params.hidden === true) continue;
          const key = layerKey(track, layer);
          const type = firstType(layer, types);
          if (type?.scope !== 'section') continue;
          const ctx: LayoutContext = {
            box: placement.box,
            axis: a,
            section,
            key,
            dpr,
            config,
            tokens: input.tokens,
            text: input.text,
            settled: input.settled ?? false,
            theme: sectionTheme,
            want,
            report: sectionReport,
            span: placement.box,
            place: placement,
            labels: 'inside',
            band,
          };
          if (out === NoSectionLays) out = [];
          out.push(layOne(track, layer, type, `${key}@section`, ctx));
        }
      return out;
    };
    const sectionBox: Box = {
      x: dataX,
      y: top,
      width: Math.max(1, width - dataX - config.inset - right),
      height: Math.max(0, height - top),
    };
    let placement: Placement;
    let sectionLays: LayerLayout[] = NoSectionLays;
    if (axis) {
      const ctx: PlaceContext = {
        width: layout.name === 'text' ? textWidth : width,
        height,
        dpr,
        config,
        tokens: input.tokens,
        text: input.text,
        model: input.model,
        box: sectionBox,
        axis,
        bands: shown,
        heights,
        ...(padded ? { pads } : {}),
        layStrip,
        report: sectionReport,
        ...(input.figure ? { figure: true } : {}),
        ...(input.windows ? { windows: true } : {}),
        ...(gridY ? { y: gridY } : {}),
      };
      placement = layout.place(section, ctx, viewport, prevPlacement);
      if (layout.window) setWindow(layout, placement, section.address, viewport, followed, input, config);
      for (let i = 0; i < placement.strips.n; i++) stripBreaks.push(placement.strips.axes[i].breaks([]));
      sectionLays = laySection(placement, axis);
    } else placement = pendingPlacement(top, sectionBox, heights, pads, config);
    // the section's tracks as bands of its first strip (on the page: tracks)
    const tracks: TrackLayout[] = [];
    const strips = placement.strips;
    let by = top;
    for (let j = 0; j < nb; j++) {
      if (j) by += config.trackGap;
      const track = shown[j];
      // a grid's track is its plot (over blocks: its cells together), its well the plot's grown by its padding
      const plot = layout.name === 'grid' && strips.n ? (placement as GridPlacement) : undefined;
      const pt = plot ? plot.padTop : strips.n ? strips.padT[j] : pads[2 * j];
      const pb = plot ? plot.padBottom : strips.n ? strips.padB[j] : pads[2 * j + 1];
      const box: Box = {
        x: dataX,
        y: plot ? plot.plotTop : strips.n ? strips.bandY[j] + strips.dy[0] : by + pt,
        width: sectionBox.width,
        // a layout of two axes overlays its bands on the plot
        height: plot ? plot.plotHeight : heights[j],
      };
      const well: Box =
        pt > 0 || pb > 0 ? { x: box.x, y: box.y - pt, width: box.width, height: box.height + pt + pb } : box;
      by += pt + heights[j] + pb;
      const layers = lays[j] ?? NoLays;
      // its header's line: its states (the "n filtered" chip among them) and slots; the label's room what its states
      // beside it leave (states on a line of their own leave it whole)
      const head =
        layout === TracksLayout
          ? trackHeaderOf(input, section, sectionId(section), track, layers, sectionLays, well, config)
          : undefined;
      const room = head ? head.header.label.width : config.labelWidth;
      const labelText = elidedLabel(state, input.text, config.labelFont, track.label, room);
      tracks.push({
        track,
        box,
        well,
        labelText,
        labelWidth: input.text.measure(config.labelFont, labelText),
        // centred on its well (its padding with it), not on its data
        labelY: head ? head.labelY : well.y + well.height / 2 + labelSize * 0.35,
        layers,
        ...(head ? { header: head.header } : {}),
      });
    }
    y = top + placement.box.height;
    // the slim rail's row, under the tracks
    const slimTop = y;
    if (railOn && ruler < 0) y += SlimRailPx;
    // the column the emphasis spans with `emphasis.selectionShape: 'column'`: the tracks laid out on the page and the
    // gaps between them
    let bandEnd = -1;
    for (const t of tracks) if (t.box.y < height) bandEnd = Math.max(bandEnd, t.box.y + t.box.height);
    const first = tracks[0]?.box;
    const band =
      first && bandEnd > first.y
        ? { x: first.x, y: first.y, width: first.width, height: bandEnd - first.y }
        : undefined;
    // each band's box and each strip's band column, in the strip's frame (tracks: its first strip's are the page's)
    const boxes: Box[] = [];
    const columns: (Box | undefined)[] = [];
    for (let i = 0; i < strips.n; i++) {
      const a = strips.axes[i];
      for (let j = 0; j < nb; j++) {
        const k = i * nb + j;
        boxes.push(
          i === 0 && layout.fillBands && !layout.axesY
            ? tracks[j].box
            : { x: a.x0, y: strips.bandY[k], width: a.width, height: strips.bandH[k] },
        );
      }
      columns.push(i === 0 && layout.fillBands && !layout.axesY ? band : stripColumn(strips, i, height));
    }
    const text = tracks.length ? directionMark(section, config) : undefined;
    const mark = text
      ? {
          text,
          x: config.inset + tracks[0].labelWidth + 6,
          width: input.text.measure(config.labelFont, text),
        }
      : undefined;
    let rail: Rail | undefined;
    if (railOn && axis) {
      rail = kept?.layout?.rail ?? new Rail();
      const own = ownPieces(section, input.model, input.data) ?? axis.pieces;
      placeRail(
        rail,
        section,
        sectionId(section),
        axis,
        own,
        ruler >= 0 ? tracks[ruler].box : undefined,
        slimTop,
        viewport,
        config,
        input.figure === true,
      );
    }
    // each rows layer's rail, upright at the right edge (beside its last numbers), along its band's scrolled rows
    let rowsRails: Rail[] | undefined;
    if (input.rails && layout === TracksLayout)
      for (let j = 0; j < tracks.length; j++) {
        const t = tracks[j];
        if (input.rails[t.track.address] === undefined) continue;
        const l = t.layers.find((x) => x.type.name === 'rows');
        const lay = l?.lay as RowsLay | undefined;
        if (!lay || lay.figure || lay.total === undefined) continue;
        const r = kept?.layout?.rowsRails?.[rowsRails?.length ?? 0] ?? new Rail();
        placeRowsRail(
          r,
          lay,
          t.track.address,
          sectionId(section),
          width - config.inset - SlimRailPx,
          input.figure === true,
        );
        (rowsRails ??= []).push(r);
      }
    const sectionLayout: SectionLayout = {
      section,
      id: sectionId(section),
      layout,
      placement,
      ...(axis ? { axis } : {}),
      ...(band ? { band } : {}),
      ...(rail ? { rail } : {}),
      ...(rowsRails ? { rowsRails } : {}),
      y: sectionTop,
      height: y - sectionTop,
      tracks,
      ...(sectionHeader ? { header: sectionHeader } : {}),
      breaks: axis && band ? (stripBreaks[0] ?? NoBreaks) : NoBreaks,
      ...(mark ? { mark } : {}),
      lays,
      stripBreaks,
      boxes,
      columns,
      sectionLays,
      pageHeight: height,
      tokens: input.tokens,
      config,
    };
    sections.push(sectionLayout);
    cache.layout = sectionLayout;
    cache.section = section;
    cache.stamp = stamp;
    cache.y = sectionTop;
    cache.axis = axis;
    cache.gridY = gridY;
    // a copy: a caller may hand the same record again with new numbers (the viewer reuses its own)
    cache.viewport = viewport ? { ...viewport } : undefined;
    state.sections.set(section.address, cache);
  }
  // a section gone from the view lets its cache go
  if (state.sections.size > input.view.sections.length) {
    const here = new Set(input.view.sections.map((s) => s.address));
    for (const k of [...state.sections.keys()]) if (!here.has(k)) state.sections.delete(k);
  }
  // a layer gone from the view lets its last layout go
  for (const key of [...state.lays.keys()]) if (!used.has(key)) state.lays.delete(key);
  // open scales: settled, the union of this layout's extents is held (laid out once more with it where it moved);
  // moving, a scale seen for the first time takes this layout's extent for the next frame
  if (input.domains) {
    // the domains are another layout's: none held here
  } else if (input.settled) {
    let moved = extents.size !== state.domains.size;
    for (const [k, e] of extents) {
      const h = state.domains.get(k);
      if (!h || h[0] !== e[0] || h[1] !== e[1]) moved = true;
    }
    if (moved) {
      state.domains = extents;
      state.domainsVersion++;
      if (!again) return layoutView(input, state, true);
    }
  } else
    for (const [k, e] of extents)
      if (!state.domains.has(k)) {
        state.domains.set(k, e);
        state.domainsVersion++;
      }
  for (const i of overlay.issues) report(i);
  const shownNodes = config.nodes.holdShown && overlay.nodes.length ? overlay.holds(showing) : [];
  return {
    width,
    height,
    dpr,
    tokens: input.tokens,
    config,
    sections,
    wants,
    issues,
    shown: shownNodes,
    ...(narrow ? { narrow } : {}),
  };
}

const NoBreaks: readonly number[] = [];

/**
 * A placement's window (a ring's), out of its key: set on a kept placement too. Active while something follows the
 * section or it follows (its window then its leader's range); shown when active or `rings.window: 'always'`, in a
 * figure only with `windows`. A function of the module, not a closure a section a frame.
 */
function setWindow(
  layout: SectionLayoutType,
  placement: Placement,
  at: string,
  viewport: Viewport | undefined,
  followed: Map<string, Viewport> | undefined,
  input: LayoutInput,
  config: RenderConfig,
): void {
  const own = followed?.get(at);
  let active = own !== undefined;
  const sections = input.view.sections;
  for (let i = 0; i < sections.length && !active; i++) if (sections[i].follow?.leader === at) active = true;
  const shown = (active || config.rings.window === 'always') && (!input.figure || input.windows === true);
  layout.window!(placement, own ?? viewport, active, shown);
}
