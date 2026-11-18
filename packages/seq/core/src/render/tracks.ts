/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { fitLevel, type Want } from '@molstar/seqx-core/data';
import type { LocationPart } from '@molstar/seqx-core/query';
import { type OverlayNode, projectRanges, type StaticModel } from '@molstar/seqx-core/resolve';
import {
  expandTemplate,
  followRange,
  templateNodes,
  type ResolvedLayer,
  type ResolvedSection,
  type ResolvedTrack,
  type ResolvedView,
  type TrackTemplate,
} from '@molstar/seqx-core/view';
import { AxisTransform, type Viewport } from './axis.js';
import { DefaultRenderConfig, fontSize, type RenderConfig } from './config.js';
import {
  Anchor,
  type DrawLists,
  drawPass,
  NoColor,
  type PassName,
  Passes,
  type PassWriter,
  passDraws,
} from './draw-lists.js';
import { RenderIssues } from './issues.js';
import type { Box, Emphasis, Hit, LayerType, LayoutContext, Overlay, StoreReader } from './layer-type.js';
import { EmptyOverlay } from './layer-type.js';
import { LettersType } from './letters.js';
import { RangesType } from './ranges.js';
import type { TextMetrics } from './text.js';
import type { SeqTokens } from './tokens.js';
import { ValuesType } from './values.js';

// The tracks layout (_spec/README.md "Render", Tracks layout): sections stacked, each with its axis transform from its
// system and range (a follower's from its leader's view), tracks with heights in px under one label column, each
// track's layers in order over its box. Tracks below the page are not laid out: the draw lists are bounded by the
// screen. The layout keeps each layer's last layout, so buffers and lanes carry over frame to frame.

/** The layer types this target draws, by name: a layer draws the first of its fallback chain found here. */
export const BuiltInLayerTypes: Readonly<Record<string, LayerType<any, any>>> = {
  letters: LettersType,
  ranges: RangesType,
  values: ValuesType,
};

export interface LayerLayout {
  readonly key: string;
  readonly layer: ResolvedLayer;
  readonly type: LayerType<any, any>;
  readonly lay: unknown;
}

export interface TrackLayout {
  readonly track: ResolvedTrack;
  /** The data box. */
  readonly box: Box;
  /** The label's measured width and its baseline (drawing never measures). */
  readonly labelWidth: number;
  readonly labelY: number;
  readonly layers: readonly LayerLayout[];
}

export interface SectionLayout {
  readonly section: ResolvedSection;
  /** Undefined while the section's length is unknown (pending). */
  readonly axis?: AxisTransform;
  readonly y: number;
  readonly height: number;
  readonly tracks: readonly TrackLayout[];
}

export interface ViewLayout {
  readonly width: number;
  readonly height: number;
  readonly dpr: number;
  readonly tokens: SeqTokens;
  readonly config: RenderConfig;
  readonly sections: readonly SectionLayout[];
  /** What the layers need from the store now. */
  readonly wants: readonly Want[];
  /** The issues of this layout (a color it cannot read, a layer it cannot draw), each once. */
  readonly issues: readonly Issue[];
  /** The dynamic nodes this layout shows, with what they name (with `nodes.holdShown`): its owner holds them. */
  readonly shown: readonly OverlayNode[];
}

export interface LayoutInput {
  readonly model: StaticModel;
  readonly view: ResolvedView;
  readonly size: { readonly width: number; readonly height: number };
  /** The window of display units of every section that does not follow another; all of each by default. */
  readonly viewport?: Viewport;
  readonly dpr?: number;
  readonly tokens: SeqTokens;
  readonly text: TextMetrics;
  readonly config?: RenderConfig;
  readonly data: StoreReader;
  readonly overlay?: Overlay;
  readonly settled?: boolean;
  readonly types?: Readonly<Record<string, LayerType<any, any>>>;
}

/** What a layout keeps between frames: each layer's prepared form (per resolved layer) and its last layout. */
export class LayoutState {
  readonly preps = new WeakMap<ResolvedLayer, unknown>();
  readonly lays = new Map<string, { type: LayerType<any, any>; lay: unknown }>();
}

/** A layer's key: its address, else its track's and its index. */
export function layerKey(track: ResolvedTrack, layer: ResolvedLayer): string {
  return layer.address ?? `${track.address}:${layer.index}`;
}

/** A section's own pieces: its axis's, or the length its source gave once open; undefined while unknown. */
function ownPieces(section: ResolvedSection, model: StaticModel, data: StoreReader): number[] | undefined {
  if (section.axis.ranges) return [...section.axis.ranges];
  const system = model.systems.get(section.axis.system);
  if (!system?.source) return undefined;
  const info = data
    .coverage(system.source)
    ?.systems.find((s) => s.name === system.name || s.aliases?.includes(system.name));
  return info?.length === undefined ? undefined : [0, info.length];
}

/**
 * A follower's pieces and window from its leader's view, continuously: the images of the leader's window edges
 * (`point`) within a run, so a follower moves with its leader as the viewer's does; the hull (or arc) of the leader's
 * positions in view where an edge maps nowhere, the leader's view has two pieces, or the window wraps a cyclic origin.
 */
function followAxis(section: ResolvedSection, leader: AxisTransform, config: RenderConfig) {
  const follow = section.follow!;
  if (leader.offsets.length === 1) {
    const at = (u: number) => leader.pieces[0] + u;
    const edge = (x: number) => follow.point(x) ?? follow.point(x - 1e-9);
    const a = edge(at(leader.u0));
    const b = edge(at(leader.u1));
    const { cyclic, length } = section.axis;
    if (
      a !== undefined &&
      b !== undefined &&
      a !== b &&
      !(cyclic && length !== undefined && Math.abs(b - a) > length / 2)
    ) {
      const lo = Math.min(a, b);
      const hi = Math.max(a, b);
      const s = Math.floor(lo);
      return { pieces: [s, Math.ceil(hi)], viewport: { start: lo - s, end: hi - s } };
    }
  }
  const pieces = followRange(section, leader.visible([]), { arc: config.followArc });
  return pieces ? { pieces } : undefined;
}

export function layoutView(input: LayoutInput, state = new LayoutState()): ViewLayout {
  const config = input.config ?? DefaultRenderConfig;
  const dpr = input.dpr ?? 1;
  const types = input.types ?? BuiltInLayerTypes;
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
  const dataW = Math.max(1, width - dataX - config.inset);
  // axes first, leaders before their followers (a follower may come first in the view)
  const bySection = new Map(input.view.sections.map((s) => [s.address, s]));
  const axes = new Map<string, AxisTransform | undefined>();
  const visiting = new Set<string>();
  const axisOf = (section: ResolvedSection): AxisTransform | undefined => {
    if (axes.has(section.address)) return axes.get(section.address);
    const at = section.address;
    let made: { pieces: number[]; viewport?: Viewport } | undefined;
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
      } else {
        visiting.add(at);
        const leader = axisOf(leaderSection);
        visiting.delete(at);
        if (leader) made = followAxis(section, leader, config);
      }
    }
    if (!made) {
      const own = ownPieces(section, input.model, input.data);
      if (own) made = { pieces: own, ...(input.viewport ? { viewport: input.viewport } : {}) };
    }
    const axis = made?.pieces.length
      ? new AxisTransform(section.axis.system, made.pieces, dataX, dataW, dpr, made.viewport)
      : undefined;
    if (!axes.has(at)) axes.set(at, axis);
    return axes.get(at);
  };
  for (const section of input.view.sections) axisOf(section);
  const labelSize = fontSize(config.labelFont);
  const sections: SectionLayout[] = [];
  const used = new Set<string>();
  const showing = new Set<string>();
  let y = config.inset;
  for (const section of input.view.sections) {
    if (sections.length) y += config.sectionGap;
    const axis = axes.get(section.address);
    if (!axis)
      report(
        issueOf(RenderIssues, 'seq.render.axis.pending', 'The section’s length is not known yet', {
          address: section.address,
          hint: RenderIssues.codes['seq.render.axis.pending'].hint,
        }),
      );
    const top = y;
    const tracks: TrackLayout[] = [];
    // a template is expanded here, over the overlay, into a track per node it selects; its source's region in view is
    // wanted, so the nodes it brings arrive
    const shown: ResolvedTrack[] = [];
    for (const track of section.tracks) {
      if (!track.template) {
        shown.push(track);
        continue;
      }
      if (axis) wantTemplate(input, section, track.template, axis, overlay, wants);
      const instances = expandTemplate(input.model, overlay, section, track, config.maxTemplateTracks);
      if (instances.length >= config.maxTemplateTracks && templateNodes(track, overlay).length > instances.length)
        report(
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
      shown.push(...instances);
    }
    for (const track of shown) {
      if (tracks.length) y += config.trackGap;
      const h = trackHeight(track, config);
      const box: Box = { x: dataX, y, width: dataW, height: h };
      const layers: LayerLayout[] = [];
      // a track below the page is not laid out: nothing of it is drawn or asked for
      if (axis && y < height)
        for (const layer of track.layers) {
          if (!layer.placed) continue;
          const key = layerKey(track, layer);
          const type = layer.fallbacks.map((t) => types[t]).find((t) => t !== undefined);
          if (!type) {
            report(
              issueOf(RenderIssues, 'seq.render.layer.type', 'No type of the layer’s fallbacks is drawn here', {
                address: key,
                hint: RenderIssues.codes['seq.render.layer.type'].hint,
              }),
            );
            continue;
          }
          let prep = state.preps.get(layer);
          if (prep === undefined) state.preps.set(layer, (prep = type.prepare(layer, input.model)));
          const kept = state.lays.get(key);
          const ctx: LayoutContext = {
            box,
            axis,
            section,
            key,
            dpr,
            config,
            tokens: input.tokens,
            text: input.text,
            settled: input.settled ?? false,
            want: (w) => wants.push(w),
            report,
          };
          const lay = type.layout(prep, ctx, input.data, overlay, kept?.type === type ? kept.lay : undefined);
          state.lays.set(key, { type, lay });
          used.add(key);
          // what a shown layer draws and goes through (a node, an item system)
          if (layer.data) showing.add(layer.data.address).add(layer.data.system);
          for (const e of layer.path?.edges ?? []) showing.add(e.from).add(e.to);
          layers.push({ key, layer, type, lay });
        }
      tracks.push({
        track,
        box,
        labelWidth: input.text.measure(config.labelFont, track.label),
        labelY: box.y + box.height / 2 + labelSize * 0.35,
        layers,
      });
      y += h;
    }
    sections.push({ section, ...(axis ? { axis } : {}), y: top, height: y - top, tracks });
  }
  // a layer gone from the view lets its last layout go
  for (const key of [...state.lays.keys()]) if (!used.has(key)) state.lays.delete(key);
  for (const i of overlay.issues) report(i);
  const shownNodes = config.nodes.holdShown && overlay.nodes.length ? overlay.holds(showing) : [];
  return { width, height, dpr, tokens: input.tokens, config, sections, wants, issues, shown: shownNodes };
}

/**
 * A template source's region in view, on the section's system as the source names it, at the level that fits; and for
 * a template below its source, each nested source's systems' regions in view (through the overlay's path from the
 * section), so a nested source that sends its nodes by region is asked for them.
 */
function wantTemplate(
  input: LayoutInput,
  section: ResolvedSection,
  template: TrackTemplate,
  axis: AxisTransform,
  overlay: Overlay,
  out: Want[],
) {
  const pieces = axis.visible([]);
  const coverage = input.data.coverage(template.source);
  const name = input.model.systems.get(section.axis.system)?.name;
  const info = coverage?.systems.find((s) => s.name === name || (name !== undefined && s.aliases?.includes(name)));
  if (coverage && info) {
    const level = fitLevel(coverage.levels, 1 / axis.devicePx);
    for (let i = 0; i + 1 < pieces.length; i += 2)
      out.push({
        source: template.source,
        system: info.name,
        start: Math.floor(pieces[i]),
        end: Math.ceil(pieces[i + 1]),
        level,
      });
  }
  if (template.depth < 1) return;
  for (const nested of overlay.nested(new Map([[template.source, template.depth]]))) {
    const cov = input.data.coverage(nested.address);
    for (const s of cov?.systems ?? []) {
      const target = overlay.systemNamed(s.name, nested.parent);
      if (!target) continue;
      const path =
        target.address === section.axis.system ? undefined : overlay.path(section.axis.system, target.address).path;
      if (!path && target.address !== section.axis.system) continue;
      const ranges = HalfOpenRanges.normalize(path ? projectRanges(path, pieces) : pieces);
      const ratio = path?.step === 1 ? 3 : path?.step === -1 ? 1 / 3 : 1;
      const level = fitLevel(cov!.levels, ratio / axis.devicePx);
      for (let i = 0; i + 1 < ranges.length; i += 2)
        out.push({
          source: nested.address,
          system: s.name,
          start: Math.floor(ranges[i]),
          end: Math.ceil(ranges[i + 1]),
          level,
        });
    }
  }
}

function trackHeight(track: ResolvedTrack, config: RenderConfig): number {
  if (track.height !== undefined) return track.height;
  let h = 0;
  for (const l of track.layers) h = Math.max(h, config.heights[l.type as keyof RenderConfig['heights']] ?? 0);
  return h || config.heights.ranges;
}

/** Draws the layout's passes (all, or those named, e.g. the emphasis alone) into the lists, cleared first. */
export function drawView(
  layout: ViewLayout,
  lists: DrawLists,
  emphasis?: Emphasis,
  passes: readonly PassName[] = Passes,
): void {
  lists.reset(passes);
  lists.width = layout.width;
  lists.height = layout.height;
  // indexed loops: an iterator a loop is an allocation where the compiler keeps it (nothing allocated once warm)
  for (let pi = 0; pi < passes.length; pi++) {
    const pass = passes[pi];
    drawPass(ViewDraws, layout, lists.writer(pass), emphasis);
  }
}

/** The view's own draw of each pass, through `drawPass` like a layer type's. */
const ViewDraws = passDraws<ViewLayout, Emphasis | undefined>({
  background: drawBackground,
  data: drawTracks,
  highlight: drawTracks,
  selection: drawTracks,
  hover: drawTracks,
  cursor: drawTracks,
  labels: drawTracks,
});

/** The background pass: the page, and each track's box. */
function drawBackground(layout: ViewLayout, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { tokens } = layout;
  out.fill(0, 0, layout.width, layout.height, tokens.colors['msui-bg']);
  for (let si = 0; si < layout.sections.length; si++) {
    const s = layout.sections[si];
    for (let ti = 0; ti < s.tracks.length; ti++) {
      const { box } = s.tracks[ti];
      if (box.y >= layout.height) continue;
      out.fill(box.x, box.y, box.width, box.height, tokens.colors['seq-track']);
    }
  }
}

/** Any other pass (the writer's): each track's own marks of it (labels, selection, cursor), then its layers'. */
function drawTracks(layout: ViewLayout, out: PassWriter, emphasis: Emphasis | undefined) {
  const { tokens, config } = layout;
  const pass = out.pass.name;
  for (let si = 0; si < layout.sections.length; si++) {
    const s = layout.sections[si];
    for (let ti = 0; ti < s.tracks.length; ti++) {
      const t = s.tracks[ti];
      if (t.box.y >= layout.height) continue;
      const { box } = t;
      switch (pass) {
        case 'labels':
          drawLabel(out, t, tokens, config);
          break;
        case 'selection':
          if (s.axis && emphasis?.selection) drawSelection(out, s.axis, box, emphasis, tokens, config);
          break;
        case 'cursor':
          if (s.axis && emphasis?.cursor?.system === s.axis.system)
            drawCursor(out, s.axis, box, emphasis.cursor, tokens, config);
          break;
      }
      if (!t.layers.length) continue;
      const clip = pass === 'data';
      if (clip) out.clip(box.x, box.y, box.width, box.height);
      for (let li = 0; li < t.layers.length; li++) {
        const l = t.layers[li];
        l.type.draw(l.lay, pass, out, emphasis);
      }
      if (clip) out.unclip();
    }
  }
}

// The label and the cursor take the layout, not its numbers: a number passed to a call that is not inlined is boxed
// (nothing allocated once warm).

/** A track's label in the label column, faded past it. */
function drawLabel(out: PassWriter, t: TrackLayout, tokens: SeqTokens, config: RenderConfig) {
  const g = out.args;
  g[0] = config.inset;
  g[1] = t.labelY;
  g[2] = t.labelWidth;
  g[3] = config.labelWidth;
  out.textArgs(t.track.label, Anchor.start, config.labelFont, tokens.colors['msui-text-muted']);
}

/** The cursor's outline around its cell, thicker than the hover's. */
function drawCursor(
  out: PassWriter,
  axis: AxisTransform,
  box: Box,
  cursor: { readonly position: number },
  tokens: SeqTokens,
  config: RenderConfig,
) {
  if (!axis.cell(cursor.position)) return;
  const w = config.outline.cursor;
  // `snap` inline, the rect through `args`: an emphasis redraw allocates nothing once warm
  const dpr = axis.dpr;
  const x0 = Math.round(axis.sx0 * dpr) / dpr;
  const x1 = Math.round(axis.sx1 * dpr) / dpr;
  const g = out.args;
  g[0] = x0 + w / 2;
  g[1] = box.y + w / 2;
  g[2] = Math.max(x1 - x0 - w, 0.5);
  g[3] = box.height - w;
  g[4] = NoColor;
  g[5] = tokens.colors['msui-accent'];
  g[6] = w;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/**
 * The selection on a section's system: an accent tint over the track and an accent outline around each range. Its part
 * is found with an indexed loop, its ranges placed with `spanIn`, snapped and tinted inline (`snap`, `withAlpha`), its
 * rects drawn through `args`: an emphasis redraw allocates nothing once warm.
 */
function drawSelection(
  out: PassWriter,
  axis: AxisTransform,
  box: Box,
  emphasis: Emphasis,
  tokens: SeqTokens,
  config: RenderConfig,
) {
  const sel = emphasis.selection;
  if (!sel) return;
  let part: LocationPart | undefined;
  for (let i = 0; i < sel.parts.length && !part; i++) if (sel.parts[i].system === axis.system) part = sel.parts[i];
  if (!part) return;
  const accent = tokens.colors['msui-accent'];
  const tint = ((accent & 0xffffff00) | Math.round(0.18 * 255)) >>> 0;
  const w = config.outline.selection;
  const dpr = axis.dpr;
  const n = axis.offsets.length;
  const g = out.args;
  for (let i = 0; i + 1 < part.ranges.length; i += 2)
    for (let k = 0; k < n; k++) {
      if (!axis.spanIn(k, part.ranges[i], part.ranges[i + 1])) continue;
      const x0 = Math.round(axis.sx0 * dpr) / dpr;
      const x1 = Math.max(Math.round(axis.sx1 * dpr) / dpr, x0 + 1 / dpr);
      g[0] = x0;
      g[1] = box.y;
      g[2] = x1 - x0;
      g[3] = box.height;
      g[4] = tint;
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
      g[0] = x0 + w / 2;
      g[1] = box.y + w / 2;
      g[2] = Math.max(x1 - x0 - w, 0.5);
      g[3] = box.height - w;
      g[4] = NoColor;
      g[5] = accent;
      g[6] = w;
      out.rectArgs();
    }
}

/** What is under a point: the topmost layer's hit in the track there. */
export function hitView(layout: ViewLayout, x: number, y: number): Hit | undefined {
  for (const s of layout.sections)
    for (const t of s.tracks) {
      const { box } = t;
      if (y < box.y || y >= box.y + box.height || x < box.x || x >= box.x + box.width) continue;
      // the topmost layer with an element under the point, else the topmost with a position
      let position: Hit | undefined;
      for (let i = t.layers.length - 1; i >= 0; i--) {
        const l = t.layers[i];
        const h = l.type.hit(l.lay, x, y);
        if (h?.element) return h;
        position ??= h;
      }
      return position;
    }
  return undefined;
}
