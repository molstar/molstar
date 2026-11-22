/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node, Value } from '@molstar/core/repr/tree';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedNavigator, ResolvedSection, ResolvedView } from '@molstar/seqx-core/view';
import type { RenderConfig } from './config.js';

// A navigator's picture (_spec/README.md "Render", Navigators, Pictures): a registry entry (`navigator-picture`) that
// says what a navigator draws of its system whole as spec tracks, laid out by the navigator's own renderer through the
// layout interface. The four built-ins (`axis`, `tracks`, `bands`, `rows`) are entries written to the same contract;
// a host's or an extension's is named `<id>:<name>`. The window is the navigator's, whatever the picture.

/** What a picture is given: the navigator, the section it binds first (its node too), its height, the config. */
export interface NavigatorPictureContext {
  readonly model: StaticModel;
  readonly view: ResolvedView;
  /** The navigator, its params as written (`navigator.params`: a namespaced picture's own, plain). */
  readonly navigator: ResolvedNavigator;
  /** Its first bound section, and that section's node in the spec (its tracks as written). */
  readonly bound: ResolvedSection;
  readonly boundNode: Node | undefined;
  /** The height its tracks may fill (CSS px; the zoom strip's track for `zoom`). */
  readonly height: number;
  /** True: the zoom strip's picture is asked for (undefined back: the axis there). */
  readonly zoom: boolean;
  readonly config: RenderConfig;
  /** The bound section's tracks the navigator leaves out (what no summary serves), by address. */
  readonly left: ReadonlySet<string>;
}

/** A picture: its tracks, in order; `rows` the id of the track copying the bound section's rows track (its rows fitted). */
export interface NavigatorPictureTracks {
  readonly tracks: readonly Node[];
  readonly rows?: string;
}

/** A navigator's picture (`navigator.content` names it). */
export interface NavigatorPicture {
  /** A built-in's bare name, or `<id>:<name>`. */
  readonly name: string;
  /** The axes it can cover: a navigator over `rows` whose picture cannot draws the `rows` picture. */
  readonly axes: readonly ('positions' | 'rows' | 'both')[];
  /** Its default thickness (CSS px), or a miniature's box `[width, height]`. */
  readonly size: number | readonly [number, number];
  /** Its tracks fill the navigator's height: the picture is made again for another height. */
  readonly fills?: boolean;
  picture(ctx: NavigatorPictureContext): NavigatorPictureTracks | undefined;
  /** The hover readout at a position of the navigator's system; the position's label when absent or undefined. */
  readout?(
    position: number,
    ctx: { readonly model: StaticModel; readonly navigator: ResolvedNavigator },
  ): string | undefined;
}

/** A navigator's default thickness by its form (CSS px), and a miniature's box. */
export const NavigatorSizes = { axis: 20, tracks: 48, bands: 24, rows: 96, both: [160, 96] } as const;

const isRuler = (t: Node) => (t.children ?? []).every((l) => l.kind !== 'layer' || l.params?.type === 'ruler');

/** Whether a track has a rows layer. */
export const hasRows = (t: Node) => (t.children ?? []).some((l) => l.kind === 'layer' && l.params?.type === 'rows');

/** One track of one layer filling `height` (no label). */
const fill = (id: string, height: number, layer: Readonly<Record<string, Value>>): Node => ({
  kind: 'track',
  id,
  params: { label: '', height: Math.max(1, height) },
  children: [{ kind: 'layer', id, params: layer }],
});

/**
 * A track copied small: its height scaled (`navigator.trackScale`, at least `trackMinPx`), no label, a rows layer
 * with no header (its tree, numbers and row columns are the view's: a picture's narrow column would hide the columns
 * and draw their chip over its top) filling what is left.
 */
function smallTrack(t: Node, config: RenderConfig): Node {
  const h = t.params?.height;
  const params: Record<string, Value> = { ...t.params, label: '' };
  if (typeof h === 'number') params.height = Math.max(config.navigator.trackMinPx, h * config.navigator.trackScale);
  const children = (t.children ?? []).map((l) => (l.kind === 'layer' && l.params?.type === 'rows' ? smallRows(l) : l));
  return { ...t, params, children };
}

/** A rows layer with no header: no tree, no numbers, no row columns, filling its track. */
function smallRows(l: Node): Node {
  const params: Record<string, Value> = { ...l.params, tree_width: 0, numbers: 'none', height: 'fill' };
  delete params.columns;
  return { ...l, params };
}

/** A small track's height: its own, else its first layer type's default, both at the picture's scale. */
export function smallHeight(t: Node, config: RenderConfig): number {
  const h = t.params?.height;
  if (typeof h === 'number') return h;
  const type = t.children?.find((l) => l.kind === 'layer')?.params?.type;
  const heights = config.heights as Readonly<Record<string, number>>;
  const base = (typeof type === 'string' ? heights[type] : undefined) ?? config.heights.ranges;
  return Math.max(config.navigator.trackMinPx, base * config.navigator.trackScale);
}

/** The thickness of an axis picture: the navigator's `size` (a number), else the default. */
function axisThickness(nav: ResolvedNavigator): number {
  return typeof nav.size === 'number' ? nav.size : NavigatorSizes.axis;
}

/** The axis: a ruler filling the navigator (its ticks and numbers at the picture's scale); the fixed fallback. */
export const AxisPicture: NavigatorPicture = {
  name: 'axis',
  axes: ['positions'],
  size: NavigatorSizes.axis,
  // a zoom strip over a genome's blocks counts in bases (`20 Mb`), as a genome's rulers do
  picture: (c) => ({
    tracks: [
      // a picture not made for the navigator's height (`fills` absent: height 0) takes its thickness, never a track of
      // one pixel, whose ticks and labels fall outside it (a blank well); over blocks the label row is the section's
      fill(
        'axis',
        c.height > 0 || c.navigator.blocks ? c.height : axisThickness(c.navigator),
        c.zoom && c.navigator.blocks ? { type: 'ruler', units: 'bases' } : { type: 'ruler' },
      ),
    ],
  }),
};

/** A chromosome's bands (an ideogram), filling the navigator and its zoom strip; none without `bands`. */
export const BandsPicture: NavigatorPicture = {
  name: 'bands',
  axes: ['positions'],
  size: NavigatorSizes.bands,
  fills: true,
  picture: (c) => {
    const data = c.navigator.bands;
    if (!data) return c.zoom ? undefined : { tracks: [] };
    return { tracks: [fill('bands', c.height, { type: 'bands', data })] };
  },
};

/** The bound section's tracks small (`tracks`, else all but rulers), or its first rows track (`rows`). */
function sectionTracks(c: NavigatorPictureContext, rowsOnly: boolean): NavigatorPictureTracks | undefined {
  if (c.zoom) return undefined;
  const nav = c.navigator;
  const tracks: Node[] = [];
  let rows: string | undefined;
  for (const t of c.boundNode?.children ?? []) {
    if (t.kind !== 'track' || t.id === undefined || t.params?.hidden === true) continue;
    const address = `${c.bound.address}/${t.id}`;
    if (c.left.has(address)) continue;
    const wanted = rowsOnly ? hasRows(t) && !rows : nav.tracks ? nav.tracks.includes(address) : !isRuler(t);
    if (!wanted) continue;
    if (hasRows(t) && !rows) rows = t.id;
    tracks.push(smallTrack(t, c.config));
  }
  return { tracks, ...(rows ? { rows } : {}) };
}

/** The section's tracks small. */
export const TracksPicture: NavigatorPicture = {
  name: 'tracks',
  axes: ['positions', 'both'],
  size: NavigatorSizes.tracks,
  picture: (c) => sectionTracks(c, false),
};

/** The section's rows summarized, filling the navigator (an overview's whole column). */
export const RowsPicture: NavigatorPicture = {
  name: 'rows',
  axes: ['positions', 'rows', 'both'],
  size: NavigatorSizes.rows,
  fills: true,
  picture: (c) => sectionTracks(c, true),
};

/**
 * The navigator pictures the core ships, by name: what the default compositions list (here, not in `built-ins`: a
 * composition draws them only with its navigators, which `renderSvg` loads when an export includes them).
 */
export const BuiltInNavigatorPictures: Readonly<Record<string, NavigatorPicture>> = {
  axis: AxisPicture,
  tracks: TracksPicture,
  bands: BandsPicture,
  rows: RowsPicture,
};

/**
 * The picture a navigator draws: the one its `content` names, or the `rows` picture where it covers `rows` and its
 * picture cannot; the axis where the table lacks it (the fixed fallback: every system has an axis).
 */
export function navigatorPicture(
  nav: ResolvedNavigator,
  pictures: Readonly<Record<string, NavigatorPicture>>,
): NavigatorPicture {
  let p = pictures[nav.content];
  if (nav.axes === 'rows' && !p?.axes.includes('rows')) p = pictures.rows;
  return p ?? AxisPicture;
}
