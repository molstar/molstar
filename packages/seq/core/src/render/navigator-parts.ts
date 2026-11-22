/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Node, type Value } from '@molstar/core/repr/tree';
import { type StaticModel } from '@molstar/seqx-core/resolve';
import {
  blockSpan,
  type LayerKinds,
  nearestOffsetOf,
  offsetOf,
  type ResolvedNavigator,
  type ResolvedSection,
  type ResolvedView,
  resolveViewNode,
} from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { type RenderConfig, renderConfig } from './config.js';
import { type PassWriter } from './draw-lists.js';
import type { WindowPart } from './layout/types.js';
import {
  AxisPicture,
  hasRows,
  type NavigatorPicture,
  NavigatorSizes,
  navigatorPicture,
  smallHeight,
} from './navigator-picture.js';
import type { RowsLay } from './rows.js';
import { type SectionLayout, type ViewLayout } from './tracks.js';
import { NavWindowFlags, navigatorPart } from './window.js';

// A navigator's parts (_spec/README.md "Render", Navigators): its picture made as a view of its own, its windows and
// their text (by block over a genome), and the helpers its renderer reads; split from `navigator.ts` by what changes
// together.

export const ZoomId = 'zoom';

/**
 * A navigator's box for a view of `width` × `height`: its `size`, else its picture's (a thickness along or across the
 * view, or a miniature's box); over both axes a picture's thickness gives way to the miniature's default box.
 */
export function navigatorBox(
  nav: ResolvedNavigator,
  width: number,
  height: number,
  pictures: Readonly<Record<string, NavigatorPicture>>,
): { readonly width: number; readonly height: number } {
  const own = navigatorPicture(nav, pictures).size;
  const size = nav.size ?? (nav.axes === 'both' && !Array.isArray(own) ? NavigatorSizes.both : own);
  if (Array.isArray(size)) return { width: size[0], height: size[1] };
  const t = size as number;
  return nav.place === 'left' || nav.place === 'right' ? { width: t, height } : { width, height: t };
}

/** What a navigator's picture is: its view, the section that draws it, the bound section's rows track there. */
export interface Picture {
  readonly view: ResolvedView;
  readonly section: ResolvedSection;
  /** The bound section's first rows track and its copy here (their addresses), when the picture draws rows. */
  readonly rows?: { readonly view: string; readonly picture: string };
}

export const PictureId = 'picture';

export function nodeAt(root: Node, address: string): Node | undefined {
  let n: Node | undefined = root;
  for (const s of address.split('/').filter(Boolean)) n = n?.children?.find((c) => c.id === s);
  return n;
}

/**
 * The navigator's picture (_spec/README.md "Render", Navigators): a view of one section over its system whole, in the
 * bound section's display order, holding what its `content` draws; `left` the tracks it leaves out (what no summary
 * serves). Resolved like any view, so its layers are the spec's layers, read and colored as the view's.
 */
export function pictureOf(
  model: StaticModel,
  view: ResolvedView,
  nav: ResolvedNavigator,
  entry: NavigatorPicture,
  kinds: LayerKinds | undefined,
  config: RenderConfig,
  left: ReadonlySet<string>,
  height: number,
  direction: 'display' | 'genomic',
  zoom = 0,
): Picture | undefined {
  const bound = view.sections.find((s) => s.address === nav.of[0]);
  if (!bound || nav.system === undefined) return undefined;
  const node = nodeAt(model.spec, bound.address);
  const same = !nav.blocks && nav.sectionAxis && nav.system === bound.axis.system;
  const dir = direction === 'genomic' ? 1 : bound.axis.dir;
  const params: Record<string, Value> = { of: nav.system, direction: dir === -1 ? 'reversed' : 'genomic' };
  let blocks = same ? (node?.children ?? []).filter((c) => c.kind === 'block') : [];
  if (nav.blocks) {
    // a navigator over blocks (_spec/README.md "Render", Navigators): a section of the bound section's own blocks (its
    // `of`, `classes`, `include` and block children: the same block list, identity in units), or of the assemblies
    // `system` names, their lattice whole (a bound section of them lends its blocks' directions and scales, never its
    // `include`)
    const of = nav.assemblies ? (nav.assemblies.length === 1 ? nav.assemblies[0] : nav.assemblies) : node?.params?.of;
    const own = !nav.assemblies || JSON.stringify(node?.params?.of) === JSON.stringify(of);
    if (of !== undefined) params.of = of as Value;
    else delete params.of;
    for (const k of nav.assemblies ? ['classes'] : ['classes', 'include', 'copies'])
      if (own && node?.params?.[k] !== undefined) params[k] = node.params[k] as Value;
    blocks = own ? (node?.children ?? []).filter((c) => c.kind === 'block') : [];
  }
  const ctx = { model, view, navigator: nav, bound, boundNode: node, height, zoom: false, config, left };
  // a picture that draws nothing here is the axis (the fixed fallback)
  const drawn = entry.picture(ctx) ?? AxisPicture.picture(ctx)!;
  const tracks = [...drawn.tracks];
  const rows: Picture['rows'] = drawn.rows
    ? { view: `${bound.address}/${drawn.rows}`, picture: `${nav.address}/${PictureId}/${drawn.rows}` }
    : undefined;
  // the rows fill what the other tracks leave of the navigator's height (an overview's whole column)
  const rowsAt = tracks.findIndex(hasRows);
  if (rowsAt >= 0 && height > 0) {
    let rest = height;
    for (let i = 0; i < tracks.length; i++) if (i !== rowsAt) rest -= smallHeight(tracks[i], config) + 1;
    const t = tracks[rowsAt];
    tracks[rowsAt] = { ...t, params: { ...t.params, height: Math.max(config.navigator.trackMinPx, rest) } };
  }
  const sections: Node[] = [{ kind: 'section', id: PictureId, params, children: [...blocks, ...tracks] }];
  // the zoom strip: the picture's own over the reachable span alone (its viewport), else the axis
  if (zoom > 0) {
    const z = { ...ctx, height: zoom, zoom: true };
    const strip = entry.picture(z) ?? AxisPicture.picture(z)!;
    // over blocks the strip has no label row (its context names the blocks, its mark the span): its height the track's
    const zp = nav.blocks ? { ...params, labels: 'none' } : params;
    sections.push({ kind: 'section', id: ZoomId, params: zp, children: [...blocks, ...strip.tracks] });
  }
  const viewNode: Node = { kind: 'view', id: nav.id, params: { layout: 'tracks' }, children: sections };
  // with the composition's layer kinds: a picture's track of an extension's type is drawn as that type
  const resolved = resolveViewNode(model, viewNode, nav.address, kinds);
  // the same system: the bound section's own axis (its pieces, collapsed introns and direction), never a guess
  const own = (s: ResolvedSection) =>
    same && !bound.blocks ? { ...s, axis: dir === bound.axis.dir ? bound.axis : genomicAxis(bound.axis) } : s;
  const out = resolved.sections.map(own);
  const section = out[0];
  if (!section || (zoom > 0 && !out[1])) return undefined;
  return { view: { ...resolved, sections: out, navigators: [] }, section, ...(rows ? { rows } : {}) };
}

/** An axis read up left to right (`navigator.direction: genomic`): a reversed axis's pieces in genomic order. */
export function genomicAxis(axis: ResolvedSection['axis']): ResolvedSection['axis'] {
  const r = axis.ranges;
  if (!r) return { ...axis, dir: 1 };
  const out: number[] = [];
  for (let k = r.length - 2; k >= 0; k -= 2) out.push(r[k], r[k + 1]);
  return { ...axis, dir: 1, ranges: out };
}

/**
 * The config a navigator's picture lays out with: no inset, label column or direction mark; its tracks and lanes small;
 * `gap` between its sections (a zoom strip's connectors).
 */
export function pictureConfig(base: RenderConfig, gap = 0): RenderConfig {
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
    labelColumn: { measure: false },
    trackGap: 1,
    sectionGap: gap,
    axis: { directionMark: 'none' },
    heights,
    // a lane at the tracks' scale: a small track keeps its lanes (its items drawn, not only their density)
    laneHeight: Math.max(min, base.laneHeight * s),
    laneGap: base.laneGap * s,
    lanePadding: base.lanePadding * s,
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
  /** Over blocks: the systems its ends are on (`p0` and `p1` positions on them). */
  s0 = '';
  s1 = '';
  r0 = 0;
  r1 = 0;
  rows = 0;
}

/** What a window says (an export's title, a slider's value): its positions on the navigator's system, its rows. */
export function windowText(model: StaticModel, nav: ResolvedNavigator, w: NavigatorWindow): string {
  if (nav.blocks) return blocksText(model, w);
  const name = (nav.system && model.systems.get(nav.system)?.name) ?? '';
  const parts: string[] = [];
  if (w.axes & NavWindowFlags.x) parts.push(`${spanText(model, nav, w.p0, w.p1)} of ${name}`);
  if (w.axes & NavWindowFlags.y) parts.push(`rows ${countText(w.r0 + 1)}–${countText(w.r1)} of ${countText(w.rows)}`);
  return parts.join(', ');
}

const countText = (n: number) => Math.round(n).toLocaleString('en-US');

/**
 * A window over blocks as text (_spec/README.md "Render", Navigators): on one block its name and positions and the
 * assembly (`chr17: 7,668,402–7,687,550 of GRCh38`), across blocks both ends (`chr16:90,100,000 – chr17:1,200,000`).
 */
function blocksText(model: StaticModel, w: NavigatorWindow): string {
  const a = model.systems.get(w.s0)?.name ?? '';
  const b = model.systems.get(w.s1)?.name ?? '';
  if (w.s0 !== w.s1) return `${a}:${countText(w.p0 + 1)} – ${b}:${countText(w.p1)}`;
  const of = assemblyName(model, w.s0);
  return `${a}: ${countText(Math.min(w.p0, w.p1) + 1)}–${countText(Math.max(w.p0, w.p1))}${of ? ` of ${of}` : ''}`;
}

/** The name of the assembly a derived sequence is of (its first name), or undefined. */
function assemblyName(model: StaticModel, system: string): string | undefined {
  const at = system.lastIndexOf('/(sequence ');
  return at > 0 ? model.assemblies.get(system.slice(0, at))?.names[0] : undefined;
}

/** A position over blocks as a readout says it (`chr17: 7,668,402`): the block's system's name, 1-based. */
export function blockPlaceText(model: StaticModel, system: string, position: number): string {
  return `${model.systems.get(system)?.name ?? ''}: ${countText(Math.floor(position) + 1)}`;
}

/**
 * Positions `[p0, p1)` on the navigator's system as text, 1-based (`chr17:7,666,403–7,688,402`); a reversed axis's
 * ends in either order, across a cyclic origin (`p1` before `p0`) the arc over it.
 */
export function spanText(model: StaticModel, nav: ResolvedNavigator, p0: number, p1: number): string {
  // a window across a cyclic origin ends before it starts: the arc over the origin, never its complement
  if (p1 < p0 && nav.system && model.systems.get(nav.system)?.cyclic) return arcText(model, nav, p0, p1);
  return arcText(model, nav, Math.min(p0, p1), Math.max(p0, p1));
}

/** An arc from `p0` on to `p1` on the navigator's system as text, 1-based (over a cyclic origin `p1` is before `p0`). */
export function arcText(model: StaticModel, nav: ResolvedNavigator, p0: number, p1: number): string {
  const label = (nav.system && model.systems.get(nav.system)?.name) ?? '';
  return `${label}:${countText(p0 + 1)}–${countText(p1)}`;
}

// the flag as a number read once: a draw reads no property of an imported object
const FlagX: number = NavWindowFlags.x;

/** Whether a window is at the fit: it moves along x alone and covers all it can reach, so it is no target (a sweep's row). */
export function windowAtFit(w: NavigatorWindow): boolean {
  return w.shown && w.axes === FlagX && !w.wrap && w.dx0 <= w.lo + 0.5 && w.dx1 >= w.hi - 0.5;
}

/**
 * Whether a positions navigator has a zoom strip (navigators.md, "The window's look", 1): over blocks by `mode`
 * (`always`; `narrow`: the window's width `px` under `minPx`, held `hysteresis` px past it while `on`; `never`), else by
 * the one-system reachable span `span` px narrower than `min`.
 */
export function wantsStrip(
  blocks: boolean,
  strip: { readonly minPx: number; readonly overBlocks?: 'always' | 'narrow' | 'never' },
  px: number,
  span: number,
  on: boolean,
  hysteresis: number,
): boolean {
  const min = strip.minPx;
  if (!blocks) return min > 0 && span < min;
  const mode = strip.overBlocks ?? 'always';
  return mode === 'always' || (mode === 'narrow' && (px < min || (on && px <= min + hysteresis)));
}

/**
 * The part along x of a window that may wrap a cyclic origin (two pieces: from its start to the axis's end, from the
 * axis's start to its end): where a piece meets the axis's end it is body, so its grips are at its true ends only.
 */
export function partOf(
  x: number,
  w: NavigatorWindow,
  c: RenderConfig['navigator']['window'],
  touch: boolean,
): WindowPart | undefined {
  if (windowAtFit(w)) return x >= w.lo && x < w.hi ? 'outside' : undefined;
  if (!w.wrap) return navigatorPart(x, w.dx0, w.dx1, w.lo, w.hi, c, touch);
  const a = navigatorPart(x, w.x0, w.x1, w.lo, w.hi, c, touch);
  if (a === 'end') return 'body';
  if (a !== 'outside') return a;
  const b = navigatorPart(x, w.w0, w.w1, w.lo, w.hi, c, touch);
  return b === 'start' ? 'body' : b;
}

export type Block = NonNullable<ResolvedSection['blocks']>[number];

/**
 * The block holding a section unit (the last starting at or before it); on a block's edge as a `start` the next one, as
 * an `end` the one before (a window's ends), by a binary search over the offsets: nothing made.
 */
export function blockIndexAt(blocks: readonly Block[], u: number, side: 'start' | 'end'): number {
  let lo = 0;
  let hi = blocks.length - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if ((blocks[mid].offset ?? 0) <= u) lo = mid;
    else hi = mid - 1;
  }
  if (side === 'start' && lo + 1 < blocks.length && u >= (blocks[lo].offset ?? 0) + blockSpan(blocks[lo])) lo++;
  if (side === 'end' && lo > 0 && u <= (blocks[lo].offset ?? 0)) lo--;
  return lo;
}

/** A position's section unit on a block (its direction and scale), its nearest edge where the block does not show it. */
export function blockUnitAtPosition(block: Block, position: number, offset: number): number | undefined {
  const o = offsetOf(block.axis, position) ?? nearestOffsetOf(block.axis, position);
  return o === undefined ? undefined : offset + o * (block.form === 'whole' ? 1 : (block.scale ?? 1));
}

export function trackHasRows(t: ResolvedSection['tracks'][number]): boolean {
  return t.layers.some((l) => l.type === 'rows');
}

export function sectionIn(view: ViewLayout, address: string): SectionLayout | undefined {
  const s = view.sections;
  for (let i = 0; i < s.length; i++) if (s[i].section.address === address) return s[i];
  return undefined;
}

/** The first rows layer's lay of a section's first strip. */
export function rowsLayOf(s: SectionLayout | undefined): RowsLay | undefined {
  if (!s) return undefined;
  const lays = s.lays;
  for (let i = 0; i < lays.length; i++) {
    const l = lays[i];
    for (let j = 0; j < l.length; j++) if (l[j].type.name === 'rows') return l[j].lay as RowsLay;
  }
  return undefined;
}

/** The folds the view's rows band shows (its session's, or its layer's `collapsed`), by identity. */
export function foldsOf(view: ViewLayout | undefined, section: string | undefined) {
  const s = view && section ? sectionIn(view, section) : undefined;
  const lay = rowsLayOf(s);
  return lay?.foldsApply ? lay.folds : undefined;
}

/** Ticks over position ranges at the picture's foot, each at least a device pixel. */
export function ticks(
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
  /** Over blocks: the system `position` is on, and the picture's block there (its index). */
  readonly system?: string;
  readonly block?: number;
  /** The square: the cell's y block (its index in the lattice) and its system. */
  readonly ySystem?: string;
  readonly yBlock?: number;
  /**
   * On the context of a navigator with a zoom strip (no target: no gesture), over the reachable span's mark (`mark`) or
   * elsewhere (`outside`); its part is `outside`.
   */
  readonly context?: 'mark' | 'outside';
}
