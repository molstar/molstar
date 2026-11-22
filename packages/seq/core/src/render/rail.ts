/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  axisBreaks,
  blockElementAt,
  nearestOffsetOf,
  piecesTotal,
  positionAt,
  type ResolvedNavigator,
  type ResolvedSection,
  type ResolvedView,
} from '@molstar/seqx-core/view';
import type { AxisTransform, Viewport } from './axis.js';
import type { RenderConfig } from './config.js';
import { Round } from './draw-lists.js';
import type { PassWriter } from './draw-lists.js';
import type { Box, Emphasis, Hit, HitTarget } from './layer-type.js';
import { collapsedSection } from './headers.js';
import type { WindowPart } from './layout/types.js';
import type { RowsLay } from './rows.js';
import type { SeqTokens } from './tokens.js';
import { floorExtent, navigatorPart } from './window.js';

// The ruler rail (_spec/README.md "Render", Navigators, The ruler rail): a linear section's default navigator, in a
// row just below its ruler's bottom border (the ruler's well grown by `SlimRailPx`; a slim row of its own under the
// tracks where it has no ruler shown), the whole axis across the data box, the section's window a thicker segment over
// it with a grip (a tick) at each end. Its grab zone is the ruler's height up from its row's foot (on an alignment's
// columns its lower `grabShare`: above it a drag selects columns), the slim row whole. Drawn grown (`grownPx`) while the pointer is over its grab zone or its
// drag runs; idle (the view cannot move) a quiet hairline with no grips and no target, so a press falls through to the
// ruler; in a figure (an export that includes navigators) a plain segment, no grips, no idle state. Collapsed breaks a
// notch each; a window over a cyclic origin two pieces, a grip at each true end.
//
// The rows rail is the same rail turned upright (one geometry: its along axis is the page's y, its line and grab zone
// across it on x): a rows layer's default navigator, at the right edge of its band beside the right gutter's numbers,
// `SlimRailPx` wide, the whole of its scrolled rows along the band's scrolled area, its window the rows in view. Idle
// only while every row is shown and the rows cannot zoom in (navigators' idle rule: a rail zoomed out to the whole band
// keeps its grips, so a grip drags it smaller again).

/**
 * The slim rail row's height, CSS px, fixed: the least a 2 px rail and its grips need to be seen and grabbed, at a
 * row's least cost (the form's size, not a reader's choice).
 */
export const SlimRailPx = 12;

/**
 * The rows' zoom bounds, px a row (fixed: 1/64 shows 10,000 rows in about 160 px, a band's whole family; 40 is twice a
 * letter row, past which a row reads no better): what the rows' zoom keys and the rows rail's grips stop at, and what
 * the rows rail's idle rule reads.
 */
export const MinRowScale = 1 / 64;
export const MaxRowScale = 40;

/**
 * A section's rail, or a rows band's (upright): where it is drawn and hit, and the window it shows (in the section's
 * own display units; a rows rail's in display rows below the pinned band).
 */
export class Rail {
  /** The section's address, id and system (a rows rail's: its rows system). */
  section = '';
  id = '';
  system = '';
  /** What says whether it can move (`railsActive`): the section's address, a rows rail's track's. */
  key = '';
  /** Upright: a rows rail (its along axis the page's y, its line and grab zone across it on x). */
  vertical = false;
  /**
   * Every row is shown and would be at the rows' zoom-in limit too (a rows rail: the window has nowhere to go): idle
   * whatever the view can do.
   */
  whole = false;
  /**
   * Along the rail: where its axis starts and its length (a section's data box's x and width; a rows rail's scrolled
   * area's y and height). Across it: the line, and the grab zone's two sides (a rows rail's: x).
   */
  x = 0;
  width = 0;
  y = 0;
  top = 0;
  bottom = 0;
  /** The ends' inset along the rail, px (the well's corner radius): already in `x` and `width`. */
  inset = 0;
  /** A slim row of its own (no ruler shown). */
  slim = false;
  /** Pinned to the bottom edge of the view's visible frame (`pinned`), outside the section's rows. */
  pinned = false;
  /** The section's display units and its window over them (end past total over a cyclic origin). */
  total = 0;
  start = 0;
  end = 0;
  /** The view can move (else drawn idle, no target): set each frame by the renderer, outside the key. */
  active = false;
  /** In a figure (an export): a plain segment, never idle. */
  figure = false;
  /** Along the rail, the last place hit outside the window (at the fit its body): its hover's position mark. */
  markAt = Number.NaN;
  /** Each collapsed break's display unit. */
  breaks: number[] = [];
  /** A section of blocks' blocks: a unit's place is a block's system and position (never the section's units). */
  blocks: ResolvedSection['blocks'] = undefined;
  /** The section's own pieces and direction: a unit's position. */
  readonly shape: { ranges: readonly number[]; dir: 1 | -1; cyclic: boolean } = { ranges: [], dir: 1, cyclic: false };
  /** The hit it gives, written in place. */
  readonly hit: { -readonly [K in keyof Hit]: Hit[K] } = {
    layer: '',
    section: '',
    system: '',
    target: 'window' as HitTarget,
    part: 'body' as WindowPart,
    navigator: 'rail',
    position: 0,
  };
}

/**
 * Sets the rail's axis along it: the data box's span, inset at both ends by the well's corner radius (never past its
 * rounded corner), so every reader of `x` and `width` (the drag, the hit, the draw)
 * maps the same pixels.
 */
function setAlong(rail: Rail, x: number, width: number, config: RenderConfig): void {
  const r = Math.max(0, Math.min(config.squircleMax, width / 4));
  rail.inset = r;
  rail.x = x + r;
  rail.width = width - 2 * r;
}

/** The rail's thickness: a line a hint, a figure's window a segment twice it (a usable rail's at rest is a line too). */
const Line = 1;
const Rest = 2;

/**
 * The display unit at a point on a rail (the whole axis along it): what a drag on it reads; `x` along the rail (a rows
 * rail's: the page's y).
 */
export function railUnitAt(rail: Rail, x: number): number {
  return ((x - rail.x) / rail.width) * rail.total;
}

/**
 * Places a section's rail: on its ruler's bottom border (`ruler`, the first shown track drawing one; on an
 * alignment's columns the grab zone its lower `grabShare`), else across the slim row at `slimTop`.
 */
export function placeRail(
  rail: Rail,
  section: ResolvedSection,
  sectionId: string,
  axis: AxisTransform,
  own: readonly number[],
  ruler: Box | undefined,
  slimTop: number,
  viewport: Viewport | undefined,
  config: RenderConfig,
  figure = false,
  pinned = false,
): void {
  rail.pinned = pinned;
  rail.section = section.address;
  rail.key = section.address;
  rail.vertical = false;
  rail.whole = false;
  rail.hit.layer = section.address;
  rail.hit.navigator = 'rail';
  rail.figure = figure;
  rail.id = sectionId;
  rail.system = section.axis.system;
  rail.blocks = section.blocks;
  setAlong(rail, axis.x0, axis.width, config);
  rail.slim = !ruler;
  if (ruler) {
    // in the row just below the ruler's bottom border (its well grown by it), clear of its ticks and labels; the grab
    // zone the ruler's height as before, down to that row's foot
    const columns = section.tracks.some((t) => t.layers.some((l) => l.type === 'rows' || l.type === 'logo'));
    const border = ruler.y + ruler.height;
    // a figure (an export that includes its rails) has no row for it: on the border
    const row = figure ? 0 : SlimRailPx;
    rail.y = border + row / 2;
    rail.bottom = border + row;
    rail.top = rail.bottom - (columns ? ruler.height * config.navigator.rail.grabShare : ruler.height);
  } else {
    rail.top = slimTop;
    rail.bottom = slimTop + SlimRailPx;
    rail.y = slimTop + SlimRailPx / 2;
  }
  const total = piecesTotal(own);
  rail.total = total;
  if (section.follow) {
    // a follower's window on its own axis: the positions it shows, as its own offsets
    const p0 = positionAt(axis, axis.u0, 'start');
    const p1 = positionAt(axis, axis.u1, 'end');
    const shape = { ranges: own, dir: section.axis.dir, cyclic: section.axis.cyclic };
    const a = p0 === undefined ? undefined : nearestOffsetOf(shape, p0);
    const b = p1 === undefined ? undefined : nearestOffsetOf(shape, p1);
    const from = a ?? 0;
    const to = b ?? total;
    // over a cyclic origin the end is before the start: the window runs on, its end past the total
    if (section.axis.cyclic && to < from) {
      rail.start = from;
      rail.end = to + total;
    } else {
      rail.start = Math.min(from, to);
      rail.end = Math.max(from, to);
    }
  } else if (viewport && section.axis.cyclic && viewport.end > total) {
    rail.start = viewport.start;
    rail.end = viewport.end;
  } else {
    rail.start = Math.max(0, Math.min(total, axis.u0));
    rail.end = Math.max(rail.start, Math.min(total, axis.u1));
  }
  rail.shape.ranges = own;
  rail.shape.dir = section.axis.dir;
  rail.shape.cyclic = section.axis.cyclic;
  const b = axisBreaks({ ranges: own, dir: section.axis.dir, cyclic: section.axis.cyclic });
  rail.breaks.length = 0;
  for (let i = 0; i < b.length; i += 3) rail.breaks.push(b[i]);
}

/**
 * Places a rows layer's rail (`track` its address): upright at `x` (its left side, `SlimRailPx` wide), along the band's
 * scrolled area (below the pinned rows), its window the rows in view, in display rows below the pinned band.
 */
export function placeRowsRail(
  rail: Rail,
  lay: RowsLay,
  track: string,
  sectionId: string,
  x: number,
  config: RenderConfig,
  figure = false,
) {
  const d = lay.display;
  const rows = Math.max(0, d.count - d.pinned);
  const area = Math.max(0, lay.box.height - lay.pinnedHeight);
  const start = Math.max(0, Math.min(rows, lay.top));
  placeUpright(
    rail,
    track,
    sectionId,
    lay.section,
    lay.rowsSystem,
    RowsRail,
    figure,
    lay.box.y + lay.pinnedHeight,
    area,
    x,
    rows,
    start,
    Math.max(start, Math.min(rows, lay.top + area / Math.max(1e-9, lay.rowScale))),
    config,
  );
  // idle only when the window can neither move nor shrink: every row shown, and shown still at the zoom-in limit (as a
  // ruler rail's: the whole axis in view and the axis fitting at the most px a position)
  rail.whole = rail.start <= 0 && rail.end >= rows && rows * MaxRowScale <= area + 1e-9;
}

/** An upright rail (a rows band's, a scrolled lanes band's): `total` units along `area` px from `y`, at x, its window. */
function placeUpright(
  rail: Rail,
  track: string,
  sectionId: string,
  section: string,
  system: string,
  navigator: string,
  figure: boolean,
  y: number,
  area: number,
  x: number,
  total: number,
  start: number,
  end: number,
  config: RenderConfig,
) {
  rail.section = section;
  rail.key = track;
  rail.vertical = true;
  rail.pinned = false;
  rail.hit.layer = track;
  rail.hit.navigator = navigator;
  rail.figure = figure;
  rail.id = sectionId;
  rail.system = system;
  setAlong(rail, y, area, config);
  rail.slim = true;
  rail.top = x;
  rail.bottom = x + SlimRailPx;
  rail.y = x + SlimRailPx / 2;
  rail.total = total;
  rail.start = start;
  rail.end = end;
  rail.whole = false;
  rail.shape.ranges = [0, total];
  rail.shape.dir = 1;
  rail.shape.cyclic = false;
  rail.breaks.length = 0;
}

/**
 * A scrolled lanes band's rail (_spec/README.md "Render", Track heights): the rows rail's geometry upright at the band's
 * right edge, the whole of its lanes along the band, its window the lanes in view; pan only (no grips, its window's ends
 * are body), drawn as the rows rail is (a hairline at rest, grown while hovered or dragged). The wheel moves it, so
 * do a drag, a press outside the window and the Page keys (the rows rail's gesture, in lanes).
 */
export function placeLanesRail(
  rail: Rail,
  box: Box,
  track: string,
  sectionId: string,
  lanes: number,
  room: number,
  top: number,
  config: RenderConfig,
  inset = 0,
) {
  placeUpright(
    rail,
    track,
    sectionId,
    sectionId,
    '',
    LanesRail,
    false,
    box.y,
    box.height,
    box.x + box.width - SlimRailPx - inset,
    lanes,
    top,
    Math.min(lanes, top + room),
    config,
  );
  // it exists only while its lanes overflow: always a target
  rail.active = true;
}

/** A lanes rail's `navigator`. */
export const LanesRail = 'lanes-rail';

/**
 * A section's rail where it scrolls inside its box (_spec/README.md "Render", Section scroll): the one rail upright at
 * `x` (the right gutter's rail slot) along its scrolled area (`y` to `end`, px: below a fixed head), the scrolled
 * content (`scrolled` px) along it, its window the px in view; pan only, always a target (it exists only while the
 * content overflows).
 */
export function placeSectionRail(
  rail: Rail,
  x: number,
  y: number,
  end: number,
  section: string,
  sectionId: string,
  scrolled: number,
  top: number,
  config: RenderConfig,
) {
  const area = Math.max(0, end - y);
  const total = Math.max(area, scrolled);
  placeUpright(
    rail,
    section,
    sectionId,
    section,
    '',
    SectionRail,
    false,
    y,
    area,
    x,
    total,
    top,
    Math.min(total, top + area),
    config,
  );
  rail.active = true;
}

/** A section rail's `navigator`. */
export const SectionRail = 'section-rail';

/** Whether a rail only pans (a lanes band's, a section's): no grips, its window's ends are body. */
export function panOnly(navigator: string | undefined): boolean {
  return navigator === LanesRail || navigator === SectionRail;
}

/** A rows rail's hit's `navigator`. */
export const RowsRail = 'rows-rail';

const Ends = new Float64Array(4);

/** The window's ends on the rail, px (the first piece over a cyclic origin: from its start to the axis's end). */
function windowEnds(rail: Rail, config: RenderConfig): void {
  const k = rail.width / Math.max(1e-9, rail.total);
  const a = rail.x + rail.start * k;
  const b = rail.x + Math.min(rail.end, rail.total) * k;
  Ends[0] = a;
  Ends[1] = b;
  Ends[2] = rail.x;
  Ends[3] = rail.x + rail.width;
  floorExtent(Ends, 0, config.navigator.window);
}

/**
 * The rect in `out.args` (along, across, length, thickness, fill; upright: turned to x, y, width, height), written
 * unless empty (numbers through the array: never boxed).
 */
function emit(out: PassWriter, vertical: boolean): void {
  const g = out.args;
  if (!(g[2] > 0) || !(g[3] > 0)) return;
  if (vertical) {
    const a = g[0];
    const l = g[2];
    g[0] = g[1];
    g[1] = a;
    g[2] = g[3];
    g[3] = l;
  }
  g[5] = 0;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/**
 * The window pass (_spec/README.md "Render", Navigators, The ruler rail): the hairline the axis's length, its breaks'
 * notches, the window a heavier segment on it with a round knob at each end where the box zooms (round caps on a pan-only
 * rail); idle (the view cannot move) a quiet hairline, no knobs; a figure's a plain segment, no knobs.
 */
export function drawRail(
  rail: Rail,
  out: PassWriter,
  emphasis: Emphasis | undefined,
  tokens: SeqTokens,
  config: RenderConfig,
) {
  if (rail.pinned) {
    // pinned over the frame's foot: its own ground, so what the frame clips never shows through it
    const g = out.args;
    g[0] = rail.x - rail.inset;
    g[1] = rail.top;
    g[2] = rail.width + 2 * rail.inset;
    g[3] = rail.bottom - rail.top;
    g[4] = tokens.colors['msui-bg'];
    emit(out, false);
  }
  const h = emphasis?.hover;
  const d = emphasis?.drag;
  const own = rail.hit.navigator;
  const lit = h?.target === 'window' && h.navigator === own && h.layer === rail.key ? h.part : undefined;
  const dragged = d !== undefined && d.navigator === own && d.layer === rail.key;
  // a place hovered outside the window (at the fit, its whole body): its position marked, the rail never grown or lit
  const place = lit === 'outside' && !dragged && !Number.isNaN(rail.markAt);
  const used = rail.active && ((lit !== undefined && !place) || dragged);
  paint(rail, out, tokens, config, used, dragged, place, dragged ? d.part : lit);
}

/** A color with its alpha set (0 to 1). */
function withAlpha(color: number, a: number): number {
  return ((color & 0xffffff00) | Math.round(a * 255)) >>> 0;
}

/** The capsule's numbers (along, length, thickness, fill), set by `pill`'s callers: a call boxes the numbers it is passed. */
const Pl = new Float64Array(4);

/** A disc's numbers (along, diameter, fill), set by its callers: a call boxes the numbers it is passed. */
const Dk = new Float64Array(3);
const DiscPoints = new Float64Array(2 * 24);

/**
 * A disc centred on the rail's line at `Dk[0]` along it (a path of 24 points: a squircle rect is no circle), its
 * diameter `Dk[1]`, its fill `Dk[2]`.
 */
function discOut(out: PassWriter, rail: Rail): void {
  const r = Dk[1] / 2;
  if (!(r > 0)) return;
  for (let i = 0; i < 24; i++) {
    const a = (i / 24) * 2 * Math.PI;
    const u = Dk[0] + r * Math.cos(a);
    const v = rail.y + r * Math.sin(a);
    DiscPoints[2 * i] = rail.vertical ? v : u;
    DiscPoints[2 * i + 1] = rail.vertical ? u : v;
  }
  const g = out.args;
  g[4] = Dk[2];
  g[5] = 0;
  g[6] = 0;
  g[7] = 0;
  out.pathArgs(DiscPoints, 24, true);
}

/** A capsule: `len` along the rail from `a`, `t` thick across, centred on the rail's line, round caps of half its thickness. */
function pillOut(out: PassWriter, rail: Rail): void {
  const a = Pl[0];
  const len = Pl[1];
  const t = Pl[2];
  const fill = Pl[3];
  if (!(len > 0) || !(t > 0)) return;
  const g = out.args;
  const c = rail.y - t / 2;
  g[0] = rail.vertical ? c : a;
  g[1] = rail.vertical ? a : c;
  g[2] = rail.vertical ? t : len;
  g[3] = rail.vertical ? len : t;
  g[4] = fill;
  g[5] = 0;
  g[6] = 0;
  g[7] = Math.min(t, len) / 2;
  g[8] = Round.both;
  out.rectArgs();
}

/**
 * The rail's picture: the hairline `msui-line`, the window a heavier segment in a mid ink (`msui-text` at 55%; the
 * `accent` tint `msui-accent` at 55%), both in full ink (`seq-window`) and `grownPx` while used; a round knob at each end
 * of the window where the box zooms (the segment's colour, a 1 px ring in `msui-bg`, growing about its centre while
 * used), round caps on the segment where it only pans. At the fit the segment is the hairline's weight and colour with
 * small knobs.
 */
function paint(
  rail: Rail,
  out: PassWriter,
  tokens: SeqTokens,
  config: RenderConfig,
  used: boolean,
  dragged: boolean,
  place: boolean,
  part: WindowPart | undefined,
): void {
  const line = tokens.colors['msui-line'];
  const full = tokens.colors['seq-window'];
  const mid = withAlpha(tokens.colors[config.navigator.rail.tint === 'accent' ? 'msui-accent' : 'msui-text'], 0.55);
  const bg = tokens.colors['msui-bg'];
  Pl[0] = rail.x;
  Pl[1] = rail.width;
  Pl[2] = Line;
  Pl[3] = line;
  pillOut(out, rail);
  if (!(rail.total > 0)) return;
  for (let i = 0; i < rail.breaks.length; i++) {
    Pl[0] = rail.x + (rail.breaks[i] / rail.total) * rail.width - 0.5;
    Pl[1] = 1;
    Pl[2] = 4;
    Pl[3] = line;
    pillOut(out, rail);
  }
  const fit = rail.start <= 0 && rail.end >= rail.total;
  const quiet = fit && !used && !rail.figure;
  const ink = used ? full : rail.figure ? mid : rail.active ? mid : withAlpha(line, 0.8);
  const t = used ? config.navigator.rail.grownPx : quiet || (!rail.active && !rail.figure) ? Line : Rest;
  windowEnds(rail, config);
  const lo = rail.x;
  const hi = rail.x + rail.width;
  const a = Math.max(lo, Ends[0]);
  const b = Math.min(hi, Ends[1]);
  const segment = quiet && rail.active ? line : ink;
  Pl[0] = a;
  Pl[1] = b - a;
  Pl[2] = t;
  Pl[3] = segment;
  pillOut(out, rail);
  if (rail.end > rail.total) {
    Pl[0] = lo;
    Pl[1] = ((rail.end - rail.total) / rail.total) * rail.width;
    Pl[2] = t;
    Pl[3] = segment;
    pillOut(out, rail);
  }
  if (!rail.active) return;
  if (place) {
    Pl[0] = Math.max(lo, Math.min(hi, rail.markAt)) - 0.75;
    Pl[1] = 1.5;
    Pl[2] = 7;
    Pl[3] = full;
    pillOut(out, rail);
  }
  if (panOnly(rail.hit.navigator)) return;
  // the knobs, circles centred on the window's true ends (a wrapped window's second piece ends before the axis's end):
  // 5 px at rest (4 at the fit), 10 hovered, 12 dragged, the hovered end one more, each about its centre; they may
  // overhang the slim row (nothing else is drawn on the rail's line, and the window pass is not clipped)
  const wraps = rail.end > rail.total;
  for (let e = 0; e < 2; e++) {
    const end = part === (e === 0 ? 'start' : 'end');
    const dia = (quiet ? 4 : used ? (dragged ? 12 : 10) : 5) + (used && end ? 1 : 0);
    const c = e === 0 ? a : wraps ? lo + ((rail.end - rail.total) / rail.total) * rail.width : b;
    Dk[0] = c;
    Dk[1] = dia + 2;
    Dk[2] = bg;
    discOut(out, rail);
    Dk[1] = dia;
    Dk[2] = segment;
    discOut(out, rail);
  }
}

/**
 * What a point on the rail's grab zone is: its window's part, or what is outside it (`touch`: the wider hit);
 * undefined idle or elsewhere. Over a cyclic origin the window's two pieces meet the axis's ends as body: its grips are
 * at its true ends only.
 */
export function hitRail(rail: Rail, px: number, py: number, config: RenderConfig, touch = false): Hit | undefined {
  // along the rail and across it (upright: the page's y and x); the grab zone reaches the inset (the old extent), a press
  // there is at the axis's end (the unit is clamped)
  const x = rail.vertical ? py : px;
  const y = rail.vertical ? px : py;
  if (!rail.active || !(rail.total > 0) || y < rail.top || y >= rail.bottom) return undefined;
  const c = config.navigator.window;
  windowEnds(rail, config);
  let part = navigatorPart(x, Ends[0], Ends[1], rail.x - rail.inset, rail.x + rail.width + rail.inset, c, touch);
  if (rail.end > rail.total) {
    // the first piece runs to the axis's end: its end there is body
    if (part === 'end') part = 'body';
    if (part === 'outside') {
      const b = rail.x + ((rail.end - rail.total) / rail.total) * rail.width;
      part = navigatorPart(x, rail.x, b, rail.x - rail.inset, rail.x + rail.width + rail.inset, c, touch);
      if (part === 'start') part = 'body';
    }
  }
  if (!part) return undefined;
  // a lanes rail only pans: its window's ends are body
  if (panOnly(rail.hit.navigator) && (part === 'start' || part === 'end')) part = 'body';
  // a window over the whole axis (the fit) has no body to drag, but the rail stays a target: its body is a place as
  // outside a window is (the hover names the place and marks it, never the whole rail lit; a sweep zooms to the region,
  // a click as elsewhere); its grips still resize it
  if (part === 'body' && rail.start <= 0 && rail.end >= rail.total) part = 'outside';
  rail.markAt = part === 'outside' ? x : Number.NaN;
  const h = rail.hit;
  h.section = rail.id;
  h.system = rail.system;
  h.part = part;
  h.block = undefined;
  const u = Math.max(0, Math.min(rail.total - 1e-9, railUnitAt(rail, x)));
  const blocks = rail.blocks;
  if (blocks?.length && blocks[0].offset !== undefined) {
    // a section of blocks: the block under the unit, its system and its position there (its direction and scale)
    let b = 0;
    while (b + 1 < blocks.length && (blocks[b + 1].offset ?? Infinity) <= u) b++;
    h.block = b;
    h.system = blocks[b].axis.system;
    h.position = blockElementAt(blocks[b], u) ?? 0;
    return h;
  }
  h.position = Math.floor(positionAt(rail.shape, u, 'start') ?? u);
  return h;
}

// Which sections and tracks draw a rail (here, not with the navigators: a viewer that lists no navigator draws rails)
// (a section's layout by name: `sectionLayoutName` of the tracks layout, which imports this module).

/** When a linear section draws its positions rail: beside a declared navigator too, replaced by one, or never. */
export type PositionsRailMode = 'always' | 'unlessNavigator' | 'none';

/** Where the positions rail sits (`navigator.railPlace`): below the ruler, under the last track, by a navigator. */
export type RailPlace = 'ruler' | 'bottom' | 'navigator' | 'pinned';

/**
 * The row a section's rail is drawn in: below its ruler, under its last track, above its first track, or pinned to the
 * bottom edge of the view's visible frame.
 */
export type RailRow = 'ruler' | 'bottom' | 'top' | 'pinned';

/**
 * The sections that draw a ruler rail by default (_spec/README.md "Render", Navigators): every linear section (tracks)
 * with `always` (the default: a declared navigator is the picture, the rail the precise control); with
 * `unlessNavigator` those no shown navigator binds along its positions (the view's declared ones and a host's
 * (`hosted`) replace it; a hidden one brings it back); with `none` no section.
 */
export function railSections(
  view: ResolvedView,
  hosted: readonly ResolvedNavigator[] = [],
  mode: PositionsRailMode = 'always',
): string[] {
  if (mode === 'none') return [];
  const bound = new Set<string>();
  if (mode === 'unlessNavigator')
    for (const list of [view.navigators, hosted])
      for (const n of list) if (!n.hidden && n.axes !== 'rows') for (const a of n.of) bound.add(a);
  return view.sections
    .filter((s) => (s.layout ?? view.layout) === 'tracks' && !bound.has(s.address) && !collapsedSection(s))
    .map((s) => s.address);
}

/**
 * The row each rail section's positions rail is drawn in (`place`, `navigator.railPlace`): `ruler` below its ruler's
 * border (a section with no ruler shown takes the `bottom` row, decided where it is laid out), `bottom` under its last
 * track, `navigator` next to a shown horizontal navigator (`top` or `bottom`) over its positions: under a top one
 * (above the section's first track, directly below it on the page), at a bottom one's side (under the last track),
 * else `ruler`.
 */
export function railRows(
  view: ResolvedView,
  sections: readonly string[],
  hosted: readonly ResolvedNavigator[] = [],
  place: RailPlace | readonly RailPlace[] = 'ruler',
): Record<string, readonly RailRow[]> {
  const places: readonly RailPlace[] = typeof place === 'string' ? [place] : place.length ? place : ['ruler'];
  const beside: Record<string, RailRow> = {};
  if (places.includes('navigator'))
    for (const list of [view.navigators, hosted])
      for (const n of list) {
        if (n.hidden || n.axes === 'rows' || (n.place !== 'top' && n.place !== 'bottom')) continue;
        for (const a of n.of) beside[a] = n.place;
      }
  const out: Record<string, readonly RailRow[]> = {};
  for (const a of sections) {
    const rows: RailRow[] = [];
    for (const p of places) {
      const row: RailRow = p === 'navigator' ? (beside[a] ?? 'ruler') : p;
      if (!rows.includes(row)) rows.push(row);
    }
    out[a] = rows;
  }
  return out;
}

/**
 * The rows a section's rail entry names (`FrameInput.rails`: `true`, one row, or a list), figure-safe: a figure has no
 * row above its tracks (its `top` is the `bottom`) and no pinned one; duplicates dropped.
 */
export function railRowList(entry: unknown, figure: boolean): readonly RailRow[] {
  const raw: readonly unknown[] = Array.isArray(entry) ? entry : [entry];
  const out: RailRow[] = [];
  for (const r of raw) {
    let row: RailRow = r === 'top' || r === 'bottom' || r === 'pinned' ? r : 'ruler';
    if (figure && row === 'top') row = 'bottom';
    if (figure && row === 'pinned') continue;
    if (!out.includes(row)) out.push(row);
  }
  return out;
}

/**
 * The rows tracks that draw a rows rail by default (_spec/README.md "Render", Navigators, The rows rail): each linear
 * section's tracks with a rows layer shown; beside a shown navigator over the section's rows (`always`, the default: the
 * overview is a picture, the rail the precise control), replaced by it (`unlessOverview`: a hidden one brings the rail
 * back), or none (`none`).
 */
export function rowsRailTracks(
  view: ResolvedView,
  hosted: readonly ResolvedNavigator[] = [],
  mode: RowsRailMode = 'always',
): string[] {
  if (mode === 'none') return [];
  const bound = new Set<string>();
  if (mode === 'unlessOverview')
    for (const list of [view.navigators, hosted])
      for (const n of list) if (!n.hidden && n.axes !== 'positions') for (const a of n.of) bound.add(a);
  const out: string[] = [];
  for (const s of view.sections) {
    if ((s.layout ?? view.layout) !== 'tracks' || bound.has(s.address)) continue;
    for (const t of s.tracks)
      if (t.layers.some((l) => l.type === 'rows' && l.placed && l.params.hidden !== true)) out.push(t.address);
  }
  return out;
}

/** When a rows band draws its rows rail: beside a navigator over its rows, unless one binds them, or never. */
export type RowsRailMode = 'always' | 'unlessOverview' | 'none';
