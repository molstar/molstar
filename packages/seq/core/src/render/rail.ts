/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { axisBreaks, nearestOffsetOf, piecesTotal, positionAt, type ResolvedSection } from '@molstar/seqx-core/view';
import type { AxisTransform, Viewport } from './axis.js';
import type { RenderConfig } from './config.js';
import type { PassWriter } from './draw-lists.js';
import type { Box, Emphasis, Hit, HitTarget } from './layer-type.js';
import type { WindowPart } from './layout/types.js';
import type { SeqTokens } from './tokens.js';
import { floorExtent, navigatorPart } from './window.js';

// The ruler rail (_spec/README.md "Render", Navigators, The ruler rail): a linear section's default navigator, with no
// row of its own: drawn on its ruler's bottom border (a slim row of its own, `SlimRailPx`, under the tracks where it
// has no ruler shown), the whole axis across the data box, the section's window a thicker segment over it with a grip
// (a tick) at each end. Its grab zone is the ruler's row (on an alignment's columns its lower `grabShare`: above it a
// drag selects columns), the slim row whole. Drawn grown (`grownPx`) while the pointer is over its grab zone or its
// drag runs; idle (the view cannot move) a quiet hairline with no grips and no target, so a press falls through to the
// ruler; in a figure (an export that includes navigators) a plain segment, no grips, no idle state. Collapsed breaks a
// notch each; a window over a cyclic origin two pieces, a grip at each true end.

/**
 * The slim rail row's height, CSS px, fixed: the least a 2 px rail and its grips need to be seen and grabbed, at a
 * row's least cost (the form's size, not a reader's choice).
 */
export const SlimRailPx = 12;

/** A section's rail: where it is drawn and hit, and the window it shows (in the section's own display units). */
export class Rail {
  /** The section's address, id and system. */
  section = '';
  id = '';
  system = '';
  /** The data box's x and width; the rail's line; its grab zone. */
  x = 0;
  width = 0;
  y = 0;
  top = 0;
  bottom = 0;
  /** A slim row of its own (no ruler shown). */
  slim = false;
  /** The section's display units and its window over them (end past total over a cyclic origin). */
  total = 0;
  start = 0;
  end = 0;
  /** The view can move (else drawn idle, no target): set each frame by the renderer, outside the key. */
  active = false;
  /** In a figure (an export): a plain segment, never idle. */
  figure = false;
  /** Each collapsed break's display unit. */
  breaks: number[] = [];
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

/** The rail's thickness at rest: a line a hint, the window a segment twice it. */
const Line = 1;
const Rest = 2;
/** A grip's height over the segment's. */
const GripRise = 3;
/** An idle window's opacity: quiet. */
const IdleAlpha = 0.4;

/** The display unit at x on a rail (the whole axis across the data box): what a drag on it reads. */
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
): void {
  rail.section = section.address;
  rail.figure = figure;
  rail.id = sectionId;
  rail.system = section.axis.system;
  rail.x = axis.x0;
  rail.width = axis.width;
  rail.slim = !ruler;
  if (ruler) {
    const columns = section.tracks.some((t) => t.layers.some((l) => l.type === 'rows' || l.type === 'logo'));
    rail.y = ruler.y + ruler.height;
    rail.bottom = ruler.y + ruler.height;
    rail.top = columns ? rail.bottom - ruler.height * config.navigator.rail.grabShare : ruler.y;
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
    rail.start = Math.min(a ?? 0, b ?? total);
    rail.end = Math.max(a ?? 0, b ?? total);
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

/** The rect in `out.args` (x, y, width, height, fill), written unless empty (numbers through the array: never boxed). */
function emit(out: PassWriter): void {
  const g = out.args;
  if (!(g[2] > 0) || !(g[3] > 0)) return;
  g[5] = 0;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/**
 * The window pass: the rail's line, its breaks' notches, the window's segment and its grips (lit, grown); idle a
 * hairline at a quiet opacity, no grips (at the fit it reads as the line); a figure's a plain segment, no grips.
 */
export function drawRail(
  rail: Rail,
  out: PassWriter,
  emphasis: Emphasis | undefined,
  tokens: SeqTokens,
  config: RenderConfig,
) {
  const g = out.args;
  const line = tokens.colors['msui-line'];
  g[0] = rail.x;
  g[1] = rail.y - Line / 2;
  g[2] = rail.width;
  g[3] = Line;
  g[4] = line;
  emit(out);
  for (let i = 0; i < rail.breaks.length; i++) {
    g[0] = rail.x + (rail.breaks[i] / rail.total) * rail.width - 0.5;
    g[1] = rail.y - 2;
    g[2] = 1;
    g[3] = 4;
    g[4] = line;
    emit(out);
  }
  if (!(rail.total > 0)) return;
  const h = emphasis?.hover;
  const d = emphasis?.drag;
  const lit = h?.target === 'window' && h.navigator === 'rail' && h.section === rail.id ? h.part : undefined;
  const dragged = d?.navigator === 'rail' && d.section === rail.id;
  const plain = rail.active || rail.figure;
  const grown = rail.active && (lit !== undefined || dragged);
  const t = grown ? config.navigator.rail.grownPx : plain ? Rest : Line;
  const ink = tokens.colors['seq-window'];
  const color = plain ? ink : ((ink & 0xffffff00) | Math.round(IdleAlpha * 255)) >>> 0;
  windowEnds(rail, config);
  const wraps = rail.end > rail.total;
  // over a cyclic origin the second piece from the axis's start to the window's true end
  const b = wraps ? rail.x + ((rail.end - rail.total) / rail.total) * rail.width : Ends[1];
  for (let k = 0; k < (wraps ? 2 : 1); k++) {
    g[0] = k === 0 ? Ends[0] : rail.x;
    g[1] = rail.y - t / 2;
    g[2] = k === 0 ? Ends[1] - Ends[0] : b - rail.x;
    g[3] = t;
    g[4] = color;
    emit(out);
  }
  if (!rail.active) return;
  const grip = tokens.colors['seq-grip'];
  const part = dragged ? d.part : lit;
  const gh = t + 2 * GripRise;
  for (let k = 0; k < 2; k++) {
    const on = part === 'body' || (k === 0 ? part === 'start' : part === 'end');
    g[0] = k === 0 ? Ends[0] - 2 : b;
    g[1] = rail.y - gh / 2;
    g[2] = 2;
    g[3] = gh;
    g[4] = on ? ink : grip;
    emit(out);
  }
}

/**
 * What a point on the rail's grab zone is: its window's part, or what is outside it (`touch`: the wider hit);
 * undefined idle or elsewhere. Over a cyclic origin the window's two pieces meet the axis's ends as body: its grips are
 * at its true ends only.
 */
export function hitRail(rail: Rail, x: number, y: number, config: RenderConfig, touch = false): Hit | undefined {
  if (!rail.active || !(rail.total > 0) || y < rail.top || y >= rail.bottom) return undefined;
  const c = config.navigator.window;
  windowEnds(rail, config);
  let part = navigatorPart(x, Ends[0], Ends[1], rail.x, rail.x + rail.width, c, touch);
  if (rail.end > rail.total) {
    // the first piece runs to the axis's end: its end there is body
    if (part === 'end') part = 'body';
    if (part === 'outside') {
      const b = rail.x + ((rail.end - rail.total) / rail.total) * rail.width;
      part = navigatorPart(x, rail.x, b, rail.x, rail.x + rail.width, c, touch);
      if (part === 'start') part = 'body';
    }
  }
  if (!part) return undefined;
  const h = rail.hit;
  h.layer = rail.section;
  h.section = rail.id;
  h.system = rail.system;
  h.part = part;
  const u = Math.max(0, Math.min(rail.total - 1e-9, railUnitAt(rail, x)));
  h.position = Math.floor(positionAt(rail.shape, u, 'start') ?? u);
  return h;
}
