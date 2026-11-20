/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RenderConfig } from './config.js';
import { Anchor, grown, NoColor, Op, type PassWriter, PathCommandPoints } from './draw-lists.js';
import {
  type HeaderIcon,
  type HeaderLine,
  HeaderMenuName,
  type HeaderSlotBox,
  type SectionHeaderLine,
} from './headers.js';
import type { Emphasis, Hit } from './layer-type.js';
import type { SeqTokens } from './tokens.js';
import type { SectionLayout, TrackLayout } from './tracks-types.js';

// The headers drawn and hit (_spec/README.md "Render", Track headers): the states in the labels pass (facts about the
// picture: always drawn, and exported), the slots and the "⋯" over a backdrop in the tip pass, over the labels, for the
// header whose actions show (an emphasis redraw: no layout); each part hit by the box it was laid out with, its hit
// the one it was laid out with (a hover makes no object).

// the anchor a number read once: a draw reads no property of an imported object (_spec/README.md "Render", Draw lists)
const AnchorStart: number = Anchor.start;
const AnchorMiddle: number = Anchor.middle;

/** A header's icon, 14 px in a 20 px slot: its 24-unit box scaled. */
const IconPx = 14;

/** Whether a box holds a point. */
function inBox(b: { x: number; y: number; width: number; height: number }, x: number, y: number): boolean {
  return x >= b.x && x < b.x + b.width && y >= b.y && y < b.y + b.height;
}

/**
 * Whether a header's actions show (`header.actions.show`): `always`; on `hover` the header the viewer names
 * (`emphasis.header`), the track a hover's hit is in (its data, its slots, its states: `within`; a section's line only
 * its own hits), the cursor's header.
 */
export function actionsShown(
  address: string,
  section: string,
  emphasis: Emphasis | undefined,
  config: RenderConfig,
  within = true,
): boolean {
  const show = config.header.actions.show;
  if (show === 'always') return true;
  if (show === 'never' || !emphasis) return false;
  if (emphasis.header === address) return true;
  const c = emphasis.cursor;
  if (c && c.part === 'header' && c.track === address) return true;
  const h = emphasis.hover;
  if (!h || h.section !== section) return false;
  const l = h.layer;
  if (l === address) return true;
  if (!within) return false;
  // a layer's key under the track: its address (`<track>/<id>`) or its index (`<track>:<n>`)
  const n = address.length;
  if (l.length <= n || !l.startsWith(address)) return false;
  const c0 = l.charCodeAt(n);
  return c0 === 47 || c0 === 58;
}

/** The states of a header line in the labels pass: a chip its text in the accent, a mark an accent squircle. */
export function drawStates(out: PassWriter, line: HeaderLine, tokens: SeqTokens, config: RenderConfig, dy: number) {
  const g = out.args;
  const states = line.states;
  for (let i = 0; i < states.length; i++) {
    const s = states[i];
    if (s.mark) {
      g[0] = s.x;
      g[1] = s.y + dy;
      g[2] = s.width;
      g[3] = s.width;
      g[4] = tokens.colors['msui-accent'];
      g[5] = NoColor;
      g[6] = 0;
      g[7] = s.width / 3;
      g[8] = 3;
      out.rectArgs();
      continue;
    }
    g[0] = s.x;
    g[1] = s.y + dy;
    g[2] = s.width;
    g[3] = s.width;
    g[4] = tokens.colors['msui-accent'];
    out.textArgs(s.text, AnchorStart, config.labelFont);
  }
}

/**
 * A header line's actions in the tip pass: a backdrop in the page's color over the label's tail (the label stays where
 * it is), each slot's icon (its label's first letter without one) in the muted ink, an entry not enabled fainter, the
 * hovered slot in the accent, and the "⋯".
 */
export function drawActions(
  out: PassWriter,
  line: HeaderLine,
  hover: Hit | undefined,
  tokens: SeqTokens,
  config: RenderConfig,
) {
  const a = line.actions;
  if (!a) return;
  const g = out.args;
  g[0] = a.x - 2;
  g[1] = a.y;
  g[2] = a.width + 2;
  g[3] = a.height;
  g[4] = tokens.colors['msui-bg'];
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
  const slots = line.slots;
  for (let i = 0; i < slots.length; i++) drawSlot(out, slots[i], hover, tokens, config);
  if (line.menu) drawSlot(out, line.menu, hover, tokens, config);
}

function drawSlot(out: PassWriter, s: HeaderSlotBox, hover: Hit | undefined, tokens: SeqTokens, config: RenderConfig) {
  const g = out.args;
  // the hovered slot (the hit as laid out, or a copy of it)
  const hovered =
    hover !== undefined && hover.target === 'header-slot' && hover.state === s.name && hover.layer === s.hit.layer;
  const ink = hovered
    ? tokens.colors['msui-accent']
    : s.enabled
      ? tokens.colors['msui-text-muted']
      : tokens.colors['msui-line'];
  const cx = s.x + s.width / 2;
  const cy = s.y + s.height / 2;
  if (s.name === HeaderMenuName) {
    // three dots across the slot's middle
    for (let k = -1; k <= 1; k++) {
      g[0] = cx + k * 4 - 1;
      g[1] = cy - 1;
      g[2] = 2;
      g[3] = 2;
      g[4] = ink;
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 1;
      g[8] = 3;
      out.rectArgs();
    }
    return;
  }
  if (s.icon) {
    iconArgs(out, s.icon, cx - IconPx / 2, cy - IconPx / 2, IconPx / 24, ink);
    return;
  }
  g[0] = cx;
  g[1] = cy + 4;
  g[2] = 0;
  g[3] = s.width;
  g[4] = ink;
  out.textArgs(s.label.charAt(0) || '?', AnchorMiddle, config.labelFont);
}

/** An icon's path stroked at `(x, y)` scaled by `k` (its commands and points copied in: nothing allocated once warm). */
function iconArgs(out: PassWriter, icon: HeaderIcon, x: number, y: number, k: number, ink: number) {
  const cmds = icon.commands;
  const pts = icon.points;
  const nc = cmds.length;
  const np = pts.length >> 1;
  out.pass.op(Op.path);
  const p = out.pass.paths;
  p.begin();
  p.commands = grown(p.commands, p.commandCount + nc);
  p.points = grown(p.points, 2 * (p.pointCount + np));
  let q = 0;
  for (let i = 0; i < nc; i++) {
    const c = cmds[i];
    p.commands[p.commandCount++] = c;
    for (let m = 0; m < PathCommandPoints[c]; m++, q += 2) {
      p.points[2 * p.pointCount] = x + pts[q] * k;
      p.points[2 * p.pointCount + 1] = y + pts[q + 1] * k;
      p.pointCount++;
    }
  }
  const j = p.n++;
  p.cmdCount[j] = p.commandCount - p.cmdStart[j];
  p.fill[j] = NoColor;
  p.stroke[j] = ink;
  p.strokeWidth[j] = 1.5;
  p.dash[j] = 0;
}

/**
 * A track's header in the tracks layout's chrome: its states in the labels pass (`dy` below its first row: a stack's
 * later blocks), its actions in the tip pass when they show.
 */
export function drawTrackHeader(
  out: PassWriter,
  s: SectionLayout,
  t: TrackLayout,
  pass: string,
  emphasis: Emphasis | undefined,
  dy: number,
) {
  const line = t.header;
  if (!line) return;
  if (pass === 'labels') drawStates(out, line, s.tokens, s.config, dy);
  else if (pass === 'tip' && dy === 0 && line.actions && actionsShown(t.track.address, s.id, emphasis, s.config))
    drawActions(out, line, emphasis?.hover, s.tokens, s.config);
}

/** A section's header line: its states (labels pass) and its actions when they show (tip pass). */
export function drawSectionHeader(out: PassWriter, s: SectionLayout, pass: string, emphasis: Emphasis | undefined) {
  const line = s.header;
  if (!line) return;
  if (pass === 'labels') drawStates(out, line, s.tokens, s.config, 0);
  else if (pass === 'tip' && line.actions && actionsShown(s.section.address, s.id, emphasis, s.config, false))
    drawActions(out, line, emphasis?.hover, s.tokens, s.config);
}

/** Whether a section draws any header actions (the tip pass walks it). */
const Actions = new WeakMap<SectionLayout, boolean>();
export function hasHeaderActions(s: SectionLayout): boolean {
  let has = Actions.get(s);
  if (has === undefined) {
    has = s.header?.actions !== undefined;
    for (const t of s.tracks) if (t.header?.actions) has = true;
    Actions.set(s, has);
  }
  return has;
}

/**
 * What a header part is under a point (_spec/README.md "Render", Track headers: Hit): a slot or the "⋯"
 * (`header-slot`), a state (`header-state`), a section's line (`section-header`); the slots are hit whether drawn or
 * not (laid out in fixed columns: a pointer reaching for one shows them). The hit is the one laid out (none made).
 */
export function hitHeaders(s: SectionLayout, x: number, y: number): Hit | undefined {
  const config = s.config;
  if (x < config.inset || x >= config.inset + config.labelWidth) return undefined;
  const sh = s.header;
  if (sh && y >= sh.line.y && y < sh.line.y + sh.line.height) return hitLine(sh, x, y) ?? sh.hit;
  for (let i = 0; i < s.tracks.length; i++) {
    const line = s.tracks[i].header;
    if (!line || y < line.line.y || y >= line.line.y + line.line.height) continue;
    const h = hitLine(line, x, y);
    if (h) return h;
  }
  return undefined;
}

function hitLine(line: HeaderLine | SectionHeaderLine, x: number, y: number): Hit | undefined {
  const slots = line.slots;
  for (let i = 0; i < slots.length; i++) if (inBox(slots[i], x, y)) return slots[i].hit;
  if (line.menu && inBox(line.menu, x, y)) return line.menu.hit;
  const states = line.states;
  for (let i = 0; i < states.length; i++) if (inBox(states[i].box, x, y)) return states[i].hit;
  return undefined;
}
