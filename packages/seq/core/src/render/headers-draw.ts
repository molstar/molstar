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
  HeaderMenuPart,
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
 * Whether a header's actions show (`header.actions.show`, _spec/README.md "Render", Track headers): `always`; on
 * `hover` only while the pointer is over the track's header (its label line, its slots, its states), the keyboard
 * cursor is on its header, or its menu is open, never on its data or a rail. A slot is hit only while they show.
 */
export function actionsShown(
  address: string,
  section: string,
  emphasis: Emphasis | undefined,
  config: RenderConfig,
): boolean {
  const show = config.header.actions.show;
  if (show === 'always') return true;
  if (show === 'never' || !emphasis) return false;
  const c = emphasis.cursor;
  if (c && c.part === 'header' && c.track === address) return true;
  // an open menu keeps its header's actions (the pointer on its way to the menu leaves the line)
  if (emphasis.menu === address) return true;
  const h = emphasis.hover;
  if (!h || h.section !== section || h.layer !== address) return false;
  const t = h.target;
  return t === 'header-label' || t === 'header-slot' || t === 'header-state' || t === 'section-header';
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
  track?: TrackLayout,
) {
  const a = line.actions;
  if (!a) return;
  const g = out.args;
  // a label the backdrop cuts ends in an ellipsis on it (the label stays where it is)
  // the track's label read here: a number passed in (a field's, or a sum) is boxed on each call
  const cut = track !== undefined && config.inset + track.labelWidth > a.x - 2;
  const grow = cut ? (line.ellipsis ?? 0) : 0;
  g[0] = a.x - 2 - grow;
  g[1] = a.y;
  g[2] = a.width + 2 + grow;
  g[3] = a.height;
  g[4] = tokens.colors['msui-bg'];
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
  if (cut) {
    g[1] = track!.labelY;
    g[2] = grow;
    g[3] = grow;
    g[4] = tokens.colors['msui-text-muted'];
    out.textArgs('…', AnchorStart, config.labelFont);
  }
  const slots = line.slots;
  for (let i = 0; i < slots.length; i++) drawSlot(out, slots[i], hover, tokens, config);
  if (line.menu) drawSlot(out, line.menu, hover, tokens, config);
}

function drawSlot(out: PassWriter, s: HeaderSlotBox, hover: Hit | undefined, tokens: SeqTokens, config: RenderConfig) {
  const g = out.args;
  // the hovered slot (the hit as laid out, or a copy of it)
  const hovered =
    hover !== undefined && hover.target === 'header-slot' && hover.state === s.part.name && hover.layer === s.hit.layer;
  const pressed = s.part.pressed === true;
  const ink = hovered || pressed
    ? tokens.colors['msui-accent']
    : s.part.enabled
      ? tokens.colors['msui-text-muted']
      : tokens.colors['msui-line'];
  const cx = s.x + s.width / 2;
  const cy = s.y + s.height / 2;
  if (pressed) {
    // pressed: an underline under the icon, as a rows layer's pressed toggle
    g[0] = cx - IconPx / 2;
    g[1] = s.y + s.height - 2;
    g[2] = IconPx;
    g[3] = 1.5;
    g[4] = ink;
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
  if (s.part === HeaderMenuPart) {
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
  if (s.part.icon) {
    iconArgs(out, s.part.icon, cx - IconPx / 2, cy - IconPx / 2, IconPx / 24, ink);
    return;
  }
  g[0] = cx;
  g[1] = cy + 4;
  g[2] = 0;
  g[3] = s.width;
  g[4] = ink;
  out.textArgs(s.part.label.charAt(0) || '?', AnchorMiddle, config.labelFont);
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
  // a host's DOM parts draw its parts over their boxes: the chrome draws none
  if (!line || line.hosted) return;
  // a stack's later rows: their labels only, the states drawn and hit once (on its first)
  if (pass === 'labels') {
    if (dy === 0) drawStates(out, line, s.tokens, s.config, 0);
  } else if (pass === 'tip' && dy === 0 && line.actions) {
    if (actionsShown(t.track.address, s.id, emphasis, s.config))
      drawActions(out, line, emphasis?.hover, s.tokens, s.config, t);
    else drawPressed(out, line, s.tokens, s.config);
  }
}

/** A section's header line: its states (labels pass) and its actions when they show (tip pass). */
export function drawSectionHeader(out: PassWriter, s: SectionLayout, pass: string, emphasis: Emphasis | undefined) {
  const line = s.header;
  if (!line || line.hosted) return;
  if (pass === 'labels') drawStates(out, line, s.tokens, s.config, 0);
  else if (pass === 'tip' && line.actions) {
    if (actionsShown(s.section.address, s.id, emphasis, s.config))
      drawActions(out, line, emphasis?.hover, s.tokens, s.config);
    else drawPressed(out, line, s.tokens, s.config);
  }
}

/** A line's pressed slots at rest (the others hidden with the hover): where they are, as the DOM parts keep them. */
function drawPressed(out: PassWriter, line: HeaderLine | SectionHeaderLine, tokens: SeqTokens, config: RenderConfig) {
  const slots = line.slots;
  for (let i = 0; i < slots.length; i++) if (slots[i].part.pressed) drawSlot(out, slots[i], undefined, tokens, config);
}

/** Whether a section draws any header actions (the tip pass walks it). */
const Actions = new WeakMap<SectionLayout, boolean>();
export function hasHeaderActions(s: SectionLayout): boolean {
  let has = Actions.get(s);
  if (has === undefined) {
    has = s.header?.actions !== undefined && !s.header.hosted;
    for (const t of s.tracks) if (t.header?.actions && !t.header.hosted) has = true;
    Actions.set(s, has);
  }
  return has;
}

/** The emphasis a hit test reads (the hover alone), one record (a hit makes no object). */
const HoverOnly: { hover?: Hit } = {};

/**
 * What a header part is under a point (_spec/README.md "Render", Track headers: Hit): a slot or the "⋯"
 * (`header-slot`), a state (`header-state`), elsewhere on a track's line `header-label`, on a section's line
 * `section-header`. A mouse's point is on its own line, so its slots hit whatever the last hover was; a finger's only
 * while its line's actions show under `hover` (a tap on the label reveals them, the next tap runs one). The hit is the
 * one laid out (none made).
 */
export function hitHeaders(
  s: SectionLayout,
  x: number,
  y: number,
  hover: Hit | undefined,
  touch = false,
): Hit | undefined {
  const config = s.config;
  if (x < config.inset || x >= config.inset + config.labelWidth) return undefined;
  HoverOnly.hover = hover;
  const mouse = !touch && config.header.actions.show !== 'never';
  const sh = s.header;
  if (sh && y >= sh.line.y && y < sh.line.y + sh.line.height)
    return hitLine(sh, x, y, mouse || actionsShown(s.section.address, s.id, HoverOnly, config)) ?? sh.hit;
  for (let i = 0; i < s.tracks.length; i++) {
    const t = s.tracks[i];
    const line = t.header;
    if (!line || y < line.line.y || y >= line.line.y + line.line.height) continue;
    const h = hitLine(line, x, y, mouse || actionsShown(t.track.address, s.id, HoverOnly, config));
    if (h) return h;
    if (line.hit) return line.hit;
  }
  return undefined;
}

function hitLine(line: HeaderLine | SectionHeaderLine, x: number, y: number, shown: boolean): Hit | undefined {
  const slots = line.slots;
  // a pressed slot is drawn at rest, so it is hit there
  for (let i = 0; i < slots.length; i++) if ((shown || slots[i].part.pressed) && inBox(slots[i], x, y)) return slots[i].hit;
  if (shown && line.menu && inBox(line.menu, x, y)) return line.menu.hit;
  const states = line.states;
  for (let i = 0; i < states.length; i++) if (inBox(states[i].box, x, y)) return states[i].hit;
  return undefined;
}
