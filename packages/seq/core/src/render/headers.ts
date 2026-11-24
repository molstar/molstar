/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Node, resolveAddress } from '@molstar/core/repr/tree';
import { idRunOf } from '@molstar/seqx-core/lang/check';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import {
  type ResolvedLayer,
  type ResolvedSection,
  type ResolvedTrack,
  type ResolvedTrackGroup,
  tracksOf,
} from '@molstar/seqx-core/view';
import { fontSize, type RenderConfig } from './config.js';
import { groupCount } from './groups-hooks.js';
import type { Box, Hit } from './layer-type.js';
import type { TextMetrics } from './text.js';

// The headers' line (_spec/README.md "Render", Track headers): what a track's header and a section's header hold in the
// label column and where. The parts are given (`LayoutInput.headers`: a host's entries made into slots and states once
// a layout; the built-in states when absent); the layout places them in fixed columns, each with the hit it gives.
// Drawn and hit by the tracks layout's chrome (`headers-draw.ts`).

/**
 * An icon as path data in a 24-unit box (Lucide's): `PathCommand` codes (move, line, cubic, close) and their points'
 * x, y pairs; stroked in the label ink, scaled to the slot's icon.
 */
export interface HeaderIcon {
  readonly commands: ArrayLike<number>;
  readonly points: ArrayLike<number>;
}

/** An action entry's slot as given: its entry's name, its label, its icon (none: the label's first letter), enabled. */
export interface HeaderSlotPart {
  readonly name: string;
  readonly label: string;
  readonly icon?: HeaderIcon;
  readonly enabled: boolean;
  /** The action is on now (the track the structure is colored by): drawn pressed, also at rest. */
  readonly pressed?: boolean;
}

/** A state entry's chip as given: its entry's name, its text ("12 filtered"), what it counts and what it is about. */
export interface HeaderStatePart {
  readonly name: string;
  readonly text: string;
  readonly count?: number;
  /** The layers it is about when not the header's own target (the filtered ones: what "n filtered" clears). */
  readonly addresses?: readonly string[];
  /**
   * A state beside the label that never cuts it (a state that appears only with a hide: "1 layer hidden"): a chip where
   * it fits beside the whole label, else a mark where that fits, else none (its words are the legend's and the menus').
   */
  readonly yields?: true;
}

/**
 * A header's parts: its action slots in slot order, its states, and whether it has a menu ("⋯": any entry applies);
 * `hosted`: a host's DOM parts draw them over the laid-out boxes, so the chrome draws none of them (a figure draws its
 * states still).
 */
export interface HeaderParts {
  readonly slots?: readonly HeaderSlotPart[];
  readonly states?: readonly HeaderStatePart[];
  readonly menu?: boolean;
  readonly hosted?: boolean;
}

/**
 * What the layout knows of a header (a track's; a section's when `track` is absent), DOM-free: the layers a filter
 * thinned (by key, with the rows it left out), the hidden layers (their keys), the section's hidden tracks.
 */
export interface HeaderInputFacts {
  readonly section: ResolvedSection;
  readonly track?: ResolvedTrack;
  /** A folded group's row (its first band): the header acts on the group, its line's address the group's. */
  readonly group?: ResolvedTrackGroup;
  readonly filtered: readonly { readonly layer: string; readonly count: number }[];
  readonly hiddenLayers: readonly string[];
  readonly hiddenTracks: readonly HiddenTrack[];
}

/** A track the spec hides (`hidden: true`): left out of the view; the application lists it (`viewer.hiddenTracks`). */
export interface HiddenTrack {
  readonly address: string;
  readonly label: string;
}

/** A header's parts made from its facts, once a layout a header (a host's entries; the built-in states by default). */
export type HeaderPartsOf = (facts: HeaderInputFacts) => HeaderParts | undefined;

/** The built-in states' entry names. */
export const HeaderStateNames = {
  filtered: 'seq:filtered',
  hiddenLayers: 'seq:hidden-layers',
} as const;

/** The "⋯"'s name in its hit (`header-slot`). */
export const HeaderMenuName = 'seq:menu';

/**
 * The built-in states (`LayoutInput.headers` absent): a track's "n filtered" (the rows its layers' filters left out, every
 * thinned layer's) and "n layers hidden". A section has none: its hidden tracks are the application's, never drawn by the
 * view (_spec/README.md, Headers). No slots and no menu: those are a host's entries.
 */
export const builtInHeaderStates: HeaderPartsOf = (f) => {
  let states: HeaderStatePart[] | undefined;
  if (f.track) {
    let rows = 0;
    for (const t of f.filtered) rows += t.count;
    if (f.filtered.length)
      (states ??= []).push({
        name: HeaderStateNames.filtered,
        text: `${rows} filtered`,
        count: rows,
        addresses: f.filtered.map((t) => t.layer),
      });
    const n = f.hiddenLayers.length;
    if (n)
      (states ??= []).push({
        name: HeaderStateNames.hiddenLayers,
        text: n === 1 ? '1 layer hidden' : `${n} layers hidden`,
        count: n,
        yields: true,
      });
  }
  return states ? { states } : undefined;
};

/** A slot as laid out: its part, its box (its column on the header's line) and the hit it gives. */
export interface HeaderSlotBox extends Box {
  readonly part: HeaderSlotPart;
  readonly hit: Hit;
}

/** The "⋯"'s part. */
export const HeaderMenuPart: HeaderSlotPart = { name: HeaderMenuName, label: 'Actions', enabled: true };

/**
 * A state as laid out: a chip (its text from `x` at the baseline `y`, `width` its measured width) or, where no chip
 * fits (a short track, a column too narrow), a mark (an accent squircle `width` px square from `x`, `y` its top);
 * `box` what is hit, `hit` what it gives.
 */
export interface HeaderStateBox {
  readonly name: string;
  readonly text: string;
  readonly count?: number;
  readonly addresses?: readonly string[];
  readonly mark: boolean;
  readonly x: number;
  readonly y: number;
  readonly width: number;
  readonly box: Box;
  readonly hit: Hit;
}

/**
 * A header's line (_spec/README.md "Render", Track headers), in the canvas's CSS px as the hits are: its box in the
 * label column, its label's box (the room the states leave it), its slots each in its fixed column, the "⋯" when it has
 * a menu, its states; `actions` the box the slots and the "⋯" take (their backdrop), absent when it draws none.
 */
export interface HeaderLine {
  /** A folded group's row: the group's address (its line acts on the group, not on the band it sits on). */
  readonly group?: string;
  readonly line: Box;
  readonly label: Box;
  readonly slots: readonly HeaderSlotBox[];
  readonly menu?: HeaderSlotBox;
  readonly states: readonly HeaderStateBox[];
  readonly actions?: Box;
  /** The width of an ellipsis in the label font: what the actions' backdrop makes room for when it cuts a label. */
  readonly ellipsis?: number;
  /** Shorter than `header.minHeightPx`: no slots, its states as marks. */
  readonly short: boolean;
  /** A host's DOM parts draw its slots, its "⋯" and its states over their boxes: the chrome draws none of them. */
  readonly hosted?: boolean;
  /**
   * What the line gives outside its slots and states: a track's `header-label`, a section's `section-header`; none on
   * a rows track (its label column is its tree and rows).
   */
  readonly hit?: Hit;
}

/** A section's header line (`section.header`): a line above its first track, its own hit `section-header`. */
export interface SectionHeaderLine extends HeaderLine {
  readonly hit: Hit;
}

const NoSlots: readonly HeaderSlotBox[] = [];
const NoStates: readonly HeaderStateBox[] = [];
const ChipGap = 6;
const MarkPx = 6;
const MarkGap = 4;

/**
 * The action columns a view lays out with (_spec/README.md "Render", Track headers): `header.slots` columns and the
 * "⋯", `header.slotPx` each, ending at the label column's right end, as many slots as leave three columns' width of
 * label; none with `header.actions.show: never`. One set for every header of the view, counted from the "⋯" leftwards
 * (the first slot next to it): every track's eye in one column, and a line reserves only the slots it has.
 */
export function actionColumns(config: RenderConfig): ActionColumns {
  let c = Columns.get(config);
  if (c) return c;
  const h = config.header;
  const end = config.inset + config.labelWidth;
  const n =
    h.actions.show === 'never' || !(h.slotPx > 0)
      ? -1
      : Math.max(0, Math.min(h.slots, Math.floor(config.labelWidth / h.slotPx) - 3));
  c = { x: n < 0 ? end : end - (n + 1) * h.slotPx, n, size: fontSize(config.labelFont) };
  Columns.set(config, c);
  return c;
}

/** The action columns' x, how many slots (−1: no actions), and the label font's size a header reads with them. */
export interface ActionColumns {
  readonly x: number;
  readonly n: number;
  readonly size: number;
}

/** By config: a layout reads them once a header (nothing made). */
const Columns = new WeakMap<RenderConfig, ActionColumns>();
const labelSizeOf = (config: RenderConfig) => actionColumns(config).size;

/** Whether a header lays out actions: a slot or a menu given, not short, not a figure, `header.actions.show` not `never`. */
function hasActions(parts: HeaderParts | undefined, short: boolean, figure: boolean, config: RenderConfig): boolean {
  return (
    !short &&
    !figure &&
    actionColumns(config).n >= 0 &&
    !!(parts && ((parts.slots && parts.slots.length) || parts.menu))
  );
}

/** What a header line is laid out from: its hits' target, its line, whether it is short, its states' own line. */
interface LineInput {
  layer: string;
  section: string;
  /** The label line's centre. */
  cy: number;
  /** The line's top and height (the label column's part of it). */
  top: number;
  height: number;
  short: boolean;
  /** States on a line of their own under the label (a track of two lines' height): their baseline; else NaN. */
  ownLine: number;
  figure: boolean;
  /** A track's line hit as `header-label` outside its parts (not a rows track's, not a section's). */
  labelHit: boolean;
  /** The label's own width (its text and indent): the states beside it leave it whole; 0 where there is none. */
  labelW: number;
}

/** The one input a layout's header lines are laid out from, set field by field (a header makes none). */
const At: LineInput = {
  layer: '',
  section: '',
  cy: 0,
  top: 0,
  height: 0,
  short: false,
  ownLine: Number.NaN,
  figure: false,
  labelHit: false,
  labelW: 0,
};

/**
 * A header's line laid out from its parts: the slots in their columns (the first `n` given, the first next to the "⋯"),
 * the "⋯" at the column's end, the states as chips ending where the line's own slots begin (the "⋯" on a line with
 * none; from the column's start on a line of their own), a state no chip fits a mark. A figure has no slots and no
 * "⋯", and no states with `header.states.export: none`.
 */
function layoutLine(
  parts: HeaderParts | undefined,
  at: LineInput,
  config: RenderConfig,
  text: TextMetrics,
): HeaderLine {
  const h = config.header;
  const cols = actionColumns(config);
  const inset = config.inset;
  const end = inset + config.labelWidth;
  const size = labelSizeOf(config);
  const slotPx = h.slotPx;
  const sy = at.cy - slotPx / 2;
  let slots = NoSlots;
  let menu: HeaderSlotBox | undefined;
  let actions: Box | undefined;
  if (parts && hasActions(parts, at.short, at.figure, config)) {
    const given = parts.slots ?? [];
    const n = Math.min(cols.n, given.length);
    const out: HeaderSlotBox[] = [];
    for (let i = 0; i < n; i++) {
      const p = given[i];
      out.push({
        part: p,
        x: end - (i + 2) * slotPx,
        y: sy,
        width: slotPx,
        height: slotPx,
        hit: { layer: at.layer, section: at.section, system: '', target: 'header-slot', state: p.name },
      });
    }
    if (n) slots = out;
    if (parts.menu)
      menu = {
        part: HeaderMenuPart,
        x: end - slotPx,
        y: sy,
        width: slotPx,
        height: slotPx,
        hit: { layer: at.layer, section: at.section, system: '', target: 'header-slot', state: HeaderMenuName },
      };
    // the line's own slots and the "⋯" column (a slotted line always has a menu: every entry is in it)
    const x = end - (n + 1) * slotPx;
    actions = { x, y: sy, width: end - x, height: slotPx };
  }
  const given = at.figure && h.states.export === 'none' ? undefined : parts?.states;
  let states = NoStates;
  // the label's room: the column, less the states beside it
  let labelEnd = end;
  if (given?.length) {
    const out: (HeaderStateBox & { readonly at: number })[] = [];
    // on the label's line the states end before the line's own actions (a gap before them) when it has any, else at
    // the column's end; on a line of their own at the column's end (the actions are on the label's line)
    const chipsEnd = actions ? actions.x - ChipGap : end;
    // states `[from, to)` as marks, right to left from `right`, centred on `cy`: the leftmost's x; a mark with no room
    // left in the column is not laid out (its words are the card's)
    const marks = (from: number, to: number, right: number, cy: number, floor = inset): number => {
      let x = right;
      let left = 0;
      for (let i = to - 1; i >= from; i--) {
        // a state that yields to the label has no mark where the label would be cut; the others cut it as before
        if (x - MarkPx < (given[i].yields ? floor : inset)) {
          if (given[i].yields) {
            left++;
            continue;
          }
          break;
        }
        x -= MarkPx;
        out.push(stateBox(given[i], i, at, true, x, cy - MarkPx / 2, MarkPx, size));
        x -= MarkGap;
      }
      // none placed and every one yielded: the label keeps the column
      return x === right && left === to - from ? end + ChipGap : x + MarkGap;
    };
    // beside the label a chip keeps it whole (a chip with no room is a mark), and a state that yields to it appears only
    // where it leaves it whole (`HeaderStatePart.yields`); the label's own room is its column, under the actions
    const keep = at.labelW > 0 ? inset + at.labelW + ChipGap : inset;
    if (at.short) labelEnd = marks(0, given.length, chipsEnd, at.cy, keep) - ChipGap;
    else if (!Number.isNaN(at.ownLine)) {
      // a line of their own under the label, from the column's start, under the action columns when clear of them;
      // all chips when every one fits, else all marks (a mix reads as two kinds of state: "a narrow column draws its
      // states as marks too", _spec/README.md "Render", Track headers)
      const own = actions && at.ownLine - size + 1 < actions.y + actions.height ? chipsEnd : end;
      const widths = given.map((g) => text.measure(config.labelFont, g.text));
      let total = 0;
      for (const w of widths) total += w + ChipGap;
      if (inset + total - ChipGap <= own) {
        let x = inset;
        for (let i = 0; i < given.length; i++) {
          out.push(stateBox(given[i], i, at, false, x, at.ownLine, widths[i], size));
          x += widths[i] + ChipGap;
        }
      } else marks(0, given.length, own, at.ownLine - size * 0.35);
    } else {
      // beside the label, the last state nearest the columns; the leading ones with no room left marks before them
      let x = chipsEnd;
      let cut = 0;
      for (let i = given.length - 1; i >= 0; i--) {
        const w = text.measure(config.labelFont, given[i].text);
        if (x - w < keep) {
          cut = i + 1;
          break;
        }
        x -= w;
        out.push(stateBox(given[i], i, at, false, x, at.cy + size * 0.35, w, size));
        x -= ChipGap;
      }
      labelEnd = cut ? marks(0, cut, x, at.cy, keep) - ChipGap : x;
    }
    // in the given order
    out.sort((a, b) => a.at - b.at);
    states = out;
  }
  const out: { -readonly [K in keyof HeaderLine]: HeaderLine[K] } = {
    line: { x: inset, y: at.top, width: config.labelWidth, height: at.height },
    label: { x: inset, y: at.cy - size / 2 - 2, width: Math.max(0, labelEnd - inset), height: size + 4 },
    slots,
    states,
    short: at.short,
  };
  if (menu) out.menu = menu;
  if (actions) {
    out.actions = actions;
    out.ellipsis = text.measure(config.labelFont, '…');
  }
  if (parts?.hosted && !at.figure) out.hosted = true;
  if (at.labelHit) out.hit = { layer: at.layer, section: at.section, system: '', target: 'header-label' };
  return out;
}

function stateBox(
  p: HeaderStatePart,
  index: number,
  at: LineInput,
  mark: boolean,
  x: number,
  y: number,
  width: number,
  size: number,
): HeaderStateBox & { readonly at: number } {
  // a chip's box: its text with 2 px around it, a line high; a mark's: its squircle with 2 px around it
  const box: Box = mark
    ? { x: x - 2, y: y - 2, width: width + 4, height: width + 4 }
    : { x: x - 2, y: y - size + 1, width: width + 4, height: size + 4 };
  return {
    at: index,
    name: p.name,
    text: p.text,
    ...(p.count !== undefined ? { count: p.count } : {}),
    ...(p.addresses !== undefined ? { addresses: p.addresses } : {}),
    mark,
    x,
    y,
    width,
    box,
    hit: {
      layer: at.layer,
      section: at.section,
      system: '',
      target: 'header-state',
      state: p.name,
      ...(p.count !== undefined ? { count: p.count } : {}),
    },
  };
}

/**
 * A track's header line (_spec/README.md "Render", Track headers) over its well, its label centred on it: its states on
 * a line of their own when it is two lines high (the label then whole, the two lines centred together), else beside
 * the label; returns the line and the label's baseline.
 */
export function trackHeaderLine(
  parts: HeaderParts | undefined,
  track: string,
  section: string,
  well: Box,
  config: RenderConfig,
  text: TextMetrics,
  figure: boolean,
  labelHit = true,
  labelW = 0,
): { readonly header: HeaderLine; readonly labelY: number } {
  const size = labelSizeOf(config);
  const short = well.height < config.header.minHeightPx;
  // a line of their own when two lines fit (a chip clear of the slots above it when the header has actions)
  const two = !short && parts?.states?.length && well.height >= 2 * (size + 6);
  const dy = two
    ? hasActions(parts, short, figure, config)
      ? Math.max(size + 4, Math.ceil(config.header.slotPx / 2 + 0.65 * size - 1))
      : size + 4
    : 0;
  const centre = well.y + well.height / 2;
  const labelY = centre + size * 0.35 - dy / 2;
  const at = At;
  at.layer = track;
  at.section = section;
  at.cy = labelY - size * 0.35;
  at.top = well.y;
  at.height = well.height;
  at.short = short;
  at.ownLine = dy ? labelY + dy : Number.NaN;
  at.figure = figure;
  at.labelHit = labelHit;
  at.labelW = labelW;
  const out = TrackLine;
  out.header = layoutLine(parts, at, config, text);
  out.labelY = labelY;
  return out;
}

/** `trackHeaderLine`'s answer, one record (read at once by its caller). */
const TrackLine: { header: HeaderLine; labelY: number } = { header: undefined as never, labelY: 0 };

/**
 * Whether a section draws its header line (`section.header`): `auto` when a section entry gives it a slot (a state or a
 * menu-only entry never draws one), `always`, `none`.
 */
export function drawsSectionHeader(parts: HeaderParts | undefined, config: RenderConfig): boolean {
  const mode = config.section.header;
  if (mode === 'none') return false;
  if (mode === 'always') return true;
  return !!parts?.slots?.length;
}

/** A section's header line, `header.sectionPx` high from `top`: its states and its slots. */
export function sectionHeaderLine(
  parts: HeaderParts | undefined,
  section: ResolvedSection,
  id: string,
  top: number,
  config: RenderConfig,
  text: TextMetrics,
  figure: boolean,
): SectionHeaderLine {
  const h = config.header.sectionPx;
  const at = At;
  at.layer = section.address;
  at.section = id;
  at.cy = top + h / 2;
  at.top = top;
  at.height = h;
  at.short = false;
  at.ownLine = Number.NaN;
  at.figure = figure;
  at.labelHit = false;
  at.labelW = 0;
  const line = layoutLine(parts, at, config, text);
  return {
    ...line,
    hit: { layer: section.address, section: id, system: '', target: 'section-header' },
  };
}

/**
 * Whether a section has collapsed: no track left in the view (every one hidden, the application lists them), no blocks
 * and no grid axis. It takes no room (no gap, line or rail), and the navigator rails, the cursor and the focus skip it.
 */
export function collapsedSection(section: ResolvedSection): boolean {
  return !section.tracks.length && !section.blocks && !section.y;
}

/**
 * The name a hidden track is listed by (_spec/README.md, Hidden tracks): its label (its id without one), else, for an empty label, its first layer's `label` or
 * `title`, else its layer's type as a word ("Ruler", "Rows"), else its id. A header's label may be empty (a ruler's, the
 * rows'): the list has no header to read it from.
 */
export function hiddenTrackName(track: Node): string {
  const text = (v: unknown) => (typeof v === 'string' && v.trim() ? v : undefined);
  const own = track.params?.label;
  // no label: the header reads the track's id; an empty one is the author's choice of none
  if (own === undefined) return track.id ?? 'Track';
  if (text(own)) return own as string;
  for (const l of track.children ?? []) {
    if (l.kind !== 'layer') continue;
    const named = text(l.params?.label) ?? text(l.params?.title);
    if (named) return named;
    const type = text(l.params?.type);
    if (type) return type[0].toUpperCase() + type.slice(1);
  }
  return track.id ?? 'Track';
}

/**
 * A layer's name, the one every list of layers says (a track header's Layers and "Hidden layers" menus, a legend's eye):
 * its `label`, else its `title`, else its track's or bridge's label where it is the only layer there (a letters layer
 * says its kind: "Sequence letters", never doubled when the label ends with it ("letters"), for its key is a scheme and the track's name alone reads as one category's), else
 * its id (its type without one). Why one: two consumers naming one layer apart read as two layers.
 */
export function layerName(
  layer: ResolvedLayer,
  owner?: { readonly label?: string; readonly layers: readonly unknown[] },
): string {
  const text = (v: unknown) => (typeof v === 'string' && v.trim() ? v : undefined);
  const label = text(owner?.label);
  const owned = owner?.layers.length === 1 ? label : undefined;
  return (
    text(layer.params.label) ??
    text(layer.params.title) ??
    // a letters layer says its kind wherever it is ("Sequence letters" with a second layer beside it too)
    (layer.type === 'letters' && label !== undefined ? (/\bletters$/i.test(label.trim()) ? label : `${label} letters`) : owned) ??
    layer.address?.slice(layer.address.lastIndexOf('/') + 1) ??
    layer.type[0].toUpperCase() + layer.type.slice(1)
  );
}

const Hidden = new WeakMap<ResolvedSection, readonly HiddenTrack[]>();
const NoHidden: readonly HiddenTrack[] = [];

/** A section's hidden tracks in spec order (each once), read from the spec its view was resolved from; kept a section. */
export function hiddenTracksOf(model: StaticModel, section: ResolvedSection): readonly HiddenTrack[] {
  let out = Hidden.get(section);
  if (out) return out;
  let found: HiddenTrack[] | undefined;
  let group: string | undefined;
  for (const { node: t, address, group: g } of tracksOf(
    resolveAddress(model.spec, idRunOf(section.address)),
    section.address,
  )) {
    // a hidden group: one entry, by its label (its `track.show` shows the group)
    if (g && g.node.params?.hidden === true) {
      if (group === g.address) continue;
      group = g.address;
      const label = g.node.params?.label;
      // its tracks counted, as the host's list says it ("Variants (3 tracks)")
      const n = (g.node.children ?? []).filter((x) => x.kind === 'track').length;
      const name = typeof label === 'string' ? label : (g.node.id ?? '');
      (found ??= []).push({ address: g.address, label: `${name} (${groupCount(n)})` });
      continue;
    }
    if (t.params?.hidden !== true) continue;
    (found ??= []).push({ address, label: hiddenTrackName(t) });
  }
  out = found ?? NoHidden;
  Hidden.set(section, out);
  return out;
}
