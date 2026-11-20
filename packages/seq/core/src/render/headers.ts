/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Node, parseAddress } from '@molstar/core/repr/tree';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedSection, ResolvedTrack } from '@molstar/seqx-core/view';
import { fontSize, type RenderConfig } from './config.js';
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
}

/** A state entry's chip as given: its entry's name, its text ("12 filtered"), what it counts and what it is about. */
export interface HeaderStatePart {
  readonly name: string;
  readonly text: string;
  readonly count?: number;
  /** What the state is about when not the header's own target (the filtered layer). */
  readonly address?: string;
}

/** A header's parts: its action slots in slot order, its states, and whether it has a menu ("⋯": any entry applies). */
export interface HeaderParts {
  readonly slots?: readonly HeaderSlotPart[];
  readonly states?: readonly HeaderStatePart[];
  readonly menu?: boolean;
}

/**
 * What the layout knows of a header (a track's; a section's when `track` is absent), DOM-free: the layers a filter
 * thinned (by key, with the rows it left out), the hidden layers (their keys), the section's hidden tracks.
 */
export interface HeaderFacts {
  readonly section: ResolvedSection;
  readonly track?: ResolvedTrack;
  readonly filtered: readonly { readonly layer: string; readonly count: number }[];
  readonly hiddenLayers: readonly string[];
  readonly hiddenTracks: readonly HiddenTrack[];
}

/** A track the spec hides (`hidden: true`): left out of the view, counted in its section's header. */
export interface HiddenTrack {
  readonly address: string;
  readonly label: string;
}

/** A header's parts made from its facts, once a layout a header (a host's entries; the built-in states by default). */
export type HeaderPartsOf = (facts: HeaderFacts) => HeaderParts | undefined;

/** The built-in states' entry names. */
export const HeaderStateNames = {
  filtered: 'seq:filtered',
  hiddenLayers: 'seq:hidden-layers',
  hiddenTracks: 'seq:hidden-tracks',
} as const;

/** The "⋯"'s name in its hit (`header-slot`). */
export const HeaderMenuName = 'seq:menu';

/**
 * The built-in states (`LayoutInput.headers` absent): a track's "n filtered" (its first thinned layer) and "n layers
 * hidden", a section's "n hidden" (its hidden tracks). No slots and no menu: those are a host's entries.
 */
export const builtInHeaderStates: HeaderPartsOf = (f) => {
  let states: HeaderStatePart[] | undefined;
  if (f.track) {
    const t = f.filtered[0];
    if (t)
      (states ??= []).push({
        name: HeaderStateNames.filtered,
        text: `${t.count} filtered`,
        count: t.count,
        address: t.layer,
      });
    const n = f.hiddenLayers.length;
    if (n)
      (states ??= []).push({
        name: HeaderStateNames.hiddenLayers,
        text: n === 1 ? '1 layer hidden' : `${n} layers hidden`,
        count: n,
      });
  } else if (f.hiddenTracks.length) {
    const n = f.hiddenTracks.length;
    states = [{ name: HeaderStateNames.hiddenTracks, text: `${n} hidden`, count: n }];
  }
  return states ? { states } : undefined;
};

/** A slot as laid out: its part, its box (its column on the header's line) and the hit it gives. */
export interface HeaderSlotBox extends Box {
  readonly name: string;
  readonly label: string;
  readonly icon?: HeaderIcon;
  readonly enabled: boolean;
  readonly hit: Hit;
}

/**
 * A state as laid out: a chip (its text from `x` at the baseline `y`, `width` its measured width) or, where no chip
 * fits (a short track, a column too narrow), a mark (an accent squircle `width` px square from `x`, `y` its top);
 * `box` what is hit, `hit` what it gives.
 */
export interface HeaderStateBox {
  readonly name: string;
  readonly text: string;
  readonly count?: number;
  readonly address?: string;
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
  readonly line: Box;
  readonly label: Box;
  readonly slots: readonly HeaderSlotBox[];
  readonly menu?: HeaderSlotBox;
  readonly states: readonly HeaderStateBox[];
  readonly actions?: Box;
  /** Shorter than `header.minHeightPx`: no slots, its states as marks. */
  readonly short: boolean;
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
 * "⋯", `header.slotPx` each, ending at the label column's right end, as many slots as leave two columns' width of
 * label; none with `header.actions.show: never`. One set for every header of the view: every track's eye in one column.
 */
export function actionColumns(config: RenderConfig): { readonly x: number; readonly n: number } {
  const h = config.header;
  const end = config.inset + config.labelWidth;
  if (h.actions.show === 'never' || !(h.slotPx > 0)) return { x: end, n: -1 };
  const n = Math.max(0, Math.min(h.slots, Math.floor(config.labelWidth / h.slotPx) - 3));
  return { x: end - (n + 1) * h.slotPx, n };
}

/** What a header line is laid out from: its hits' target, its line, whether it is short, its states' own line. */
interface LineInput {
  readonly layer: string;
  readonly section: string;
  /** The label line's centre. */
  readonly cy: number;
  readonly line: Box;
  readonly short: boolean;
  /** States on a line of their own under the label (a track of two lines' height): their baseline. */
  readonly ownLine?: number;
  readonly figure: boolean;
}

/**
 * A header's line laid out from its parts: the slots in their columns (the first `n` given), the "⋯" at the column's
 * end, the states as chips ending where the action columns begin (from the column's start on a line of their own), a
 * state no chip fits a mark. A figure has no slots and no "⋯", and no states with `header.states.export: none`.
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
  const size = fontSize(config.labelFont);
  const slotPx = h.slotPx;
  const sy = at.cy - slotPx / 2;
  let slots = NoSlots;
  let menu: HeaderSlotBox | undefined;
  let actions: Box | undefined;
  if (!at.short && !at.figure && cols.n >= 0 && (parts?.slots?.length || parts?.menu)) {
    const given = parts.slots ?? [];
    const n = Math.min(cols.n, given.length);
    const out: HeaderSlotBox[] = [];
    for (let i = 0; i < n; i++) {
      const p = given[i];
      out.push({
        name: p.name,
        label: p.label,
        ...(p.icon ? { icon: p.icon } : {}),
        enabled: p.enabled,
        x: cols.x + i * slotPx,
        y: sy,
        width: slotPx,
        height: slotPx,
        hit: { layer: at.layer, section: at.section, system: '', target: 'header-slot', state: p.name },
      });
    }
    if (n) slots = out;
    if (parts.menu)
      menu = {
        name: HeaderMenuName,
        label: 'Actions',
        enabled: true,
        x: end - slotPx,
        y: sy,
        width: slotPx,
        height: slotPx,
        hit: { layer: at.layer, section: at.section, system: '', target: 'header-slot', state: HeaderMenuName },
      };
    actions = { x: cols.x, y: sy, width: end - cols.x, height: slotPx };
  }
  const given = at.figure && h.states.export === 'none' ? undefined : parts?.states;
  let states = NoStates;
  // the label's room: the column, less the states beside it
  let labelEnd = end;
  if (given?.length) {
    const out: (HeaderStateBox & { readonly at: number })[] = [];
    // on the label's line the states end before the action columns (a gap before them) when the header has actions,
    // else at the column's end; on a line of their own at the column's end (the actions are on the label's line)
    const chipsEnd = actions ? cols.x - ChipGap : end;
    // states `[from, to)` as marks, right to left from `right`, centred on `cy`: the leftmost's x
    const marks = (from: number, to: number, right: number, cy: number): number => {
      let x = right;
      for (let i = to - 1; i >= from; i--) {
        x -= MarkPx;
        out.push(stateBox(given[i], i, at, true, x, cy - MarkPx / 2, MarkPx, size));
        x -= MarkGap;
      }
      return x + MarkGap;
    };
    if (at.short) labelEnd = marks(0, given.length, chipsEnd, at.cy) - ChipGap;
    else if (at.ownLine !== undefined) {
      // a line of their own under the label, from the column's start; what does not fit in the column a mark (the
      // first always a chip)
      let x = inset;
      let i = 0;
      for (; i < given.length; i++) {
        const w = text.measure(config.labelFont, given[i].text);
        if (i > 0 && x + w > end) break;
        out.push(stateBox(given[i], i, at, false, x, at.ownLine, w, size));
        x += w + ChipGap;
      }
      if (i < given.length) marks(i, given.length, end, at.ownLine - size * 0.35);
    } else {
      // beside the label, the last state nearest the columns; the leading ones with no room left marks before them
      let x = chipsEnd;
      let cut = 0;
      for (let i = given.length - 1; i >= 0; i--) {
        const w = text.measure(config.labelFont, given[i].text);
        if (x - w < inset) {
          cut = i + 1;
          break;
        }
        x -= w;
        out.push(stateBox(given[i], i, at, false, x, at.cy + size * 0.35, w, size));
        x -= ChipGap;
      }
      labelEnd = cut ? marks(0, cut, x + ChipGap, at.cy) - ChipGap : x;
    }
    // in the given order
    out.sort((a, b) => a.at - b.at);
    states = out;
  }
  return {
    line: at.line,
    label: { x: inset, y: at.cy - size / 2 - 2, width: Math.max(0, labelEnd - inset), height: size + 4 },
    slots,
    ...(menu ? { menu } : {}),
    states,
    ...(actions ? { actions } : {}),
    short: at.short,
  };
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
    ...(p.address !== undefined ? { address: p.address } : {}),
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

/** The chip on a line of its own under the label: its baseline below the label's, when the track is two lines high. */
export function ownLineDy(config: RenderConfig, height: number): number {
  const size = fontSize(config.labelFont);
  return height >= 2 * (size + 6) ? size + 4 : 0;
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
): { readonly header: HeaderLine; readonly labelY: number } {
  const size = fontSize(config.labelFont);
  const short = well.height < config.header.minHeightPx;
  const dy = !short && parts?.states?.length ? ownLineDy(config, well.height) : 0;
  const centre = well.y + well.height / 2;
  const labelY = centre + size * 0.35 - dy / 2;
  const header = layoutLine(
    parts,
    {
      layer: track,
      section,
      cy: labelY - size * 0.35,
      line: { x: config.inset, y: well.y, width: config.labelWidth, height: well.height },
      short,
      ...(dy ? { ownLine: labelY + dy } : {}),
      figure,
    },
    config,
    text,
  );
  return { header, labelY };
}

/**
 * Whether a section draws its header line (`section.header`): `auto` when it has a state or a slot, `always`, `none`.
 */
export function drawsSectionHeader(parts: HeaderParts | undefined, config: RenderConfig): boolean {
  const mode = config.section.header;
  if (mode === 'none') return false;
  if (mode === 'always') return true;
  return !!(parts?.states?.length || parts?.slots?.length);
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
  const line = layoutLine(
    parts,
    {
      layer: section.address,
      section: id,
      cy: top + h / 2,
      line: { x: config.inset, y: top, width: config.labelWidth, height: h },
      short: false,
      figure,
    },
    config,
    text,
  );
  return {
    ...line,
    hit: { layer: section.address, section: id, system: '', target: 'section-header' },
  };
}

/** A track's layers hidden by the spec (`hidden: true`), by key. */
export function hiddenLayersOf(
  track: ResolvedTrack,
  key: (layer: ResolvedTrack['layers'][number]) => string,
): string[] {
  let out: string[] | undefined;
  for (const l of track.layers) if (l.params.hidden === true) (out ??= []).push(key(l));
  return out ?? (NoKeys as string[]);
}
const NoKeys: readonly string[] = [];

/** The spec's node at an address (a path of ids), undefined when there is none. */
function nodeAt(root: Node, address: string): Node | undefined {
  const p = parseAddress(address).pattern;
  if (!p) return undefined;
  let n: Node | undefined = root;
  for (const s of p.segments) {
    if (typeof s !== 'string') return undefined;
    n = n?.children?.find((c) => c.id === s);
  }
  return n;
}

const Hidden = new WeakMap<ResolvedSection, readonly HiddenTrack[]>();
const NoHidden: readonly HiddenTrack[] = [];

/** A section's hidden tracks in spec order (each once), read from the spec its view was resolved from; kept a section. */
export function hiddenTracksOf(model: StaticModel, section: ResolvedSection): readonly HiddenTrack[] {
  let out = Hidden.get(section);
  if (out) return out;
  let found: HiddenTrack[] | undefined;
  for (const t of nodeAt(model.spec, section.address)?.children ?? []) {
    if (t.kind !== 'track' || t.id === undefined || t.params?.hidden !== true) continue;
    const label = t.params?.label;
    (found ??= []).push({ address: `${section.address}/${t.id}`, label: typeof label === 'string' ? label : t.id });
  }
  out = found ?? NoHidden;
  Hidden.set(section, out);
  return out;
}
