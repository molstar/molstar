/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import type { ResolvedLayer, ResolvedSection, ResolvedTrack, ResolvedTrackGroup } from '@molstar/seqx-core/view';
import type { RenderConfig } from './config.js';
import { RenderIssues } from './issues.js';
import type { Box, Hit, LayerType } from './layer-type.js';
import { ChevronPx, CountGapPx, groupCount } from './groups-hooks.js';
import { lanesHeight } from './sizing.js';

// A section's track groups in the tracks layout (_spec/README.md "Render", Track groups): the bands a folded group
// draws in place of its tracks (its summary: each member layer in its type's summary form, or its authored summary
// track, or its header alone), the room a header line and a group's gaps take before a band, and the header lines.
// Applied once a layout, never a frame's work: a fold changes the bands, a warm frame reuses them.

/** A layer type's summary form (`LayerType.collapse`): one shared lane, a thin strip a member, or one row. */
export type CollapseForm = 'lane' | 'strip' | 'row';

/**
 * The built-in types' summary forms (_spec/README.md "Render", Layer types: `collapse`), read where a type does not say
 * its own: the ranges family and pairs one lane, values a strip, rows their consensus row; letters, logos, rulers,
 * matrices, grids, frames and bands none (left out, counted).
 */
const BuiltInCollapse: Readonly<Record<string, CollapseForm>> = {
  ranges: 'lane',
  secondary: 'lane',
  variants: 'lane',
  pairs: 'lane',
  values: 'strip',
  rows: 'row',
};

/** A type's summary form, its own or the built-ins'; undefined: left out of a summary. */
export function collapseOf(type: LayerType<any, any> | undefined): CollapseForm | undefined {
  if (!type) return undefined;
  return type.collapse ?? BuiltInCollapse[type.name];
}

/** A group's header as the layout places it: its line (open) or its folded row's label column, and its count. */
export interface GroupHeader {
  readonly group: ResolvedTrackGroup;
  readonly folded: boolean;
  /** Its chevron drawn (chrome: an export's figure draws none, its label at the column's start). */
  readonly chevron: boolean;
  /** The header's box in the label column: the line above an open group's first track, a folded row's column. */
  readonly box: Box;
  /** Its label as drawn (cut with an ellipsis before the count), measured, and the baseline of its line. */
  readonly label: string;
  readonly labelWidth: number;
  readonly baseline: number;
  /** The quiet count beside a folded group's label ("3 tracks"), and what its summary leaves out ("+2 not shown"). */
  readonly count?: string;
  readonly countWidth: number;
  /** Its hit (`group-header`, the group's address in `layer`): made once a layout, a hover allocates nothing. */
  readonly hit: Hit;
  /** Its bands in the section's shown list, `[first, last]`. */
  readonly first: number;
  readonly last: number;
}

/** What folding a section's groups gives the layout: its bands, the room before each, the groups' bands. */
export interface FoldedBands {
  readonly shown: ResolvedTrack[];
  readonly trackOf: number[];
  readonly sourceOf: ResolvedTrack[];
  /** Px before each band past the layout's `trackGap` (a header line, a group's wider gap; a strip's narrower one). */
  readonly leads?: number[];
  /** Each group drawn, with its bands (the header boxes are the layout's, `groupHeaderBoxes`). */
  readonly groups?: readonly {
    group: ResolvedTrackGroup;
    folded: boolean;
    first: number;
    last: number;
    count?: string;
    /** A folded row of strips alone: the px above its first strip, the header's line. */
    pad?: number;
  }[];
  /** A summary layer's original (its type's prep is the original's: a toggle prepares nothing again). */
  readonly originals?: ReadonlyMap<ResolvedLayer, ResolvedLayer>;
}

const NoIds: readonly string[] = [];

/** The groups folded now: the viewport's set (session state), else each group's `collapsed`. */
export function foldedIds(section: ResolvedSection, own: readonly string[] | undefined): readonly string[] {
  if (own) return own;
  const out: string[] = [];
  for (const g of section.trackGroups ?? []) if (g.collapsed) out.push(g.address.slice(g.address.lastIndexOf('/') + 1));
  return out.length ? out : NoIds;
}

/** A summary layer: the member's layer keyed `<layer>#summary` and its params in the summary form (no labels). */
const Summaries = new WeakMap<ResolvedLayer, ResolvedLayer>();
function summaryLayer(l: ResolvedLayer, form: CollapseForm): ResolvedLayer {
  let s = Summaries.get(l);
  if (!s) {
    const params = form === 'lane' ? { ...l.params, lanes: 'one' } : l.params;
    s = { ...l, params, ...(l.address !== undefined ? { address: `${l.address}#summary` } : {}) };
    Summaries.set(l, s);
  }
  return s;
}

/**
 * The synthetic tracks of a folded group, kept a group and a shape so a frame's key holds: by the resolved group (one
 * view's, so one viewer's: two viewers of one spec keep their own), then the band's address.
 */
const Bands = new WeakMap<ResolvedTrackGroup, Map<string, ResolvedTrack>>();
function band(
  group: ResolvedTrackGroup,
  address: string,
  label: string,
  height: number,
  layers: readonly ResolvedLayer[],
): ResolvedTrack {
  let own = Bands.get(group);
  if (!own) Bands.set(group, (own = new Map()));
  const k = own.get(address);
  if (
    k &&
    k.height === height &&
    k.layers.length === layers.length &&
    k.layers.every((l, i) => l === layers[i]) &&
    k.label === label
  )
    return k;
  const t: ResolvedTrack = { address, label, height, layers };
  own.set(address, t);
  return t;
}

/**
 * Folds a section's groups over its shown tracks (templates expanded): in `tracks`, a folded group's tracks give way
 * to its summary bands and an open group's authored summary track is left out; elsewhere (`text`, `rings`, `grid`)
 * every group open with no header and its summary track left out, said once when some are folded
 * (`seq.render.groups.layout`). Undefined when the section has no groups (nothing to do).
 */
export function foldGroups(
  section: ResolvedSection,
  shown: readonly ResolvedTrack[],
  trackOf: readonly number[],
  sourceOf: readonly ResolvedTrack[],
  own: readonly string[] | undefined,
  tracksLayout: boolean,
  config: RenderConfig,
  typeOf: (layer: ResolvedLayer) => LayerType<any, any> | undefined,
  report: (issue: Issue) => void,
): FoldedBands | undefined {
  if (!section.trackGroups) return undefined;
  // a pan lays the section out again with the same tracks and folds: the folding kept (a layout allocates none for it)
  const m = FoldMemos.get(section);
  if (
    tracksLayout &&
    m &&
    m.own === own &&
    m.config === config &&
    sameList(m.shown, shown) &&
    sameList(m.trackOf, trackOf) &&
    sameList(m.sourceOf, sourceOf)
  )
    return m.out;
  const out = fold(section, shown, trackOf, sourceOf, own, tracksLayout, config, typeOf, report);
  if (tracksLayout && out)
    FoldMemos.set(section, { own, config, shown: [...shown], trackOf: [...trackOf], sourceOf: [...sourceOf], out });
  return out;
}

/** A section's last folding, by the resolved section (a new view resolves new sections: none kept past it). */
const FoldMemos = new WeakMap<
  ResolvedSection,
  {
    readonly own: readonly string[] | undefined;
    readonly config: RenderConfig;
    readonly shown: readonly ResolvedTrack[];
    readonly trackOf: readonly number[];
    readonly sourceOf: readonly ResolvedTrack[];
    readonly out: FoldedBands;
  }
>();

function sameList<T>(a: readonly T[], b: readonly T[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

function fold(
  section: ResolvedSection,
  shown: readonly ResolvedTrack[],
  trackOf: readonly number[],
  sourceOf: readonly ResolvedTrack[],
  own: readonly string[] | undefined,
  tracksLayout: boolean,
  config: RenderConfig,
  typeOf: (layer: ResolvedLayer) => LayerType<any, any> | undefined,
  report: (issue: Issue) => void,
): FoldedBands | undefined {
  const groups = section.trackGroups;
  if (!groups) return undefined;
  const folded = foldedIds(section, own);
  const isFolded = (g: string) => folded.includes(g.slice(g.lastIndexOf('/') + 1));
  const out: FoldedBands & {
    shown: ResolvedTrack[];
    trackOf: number[];
    sourceOf: ResolvedTrack[];
    leads: number[];
    groups: {
      group: ResolvedTrackGroup;
      folded: boolean;
      first: number;
      last: number;
      count?: string;
      pad?: number;
    }[];
    originals: Map<ResolvedLayer, ResolvedLayer>;
  } = { shown: [], trackOf: [], sourceOf: [], leads: [], groups: [], originals: new Map() };
  if (!tracksLayout) {
    if (folded.length)
      report(
        issueOf(RenderIssues, 'seq.render.groups.layout', `Groups are drawn open in a ${section.layout} section`, {
          address: section.address,
          hint: RenderIssues.codes['seq.render.groups.layout'].hint,
        }),
      );
    for (let j = 0; j < shown.length; j++) {
      if (shown[j].summaryOf) continue;
      out.shown.push(shown[j]);
      out.trackOf.push(trackOf[j]);
      out.sourceOf.push(sourceOf[j]);
    }
    return { shown: out.shown, trackOf: out.trackOf, sourceOf: out.sourceOf };
  }
  const wider = config.groupGap - config.trackGap;
  let prevGroup: string | undefined;
  /** The next band's lead: past a group's end, its wider gap. */
  const push = (t: ResolvedTrack, ti: number, source: ResolvedTrack, lead: number) => {
    out.shown.push(t);
    out.trackOf.push(ti);
    out.sourceOf.push(source);
    out.leads.push(lead);
  };
  const done = new Set<string>();
  for (let j = 0; j < shown.length; j++) {
    const t = shown[j];
    const g = t.group;
    const after = out.shown.length > 0 && prevGroup !== g ? wider : 0;
    if (!g) {
      push(t, trackOf[j], sourceOf[j], after);
      prevGroup = undefined;
      continue;
    }
    const group = groups.find((x) => x.address === g);
    if (!group) continue;
    if (!isFolded(g)) {
      if (t.summaryOf) continue;
      const first = prevGroup !== g;
      push(t, trackOf[j], sourceOf[j], first ? after + config.header.groupPx : 0);
      if (first) out.groups.push({ group, folded: false, first: out.shown.length - 1, last: out.shown.length - 1 });
      else out.groups[out.groups.length - 1].last = out.shown.length - 1;
      prevGroup = g;
      continue;
    }
    prevGroup = g;
    if (done.has(g)) continue;
    done.add(g);
    const start = out.shown.length;
    let pad = 0;
    const members = shown.filter((m) => m.group === g && !m.summaryOf);
    let count = groupCount(members.length);
    if (group.summary) {
      // the authored summary: that track, drawn as any track, only while folded (EC-103)
      const k = shown.findIndex((m) => m.address === group.summary);
      if (k >= 0) push(shown[k], trackOf[k], sourceOf[k], after);
    } else if (group.none) {
      push(band(group, `${g}#header`, '', config.header.groupPx, []), -1, t, after);
    } else {
      // the automatic summary: one lane of every lane-form layer, a strip a values member, at most `maxStrips`
      const lane: ResolvedLayer[] = [];
      const strips: ResolvedLayer[] = [];
      let left = 0;
      for (const m of members)
        for (const l of m.layers) {
          if (!l.placed || l.params.hidden === true) continue;
          const form = collapseOf(typeOf(l));
          if (!form || form === 'row') continue;
          const s = summaryLayer(l, form);
          out.originals.set(s, l);
          if (form === 'lane') lane.push(s);
          else if (strips.length < config.group.summary.maxStrips) strips.push(s);
          else left++;
        }
      if (left) count = groupCount(members.length, left);
      // one lane when a member draws in one (or nothing would be drawn: the empty lane), the strips under it; strips
      // alone stand with no lane above them, the row padded at its top to the header's line
      const laned = lane.length > 0 || strips.length === 0;
      const stripsPx = strips.length * (config.group.summary.stripPx + 1) - 1;
      pad = laned ? 0 : Math.max(0, config.header.groupPx - stripsPx);
      if (laned) push(band(group, `${g}#summary`, '', lanesHeight(1, config), lane), -1, t, after);
      strips.forEach((s, k) =>
        push(
          band(group, `${g}#strip${k}`, '', config.group.summary.stripPx, [s]),
          -1,
          t,
          // a strip under the lane (or the one before): a 1 px gap
          !laned && k === 0 ? after + pad : 1 - config.trackGap,
        ),
      );
    }
    out.groups.push({ group, folded: true, first: start, last: out.shown.length - 1, count, ...(pad ? { pad } : {}) });
  }
  return out;
}

/**
 * The group headers of a laid out section: each group's box (an open group's line above its first well, a folded
 * group's rows in the label column), its label cut with an ellipsis before its count, measured once a layout.
 */
export function groupHeaders(
  groups: NonNullable<FoldedBands['groups']>,
  tracks: readonly { readonly well: Box }[],
  section: ResolvedSection,
  config: RenderConfig,
  measure: (text: string) => number,
  elide: (text: string, room: number) => string,
  scroll?: { readonly y: number; readonly end: number; readonly fixed: number },
  figure = false,
): GroupHeader[] {
  // a pan with the same bands where they were: the headers kept (a layout allocates none for them)
  const m = HeaderMemos.get(groups);
  if (
    m &&
    m.config === config &&
    m.figure === figure &&
    m.scrollY === (scroll?.y ?? -1) &&
    m.scrollEnd === (scroll?.end ?? -1) &&
    sameWells(m.wells, tracks)
  )
    return m.out;
  const out = headersOf(groups, tracks, section, config, measure, elide, scroll, figure);
  const wells = new Float64Array(2 * tracks.length);
  for (let i = 0; i < tracks.length; i++) {
    wells[2 * i] = tracks[i].well.y;
    wells[2 * i + 1] = tracks[i].well.height;
  }
  HeaderMemos.set(groups, { config, figure, scrollY: scroll?.y ?? -1, scrollEnd: scroll?.end ?? -1, wells, out });
  return out;
}

/** The last headers of a folding (its `groups`, kept by `foldGroups`), with what they were laid out from. */
const HeaderMemos = new WeakMap<
  object,
  {
    readonly config: RenderConfig;
    readonly figure: boolean;
    readonly scrollY: number;
    readonly scrollEnd: number;
    readonly wells: Float64Array;
    readonly out: GroupHeader[];
  }
>();

function sameWells(wells: Float64Array, tracks: readonly { readonly well: Box }[]): boolean {
  if (wells.length !== 2 * tracks.length) return false;
  for (let i = 0; i < tracks.length; i++)
    if (wells[2 * i] !== tracks[i].well.y || wells[2 * i + 1] !== tracks[i].well.height) return false;
  return true;
}

function headersOf(
  groups: NonNullable<FoldedBands['groups']>,
  tracks: readonly { readonly well: Box }[],
  section: ResolvedSection,
  config: RenderConfig,
  measure: (text: string) => number,
  elide: (text: string, room: number) => string,
  scroll: { readonly y: number; readonly end: number; readonly fixed: number } | undefined,
  figure: boolean,
): GroupHeader[] {
  const out: GroupHeader[] = [];
  const chevron = figure ? 0 : ChevronPx;
  const size = Number.parseFloat(config.labelFont) || 12;
  for (const g of groups) {
    const first = tracks[g.first]?.well;
    const last = tracks[g.last]?.well;
    if (!first || !last) continue;
    const x = config.inset;
    const width = config.labelWidth;
    const pad = g.pad ?? 0;
    const box: Box = g.folded
      ? { x, y: first.y - pad, width, height: Math.max(1, last.y + last.height - first.y + pad) }
      : { x, y: first.y - config.header.groupPx, width, height: config.header.groupPx };
    // a scrolled section's header wholly in its scrolled area, as a track's header line (or in the fixed head)
    if (scroll && g.first >= scroll.fixed && (box.y < scroll.y || box.y + box.height > scroll.end)) continue;
    const line = g.folded ? Math.min(box.height, Math.max(config.header.groupPx, first.height)) : box.height;
    const countWidth = g.count ? measure(g.count) : 0;
    // the label cut before the count, the count before the chevron: the chevron never dropped
    const fits = countWidth > 0 && chevron + countWidth + CountGapPx <= width;
    const label = elide(g.group.label, Math.max(0, width - chevron - (fits ? countWidth + CountGapPx : 0)));
    out.push({
      group: g.group,
      folded: g.folded,
      chevron: !figure,
      box,
      label,
      labelWidth: measure(label),
      baseline: box.y + line / 2 + size * 0.35,
      ...(fits && g.count ? { count: g.count } : {}),
      countWidth: fits ? countWidth : 0,
      first: g.first,
      last: g.last,
      hit: { layer: g.group.address, section: section.address, system: section.axis.system, target: 'group-header' },
    });
  }
  return out;
}
