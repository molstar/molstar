/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Framework-free helpers of the default UIs (seqspec:interaction#the-default-ui, seqspec:interaction#the-dom-element):
// React's parts and the DOM parts call the same code for a section's header (its label and switches), its hidden
// tracks, the slots a chrome leaves, a track's menu, the legend, the issue badge and list, a label as text, and the
// plain `Description`. The viewer has no markdown. Edge-rate: each runs when the document, the assembled view or a UI
// state changes, never per pointer move.

import {
  type DocNode,
  dataDomain,
  Layout,
  resolveStyle,
  SectionAxis,
  sectionHighlights,
  type ThemeRef,
} from '@molstar/seq-core/format';
import {
  type Numbering,
  Numbering as Numberings,
  OrderedSet,
  type PositionLabel,
  type Severity,
} from '@molstar/seq-core/model';
import { type ActionEntry, nodePath, trackDataRefs } from './actions.js';
import { type ColorSchemePicker, colorSchemePickers } from './color-schemes.js';
import { SeqCommands } from './commands.js';
import { type LegendMode, SeqViewerConfig } from './config.js';
import type { SeqIssue, SeqViewerContext } from './context.js';
import { severityRank } from './graph.js';
import { HitRegion } from './hit.js';
import { itemActionsOf } from './hover-card.js';
import type { HoverLabel } from './labels.js';
import { EntryKind } from './layout.js';
import { type LegendItem, type LegendRamp, type LegendSize, trackLegends } from './legend.js';
import type { SeqLociAny } from './lineage.js';
import { sameData } from './same.js';
import { legendNumber } from './themes.js';

/** An entry of a section's `positions`: a data ref, or its rows after a layer's filter (seqspec:format#sections). */
export type PositionEntry =
  | string
  | { readonly data: string; readonly filter: { readonly values: readonly unknown[] } };

export interface MenuItem {
  readonly key: string;
  readonly label: string;
  readonly checked?: boolean;
  readonly danger?: boolean;
  run(): void;
}
export interface SectionSwitch {
  /** The action it runs with `{ section, value }`. */
  readonly action: string;
  readonly name: string;
  readonly value: string | null;
  readonly options: readonly { readonly value: string | null; readonly label: string }[];
}
export interface SectionInfo {
  readonly label: string | null;
  readonly switches: readonly SectionSwitch[];
  /** The `positions` filter's label, or null: every position shown. */
  readonly filter: string | null;
}

const nodeOf = (ctx: SeqViewerContext, ref: string) => nodePath(ctx.document.shown.value, ref)?.at(-1) ?? null;
const labelOf = (n: DocNode | null, ref: string) => String(n?.params?.label || ref);
const choices = (...values: string[]) => values.map((value) => ({ value, label: value }));

const la: PositionLabel = { num: 0, suffix: -1 as never, name: -1 as never };
const lb: PositionLabel = { ...la };
const segs = (n: Numbering) =>
  n.method === 'segments'
    ? Array.from(n.indexStart, (s, k) => `${s},${n.labelStart[k]},${n.segmentLength[k]}`).join(';')
    : n.method === 'offset' && !n.skipZero
      ? `0,${n.first},${n.length}`
      : null;
const alike = new WeakMap<Numbering, WeakMap<Numbering, boolean>>();
/**
 * Whether two numberings label every position alike (`label` and `ordinal` on most structures, often `auth` too):
 * offsets and segments compared as segments, others position by position; once per pair.
 */
function sameLabels(a: Numbering, b: Numbering): boolean {
  if (a === b) return true;
  let known = alike.get(a);
  const hit = known?.get(b);
  if (hit !== undefined) return hit;
  let same = a.length === b.length;
  const [sa, sb] = [segs(a), segs(b)];
  if (same && sa !== null && sb !== null) same = sa === sb;
  else
    for (let i = 0; same && i < a.length; i++) {
      const [x, y] = [Numberings.labelOf(a, i as never, la), Numberings.labelOf(b, i as never, lb)];
      same = x === y && la.num === lb.num && la.name === lb.name;
      if (same && la.suffix !== lb.suffix)
        same = a.method === 'labels' && b.method === 'labels' && a.strings.get(la.suffix) === b.strings.get(lb.suffix);
    }
  if (!known) {
    known = new WeakMap();
    alike.set(a, known);
  }
  known.set(b, same);
  return same;
}

/** What `SectionHeader` shows: the shown document's section params, the set's numberings, presets, `ui.style`. */
export function sectionInfo(ctx: SeqViewerContext, ref: string): SectionInfo {
  const p = nodeOf(ctx, ref)?.params ?? {};
  const view = ctx.graph.assembled.value?.view;
  // The section's number in the view read here, not the layout's: that one follows its own updates, no input of the
  // callers'.
  const s = view ? sectionNumber(view, ref) : -1;
  const switches: SectionSwitch[] = [];
  const add = (action: string, name: string, value: unknown, options: SectionSwitch['options']) =>
    switches.push({ action, name, value: value == null ? null : String(value), options });
  if (view && s >= 0) {
    const { resolved: rv, set } = view;
    // Per unit, its numberings by name (its own first, then its molecule's).
    const byUnit: Map<string, Numbering>[] = [];
    let [ligand, copies, tracks] = [false, false, 0];
    for (let b = rv.sections.blockStart[s]; b < rv.sections.blockStart[s] + rv.sections.blockCount[s]; b++) {
      ligand ||= rv.blocks.ligand[b] === 1;
      const members = rv.blocks.members[b];
      const units = members ? Array.from(members) : rv.blocks.unit[b] >= 0 ? [rv.blocks.unit[b]] : [];
      copies ||= units.length > 1;
      for (const u of units) {
        copies ||= set.groups.members[set.units.group[u]].length > 1;
        const mol = set.molecules[set.units.molecule[u]];
        const own = new Map<string, Numbering>();
        for (const n of [...(set.units.numberings[u]?.items ?? []), ...(mol?.numberings.items ?? [])])
          if (!own.has(set.strings.get(n.name))) own.set(set.strings.get(n.name), n);
        byUnit.push(own);
      }
    }
    // A numbering that labels every unit's positions as one already offered would switch to nothing new: left out,
    // unless the section names it (seqspec:interaction#the-default-ui).
    const names: string[] = [];
    for (const name of new Set(byUnit.flatMap((m) => [...m.keys()]))) {
      const same = (o: string) =>
        byUnit.every((m) => {
          const [x, y] = [m.get(name), m.get(o)];
          return x && y ? sameLabels(x, y) : !x && !y;
        });
      if (name === p.numbering || !names.some(same)) names.push(name);
    }
    for (let t = 0; t < rv.tracks.count; t++) if (rv.tracks.section[t] === s) tracks++;
    if (names.length > 1)
      add('section.numbering', 'Numbering', p.numbering, [{ value: null, label: 'default' }, ...choices(...names)]);
    // A switch whose own value removes what offers it (the heme drawn as atoms is no ligand block; each copy its own
    // block) stays offered while the section sets it: both values stay reachable.
    if (copies || p.copies != null) add('section.copies', 'Copies', p.copies ?? 'one', choices('one', 'each'));
    if (ligand || p.ligand_form != null)
      add('section.ligandForm', 'Ligand', p.ligand_form ?? 'position', choices('position', 'atoms'));
    if (rv.layout === Layout.tracks && rv.sections.axis[s] === SectionAxis.blocks && tracks > 1)
      add('section.rowSizing', 'Rows', p.row_sizing ?? 'uniform', choices('uniform', 'natural'));
    const geometry = rv.layout === Layout.text ? 'wrapped' : 'linear';
    const presets = ctx.registry.list('presets').filter((x) => x.geometry === geometry && !x.name.includes(':'));
    if (presets.length > 1)
      add('section.style', 'Style', ctx.ui.style.value.get(ref) ?? resolveStyle(p.style).preset, [
        { value: null, label: 'auto' },
        ...choices(...presets.map((x) => x.name)),
      ]);
  }
  // An entry is a data ref, or `{ data, filter }`: named by the filter's values, else the data's label.
  const positions = p.positions == null ? null : ([p.positions].flat() as PositionEntry[]);
  const named = (e: PositionEntry) =>
    typeof e === 'string' ? labelOf(nodeOf(ctx, e), e) : e.filter.values.length ? e.filter.values.join(', ') : e.data;
  return {
    label: p.label == null ? null : String(p.label),
    switches,
    filter: positions?.map(named).join(', ') ?? null,
  };
}

function sectionNumber(view: ShownView, ref: string): number {
  const sections = view.resolved.sections;
  for (let i = 0; i < sections.count; i++) if (view.set.strings.get(sections.ref[i]) === ref) return i;
  return -1;
}

type ShownView = NonNullable<NonNullable<SeqViewerContext['graph']['assembled']['value']>['view']>;
/** What a section's header and switches show, computed once per change of their inputs (`sectionChrome`). */
export interface SectionChrome extends SectionInfo {
  /** The section's hidden tracks (`HiddenChips`). */
  readonly hidden: readonly { ref: string; label: string }[];
}
interface ChromeEntry {
  readonly shown: unknown;
  readonly assembled: unknown;
  readonly style: unknown;
  readonly registry: unknown;
  readonly value: SectionChrome;
}
const chromes = new WeakMap<object, Map<string, ChromeEntry>>();
/**
 * A section's chrome (`sectionInfo` and `hiddenTracks`) for both UIs: the header, its switches and its hidden chips
 * read this one value, so the numbering and presets work is done once per change of the document, the assembled view,
 * `ui.style` or the registry (the inputs it reads), not once per part and per render. A part whose inputs did not change
 * gets the same object back; one whose parts did not change (by `sameData`) gets those same objects.
 */
export function sectionChrome(ctx: SeqViewerContext, section: string): SectionChrome {
  const [shown, assembled, style, registry] = [
    ctx.document.shown.value,
    ctx.graph.assembled.value,
    ctx.ui.style.value,
    ctx.registry.version.value,
  ];
  const by = chromes.get(ctx) ?? new Map<string, ChromeEntry>();
  chromes.set(ctx, by);
  const known = by.get(section);
  if (
    known &&
    known.shown === shown &&
    known.assembled === assembled &&
    known.style === style &&
    known.registry === registry
  )
    return known.value;
  const next = { ...sectionInfo(ctx, section), hidden: hiddenTracks(ctx, section) };
  const prev = known?.value;
  let value: SectionChrome = next;
  if (prev) {
    const [switches, hidden] = [
      sameData(prev.switches, next.switches) ? prev.switches : next.switches,
      sameData(prev.hidden, next.hidden) ? prev.hidden : next.hidden,
    ];
    value =
      prev.label === next.label && prev.filter === next.filter && switches === prev.switches && hidden === prev.hidden
        ? prev
        : { label: next.label, filter: next.filter, switches, hidden };
  }
  by.set(section, { shown, assembled, style, registry, value });
  return value;
}
/** The states `sectionChrome` reads, for a part's subscriptions. */
export const sectionChromeStates = (ctx: SeqViewerContext) =>
  [ctx.document.shown, ctx.graph.assembled, ctx.ui.style, ctx.registry.version] as const;

/** The section's hidden tracks (`HiddenChips`). */
export function hiddenTracks(ctx: SeqViewerContext, section: string): { ref: string; label: string }[] {
  const out: { ref: string; label: string }[] = [];
  const walk = (n: DocNode) => {
    if (n.kind === 'track' && n.params?.hidden === true) out.push({ ref: n.ref ?? '', label: labelOf(n, n.ref ?? '') });
    for (const c of n.children ?? []) walk(c);
  };
  const node = nodeOf(ctx, section);
  if (node) walk(node);
  return out;
}

export const trackLabel = (ctx: SeqViewerContext, track: string) => labelOf(nodeOf(ctx, track), track);

/**
 * The default compositions (seqspec:interaction#slots-and-chrome): the parts the default layout places under each
 * chrome, where `full` places them. `RowControl` has no default part: it is placed only when a host gives one. The parts
 * other parts place (`HiddenChips` and `SectionSwitches` in a section's header) and the primitives are not here.
 */
export const ChromeCompositions = {
  full: [
    'ViewTabs',
    'Search',
    'Busy',
    'Export',
    'Issues',
    'SnapshotPlayer',
    'Description',
    'EmptyState',
    'ErrorState',
    'SectionHeader',
    'TrackActions',
    'RowControl',
    'Legend',
    'Tooltip',
    'HoverCard',
  ],
  clean: [
    'Search',
    'Busy',
    'SnapshotPlayer',
    'Description',
    'EmptyState',
    'ErrorState',
    'RowControl',
    'Tooltip',
    'HoverCard',
  ],
} as const;
/** The parts `full` places and `clean` does not: a host's slot (or `"default"`) places one under `clean` too. */
export const CleanParts = ['ViewTabs', 'SectionHeader', 'TrackActions', 'Legend', 'Export', 'Issues'] as const;
const PLACEABLE: readonly string[] = ChromeCompositions.full;
const PRIMITIVES = ['Hint', 'ScrollArea', 'Menu'];

const reported = new WeakMap<SeqViewerContext, Set<string>>();
/** React's part: a component function, or a memo or forward-ref object. */
export const isReactPart = (v: unknown) => typeof v === 'function' || (typeof v === 'object' && v !== null);
/** The DOM's part: a function it calls. */
export const isDomPart = (v: unknown) => typeof v === 'function';

/**
 * The default UIs' slots: each part's default, else the host's: a part, "default", or "none" (null: left out wherever it
 * is placed; a primitive keeps its default). A key set to `undefined` (a part not loaded yet) is not set. A value that is
 * not a part of the UI's (`isPart`: React's or the DOM's) falls back to the default and reports `seqv.ui.slot`, once per
 * slot. The chrome never changes a slot: a
 * part placed by a host's own composition is the same part under either chrome (`placedParts` is the composition).
 */
export function resolveSlots<T extends object>(
  ctx: SeqViewerContext,
  defaults: T,
  host: { readonly [K in keyof T]?: T[K] | 'none' | 'default' } | undefined,
  isPart: (v: unknown) => boolean,
): T {
  const out = { ...defaults } as Record<string, unknown>;
  for (const [k, v] of Object.entries(host ?? {})) {
    if (v === undefined) continue;
    if (v === 'default' || (v === 'none' && PRIMITIVES.includes(k))) out[k] = defaults[k as keyof T];
    else if (v === 'none') out[k] = null;
    else if (isPart(v)) out[k] = v;
    else {
      out[k] = defaults[k as keyof T];
      const once = reported.get(ctx) ?? new Set<string>();
      reported.set(ctx, once);
      if (!once.has(k)) {
        once.add(k);
        ctx.report('seqv.ui.slot', `the UI slot "${k}" is not a part: its default is used`);
      }
    }
  }
  return out as T;
}

/**
 * The parts the default layout places: the chrome's composition (this UI's option, else the context's `UI.Chrome`),
 * without the legend when it is hidden (this UI's `legend`, else `UI.Legend`) and the track menu under
 * `UI.TrackActions: none`, then each part a host's slot names (a part or `"default"`) placed where `full` places it, and
 * none of those set to `"none"`.
 */
export function placedParts(
  ctx: SeqViewerContext,
  chrome: 'full' | 'clean' | undefined,
  host: { readonly [k: string]: unknown } | undefined,
  legend?: LegendMode,
): ReadonlySet<string> {
  const out = new Set<string>(ChromeCompositions[chrome ?? ctx.config.get(SeqViewerConfig.UI.Chrome)]);
  if ((legend ?? ctx.config.get(SeqViewerConfig.UI.Legend)) === 'hidden') out.delete('Legend');
  if (ctx.config.get(SeqViewerConfig.UI.TrackActions) === 'none') out.delete('TrackActions');
  for (const [k, v] of Object.entries(host ?? {}))
    if (v === 'none') out.delete(k);
    else if (v !== undefined && PLACEABLE.includes(k)) out.add(k);
  return out;
}

/**
 * Where a track menu is in its "Colors" submenu (seqspec:interaction#color-schemes-at-runtime): `layer` its schemes are
 * listed for, null: the menu itself. The owner keeps it and answers `open` by showing the menu again with the other
 * list (both UIs).
 */
export interface ColorsDrill {
  readonly layer: string | null;
  open(layer: string | null): void;
}

/** The track's index in the assembled view, or -1. Menu-rate. */
function trackIndex(ctx: SeqViewerContext, track: string): number {
  const view = ctx.graph.assembled.value?.view;
  const id = view?.set.strings.find(track) ?? -1;
  if (!view || id < 0) return -1;
  return view.resolved.tracks.ref.indexOf(id);
}

/**
 * The items a track shows: its layers' data after their `filter` (a track of "Motif" over all of a feature set covers
 * those items, not the set), as item loci; null when its data names none (streamed, or not loaded).
 */
function shownLoci(ctx: SeqViewerContext, track: string): SeqLociAny | null {
  const view = ctx.graph.assembled.value?.view;
  const t = trackIndex(ctx, track);
  if (!view || t < 0) return null;
  const { resolved: rv, set } = view;
  const items = new Map<number, Set<number>>();
  for (let i = 0; i < rv.layers.count; i++) {
    if (rv.layers.track[i] !== t || set.strings.get(rv.layers.type[i]) === 'rows') continue;
    for (let d = rv.layers.dataStart[i]; d < rv.layers.dataStart[i] + rv.layers.dataCount[i]; d++) {
      const ann = set.annotations[rv.data[d]];
      if (!ann) continue;
      const kept = rv.kept[d];
      const into = items.get(ann.id) ?? new Set<number>();
      if (!kept) for (let k = 0; k < ann.itemCount; k++) into.add(k);
      else for (let r = 0; r < ann.count; r++) if (kept[r]) into.add(ann.item[r]);
      items.set(ann.id, into);
    }
  }
  const elements = [...items]
    .filter(([, s]) => s.size > 0)
    .map(([annotation, s]) => ({
      annotation: annotation as never,
      items: OrderedSet.ofSortedArray(Int32Array.from(s).sort()),
    }));
  return elements.length > 0 ? ({ kind: 'item-loci', set, elements } as SeqLociAny) : null;
}

/** Runs `layer.colorScheme` on each layer (a merged legend group's, or one), writing as `UI.ColorSchemeWrites` says. */
export function pickScheme(ctx: SeqViewerContext, layers: readonly string[], scheme: string | null) {
  const write = ctx.config.get(SeqViewerConfig.UI.ColorSchemeWrites);
  for (const layer of layers) ctx.actions.run('layer.colorScheme', { layer, scheme, write });
}

/** A picker's choices as menu items (a legend title's popover, a track's "Colors" submenu): the shown one checked. */
export function schemeItems(
  ctx: SeqViewerContext,
  layers: readonly string[],
  p: Pick<ColorSchemePicker, 'value' | 'options'>,
): MenuItem[] {
  return p.options.map((o) => ({
    key: `scheme:${o.value}`,
    label: o.label,
    checked: o.value === p.value,
    run: () => pickScheme(ctx, layers, o.value),
  }));
}

/**
 * A track header's menu (`TrackActions`): `UI.TrackActions` in order, each built-in when it applies. `colors` needs the
 * owner's `drill`: "Colors ▸" per layer that more than one scheme fits, and, while `drill.layer` is one of them, that
 * layer's schemes after a "‹ Colors" back item instead of the menu.
 */
export function trackMenu(
  ctx: SeqViewerContext,
  section: string,
  track: string,
  drill?: ColorsDrill,
  row = -1,
): MenuItem[] {
  const list = ctx.config.get(SeqViewerConfig.UI.TrackActions);
  const pickers =
    drill && ctx.config.get(SeqViewerConfig.UI.ColorSchemes) && list !== 'none' && list.includes('colors')
      ? colorSchemePickers(ctx, trackIndex(ctx, track))
      : [];
  const sub = drill?.layer == null ? undefined : pickers.find((p) => p.layer === drill.layer);
  if (drill && sub) {
    const back = { key: 'colors-back', label: '‹ Colors', run: () => drill.open(null) };
    return [back, ...schemeItems(ctx, [sub.layer], sub)];
  }
  const node = nodeOf(ctx, track);
  const refs = node ? trackDataRefs(node) : [];
  const out: MenuItem[] = [];
  // The keyboard cursor's item on this track: its item actions first, the hover card's keyboard path (Shift+F10). The
  // key names the item, so `sameMenu` sees a menu on another item as changed (its `run` holds the target).
  const found = itemActionsOf(ctx, 'keyboard');
  if (found?.target.section === section && found.target.track === track) {
    const { annotation, item, split } = found.target;
    for (const a of found.actions)
      out.push({
        key: `item:${a.name}:${annotation}:${item}:${split}`,
        label: a.title,
        run: () => ctx.actions.run(a.name, found.target, 'keyboard'),
      });
  }
  const action = (name: string, danger?: boolean, label?: string, key = name) =>
    ctx.actions.can(name, track) &&
    out.push({
      key,
      label: label ?? (ctx.registry.get('actions', name) as ActionEntry | undefined)?.title ?? name,
      danger,
      run: () => ctx.actions.run(name, track),
    });
  for (const a of list === 'none' ? [] : list) {
    if (typeof a === 'object') {
      const split = splitValue(ctx, track, row);
      if (a.when?.(ctx, section, track, split) ?? true)
        out.push({
          key: a.id,
          label: a.labelFor?.(ctx, section, track, split) ?? a.label,
          danger: a.danger,
          ...(a.checked ? { checked: a.checked(ctx, section, track, split) } : {}),
          run: () => a.run(ctx, section, track, split),
        });
    } else if (a === 'filter') action('track.filterView');
    else if (a === 'retry') action('stream.retry');
    else if (a === 'hide') action('track.hide');
    else if (a === 'remove') action('track.remove', true);
    else if (a === 'focus' && refs.length > 0)
      out.push({
        key: 'focus',
        label: 'Focus',
        run: () => {
          const loci = shownLoci(ctx, track);
          // A failed command is reported (`seqv.command.error`) by the commands; the catch only keeps it from being unhandled.
          void ctx.commands
            .dispatch(SeqCommands.Focus, { target: loci ? { loci } : { refs }, section, origin: 'api' })
            .catch(() => {});
        },
      });
    // A split track's rows (seqspec:format#layer-types): the toggle's keyboard path, through this menu.
    else if (a === 'rows' && node?.kind === 'track')
      action('group.toggle', undefined, ctx.ui.expanded.value.has(track) ? 'Collapse rows' : 'Expand rows', 'rows');
    else if (a === 'colors' && drill)
      for (const p of pickers)
        out.push({
          key: `colors:${p.layer}`,
          // Two layers with schemes: each named by its type ("Colors, letters ▸").
          label: pickers.length > 1 ? `Colors, ${layerKind(ctx, p.layer)} ▸` : 'Colors ▸',
          run: () => drill.open(p.layer),
        });
  }
  return out;
}

/**
 * What a track's menu part shows and decides "changed" by (`sameMenu`): the track's label, the items, the drilled layer
 * and `UI.TrackActions` itself. The items' `run` functions are built anew on every call, so they are not compared; the
 * only handlers that differ between two calls with the same items are a host's in that list, which is compared by
 * identity, and an item action's, whose key names its item.
 */
export interface TrackMenuModel {
  readonly label: string;
  readonly items: readonly MenuItem[];
  readonly layer: string | null;
  readonly list: unknown;
}
export function trackMenuModel(
  ctx: SeqViewerContext,
  section: string,
  track: string,
  drill: ColorsDrill,
  row = -1,
): TrackMenuModel {
  return {
    label: trackLabel(ctx, track),
    items: trackMenu(ctx, section, track, drill, row),
    layer: drill.layer,
    list: ctx.config.get(SeqViewerConfig.UI.TrackActions),
  };
}
export function sameMenu(a: TrackMenuModel, b: TrackMenuModel): boolean {
  if (a.label !== b.label || a.layer !== b.layer || !Object.is(a.list, b.list) || a.items.length !== b.items.length)
    return false;
  for (let i = 0; i < a.items.length; i++) {
    const [x, y] = [a.items[i], b.items[i]];
    if (x.key !== y.key || x.label !== y.label || x.checked !== y.checked || x.danger !== y.danger) return false;
  }
  return true;
}

/** A layer's type name ("letters", "ranges"), else its ref. */
function layerKind(ctx: SeqViewerContext, layer: string): string {
  const view = ctx.graph.assembled.value?.view;
  const i = view ? view.resolved.layers.ref.indexOf(view.set.strings.find(layer)) : -1;
  const type = view && i >= 0 ? view.resolved.layers.type[i] : -1;
  return view && type >= 0 ? view.set.strings.get(type) : layer;
}

/**
 * A legend group's color scheme picker: its title opens it. It switches the layer whose legend the group shows, on
 * every track the group stands for.
 */
export interface LegendColors {
  readonly layers: readonly string[];
  /** The scheme shown: the override, else the document's scheme name; null: the document's own color. */
  readonly value: string | null;
  /** `value`'s label, which the title's icon says; null: the document's. */
  readonly label: string | null;
  readonly options: ColorSchemePicker['options'];
}
/** A group's color range editor: the layers it sets, the ramp's ends, the data's own span, and whether it is the person's. */
export interface LegendRange {
  readonly layers: string[];
  readonly lo: number;
  readonly hi: number;
  /** The data's span (lo > hi: none); the ramp marks it, so a range wider than the data shows where the data sits. */
  readonly dataLo: number;
  readonly dataHi: number;
  /** Set at runtime (`ui.colorRanges`): the editor offers a reset. */
  readonly edited: boolean;
}
/** Sets one end of a group's range from the text typed (the other as shown); false, nothing set, for a bad number or order. */
export function editRange(ctx: SeqViewerContext, r: LegendRange, end: 'lo' | 'hi', text: string): boolean {
  const v = text.trim() === '' ? Number.NaN : Number(text);
  const range: [number, number] = end === 'lo' ? [v, r.hi] : [r.lo, v];
  let ok = true;
  for (const layer of r.layers) ok = ctx.actions.run('layer.colorRange', { layer, range }) && ok;
  return ok;
}
/** The data's span as shares of the ramp (0..1, clamped to its ends); null: the data has none to mark. */
export function dataSpan(r: LegendRange): readonly [number, number] | null {
  if (!(r.dataLo <= r.dataHi) || !(r.lo < r.hi)) return null;
  const at = (x: number) => Math.min(1, Math.max(0, (x - r.lo) / (r.hi - r.lo)));
  return [at(r.dataLo), at(r.dataHi)];
}
/** The layers' own ranges again (the document's, or the data's). */
export function resetRange(ctx: SeqViewerContext, r: LegendRange) {
  for (const layer of r.layers) ctx.actions.run('layer.colorRange', { layer, range: null });
}
/**
 * A legend group (seqspec:interaction#the-legend): the legend of one track's shown layers, every channel's, or of
 * several tracks whose legends are identical (the same items, ramps, sizes and picker), merged.
 */
export interface LegendGroup {
  /** The first track's `section/trackRef`: unique, stable while that track is shown, the key of `ui.legendOpen`. */
  readonly key: string;
  /** The first track's section (its layout gives the legend's grid: the label column, the data area). */
  readonly section: string;
  /** The first track. */
  readonly ref: string;
  /** Every track the group stands for, in order. */
  readonly refs: readonly string[];
  /** The title: the track's label (else its ref); merged, the tracks' common name and their count ("Domains (10)"). */
  readonly track: string;
  readonly items: readonly LegendItem[];
  readonly ramps: readonly LegendRamp[];
  readonly sizes: readonly LegendSize[];
  /** Null: no scheme switch (`UI.ColorSchemes` off, or fewer than two schemes fit the layer). */
  readonly colors: LegendColors | null;
  /** Null: the group has no single ramp a range could set (a scheme or palette over values). */
  readonly range: LegendRange | null;
}
/**
 * Whether a layer's color is a scheme or a normalized continuous or discrete palette: what a range can set (an
 * absolute palette's checkpoints are its range).
 */
const rangeable = (r: ThemeRef | null | undefined) => {
  if (r?.theme === 'scheme') return true;
  const palette = (r?.params as { palette?: { kind?: string; mode?: string } } | undefined)?.palette;
  return (
    r?.theme === 'field_palette' && (palette?.kind ?? 'categorical') !== 'categorical' && palette?.mode !== 'absolute'
  );
};
/** The lists of a view that is not there: one object, so a part comparing lists by identity finds no change. */
const NO_ENTRIES: readonly never[] = [];
/** The groups per assembled view: both UIs, and each of their renders, read one list until the view changes. */
const legends = new WeakMap<object, readonly LegendGroup[]>();

/**
 * A merged group's title: the labels' common leading words, else their common trailing words, else the layer type, as a
 * plural with the count: "Gene G1" … "Gene G5" → "Genes (5)", "G1.1 domain" … → "Domains (10)".
 */
export function groupTitle(labels: readonly string[], kind: string): string {
  if (labels.length === 1) return labels[0];
  const ws = labels.map((l) => l.split(/\s+/).filter(Boolean));
  const common = (at: (w: readonly string[], k: number) => string | undefined) => {
    const out: string[] = [];
    for (let k = 0; ; k++) {
      const w = at(ws[0], k);
      if (w === undefined || ws.some((x) => at(x, k) !== w)) return out;
      out.push(w);
    }
  };
  const lead = common((w, k) => w[k]);
  const trail = common((w, k) => w[w.length - 1 - k]).reverse();
  const name = (lead.length > 0 ? lead : trail.length > 0 ? trail : [kind]).join(' ');
  const cap = name.charAt(0).toUpperCase() + name.slice(1);
  return `${/[a-rt-z]$/.test(cap) ? `${cap}s` : cap} (${labels.length})`;
}

/**
 * The legend's groups from the shown layers' themes (renderer `legend`), in track order: a track with nothing to list
 * has none (its color schemes are in its menu), and tracks whose legends are identical share one, keyed by their
 * content and picker as JSON. Items equal in label, color and glyph are listed once. Computed once per assembled view
 * (a theme, a hidden track, a picked scheme or the skin's `dark` assemble a new one).
 */
export function legendGroups(ctx: SeqViewerContext): readonly LegendGroup[] {
  const view = ctx.graph.assembled.value?.view;
  if (!view) return NO_ENTRIES;
  const known = legends.get(view);
  if (known) return known;
  const { resolved: rv, set } = view;
  const out = trackLegends(ctx.registry, view);
  const picking = ctx.config.get(SeqViewerConfig.UI.ColorSchemes);
  type Merged = Omit<LegendGroup, 'track' | 'refs' | 'colors'> & {
    refs: string[];
    labels: string[];
    kind: string;
    colors: (LegendColors & { layers: string[] }) | null;
    range: LegendRange | null;
  };
  const merged = new Map<string, Merged>();
  const rangeOf = (i: number, ramp: LegendRamp): LegendRange => {
    const layer = set.strings.get(rv.layers.ref[i]);
    const anns = Array.from(rv.data.subarray(rv.layers.dataStart[i], rv.layers.dataStart[i] + rv.layers.dataCount[i]));
    const [dataLo, dataHi] = dataDomain(anns.flatMap((id) => (set.annotations[id] ? [set.annotations[id]] : [])));
    return { layers: [layer], lo: ramp.lo, hi: ramp.hi, dataLo, dataHi, edited: ctx.ui.colorRanges.value.has(layer) };
  };
  for (const [t, g] of out) {
    const ref = set.strings.get(rv.tracks.ref[t]);
    const own = new Set(g.layers.map((i) => set.strings.get(rv.layers.ref[i])));
    // The picker of the layer the legend shows (another layer's is in the track's menu).
    const p = picking ? colorSchemePickers(ctx, t).find((x) => own.has(x.layer)) : undefined;
    const content = { items: g.items, ramps: g.ramps, sizes: g.sizes };
    const id = JSON.stringify([content, p ? [p.value, p.options.map((o) => o.value)] : null]);
    const label = rv.tracks.label[t] >= 0 ? set.strings.get(rv.tracks.label[t]) : '';
    // A single ramp over a scheme's or palette's color: its layer's range can be set (the data's span marked).
    const ranged = g.ramps.length === 1 ? g.layers.find((i) => rangeable(rv.layers.themes[i][0])) : undefined;
    const rangeLayer = ranged === undefined ? '' : set.strings.get(rv.layers.ref[ranged]);
    const m = merged.get(id);
    if (m) {
      m.refs.push(ref);
      m.labels.push(label || ref);
      if (p && m.colors && !m.colors.layers.includes(p.layer)) m.colors.layers.push(p.layer);
      if (m.range && ranged !== undefined && !m.range.layers.includes(rangeLayer)) m.range.layers.push(rangeLayer);
      continue;
    }
    const section = ctx.layout.sectionRef(rv.tracks.section[t]) ?? '';
    const entry = p?.value ? ctx.registry.svs.get('colorSchemes', p.value) : undefined;
    merged.set(id, {
      key: `${section}/${ref}`,
      section,
      ref,
      refs: [ref],
      labels: [label || ref],
      kind: set.strings.get(rv.layers.type[g.layers[0]]),
      ...content,
      colors: p
        ? { layers: [p.layer], value: p.value, label: p.value ? (entry?.label ?? p.value) : null, options: p.options }
        : null,
      range: ranged === undefined ? null : rangeOf(ranged, g.ramps[0]),
    });
  }
  const groups = [...merged.values()].map(
    ({ labels, kind, ...m }): LegendGroup => ({ ...m, track: groupTitle(labels, kind) }),
  );
  legends.set(view, groups);
  return groups;
}

/** One entry of a section's `highlights` as the legend lists it: a location set kept, with its removal. */
export interface SectionHighlightEntry {
  /** `section/n`: unique, stable while the entry is there. */
  readonly key: string;
  readonly section: string;
  /** The entry's data ref (`highlight.remove`'s `data`). */
  readonly data: string;
  /** The set's label (else its ref). */
  readonly label: string;
  readonly style: 'band' | 'highlight' | 'accent';
}
const highlightLists = new WeakMap<object, readonly SectionHighlightEntry[]>();
/** The shown view's sections' `highlights` entries, in section then entry order; once per assembled view. */
export function sectionHighlightEntries(ctx: SeqViewerContext): readonly SectionHighlightEntry[] {
  const view = ctx.graph.assembled.value?.view;
  if (!view) return NO_ENTRIES;
  const known = highlightLists.get(view);
  if (known) return known;
  const { resolved: rv, set } = view;
  const out: SectionHighlightEntry[] = [];
  for (let s = 0; s < rv.sections.count; s++) {
    const section = set.strings.get(rv.sections.ref[s]);
    for (const [n, h] of sectionHighlights(rv, set, s).entries())
      out.push({
        key: `${section}/${n}`,
        section,
        data: h.data,
        label: labelOf(nodeOf(ctx, h.data), h.data),
        style: h.style,
      });
  }
  highlightLists.set(view, out);
  return out;
}

/**
 * Whether the legend is open: as the person's toggle or `view.legend` left it (`ui.legendExpanded`), else as this UI's
 * `legend` option says, else `UI.Legend`.
 */
export const legendExpanded = (ctx: SeqViewerContext, mode?: LegendMode) =>
  ctx.ui.legendExpanded.value ?? (mode ?? ctx.config.get(SeqViewerConfig.UI.Legend)) === 'expanded';

/** A group as shown: at most `UI.LegendItems` items and the rest as "+n more", unless it is open (`ui.legendOpen`). */
export function legendShown(ctx: SeqViewerContext, g: LegendGroup): { items: readonly LegendItem[]; more: number } {
  const cap = ctx.config.get(SeqViewerConfig.UI.LegendItems);
  const open = ctx.ui.legendOpen.value.has(g.key);
  return open || g.items.length <= cap
    ? { items: g.items, more: 0 }
    : { items: g.items.slice(0, cap), more: g.items.length - cap };
}

/**
 * The legend's groups as shown: at most `UI.LegendGroups`, the rest as "+n more", unless the legend is open
 * (`ui.legendOpen` holds `key`, the first group's `section/*`); `toggle`: a "+n more" or "fewer" control is shown.
 */
export function legendShownGroups(
  ctx: SeqViewerContext,
  groups: readonly LegendGroup[],
): { groups: readonly LegendGroup[]; more: number; key: string; toggle: boolean } {
  const cap = ctx.config.get(SeqViewerConfig.UI.LegendGroups);
  const key = groups.length > 0 ? `${groups[0].section}/*` : '';
  const open = ctx.ui.legendOpen.value.has(key);
  if (groups.length <= cap) return { groups, more: 0, key, toggle: false };
  return open
    ? { groups, more: 0, key, toggle: true }
    : { groups: groups.slice(0, cap), more: groups.length - cap, key, toggle: true };
}

/**
 * A group's title as both UIs show it: its name, whole; with a picker, a palette icon after it whose label (and the
 * hint's) is the scheme shown, "Colors: Zappo" ("Colors: document": the document's own color).
 */
export const legendTitle = (g: LegendGroup) => ({
  name: g.track,
  colors: g.colors ? `Colors: ${g.colors.label ?? 'document'}` : null,
});

/** A ramp's end as text (both UIs' legends), as the themes' legend samples read: `legendNumber`. */
export const rampEnd = legendNumber;

/** Issues of one ref, as both UIs list them (null: the document's own). */
export interface IssueGroup {
  readonly ref: string | null;
  readonly issues: SeqIssue[];
}
/** Issues grouped by ref, in the issue manager's order (null: the document's own). */
export function issuesByRef(issues: readonly SeqIssue[]): IssueGroup[] {
  const out: IssueGroup[] = [];
  for (const i of issues) {
    const last = out.at(-1);
    if (last && last.ref === i.ref) last.issues.push(i);
    else out.push({ ref: i.ref, issues: [i] });
  }
  return out;
}

/**
 * The issue badge as both UIs show it (seqspec:interaction#the-default-ui): errors and warnings counted as issues, the
 * badge marked by the most severe; notes alone a quiet "N notes", never issues. Null: no issues.
 */
export function issueBadge(issues: readonly SeqIssue[]): { readonly severity: Severity; readonly text: string } | null {
  let errors = 0;
  let warnings = 0;
  for (const i of issues) {
    if (i.severity === 'error') errors++;
    else if (i.severity === 'warning') warnings++;
  }
  const n = errors + warnings;
  if (n > 0) return { severity: errors > 0 ? 'error' : 'warning', text: `${n} ${n === 1 ? 'issue' : 'issues'}` };
  const notes = issues.length;
  return notes === 0 ? null : { severity: 'info', text: `${notes} ${notes === 1 ? 'note' : 'notes'}` };
}

/**
 * The issue list as both UIs show it: the `shown` most severe (errors, then warnings, then notes, the manager's order
 * within one), the errors and warnings grouped by ref, then the notes grouped by ref under their own heading; `more`:
 * how many were left out.
 */
export function issueList(
  issues: readonly SeqIssue[],
  shown = 40,
): { readonly problems: IssueGroup[]; readonly notes: IssueGroup[]; readonly more: number } {
  const n = Math.min(shown, issues.length);
  const order = issues.map((_, k) => k).sort((a, b) => severityRank(issues[a]) - severityRank(issues[b]) || a - b);
  const kept = new Uint8Array(issues.length);
  for (let k = 0; k < n; k++) kept[order[k]] = 1;
  const problems: SeqIssue[] = [];
  const notes: SeqIssue[] = [];
  for (let k = 0; k < issues.length; k++)
    if (kept[k]) (issues[k].severity === 'info' ? notes : problems).push(issues[k]);
  return { problems: issuesByRef(problems), notes: issuesByRef(notes), more: issues.length - n };
}

/** A hover label as one line (the live region's announcement): its title, lines, and a header's legend items. */
export const labelText = (l: HoverLabel | null) =>
  l
    ? [
        l.title,
        ...l.lines.map((x) => `${x.key} ${x.value}`),
        ...(l.legend?.items.map((i) => i.label) ?? []),
        l.legend?.more ? `+${l.legend.more} more` : '',
      ]
        .filter(Boolean)
        .join(', ')
    : '';

/** The current view's sections that are not hidden, in order (`ViewSurface`). */
export function shownSections(ctx: SeqViewerContext): string[] {
  const view = ctx.graph.assembled.value?.view;
  const s = view?.resolved.sections;
  const out: string[] = [];
  for (let i = 0; view && s && i < s.count; i++) if (!s.hidden[i]) out.push(view.set.strings.get(s.ref[i]));
  return out;
}

/** A surface's accessible name: the section's label, else its view's. */
export function sectionLabel(ctx: SeqViewerContext, section: string): string {
  const path = nodePath(ctx.document.shown.value, section);
  const own = path?.at(-1)?.params?.label;
  return String(own || path?.findLast((n) => n.kind === 'view')?.params?.label || section);
}

/** The current snapshot's description, else the document's. */
export function descriptionOf(ctx: SeqViewerContext): { text: string; format: 'markdown' | 'plaintext' } | null {
  const shown = ctx.document.shown.value;
  const snap = ctx.document.base.value?.snapshots?.find((s) => s.ref === ctx.snapshots.current.value);
  const from = snap?.description ? snap : shown?.metadata?.description ? shown.metadata : null;
  return from?.description ? { text: from.description, format: from.description_format ?? 'markdown' } : null;
}

/**
 * What the default `Description` shows: a `plaintext` description as paragraphs split on blank lines, text only; null
 * for a `markdown` one or none (a host's part renders markdown, with its own library).
 */
export function plainDescription(ctx: SeqViewerContext): string[] | null {
  const d = descriptionOf(ctx);
  if (d?.format !== 'plaintext') return null;
  return d.text
    .replace(/\r\n?/g, '\n')
    .split(/\n\s*\n/)
    .map((p) =>
      p
        .split('\n')
        .filter((l) => l.trim() !== '')
        .join(' '),
    )
    .filter((p) => p !== '');
}

/** A root error (`ErrorState`): nothing assembled while nodes failed; their messages, at most three. */
export function rootErrors(ctx: SeqViewerContext): string[] {
  if (!ctx.document.shown.value || ctx.graph.assembled.value?.view) return [];
  const out: string[] = [];
  for (const [ref, s] of ctx.graph.statuses.value)
    if (s.state === 'error') out.push(`${ref}: ${s.issues[0]?.message ?? 'failed'}`);
  return out.slice(0, 3);
}

/** A split track's row value at display row `row` (the menu's header line), null: not a split row. */
function splitValue(ctx: SeqViewerContext, track: string, row: number): string | null {
  const view = ctx.graph.assembled.value?.view;
  const t = row >= 0 && view ? trackIndex(ctx, track) : -1;
  const k = t >= 0 && view ? view.resolved.tracks.split[t] : -1;
  const labels =
    k >= 0
      ? (view?.prepared[k]?.state as { split?: { labels: readonly string[] } } | undefined)?.split?.labels
      : undefined;
  return labels && row < labels.length ? labels[row] : null;
}

/** The track whose header the pointer is over in this section, else the keyboard cursor's track there (`TrackActions`). */
export function headerTrack(ctx: SeqViewerContext, section: string): string | null {
  const h = ctx.hover.live.current;
  const view = ctx.graph.assembled.value?.view;
  if (view && h.region === HitRegion.header && h.track >= 0 && h.section === ctx.layout.sectionIndex(section))
    return view.set.strings.get(view.resolved.tracks.ref[h.track]);
  const c = ctx.ui.keyCursor.value;
  return c?.section === section ? c.track : null;
}

/** A track's header line on screen (`TrackActions` sits on it): a row of a rows track, else the track's entry. Reused. */
export interface HeaderLine {
  /** The layout entry; -1: none (not placed yet, or out of view). */
  entry: number;
  /** The row within a lazy run's entry (one entry for all its rows), else 0. */
  row: number;
  /** Surface px. */
  top: number;
  height: number;
}
export const headerLine = (): HeaderLine => ({ entry: -1, row: 0, top: 0, height: 0 });
function setLine(out: HeaderLine, entry: number, row: number, top: number, height: number) {
  out.entry = entry;
  out.row = row;
  out.top = top + row * height;
  out.height = height;
}

/**
 * The track's header line on screen, into `out`: the one under the pointer; else `out`'s previous line while it is
 * still the track's and in view (the pointer moved onto the menu, or the view scrolled under it); else the first in
 * view. False (`out.entry` -1): none in view. A rows track has a line per row, so this runs on every hover change, not
 * only when the track changes. Pointer-rate: allocation-free.
 */
export function trackLine(ctx: SeqViewerContext, section: string, track: string, out: HeaderLine): boolean {
  const view = ctx.graph.assembled.value?.view;
  const l = ctx.layout.get(section)?.current;
  const vp = ctx.viewports.get(section)?.live.current;
  const prev = out.entry;
  const prevRow = out.row;
  out.entry = -1;
  if (!view || !l || !vp) return false;
  const tracks = view.resolved.tracks;
  let t = tracks.count - 1;
  while (t >= 0 && view.set.strings.get(tracks.ref[t]) !== track) t--;
  if (t < 0) return false;
  const p = ctx.pointer.current;
  const over = p.inside && p.section === ctx.layout.sectionIndex(section);
  for (let e = 0; e < l.entryCount; e++) {
    if (l.entryKind[e] !== EntryKind.track || l.entryTrack[e] !== t) continue;
    const pinned = e < l.pinnedCount;
    const top = pinned ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - vp.scrollY;
    const run = l.entryRun[e];
    const h = run > 0 ? l.entryHeight[e] / run : l.entryHeight[e];
    if (!(h > 0)) continue;
    const min = pinned ? 0 : l.scrolledTop;
    // A lazy run's rows in view, else the entry alone (row 0).
    const first = run > 0 ? Math.max(0, Math.floor((min - top) / h)) : 0;
    const last = run > 0 ? Math.min(run - 1, Math.ceil((vp.height - top) / h) - 1) : 0;
    if (first > last || top + first * h + h <= min || top + first * h >= vp.height) continue;
    if (over && p.y >= top + first * h && p.y < top + (last + 1) * h) {
      setLine(out, e, Math.floor((p.y - top) / h), top, h);
      return true;
    }
    if (e === prev && prevRow >= first && prevRow <= last) {
      setLine(out, e, prevRow, top, h);
      if (!over) return true;
    } else if (out.entry < 0) setLine(out, e, first, top, h);
  }
  return out.entry >= 0;
}
