/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/seq-molstar-poc/viewer`, the viewer half (seqspec:libraries#the-mol-adapter): a track's colors as Seq* runs,
// read through the color themes the viewer prepared (so any theme, a scheme's dark variant included, reaches 3D as
// drawn); a structure shown through MVS composed from the viewer's state at that moment (the structure an item names,
// the colors of the track chosen for it, the highlights, the selection); host actions for the hover card and the track
// menu; and the link between a viewer and a 3D view over `ctx.systems`.

import { equal } from '@molstar/seq-core/format';
import { type Annotation, Annotations, Extents, KeyKind, type Loci } from '@molstar/seq-core/model';
import {
  type ActionEntry,
  type Coordinates,
  type CustomTrackAction,
  type ItemTarget,
  type LociSize,
  lociSize,
  type SeqViewerContext,
  ThemeLocation,
  type ViewModel,
} from '@molstar/seq-viewer';
import { OrderedSet } from '@molstar/core/data/int';
import { StructureElement } from '@molstar/model/model/structure';
import { type ComposeInput, complexOf, composeStructureMvs, type StructureRef } from './compose.js';
import type { StructureLink } from './link.js';
import type { MolstarDefaults } from './mvs.js';
import { complexUnits } from './seq-units.js';

const hex = (c: number) => `#${c.toString(16).padStart(6, '0')}`;
const COLOR_TYPES = new Set(['ranges', 'values', 'color', 'points', 'variants', 'emphasis', 'secondary', 'letters']);

/** An item as "Show in 3D" reads it: the viewer's `ItemTarget` (the hover card's) or its annotation and item. */
export type ItemRef = Pick<ItemTarget, 'annotation' | 'item'> & {
  readonly track?: string | null;
  readonly split?: number;
};

function annotationOf(view: ViewModel, ref: string): Annotation | null {
  for (const a of view.set.annotations) if (a && view.set.strings.get(a.key) === ref) return a;
  return null;
}
const structureField = (view: ViewModel, a: Annotation) => {
  const id = view.set.strings.find('structure');
  const f = a.fields.find((x) => x.name === id);
  return f?.kind === 'json' ? (f.data as readonly (StructureRef | null | undefined)[]) : null;
};
const sameSource = (a: StructureRef['source'], b: StructureRef['source']) => equal(a, b);
const sameLoci = (a: StructureElement.Loci | null, b: StructureElement.Loci | null) =>
  a === null || b === null ? a === b : StructureElement.Loci.areEqual(a, b);

/**
 * The structure an item stands for (its rows' JSON field `structure`) and its annotation's `custom.molstar`. Its
 * chains in front only when it was hovered on its own row (an expanded split row): else the whole structure.
 */
export function itemStructure(
  view: ViewModel,
  target: ItemRef,
): { structure: StructureRef; defaults: MolstarDefaults } | null {
  const a = annotationOf(view, target.annotation);
  const field = a && structureField(view, a);
  if (!a || !field) return null;
  for (let r = 0; r < a.count; r++) {
    if (a.item[r] !== target.item) continue;
    const s = field[r];
    if (!s?.source) continue;
    const structure = (target.split ?? -1) >= 0 ? s : { ...s, chains: undefined, authChains: undefined };
    return { structure, defaults: (a.custom.doc.molstar ?? {}) as MolstarDefaults };
  }
  return null;
}

/** A JSON field of an annotation by name; null: none. */
const jsonField = (view: ViewModel, a: Annotation, name: string) => {
  const id = view.set.strings.find(name);
  const f = a.fields.find((x) => x.name === id);
  return f?.kind === 'json' ? (f.data as readonly unknown[]) : null;
};
/** Rows that name structures (`structure`, or a list in `structures`) describe structures: they color none. */
const describesStructures = (view: ViewModel, a: Annotation) =>
  jsonField(view, a, 'structure') !== null || jsonField(view, a, 'structures') !== null;

/**
 * The structures a track's data names (its drawn rows' `structure` field, or each one of their `structures`), each
 * once, in row order, whole (no chains).
 */
export function trackStructures(view: ViewModel, track: string): StructureRef[] {
  const rv = view.resolved;
  const out: StructureRef[] = [];
  const add = (s: StructureRef | null | undefined) => {
    if (s?.source && !out.some((x) => sameSource(x.source, s.source)))
      out.push({ ...s, chains: undefined, authChains: undefined });
  };
  for (let i = 0; i < rv.layers.count; i++) {
    if (view.set.strings.get(rv.tracks.ref[rv.layers.track[i]]) !== track) continue;
    for (let d = rv.layers.dataStart[i]; d < rv.layers.dataStart[i] + rv.layers.dataCount[i]; d++) {
      const a = view.set.annotations[rv.data[d]];
      if (!a) continue;
      const one = jsonField(view, a, 'structure') as readonly (StructureRef | null)[] | null;
      const many = jsonField(view, a, 'structures') as readonly (readonly StructureRef[] | null)[] | null;
      for (let r = 0; r < a.count; r++) {
        if (rv.kept[d] && !rv.kept[d]?.[r]) continue;
        add(one?.[r]);
        for (const s of many?.[r] ?? []) add(s);
      }
    }
  }
  return out;
}

/**
 * A layer's colors as drawn: each data row's color from the layer's prepared color theme (its filter's rows only), or
 * a `letters` layer's per position, appended as runs in their system with their colors. Rows of an annotation that
 * names structures (coverage) color nothing; dense values are left out. Edge rate.
 */
/** Rows whose field `field` reads `value` (a split row's: its rows alone). */
export interface RowFilter {
  readonly field: string;
  readonly value: string;
}

export function layerColors(
  view: ViewModel,
  layer: number,
  runs: Extents,
  colors: string[],
  where: RowFilter | null = null,
): number {
  const theme = view.prepared[layer]?.themes[0];
  if (theme?.kind !== 'color') return 0;
  const rv = view.resolved;
  const set = view.set;
  const loc = ThemeLocation.create();
  const before = runs.count;
  if (set.strings.get(rv.layers.type[layer]) === 'letters') {
    const section = rv.tracks.section[rv.layers.track[layer]];
    for (let b = 0; b < rv.blocks.count; b++) {
      const u = rv.blocks.unit[b];
      if (rv.blocks.section[b] !== section || u < 0 || rv.blocks.alignment[b] >= 0) continue;
      const code = set.molecules[set.units.molecule[u]]?.code;
      let at = -1;
      for (let i = set.units.rangeStart[u]; code && i < set.units.rangeStop[u]; i++) {
        ThemeLocation.reset(loc);
        [loc.layer, loc.block, loc.unit, loc.index, loc.letter] = [layer, b, u, i, code[i]];
        const c = theme.value(loc);
        const last = runs.count - 1;
        if (c < 0) at = -1;
        else if (at >= 0 && colors[colors.length - 1] === hex(c) && runs.stop[last] === i) runs.stop[last] = i + 1;
        else {
          Extents.add(runs, { kind: KeyKind.unit === 1 ? 1 : 1, id: u, row: -1 }, i, i + 1);
          colors.push(hex(c));
          at = i;
        }
      }
    }
    return runs.count - before;
  }
  for (let d = rv.layers.dataStart[layer]; d < rv.layers.dataStart[layer] + rv.layers.dataCount[layer]; d++) {
    const a = set.annotations[rv.data[d]];
    if (!a || a.geometry === 'rows' || a.geometry === 'pairs' || describesStructures(view, a)) continue;
    // The layer's filter: only the rows it draws.
    const kept = rv.kept[d];
    for (let r = 0; r < a.count; r++) {
      if (kept && !kept[r]) continue;
      if (where && String(Annotations.rowValue(a, set, r, where.field) ?? '-') !== where.value) continue;
      ThemeLocation.reset(loc);
      loc.layer = layer;
      loc.annotation = a.id;
      loc.row = r;
      loc.item = a.item[r];
      loc.unit = a.target.unit[r];
      loc.index = a.start[r];
      loc.value = a.sparse ? a.sparse[r] : Number.NaN;
      loc.category = a.category ? a.category[r] : -1;
      const c = theme.value(loc);
      if (c < 0) continue;
      const unit = a.target.unit[r];
      const system = unit >= 0 ? { kind: 1, id: unit, row: -1 } : { kind: 0, id: a.target.molecule[r], row: -1 };
      if (system.id < 0 || a.start[r] < 0) continue;
      Extents.add(runs, system, a.start[r], Math.max(a.start[r] + 1, a.stop[r]));
      colors.push(hex(c));
    }
  }
  return runs.count - before;
}

/** Items this long or shorter (residues) are sites: a layer of only such items shows them as sticks in 3D. */
export const SITE_LENGTH = 3;
const SITE_TYPES = new Set(['points', 'variants']);

/** A layer of sites: points or variants by kind, or one whose drawn rows are all at most `SITE_LENGTH` long. */
function isSiteLayer(view: ViewModel, layer: number): boolean {
  const rv = view.resolved;
  const type = view.set.strings.get(rv.layers.type[layer]);
  if (SITE_TYPES.has(type)) return true;
  if (type !== 'ranges') return false;
  let any = false;
  for (let d = rv.layers.dataStart[layer]; d < rv.layers.dataStart[layer] + rv.layers.dataCount[layer]; d++) {
    const a = view.set.annotations[rv.data[d]];
    for (let r = 0; a && r < a.count; r++) {
      if (rv.kept[d] && !rv.kept[d]?.[r]) continue;
      if (a.stop[r] - a.start[r] > SITE_LENGTH) return false;
      any = true;
    }
  }
  return any;
}

/**
 * A track's colors (every color layer of it, in order) as runs and colors; `sticks` when every one of its color layers
 * is a site layer (by kind or by its items' length, never by name): its residues are shown as ball and stick too.
 */
export function trackColors(
  view: ViewModel,
  track: string,
  where: RowFilter | null = null,
): { runs: Extents; colors: string[]; sticks: boolean } {
  const rv = view.resolved;
  const runs = Extents.create(16);
  const colors: string[] = [];
  let sticks = true;
  let any = false;
  for (let i = 0; i < rv.layers.count; i++) {
    const t = rv.layers.track[i];
    if (t < 0 || view.set.strings.get(rv.tracks.ref[t]) !== track) continue;
    if (!COLOR_TYPES.has(view.set.strings.get(rv.layers.type[i]))) continue;
    if (layerColors(view, i, runs, colors, where) > 0) {
      any = true;
      sticks &&= isSiteLayer(view, i);
    }
  }
  return { runs, colors, sticks: any && sticks };
}

/**
 * The tracks whose data is the structure's own (rows on its complex's units, as a model's pLDDT): its default
 * coloring. By where the data sits, never by name.
 */
export function ownTracks(view: ViewModel, structure: StructureRef): string[] {
  const complex = complexOf(view.set, structure);
  if (!complex) return [];
  const units = new Set<number>();
  const t = { kind: 0, id: 0, row: -1 };
  for (const u of complexUnits(view.set, complex))
    if (view.set.byName(u.ref, t) && t.kind === KeyKind.unit) units.add(t.id);
  const rv = view.resolved;
  const out: string[] = [];
  for (let i = 0; i < rv.layers.count; i++) {
    const track = view.set.strings.get(rv.tracks.ref[rv.layers.track[i]]);
    for (let d = rv.layers.dataStart[i]; d < rv.layers.dataStart[i] + rv.layers.dataCount[i]; d++) {
      const a = view.set.annotations[rv.data[d]];
      // Rows that name structures (its chains' coverage) color nothing.
      if (!a || a.count === 0 || !units.has(a.target.unit[0]) || describesStructures(view, a)) continue;
      if (!out.includes(track)) out.push(track);
    }
  }
  return out;
}

export interface ComposeOptions {
  /** The tracks whose colors go to 3D, in order (a later one wins); unset: the structure's own (`ownTracks`). */
  readonly tracks?: readonly string[];
  /** Only their rows that pass this (a split row's value: that row alone). */
  readonly where?: RowFilter | null;
  /** The selection's color in 3D (default `#e8a33d`) and its label. */
  readonly selectionColor?: string;
  readonly selectionLabel?: string;
  /** Highlight channels' color (default `#f2c94c`); null: highlights stay out. */
  readonly highlightColor?: string | null;
  readonly services?: ComposeInput['services'];
  /** The 3D canvas's background (the host's surface, light or dark). */
  readonly background?: string;
}

function lociRuns(list: readonly (Loci | null)[]): Extents {
  const e = Extents.create(8);
  for (const l of list) if (l) Extents.addLoci(e, l);
  return e;
}

/**
 * A structure composed with the viewer's state now: the chosen tracks' colors as drawn (`tracks`, else the structure's
 * own), mapped through the document's residue maps onto its chains, residues they do not reach in the neutral base; the
 * highlight channels; the selection. Null: nothing shown yet, or the source resolves to no URL.
 */
export function composeStructureView(
  ctx: SeqViewerContext,
  structure: StructureRef,
  o: ComposeOptions & { readonly defaults?: MolstarDefaults } = {},
) {
  const view = ctx.graph.assembled.value?.view;
  if (!view) return null;
  const layers: { runs: Extents; colors: readonly string[]; sticks?: boolean }[] = [];
  for (const track of o.tracks ?? ownTracks(view, structure)) {
    const c = trackColors(view, track, o.where ?? null);
    if (c.runs.count) layers.push(c);
  }
  if (o.highlightColor !== null) {
    const lit: Loci[] = [];
    for (const e of ctx.highlight.channels.value.values()) if ('set' in e.target) lit.push(e.target as Loci);
    const runs = lociRuns(lit);
    if (runs.count) layers.push({ runs, colors: Array(runs.count).fill(o.highlightColor ?? '#f2c94c') });
  }
  const sel = ctx.selection.current.value;
  const selected = lociRuns([sel.seq, sel.aligned, sel.rows, sel.items, sel.boundaries] as (Loci | null)[]);
  return composeStructureMvs({
    set: view.set,
    structure,
    defaults: o.defaults,
    layers,
    selection: selected.count
      ? { runs: selected, label: o.selectionLabel ?? 'Selected in the sequence', color: o.selectionColor }
      : null,
    services: o.services,
    background: o.background,
  });
}

/**
 * The tracks that can color a structure: those whose colors reach its residues (a sequence's features through its
 * mappings; a structure's own data only on it), its own first. Edge rate (a menu, a chooser): composes per track.
 */
export function colorTracks(ctx: SeqViewerContext, structure: StructureRef): { track: string; label: string }[] {
  const view = ctx.graph.assembled.value?.view;
  if (!view) return [];
  const rv = view.resolved;
  const own = ownTracks(view, structure);
  const out: { track: string; label: string }[] = [];
  for (let t = 0; t < rv.tracks.count; t++) {
    const track = view.set.strings.get(rv.tracks.ref[t]);
    if (out.some((x) => x.track === track)) continue;
    const r = composeStructureView(ctx, structure, { tracks: [track], highlightColor: null });
    if (r?.tables.length) out.push({ track, label: view.set.strings.get(rv.tracks.label[t]) || track });
  }
  return [...out.filter((x) => own.includes(x.track)), ...out.filter((x) => !own.includes(x.track))];
}

/** "Show in 3D" for an item: its structure and its annotation's defaults, colored by `o.tracks` or its own. */
export function composeMvs(ctx: SeqViewerContext, target: ItemRef, o: ComposeOptions = {}) {
  const view = ctx.graph.assembled.value?.view;
  const s = view ? itemStructure(view, target) : null;
  return s ? composeStructureView(ctx, s.structure, { ...o, defaults: s.defaults }) : null;
}

/** What a host's 3D view shows: the structure, the track coloring it (null: its own default), the composed state. */
export interface Shown {
  readonly structure: StructureRef;
  readonly colorBy: string | null;
  /** Only the coloring track's rows that pass this (one split row of it); null: all its rows. */
  readonly where?: RowFilter | null;
  readonly defaults?: MolstarDefaults;
}

/** The filter for a split track's row by its value: the layer's split field and the value. */
export function splitWhere(view: ViewModel, track: string, value: string): RowFilter | null {
  const rv = view.resolved;
  for (let i = 0; i < rv.layers.count; i++) {
    if (view.set.strings.get(rv.tracks.ref[rv.layers.track[i]]) !== track) continue;
    const sp = rv.layers.params[i]?.split as string | { by: string } | null | undefined;
    const field = typeof sp === 'string' ? sp : sp?.by;
    if (field) return { field, value };
  }
  return null;
}

/** A split row's filter: the split field of the track's layer and the row's value; null: not on a split row. */
export function splitRow(view: ViewModel, track: string, split: number): RowFilter | null {
  if (split < 0) return null;
  const rv = view.resolved;
  for (let i = 0; i < rv.layers.count; i++) {
    if (view.set.strings.get(rv.tracks.ref[rv.layers.track[i]]) !== track) continue;
    const sp = rv.layers.params[i]?.split as string | { by: string } | null | undefined;
    const labels = (view.prepared[i]?.state as { split?: { labels: readonly string[] } } | undefined)?.split?.labels;
    const field = typeof sp === 'string' ? sp : sp?.by;
    if (field && labels && split < labels.length) return { field, value: labels[split] };
  }
  return null;
}

/**
 * The host's "Color 3D by …" item action (`molstar:color-3d`, for `UI.ItemActions`): on an item of a track whose
 * colors reach the shown structure, it colors the structure by that track, by the hovered split row alone when the
 * item was hovered on one ("Color 3D by Oligomerization"). Coverage items (naming structures) do not offer it.
 */
export function colorByItemAction(shown: () => Shown | null, show: (s: Shown) => void): ActionEntry<ItemRef> {
  const rowOf = (ctx: SeqViewerContext, target: ItemRef) => {
    const view = ctx.graph.assembled.value?.view;
    return view && target.track ? splitRow(view, target.track, target.split ?? -1) : null;
  };
  return {
    name: 'molstar:color-3d',
    title: 'Color 3D by this track',
    when: (ctx, target) => {
      const s = shown();
      return !!s && !!target.track && colorTracks(ctx, s.structure).some((x) => x.track === target.track);
    },
    titleFor: (ctx, target) => {
      const row = rowOf(ctx, target);
      return row ? `Color 3D by ${row.value}` : 'Color 3D by this track';
    },
    run: (ctx, target) => {
      const s = shown();
      if (s && target.track) show({ ...s, colorBy: target.track, where: rowOf(ctx, target) });
    },
  };
}

/** The hover card's button for an item against what 3D shows. */
export function showTitle(item: StructureRef, shown: Shown | null): string {
  if (!shown || !sameSource(item.source, shown.structure.source) || item.assembly !== shown.structure.assembly) {
    return 'Show in 3D';
  }
  const want = item.chains?.length ? [...item.chains].sort().join(', ') : null;
  const now = shown.structure.chains?.length ? [...shown.structure.chains].sort().join(', ') : null;
  if (want === now) return 'Shown in 3D';
  // Named by the author's chain IDs, what papers and the PDB's pages use.
  const name = [...(item.authChains ?? item.chains ?? [])].sort().join(', ');
  return want ? `Bring chain ${name} to front` : 'Show all chains';
}

/**
 * The host's "Show in 3D" item action (`molstar:show-in-3d`, for `UI.ItemActions`): offered on every item that names
 * a structure, whatever is shown; its title says what it does against what is shown (`shown`): "Show in 3D", "Shown
 * in 3D", "Bring chain D to front", "Show all chains". `show` gets what to show (the host composes and loads it).
 */
export function showInStructureAction(show: (s: Shown) => void, shown?: () => Shown | null): ActionEntry<ItemRef> {
  const structureOf = (ctx: SeqViewerContext, target: ItemRef) => {
    const view = ctx.graph.assembled.value?.view;
    return view ? itemStructure(view, target) : null;
  };
  return {
    name: 'molstar:show-in-3d',
    title: 'Show in 3D',
    when: (ctx, target) => structureOf(ctx, target) !== null,
    titleFor: (ctx, target) => {
      const s = structureOf(ctx, target);
      return s ? showTitle(s.structure, shown?.() ?? null) : 'Show in 3D';
    },
    run: (ctx, target) => {
      const s = structureOf(ctx, target);
      // Shown already: the coloring stays as chosen.
      const now = shown?.() ?? null;
      const same = !!now && !!s && sameSource(now.structure.source, s.structure.source);
      if (s)
        show({
          structure: s.structure,
          colorBy: same ? now.colorBy : null,
          where: same ? now.where : null,
          defaults: s.defaults,
        });
    },
  };
}

/**
 * Track menu items (`UI.TrackActions`): "Show <label> in 3D" for each structure a track's data names (a structure's
 * coverage row offers it; the summary, every one); "Apply to 3D" where the track's colors reach the shown
 * structure, checked on the track coloring it.
 */
export function structureTrackActions(o: {
  readonly structures: readonly { readonly label: string; readonly structure: StructureRef }[];
  readonly shown: () => Shown | null;
  readonly show: (s: Shown) => void;
}): CustomTrackAction[] {
  const names = (ctx: SeqViewerContext, track: string) => {
    const view = ctx.graph.assembled.value?.view;
    return view ? trackStructures(view, track) : [];
  };
  const coloring = (ctx: SeqViewerContext, s: Shown) => {
    const view = ctx.graph.assembled.value?.view;
    return s.colorBy ?? (view ? (ownTracks(view, s.structure)[0] ?? null) : null);
  };
  return [
    ...o.structures.map(
      (x): CustomTrackAction => ({
        id: `molstar:show:${x.label}`,
        label: `Show ${x.label} in 3D`,
        when: (ctx, _section, track) => names(ctx, track).some((s) => sameSource(s.source, x.structure.source)),
        run: () => o.show({ structure: x.structure, colorBy: null }),
      }),
    ),
    {
      id: 'molstar:color-by',
      label: 'Apply to 3D',
      // On a split row's line, that row alone; on the track's own line, the whole track.
      labelFor: (_ctx, _section, _track, split) => (split === null ? 'Apply to 3D' : `Apply ${split} to 3D`),
      when: (ctx, _section, track, split) => {
        const s = o.shown();
        if (!s) return false;
        const view = ctx.graph.assembled.value?.view;
        const where = split === null || !view ? null : splitWhere(view, track, split);
        return !!composeStructureView(ctx, s.structure, { tracks: [track], where, highlightColor: null })?.tables
          .length;
      },
      checked: (ctx, _section, track, split) => {
        const s = o.shown();
        return !!s && coloring(ctx, s) === track && (s.where?.value ?? null) === split;
      },
      run: (ctx, _section, track, split) => {
        const s = o.shown();
        const view = ctx.graph.assembled.value?.view;
        const where = split === null || !view ? null : splitWhere(view, track, split);
        if (s) o.show({ ...s, colorBy: track, where });
      },
    },
  ];
}

/** A 3D view as the link needs it: Mol*'s plugin implements it (`@molstar/seq-molstar-poc/plugin`), so could any host. */
export interface StructureView3D {
  highlight(loci: StructureElement.Loci | null): void;
  select(loci: StructureElement.Loci | null): void;
  focus(loci: StructureElement.Loci): void;
  /** Mol* hover and selection, as loci of the linked structure (null: off it). Returns the unsubscribe. */
  onHover(fn: (loci: StructureElement.Loci | null) => void): () => void;
  onSelect(fn: (loci: StructureElement.Loci | null) => void): () => void;
}

/** Runs in a link's chain systems as `Coordinates` by key, what the viewer's channels take. Edge rate. */
export function coordinatesOf(link: StructureLink, runs: Extents): Coordinates {
  const systems: string[] = [];
  const out: number[] = [];
  for (let i = 0; i < runs.count; i++) {
    const m = link.memberOf(runs, i);
    if (!m) continue;
    let k = systems.indexOf(m.key);
    if (k < 0) k = systems.push(m.key) - 1;
    out.push(k, runs.start[i], runs.stop[i]);
  }
  return { systems, runs: Int32Array.from(out) };
}

/**
 * How the link reacts to a Seq* set by its size (seqspec:libraries#the-mol-adapter), read before it resolves: within
 * `maxItems` items and `maxLength` positions, residue by residue; past either, `over`: the whole chains it touches, or
 * nothing.
 */
export interface LinkReaction {
  readonly maxItems: number;
  readonly maxLength: number;
  readonly over: 'chains' | 'skip';
}
/** The link's reactions (`linkViews`' options): its composition, as the viewer spec is the viewer's. */
export interface LinkViewsOptions {
  /** Seq*'s hover to 3D and back; false: neither. */
  readonly hover?: false | Partial<LinkReaction>;
  /** Selections both ways; false: neither. */
  readonly select?: false | Partial<LinkReaction>;
  /** A selection focuses the other view. */
  readonly focus?: boolean;
}
export const DefaultLinkViewsOptions: {
  readonly hover: LinkReaction;
  readonly select: LinkReaction;
  readonly focus: boolean;
} = {
  hover: { maxItems: 2_000, maxLength: 20_000, over: 'chains' },
  select: { maxItems: 20_000, maxLength: 200_000, over: 'chains' },
  focus: true,
};

/** The whole chains a set's runs touch, as one Mol* loci (a large set's summary). Edge rate. */
function wholeChains(link: StructureLink, runs: Extents): StructureElement.Loci {
  const members = link.toMembers(runs, Extents.create(16));
  const units = new Map<number, { unit: StructureElement.Loci['elements'][number]['unit']; lo: number; hi: number }>();
  for (let i = 0; i < members.count; i++) {
    const c = link.memberOf(members, i)?.chain;
    if (!c || c.count === 0) continue;
    const u = units.get(c.unit.id) ?? { unit: c.unit, lo: c.start[0], hi: c.end[c.count - 1] };
    [u.lo, u.hi] = [Math.min(u.lo, c.start[0]), Math.max(u.hi, c.end[c.count - 1])];
    units.set(c.unit.id, u);
  }
  const elements = [...units.values()].map((u) => ({
    unit: u.unit,
    indices: OrderedSet.ofBounds(u.lo, u.hi) as OrderedSet<StructureElement.UnitIndex>,
  }));
  return StructureElement.Loci(link.structure, elements);
}

/** A set's loci in 3D under a reaction: residue by residue within its limits, else its chains or nothing. */
function react(link: StructureLink, r: LinkReaction, size: LociSize | null, read: () => readonly (Loci | null)[]) {
  if (!size) return null;
  const within = size.items <= r.maxItems && size.length <= r.maxLength;
  if (!within && r.over === 'skip') return null;
  const runs = lociRuns(read());
  if (runs.count === 0) return null;
  return within ? link.toLoci(runs) : wholeChains(link, runs);
}

/**
 * Hover and selection both ways between a viewer and a 3D view over one link: Seq* hover lights the residues in 3D and
 * the 3D hover lights Seq*'s (owner `3d`); a selection in either selects in the other, and a 3D selection also focuses
 * the Seq* view on it, a Seq* one the 3D view. A Seq* set is sent by its size (`LinkReaction`, read before it resolves):
 * `DefaultLinkViewsOptions` unless `o` says otherwise. Replays and carries are no one's move
 * (seqspec:interaction#replays-and-carries): linking clears nothing in 3D, and a carry (origin `graph`) is sent only
 * where it changed the selection's 3D image last seen (a 3D selection's own included), never focused. Returns the unlink.
 */
export function linkViews(
  ctx: SeqViewerContext,
  view3d: StructureView3D,
  link: StructureLink,
  o: LinkViewsOptions = {},
): () => void {
  const stops: (() => void)[] = [];
  const scratch = Extents.create(16);
  if (o.hover !== false) {
    const r = { ...DefaultLinkViewsOptions.hover, ...o.hover };
    // What was last sent (undefined: not known, the 3D view hovered its own): the same again (a relayout's refresh,
    // the replay at link time) is not sent.
    let lit: StructureElement.Loci | null | undefined = null;
    const sub = ctx.hover.live.version.subscribe(() => {
      if (ctx.hover.live.current.owner === '3d') {
        lit = undefined;
        return;
      }
      const next = react(link, r, ctx.hover.size(), () => [ctx.hover.toLoci() as Loci | null]);
      if (lit !== undefined && sameLoci(next, lit)) return;
      lit = next;
      view3d.highlight(next);
    });
    stops.push(() => sub.unsubscribe());
    stops.push(
      view3d.onHover((loci) => {
        scratch.count = 0;
        ctx.hover.setLoci(loci ? coordinatesOf(link, link.fromLoci(loci, scratch)) : null, '3d');
      }),
    );
  }
  if (o.select !== false) {
    const r = { ...DefaultLinkViewsOptions.select, ...o.select };
    // What was last sent to 3D: Mol*'s event echoing it is not sent back (Seq*'s echo is told by its origin). `seen`:
    // the 3D image of the Seq* selection last seen, the one 3D made too: a carry is compared with it, so a 3D selection
    // Seq* holds only in part (an atom, a unit no system covers) is not replaced by its image on a data edit.
    let sent3d: StructureElement.Loci | null = null;
    let seen: StructureElement.Loci | null = null;
    let replay = true;
    const sub = ctx.selection.current.subscribe((sel) => {
      const first = replay;
      replay = false;
      const parts = [sel.seq, sel.aligned, sel.rows, sel.items, sel.boundaries] as (Loci | null)[];
      const size = { items: 0, length: 0 };
      for (const l of parts) {
        const z = lociSize(l);
        size.items += z?.items ?? 0;
        size.length += z?.length ?? 0;
      }
      const loci = react(link, r, size, () => parts);
      const next = loci?.elements.length ? loci : null;
      // The replay at link time sends a selection Seq* holds, never a clear; a carry only what it changed.
      const same = sameLoci(next, seen);
      seen = next;
      if (sel.origin === 'host:3d' || ((first || sel.origin === 'graph') && same)) return;
      sent3d = next;
      view3d.select(sent3d);
      if (sent3d && o.focus !== false && !first && sel.origin !== 'graph') view3d.focus(sent3d);
    });
    stops.push(() => sub.unsubscribe());
    stops.push(
      view3d.onSelect((loci) => {
        const same = loci === null ? sent3d === null : sent3d !== null && StructureElement.Loci.areEqual(loci, sent3d);
        if (same) return;
        sent3d = loci;
        scratch.count = 0;
        const coords = loci ? coordinatesOf(link, link.fromLoci(loci, scratch)) : null;
        const seq = coords?.runs.length ? coords : null;
        ctx.selection.set(seq, 'host:3d');
        if (seq && o.focus !== false) ctx.focus.fromLoci(seq, { origin: 'host:3d' });
      }),
    );
  }
  return () => {
    for (const s of stops) s();
  };
}
