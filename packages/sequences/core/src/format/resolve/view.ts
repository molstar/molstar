/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `composeView` and `resolveView`: sections and layers concatenated into the resolved view's tables, paths shared,
// `links` between blocks of different sections (memoized per set), `systems`, scroll groups checked (in `composeView`,
// so a viewer composing the stages itself gets the check). Our viewer and the Nightingale materializer both draw from
// this, so they cannot read the format differently.

import {
  type ComposedPath,
  type CoordinateMap,
  type CoordinateSystem,
  CoordinateSystemKind,
  Issue,
  Mappings,
  type NameId,
  type Numbering,
  OrderedSet,
  type SeqSet,
  type SortedArray,
  type UnitId,
} from "@molstar/seq-core/model";
import type { SeqViewSpecDoc, ThemeRef } from "../doc.js";
import { type DocIndex, indexDoc } from "../load/compile.js";
import { issueOf } from "../load/issues.js";
import type { SvsRegistries } from "../registry.js";
import type { StyleRenderProps } from "../tables/style.js";
import { type ResolvedSystems, systemsOf } from "./coordinates.js";
import { followsOf, type ResolvedFollows } from "./follows.js";
import { type ResolvedLayer, resolveLayer } from "./layer.js";
import { blockSystem, pathTo } from "./reach.js";
import { type BlockRow, childrenOf, depthOf, type ResolvedSection, resolveSection } from "./section.js";
import { forEachOf, instantiate, templateNodes } from "./template.js";

export type TableColumns<K extends string, T> = { readonly [k in K]: T };
/** A table: `count` rows, a column each. */
export type TableRows<T> = { readonly count: number } & T;
export type TableParams = TableColumns<"params", readonly Readonly<Record<string, unknown>>[]>;
export type TableLists<K extends string, T> = TableColumns<K, readonly T[]>;

/** One view's tables: what viewers draw from. */
export interface ResolvedView {
  readonly ref: NameId;
  readonly layout: number;
  readonly sections: TableRows<
    TableColumns<"ref" | "heightRows" | "style" | "scrollGroup" | "blockStart" | "blockCount", Int32Array> &
      TableColumns<"axis" | "arrangement" | "ligandForm" | "height" | "rowSizing" | "hidden", Uint8Array> &
      TableLists<"focus", { block: number; start: number; stop: number } | null> &
      TableLists<"render", StyleRenderProps | null> &
      TableParams
  >;
  readonly blocks: TableRows<
    TableColumns<
      "section" | "unit" | "alignment" | "base" | "start" | "stop" | "first" | "length" | "label" | "numberingPath",
      Int32Array
    > &
      TableColumns<"ligand", Uint8Array> &
      TableColumns<"dir", Int8Array> &
      TableLists<"layout", CoordinateMap> &
      TableLists<"members", SortedArray<UnitId> | null> &
      TableLists<"positions" | "varying", OrderedSet | null> &
      TableLists<"numbering", Numbering | null> &
      TableParams
  >;
  readonly groups: TableRows<
    TableColumns<"ref" | "section" | "parent" | "label", Int32Array> &
      TableColumns<"collapsed" | "hidden", Uint8Array> &
      TableParams
  >;
  readonly tracks: TableRows<
    TableColumns<"ref" | "section" | "group" | "label" | "laneStart" | "laneCount" | "split", Int32Array> &
      TableColumns<"height", Float64Array> &
      TableColumns<"hidden" | "expanded", Uint8Array> &
      TableParams
  >;
  readonly lanes: ResolvedLanes;
  readonly layers: TableRows<
    TableColumns<
      | "ref"
      | "track"
      | "type"
      | "declared"
      | "domainGroup"
      | "dataStart"
      | "dataCount"
      | "placementStart"
      | "placementCount",
      Int32Array
    > &
      TableColumns<"sizing", Uint8Array> &
      TableColumns<"domainMin" | "domainMax", Float64Array> &
      TableLists<"themes", readonly (ThemeRef | null)[]> &
      /** A template instance's node (seqspec:format#placement); null: a document layer. */
      TableLists<"node", string | null> &
      /** A matrix's states drawn, in order (seqspec:format#keyed-sets); null: all, in data order. */
      TableLists<"members", Int32Array | null> &
      TableParams
  >;
  readonly data: Int32Array;
  readonly kept: readonly (Uint8Array | null)[];
  readonly placements: TableRows<
    TableColumns<"layer" | "block" | "annotation" | "slot" | "member" | "lane" | "path", Int32Array>
  >;
  readonly paths: readonly ComposedPath[];
  readonly links: TableRows<TableColumns<"fromBlock" | "toBlock" | "path", Int32Array>>;
  /** Shared domains per (group, field), of the layers that share a `domain_group` (domains.ts). */
  readonly domains: TableRows<TableColumns<"group" | "field", Int32Array> & TableColumns<"min" | "max", Float64Array>>;
  /** Section pairs that follow each other (`follow`), with a route per block of either side (follows.ts). */
  readonly follows: ResolvedFollows;
  readonly systems: ResolvedSystems;
}

/**
 * A `rows` track's rows in order, by lane (read them through `Lanes`): lanes `0..row.length` stored, then each lazy
 * run's, never stored (seqspec:sources#data-by-reference: a lazy alignment's million rows are never listed).
 */
export interface ResolvedLanes {
  /** Lane ids, the stored and the lazy runs'. */
  readonly count: number;
  readonly row: Int32Array;
  readonly depth: Int32Array;
  /** Lane `start + k` is row `order[k]` (null: `k`, data order); depth by `parent` (the alignment's, shared). */
  readonly lazy: TableRows<
    TableColumns<"start" | "length", Int32Array> &
      TableLists<"order", Int32Array | null> &
      TableLists<"parent", Int32Array>
  >;
}

/** A lane's row and depth, stored or in a lazy run; -1: no such lane. Allocates nothing. */
export const Lanes = {
  rowOf(lanes: ResolvedLanes, lane: number): number {
    if (lane < lanes.row.length) return lane >= 0 ? lanes.row[lane] : -1;
    const z = lanes.lazy;
    for (let i = 0; i < z.count; i++) {
      const k = lane - z.start[i];
      const order = z.order[i];
      if (k >= 0 && k < z.length[i]) return order ? order[k] : k;
    }
    return -1;
  },
  depthOf(lanes: ResolvedLanes, lane: number): number {
    if (lane < lanes.row.length) return lane >= 0 ? lanes.depth[lane] : -1;
    const row = Lanes.rowOf(lanes, lane);
    const z = lanes.lazy;
    for (let i = 0; row >= 0 && i < z.count; i++) {
      if (lane >= z.start[i] && lane < z.start[i] + z.length[i]) return depthOf(z.parent[i], row);
    }
    return -1;
  },
};

type Row = Record<string, unknown>;
type Typed = typeof Int32Array | typeof Uint8Array | typeof Int8Array | typeof Float64Array;
/** Rows as struct-of-arrays: typed columns by `typed`, every other key of `keys` a plain list. */
function pack(rows: readonly Row[], typed: Record<string, Typed>, plain: string[] = []) {
  const out: Row = { count: rows.length };
  for (const [k, T] of Object.entries(typed)) out[k] = (T as typeof Float64Array).from(rows, (r) => r[k] as number);
  for (const k of plain) out[k] = rows.map((r) => r[k]);
  return out;
}
const cols = (T: Typed, keys: string) => Object.fromEntries(keys.split(" ").map((k) => [k, T]));

/** Sections and their layers as one view's tables, with `links` between blocks of different sections. */
export function composeView(
  view: string,
  sections: readonly ResolvedSection[],
  layers: readonly (readonly ResolvedLayer[])[],
  set: SeqSet,
  issues: Issue[] = [],
): ResolvedView {
  const paths: ComposedPath[] = [];
  const index = new Map<ComposedPath, number>();
  const pathOf = (p: ComposedPath | undefined) => {
    if (p && !index.has(p)) index.set(p, paths.push(p) - 1);
    return p ? (index.get(p) as number) : -1;
  };
  const [sec, blocks, groups, tracks, lrows, placements] = [[], [], [], [], [], []] as Row[][];
  // Lanes: every section's stored ones first, then the lazy runs (no lane listed per row of a lazy alignment).
  let stored = 0;
  for (const s of sections) stored += s.lanes.row.length;
  const lanes = { row: new Int32Array(stored), depth: new Int32Array(stored) };
  const lazy = {
    start: [] as number[],
    length: [] as number[],
    order: [] as (Int32Array | null)[],
    parent: [] as Int32Array[],
  };
  let [l0, z0] = [0, stored];
  const [data, kept] = [[] as number[], [] as (Uint8Array | null)[]];
  scrollGroups(sections, set, issues);
  sections.forEach((s, section) => {
    const [b0, g0, t0, n] = [blocks.length, groups.length, tracks.length, s.lanes.row.length];
    const laneOf = (x: number) => (x < 0 ? -1 : x < n ? x + l0 : x - n + z0);
    const focus = s.row.focus && { ...s.row.focus, block: s.row.focus.block + b0 };
    sec.push({ ...s.row, focus, blockStart: b0, blockCount: s.blocks.length });
    for (const b of s.blocks) blocks.push({ ...b, section, numberingPath: pathOf(s.paths[b.numberingPath]) });
    for (const g of s.groups) groups.push({ ...g, section, parent: g.parent < 0 ? -1 : g.parent + g0 });
    // A track's split layer by ref: a layer left out shifts no index.
    const composed = new Map<string, number>();
    lanes.row.set(s.lanes.row, l0);
    lanes.depth.set(s.lanes.depth, l0);
    for (const z of s.lanes.lazy) {
      lazy.start.push(laneOf(z.start));
      lazy.length.push(z.length);
      lazy.order.push(z.order);
      lazy.parent.push(z.parent);
    }
    // The track from the section, not the layer: a reused layer's may be stale after a track moved. By ref, once.
    const trackOf = new Map<string, number>();
    for (const x of s.layers) trackOf.set(x.ref, x.track);
    (layers[section] ?? []).forEach((l) => {
      const at = { dataStart: data.length, placementStart: placements.length };
      const track = trackOf.get(l.ref);
      composed.set(l.ref, lrows.length);
      lrows.push({
        ...l,
        ...at,
        ref: set.strings.intern(l.ref),
        track: (track ?? l.track) + t0,
        domainGroup: l.group === null ? -1 : set.strings.intern(l.group),
        dataCount: l.data.length,
        placementCount: l.placements.length,
      });
      for (const d of l.data) data.push(d);
      for (const k of l.kept) kept.push(k);
      for (const p of l.placements) {
        const lane = laneOf(p.lane);
        placements.push({ ...p, layer: lrows.length - 1, block: p.block + b0, lane, path: pathOf(l.paths[p.path]) });
      }
    });
    // A track whose every layer is pending on dynamic nodes not held is not shown until one arrives
    // (seqspec:sources#dynamic-nodes).
    const pending = new Map<number, boolean>();
    for (const l of layers[section] ?? []) {
      const tk = trackOf.get(l.ref) ?? l.track;
      pending.set(tk, (pending.get(tk) ?? true) && l.pending);
    }
    s.tracks.forEach((t, k) => {
      const laneStart = t.laneCount > 0 ? laneOf(t.laneStart) : t.laneStart + l0;
      const split = t.split < 0 ? -1 : (composed.get(s.layers[t.split]?.ref ?? "") ?? -1);
      const hidden = pending.get(k) === true ? 1 : t.hidden;
      tracks.push({ ...t, section, group: t.group < 0 ? -1 : t.group + g0, laneStart, split, hidden });
    });
    l0 += n;
    for (const z of s.lanes.lazy) z0 += z.length;
  });
  const links: Row[] = [];
  const system = blocks.map((b) => blockSystem(b as unknown as BlockRow));
  for (let i = 0; i < blocks.length; i++) {
    for (let j = i + 1; j < blocks.length; j++) {
      if (blocks[i].section === blocks[j].section) continue;
      const [a, b] = [system[i], system[j]];
      const p = linkOf(set, a, b);
      // Columns reach other systems through each of their rows (seqspec:model#projection): two alignments joined row by
      // row (a protein's and its cDNA's) are not one link, and that is no error.
      const columns = a.kind === CoordinateSystemKind.columns || b.kind === CoordinateSystemKind.columns;
      if (p === "ambiguous" && !columns)
        issues.push(Issue.error("svs.mapping.ambiguous", "", `blocks ${i} and ${j}: several paths`));
      else if (p !== undefined && p !== "ambiguous")
        links.push({ fromBlock: i, toBlock: j, path: pathOf(p ?? Mappings.compose(set, a, [])) });
    }
  }
  // The shared domains: per group the value's (the members' `domainMin/Max`, equal) and the other fields' they used.
  const domains: Row[] = [];
  const seen = new Set<string>();
  for (const l of lrows as unknown as (ResolvedLayer & { domainGroup: number })[]) {
    if (l.group === null) continue;
    const add = (field: string, min: number, max: number) => {
      if (!seen.has(`${l.group}\0${field}`) && seen.add(`${l.group}\0${field}`) && min <= max)
        domains.push({ group: l.domainGroup, field: set.strings.intern(field), min, max });
    };
    add("value", l.domainMin, l.domainMax);
    for (const f of l.shared) add(f.field, f.min, f.max);
  }
  const I = Int32Array;
  const U = Uint8Array;
  const F = Float64Array;
  const hiddenGroup = (g: number): boolean =>
    g >= 0 && (groups[g].hidden === 1 || hiddenGroup(groups[g].parent as number));
  const shownTrack = (t: number) => t >= 0 && tracks[t].hidden !== 1 && !hiddenGroup(tracks[t].group as number);
  const shownBlock = (b: number) => sec[blocks[b].section as number].hidden !== 1;
  const systems = systemsOf(set, {
    sections: sec.map((x, k) => ({ hidden: x.hidden === 1, blocks: system.filter((_, b) => blocks[b].section === k) })),
    paths,
    placements: placements.map((x) => ({
      path: x.path as number,
      shown: shownBlock(x.block as number) && shownTrack(lrows[x.layer as number].track as number),
    })),
    links: links.map((x) => ({
      path: x.path as number,
      shown: shownBlock(x.fromBlock as number) && shownBlock(x.toBlock as number),
    })),
  });
  const follows = followsOf(set, view, sections, pathOf, issues);
  return {
    ref: set.strings.intern(view),
    layout: sections[0]?.layout ?? 0,
    sections: pack(
      sec,
      {
        ...cols(I, "ref heightRows style scrollGroup blockStart blockCount"),
        ...cols(U, "axis arrangement ligandForm height rowSizing hidden"),
      },
      ["focus", "render", "params"],
    ),
    blocks: pack(
      blocks,
      {
        ...cols(I, "section unit alignment base start stop first length label numberingPath"),
        ligand: U,
        dir: Int8Array,
      },
      "layout members positions varying numbering params".split(" "),
    ),
    groups: pack(groups, { ...cols(I, "ref section parent label"), ...cols(U, "collapsed hidden") }, ["params"]),
    tracks: pack(
      tracks,
      {
        ...cols(I, "ref section group label laneStart laneCount split"),
        height: Float64Array,
        ...cols(U, "hidden expanded"),
      },
      ["params"],
    ),
    lanes: {
      ...lanes,
      count: z0,
      lazy: { count: lazy.start.length, ...lazy, start: I.from(lazy.start), length: I.from(lazy.length) },
    },
    layers: pack(
      lrows,
      {
        ...cols(I, "ref track type declared domainGroup dataStart dataCount placementStart placementCount"),
        sizing: U,
        ...cols(F, "domainMin domainMax"),
      },
      ["themes", "node", "members", "params"],
    ),
    data: I.from(data),
    kept,
    placements: pack(placements, cols(I, "layer block annotation slot member lane path")),
    paths,
    links: pack(links, cols(I, "fromBlock toBlock path")),
    domains: pack(domains, { ...cols(I, "group field"), ...cols(F, "min max") }),
    follows,
    systems,
  } as unknown as ResolvedView;
}

/** Each set's links by block coordinate system pair: a rerun of the view stage searches no path. */
const linkMemo = new WeakMap<SeqSet, Map<string, ReturnType<typeof pathTo>>>();
function linkOf(set: SeqSet, a: CoordinateSystem, b: CoordinateSystem): ReturnType<typeof pathTo> {
  let memo = linkMemo.get(set);
  if (!memo) {
    memo = new Map();
    linkMemo.set(set, memo);
  }
  const key = `${a.kind}:${a.id}:${a.row}>${b.kind}:${b.id}:${b.row}`;
  if (memo.has(key)) return memo.get(key);
  // Copies never reach each other (seqspec:format#blocks), nor do coordinate systems no mapping, unit or row joins.
  const u = CoordinateSystemKind.unit;
  const copies = a.kind === u && b.kind === u && a.id !== b.id && set.units.molecule[a.id] === set.units.molecule[b.id];
  const p = copies || componentOf(set, a) !== componentOf(set, b) ? undefined : pathTo(set, a, b);
  memo.set(key, p);
  return p;
}

/** The mapping graph's components (union-find over molecules, then alignments), once per set. */
const components = new WeakMap<SeqSet, Int32Array>();
function componentOf(set: SeqSet, s: CoordinateSystem): number {
  const n = set.molecules.length;
  const node = (x: CoordinateSystem) =>
    x.kind === CoordinateSystemKind.molecule
      ? x.id
      : x.kind === CoordinateSystemKind.unit
        ? set.units.molecule[x.id]
        : n + x.id;
  let up = components.get(set);
  if (!up) {
    const parent = Int32Array.from({ length: n + set.alignments.length }, (_, i) => i);
    const find = (i: number) => {
      let r = i;
      while (parent[r] !== r) r = parent[r];
      return r;
    };
    for (const m of set.mappings) if (m) parent[find(node(m.from))] = find(node(m.to));
    for (let i = 0; i < parent.length; i++) parent[i] = find(i);
    up = parent;
    components.set(set, up);
  }
  return up[node(s)] ?? -1;
}

/** Sections of one scroll group must have equal axes (`svs.section.scroll_group`). */
function scrollGroups(sections: readonly ResolvedSection[], set: SeqSet, issues: Issue[]) {
  const groups = new Map<number, ResolvedSection>();
  for (const s of sections) {
    const first = groups.get(s.row.scrollGroup);
    if (s.row.scrollGroup >= 0 && !first) groups.set(s.row.scrollGroup, s);
    else if (first && !equalAxes(first, s)) {
      const message = `scroll group "${set.strings.get(s.row.scrollGroup)}": its axis differs from "${first.ref}"'s`;
      issues.push(Issue.error("svs.section.scroll_group", s.path, message));
    }
  }
}

/**
 * Equal axes (seqspec:format#sections): equal units or members, pieces (range, direction) and shown positions, in
 * order.
 */
const equalAxes = (a: ResolvedSection, b: ResolvedSection) =>
  a.blocks.length === b.blocks.length &&
  a.blocks.every((x, i) => {
    const y = b.blocks[i];
    const shown =
      x.positions === y.positions ||
      (x.positions !== null && y.positions !== null && OrderedSet.areEqual(x.positions, y.positions));
    return (
      x.unit === y.unit &&
      x.alignment === y.alignment &&
      `${x.members}` === `${y.members}` &&
      x.length === y.length &&
      x.start === y.start &&
      x.stop === y.stop &&
      x.dir === y.dir &&
      shown
    );
  });

/**
 * One view: `resolveSection` per section, `resolveLayer` per layer, then `composeView`. A document
 * is normalized first; a DocIndex is taken as an index of a normalized document. `supported`: the drawing viewer's
 * layer types for the view's layout, its themes and glyphs (`supportedAll` for every registered one).
 */
export function resolveView(
  doc: SeqViewSpecDoc | DocIndex,
  set: SeqSet,
  view: string,
  reg: SvsRegistries,
  supported: ReadonlySet<string>,
): { view: ResolvedView; issues: Issue[] } {
  const index = "nodes" in doc ? doc : indexDoc(doc, reg);
  const issues: Issue[] = [];
  const v = index.nodes.get(view);
  if (v?.kind !== "view") {
    issues.push(Issue.error("svs.ref.unresolved", "", `"${view}" is not a view`));
    return { view: composeView(view, [], [], set), issues };
  }
  const sections: ResolvedSection[] = [];
  const layers: ResolvedLayer[][] = [];
  for (const c of childrenOf(index, v)) {
    if (c.kind !== "section") continue;
    const r = resolveSection(index, set, c.ref, reg);
    for (const i of r.issues) issues.push(i);
    // A template track's instances in place of its layers; an instance placing nothing left out (`reaching`).
    const s = r.section;
    const instances = new Map<number, string[]>();
    const own = new Map<string, ResolvedLayer>();
    const resolve = (ref: string, node?: string) => {
      const x = resolveLayer(s, ref, index, set, reg, supported, node);
      for (const i of x.issues) issues.push(i);
      return x.layer;
    };
    s.tracks.forEach((tr, k) => {
      const ref = set.strings.get(tr.ref);
      const fe = forEachOf(index.nodes.get(ref));
      const ls = s.layers.filter((l) => l.track === k);
      if (!fe) for (const l of ls) own.set(l.ref, resolve(l.ref));
      else
        for (const node of templateNodes(ref, index, set, undefined, reg, (where, m) =>
          issues.push(issueOf(reg, "svs.layer.data", `${index.nodes.get(ref)?.path}/params/for_each/${where}`, m)),
        )) {
          const xs = ls.map((l) => resolve(l.ref, node));
          if (fe.reaching && xs.every((x) => x.placements.length === 0)) continue;
          for (const x of xs) own.set(x.ref, x);
          instances.set(k, [...(instances.get(k) ?? []), node]);
        }
      if (fe && !instances.has(k)) instances.set(k, []);
    });
    const section = instantiate(s, instances, index, set);
    sections.push(section);
    layers.push(section.layers.flatMap((l) => own.get(l.ref) ?? []));
  }
  return { view: composeView(view, sections, layers, set, issues), issues };
}
