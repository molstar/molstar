/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `resolveSection` (seqspec:format#sections, seqspec:format#blocks, seqspec:format#placement, seqspec:format#focus):
// one section's axis (explicit, complex and default blocks, copies merged, ligand forms; each block a piece of its
// coordinate system with its layout map), `positions` and coverage, numberings through paths, focus, groups, tracks and
// lanes. Blocks of an unavailable node are left out quietly: nothing cascades. Indices are the block coordinate
// system's: a ranged unit's from its `rangeStart`.

import {
  type ComposedPath,
  type CoordinateMap,
  type CoordinateSystem,
  CoordinateSystemKind,
  Coords,
  type Issue,
  KeyKind,
  type KeyTarget,
  type MapHit,
  MoleculeKind,
  type NameId,
  NO_LABEL,
  Numbering,
  OrderedSet,
  type PositionLabel,
  type SeqSet,
  SeqSets,
  type SortedArray,
  SortedRanges,
  type UnitId,
} from "@molstar/seq-core/model";
import type { DocIndex } from "../load/compile.js";
import { issueOf } from "../load/issues.js";
import type { SvsRegistries } from "../registry.js";
import type { IndexedNode } from "../schema/tree-schema.js";
import { resolveStyle, type StyleRenderProps } from "../tables/style.js";
import { Arrangement, Height, Layout, LigandForm, RowSizing, SectionAxis } from "./codes.js";
import { layoutOf } from "./coordinates.js";
import { alignmentRows, type KeyFilterValue, keyed, type OrderKeyValue } from "./keyed.js";
import { blockSystem, columnsSystem, cover, endsOf, matcher, pathTo, reach, steps, unitSystem } from "./reach.js";

export type ResolvedParams = Readonly<Record<string, unknown>>;
export interface BlockRow {
  /** The block nodes that gave it (a complex block's units, merged copies); none for default blocks. */
  readonly nodes: readonly string[];
  readonly unit: number;
  readonly alignment: number;
  /** Merged copies; null: one unit. */
  readonly members: SortedArray<UnitId> | null;
  /** Its coordinate system's domain start: per-position arrays index `i - base`. */
  readonly base: number;
  /** The piece, `[start, stop)` as written (`range`), laid out in `dir`; `first`: the index at shown position 0. */
  readonly start: number;
  readonly stop: number;
  readonly dir: number;
  readonly first: number;
  /** The layout map: block coordinate system → shown positions. */
  readonly layout: CoordinateMap;
  readonly length: number;
  readonly label: number;
  readonly ligand: number;
  readonly positions: OrderedSet | null;
  readonly varying: OrderedSet | null;
  readonly numbering: Numbering | null;
  /** Into the section's paths, from the numbering's coordinate system to the block's; -1: on the block itself. */
  readonly numberingPath: number;
  readonly params: ResolvedParams;
}
export interface GroupRow {
  readonly ref: number;
  readonly parent: number;
  readonly label: number;
  readonly collapsed: number;
  readonly hidden: number;
  readonly params: ResolvedParams;
}
export interface TrackRow {
  readonly ref: number;
  readonly group: number;
  readonly label: number;
  readonly height: number;
  readonly hidden: number;
  readonly laneStart: number;
  readonly laneCount: number;
  /**
   * Its split layer, an index into the section's `layers` (seqspec:format#layer-types; matched by ref when composed);
   * -1: none.
   */
  readonly split: number;
  /** That layer's rows at open: `expanded`, else expanded for `variants`. */
  readonly expanded: number;
  readonly params: ResolvedParams;
}
/**
 * A lazy alignment's every row as lanes, never stored: lane `start + k` is row `order[k]` (null: `k`, data order), its
 * depth by `parent` (the alignment's, shared, not copied).
 */
export interface LazyLanes {
  readonly start: number;
  readonly length: number;
  readonly order: Int32Array | null;
  readonly parent: Int32Array;
}
/** A row's depth: the variants above it (`parent` chains, bounded by its length). */
export function depthOf(parent: Int32Array, row: number): number {
  let d = 0;
  for (let q = parent[row]; q >= 0 && d < parent.length; q = parent[q]) d++;
  return d;
}
export interface SectionRow {
  readonly ref: number;
  readonly axis: number;
  readonly arrangement: number;
  readonly ligandForm: number;
  readonly height: number;
  readonly heightRows: number;
  readonly rowSizing: number;
  readonly hidden: number;
  /** The preset `style` names (an inline preset's `preset`); -1: none. */
  readonly style: number;
  /** An inline preset's keys as render props; null: none. */
  readonly render: StyleRenderProps | null;
  readonly scrollGroup: number;
  readonly focus: { block: number; start: number; stop: number } | null;
  readonly params: ResolvedParams;
}
/** One section's rows of the resolved tables, indices local to it. */
export interface ResolvedSection {
  readonly ref: string;
  /** The section node's path, for issues about it. */
  readonly path: string;
  readonly layout: number;
  readonly row: SectionRow;
  readonly blocks: readonly BlockRow[];
  /**
   * The blocks the `positions` filter left out (it reaches none of their positions), in order: a `pairs` layer keeps
   * the pairs with an end there, drawing the other end.
   */
  readonly leftOut: readonly BlockRow[];
  readonly groups: readonly GroupRow[];
  readonly tracks: readonly TrackRow[];
  /**
   * Stored lanes first (`row`, `depth`), then each lazy run's (seqspec:sources#data-by-reference: a million rows are
   * never listed).
   */
  readonly lanes: {
    readonly row: readonly number[];
    readonly depth: readonly number[];
    readonly lazy: readonly LazyLanes[];
  };
  /**
   * Its layers in order, by ref, with their track and declared type: `composeView` takes each layer's track from here.
   */
  readonly layers: readonly { readonly ref: string; readonly track: number; readonly type: string }[];
  readonly paths: readonly ComposedPath[];
  /** What `resolveLayer` reads (the axis kind and the layout, units, members, alignments, lengths, ligand forms, block nodes, lanes, and the layers
   * of each track with lanes): equal keys give equal layers, so a viewer reuses them while `positions`, a numbering
   * or the order of tracks changes. */
  readonly axisKey: string;
}

const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const kindOf = (set: SeqSet, ref: unknown) => (typeof ref === "string" && set.byName(ref, t) ? t.kind : -1);
const head = (ref: string) => ref.split("/")[0];
export const refsOf = (v: unknown): string[] => (typeof v === "string" ? [v] : Array.isArray(v) ? v : []);
/** A node's params with its `custom`, as the resolved tables carry them. */
export const paramsOf = (n: IndexedNode): ResolvedParams =>
  n.node.custom ? { ...n.params, custom: n.node.custom } : n.params;
export const childrenOf = (doc: DocIndex, n: IndexedNode) =>
  (n.node.children ?? []).map((c) => doc.nodes.get(c.ref ?? c.kind) as IndexedNode);
/** The most positions one block lays out: the Int32 the axis and the model hold (VL-28). */
const MaxBlockPositions = 0x7fffffff;
const flag = (v: unknown) => (v === true ? 1 : 0);
const LABEL = /^(-?[0-9]+)([A-Za-z]+|\.[0-9]+)?$/;

/** A document label (an int, `"27A"`, `"111.1"` or a name) in a numbering's terms; false when it is none. */
export function toLabel(set: SeqSet, n: Numbering, v: unknown, out: PositionLabel): boolean {
  out.num = NO_LABEL;
  out.suffix = out.name = -1 as NameId;
  out.name = n.method === "names" ? set.strings.find(String(v)) : out.name;
  if (n.method === "names") return out.name >= 0;
  const m = LABEL.exec(String(v));
  if (!m) return false;
  out.num = Number(m[1]);
  if (m[2]) out.suffix = set.strings.intern(m[2]);
  return true;
}

interface Entry {
  unit: number;
  alignment: number;
  node: string | null;
  params: ResolvedParams;
}

/** The section's axis, explicit or by default. */
function entriesOf(doc: DocIndex, set: SeqSet, blocks: IndexedNode[], say: Say): Entry[] {
  const out: Entry[] = [];
  const push = (unit: number, alignment: number, b: IndexedNode | null) =>
    out.push({ unit, alignment, node: b?.ref ?? null, params: b ? paramsOf(b) : {} });
  for (const b of blocks) {
    const param = ["unit", "sequence", "ligand", "complex", "alignment"].find((k) => b.params[k] != null) as string;
    const ref = b.params[param] as string;
    const at = `${b.path}/params/${param}`;
    // A ref under a node the set does not hold (a streamed source's dynamic node not held yet) is left out quietly.
    const k = kindOf(set, ref);
    if (k < 0 && kindOf(set, head(ref)) < 0) continue;
    if (k < 0) say("svs.ref.derived_unknown", at, `"${ref}" names nothing its node yielded`);
    else if (param === "alignment" || param === "unit") {
      const unit = param === "unit";
      if (k === (unit ? KeyKind.unit : KeyKind.alignment)) push(unit ? t.id : -1, unit ? -1 : t.id, b);
    } else if (param === "complex") for (const u of set.complexes[t.id]?.units ?? []) push(u, -1, b);
    else if (
      k !== KeyKind.molecule ||
      set.molecules[t.id]?.kind !== (param === "ligand" ? MoleculeKind.ligand : MoleculeKind.polymer)
    ) {
      say("svs.ref.wrong_kind", at, `"${ref}" is not a ${param}`);
    } else if (kindOf(set, `${ref}/unit`) !== KeyKind.unit) {
      say("svs.ref.wrong_kind", at, `"${ref}" is a provider's ${param}, with no own unit: name its units`);
    } else push(t.id, -1, b);
  }
  if (blocks.length > 0) return out;
  const data = doc.order.filter((n) => n.parent === "data");
  const units: number[] = [];
  for (const d of data) {
    const k = kindOf(set, d.ref);
    if (d.kind === "unit" && k === KeyKind.unit) units.push(t.id);
    if (d.kind === "complex" && k === KeyKind.complex) for (const u of set.complexes[t.id]?.units ?? []) units.push(u);
  }
  if (units.length === 0) {
    for (const d of data) {
      if ((d.kind === "sequence" || d.kind === "ligand") && kindOf(set, `${d.ref}/unit`) === KeyKind.unit) {
        units.push(t.id);
      }
    }
  }
  const isLigand = (u: number) => set.molecules[set.units.molecule[u]]?.kind === MoleculeKind.ligand;
  for (const u of [...units.filter((u) => !isLigand(u)), ...units.filter(isLigand)]) push(u, -1, null);
  const aln = units.length
    ? undefined
    : data.find((d) => d.kind === "alignment" && kindOf(set, d.ref) === KeyKind.alignment);
  if (aln) push(-1, t.id, null);
  return out;
}

type Say = (code: string, path: string, message: string) => void;
const ordinal = (set: SeqSet) => set.strings.intern("ordinal");

/** A block's numbering: its unit's or molecule's, else one a path reaches, carried through it. */
function numberingOf(set: SeqSet, b: BlockRow, name: unknown, paths: ComposedPath[], say: Say, at: string) {
  const none = { numbering: null as Numbering | null, path: -1 };
  const u = b.unit;
  const mol = u >= 0 ? set.molecules[set.units.molecule[u]] : null;
  const fallback = () => {
    if (!mol) return none;
    const all = [...(set.units.numberings[u]?.items ?? []), ...mol.numberings.items];
    return { numbering: all.find((n) => n.isDefault) ?? mol.numberings.byName(ordinal(set)), path: -1 };
  };
  if (typeof name !== "string") return fallback();
  const id = set.strings.intern(name);
  const aln = set.alignments[b.alignment];
  const from: CoordinateSystem = aln
    ? { kind: CoordinateSystemKind.row, id: aln.id, row: aln.referenceRow }
    : unitSystem(u);
  const own = aln ? aln.numberingsOf(aln.referenceRow)?.byName(id) : SeqSets.numberingOf(set, u as UnitId, id);
  const found: { numbering: Numbering; path: ComposedPath | null }[] = [];
  const rowPath = aln ? pathTo(set, from, columnsSystem(aln.id)) : null;
  if (own) found.push({ numbering: own, path: rowPath && rowPath !== "ambiguous" ? rowPath : null });
  else if (mol?.kind === MoleculeKind.ligand) return fallback();
  let best = Number.POSITIVE_INFINITY;
  for (const m of own ? [] : set.molecules) {
    const n = m?.numberings.byName(id);
    const mid = { kind: CoordinateSystemKind.molecule, id: m?.id as number, row: -1 };
    const p = n ? pathTo(set, mid, blockSystem(b)) : undefined;
    // Labels are never carried across a codon step, either way: a residue's label does not name a nucleotide.
    if (!n || !p || p === "ambiguous" || p.codon >= 0 || steps.length > best) continue;
    if (steps.length < best) found.length = 0;
    best = steps.length;
    found.push({ numbering: n, path: p });
  }
  if (found.length === 0) {
    say(
      "svs.numbering.fallback",
      at,
      `no numbering "${name}" on this block or anything a path reaches from it; its default numbering is used`,
    );
    return fallback();
  }
  if (found.length > 1 && !found.every((f) => sameLabels(f, found[0], b.start, b.stop))) {
    say("svs.mapping.ambiguous", at, `paths to ${found.length} numberings "${name}" give different labels`);
  }
  const p = found[0].path;
  return { numbering: found[0].numbering, path: p ? paths.push(p) - 1 : -1 };
}

const la: PositionLabel = { num: 0, suffix: -1 as NameId, name: -1 as NameId };
const lb: PositionLabel = { num: 0, suffix: -1 as NameId, name: -1 as NameId };
const hit = { element: -1, phase: 0, first: -1 };
const takeFirst = (h: MapHit, o: typeof hit) => {
  if (o.first < 0) o.first = h.element;
};
/** The path's from index whose image is block index `i`, the first copy; -1: none. */
function sourceAt(p: ComposedPath, i: number): number {
  if (!Coords.sameSystem(p.map.from, p.from)) return Coords.toAt(p.map, i, 0);
  hit.first = -1;
  Coords.forEachFromAt(p.map, i, hit, takeFirst, hit);
  return hit.first;
}
/** Two numberings label the block's indices `[lo, hi)` alike, each through its path. */
function sameLabels(a: { numbering: Numbering; path: ComposedPath | null }, b: typeof a, lo: number, hi: number) {
  const at = (x: typeof a, i: number, out: PositionLabel) => {
    const j = x.path ? sourceAt(x.path, i) : i;
    return j >= 0 && Numbering.labelOf(x.numbering, j as never, out);
  };
  for (let i = lo; i < hi; i++) {
    const [ha, hb] = [at(a, i, la), at(b, i, lb)];
    if (ha !== hb || (ha && (la.num !== lb.num || la.suffix !== lb.suffix || la.name !== lb.name))) return false;
  }
  return true;
}

/** Positions where merged members' observed sets differ. */
function varyingOf(set: SeqSet, members: Int32Array): OrderedSet | null {
  const { rangeStart, rangeStop, observed } = set.units;
  const has = (u: number, i: number) => !observed[u] || SortedRanges.has(observed[u] as never, i);
  const out: number[] = [];
  for (let i = rangeStart[members[0]]; i < rangeStop[members[0]]; i++) {
    if (members.some((m) => has(m, i) !== has(members[0], i))) out.push(i);
  }
  return out.length ? OrderedSet.ofSortedArray(Int32Array.from(out)) : null;
}

/** The positions the masks mark, or null when none (the block is not reached). */
function marked(mask: Uint8Array, lo: number, hi: number): OrderedSet | null {
  const out: number[] = [];
  // A piece across a cyclic origin (hi < lo; hi === lo: all of it from lo) holds [0, hi) and [lo, L).
  const wraps = hi <= lo;
  for (let i = wraps ? 0 : lo; i < (wraps ? mask.length : hi); i++)
    if (mask[i] && (!wraps || i < hi || i >= lo)) out.push(i);
  return out.length ? OrderedSet.ofSortedArray(Int32Array.from(out)) : null;
}

/**
 * A unit over all of a cyclic molecule (seqspec:model#cyclic-systems): its pieces may wrap the origin and start at
 * another.
 */
function cyclicUnit(set: SeqSet, u: number): boolean {
  const m = set.molecules[set.units.molecule[u]];
  return !!m?.cyclic && set.units.rangeStart[u] === 0 && set.units.rangeStop[u] === m.length;
}

/** The shortest arc of a cyclic system of `length` holding the positions (`Coords.cyclicHull`): start > stop wraps. */
function arcOf(s: OrderedSet, length: number): { start: number; stop: number } {
  const runs: number[] = [];
  OrderedSet.forEach(
    s,
    (i) => {
      if (runs.length && runs[runs.length - 1] === i) runs[runs.length - 1]++;
      else runs.push(i, i + 1);
    },
    null,
  );
  const n = Coords.cyclicHull(Int32Array.from(runs), runs.length >> 1, length, runsA);
  return { start: runsA[0], stop: runsA[2 * n - 1] };
}

const runsA = new Int32Array(4);
const runsB = new Int32Array(4);
/** A block's piece and [s, t) meet; either may wrap a cyclic origin (a piece's start === stop: all of it). */
function meets(b: BlockRow, s: number, t: number): boolean {
  const L = b.base + b.length;
  const n = b.start === b.stop ? Coords.cyclicRuns(b.base, L, L, runsA) : Coords.cyclicRuns(b.start, b.stop, L, runsA);
  const m = Coords.cyclicRuns(s, t, L, runsB);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < m; j++) if (runsA[2 * i] < runsB[2 * j + 1] && runsB[2 * j] < runsA[2 * i + 1]) return true;
  return false;
}

/** A 32-bit hash of an order's contents (an axis key's; exported for its test). */
export function orderHash(order: ArrayLike<number>): string {
  let h = 2166136261;
  for (let i = 0; i < order.length; i++) h = Math.imul(h ^ order[i], 16777619);
  return (h >>> 0).toString(36);
}

/** One section's axis, blocks, groups, tracks, lanes and focus. */
export function resolveSection(
  doc: DocIndex,
  set: SeqSet,
  section: string,
  reg: SvsRegistries,
): { section: ResolvedSection; issues: Issue[] } {
  const issues: Issue[] = [];
  const say: Say = (code, path, message) => issues.push(issueOf(reg, code, path, message));
  const n = doc.nodes.get(section) as IndexedNode;
  const p = n.params;
  const view = doc.nodes.get(n.parent ?? "");
  const kids = childrenOf(doc, n);
  const str = (v: unknown) => (typeof v === "string" ? set.strings.intern(v) : -1);
  const paths: ComposedPath[] = [];

  // Blocks, copies merged at the first one's place.
  const entries = entriesOf(
    doc,
    set,
    kids.filter((c) => c.kind === "block"),
    say,
  );
  // Copies merge when they show the same piece: two pieces of one unit stay two blocks.
  const merged: Entry[][] = [];
  const byGroup = new Map<string, Entry[]>();
  const dirOf = (e: Entry) => ((e.params.direction ?? p.direction) === -1 ? -1 : 1);
  for (const e of entries) {
    const piece = `${e.params.range}|${e.params.origin}|${dirOf(e)}`;
    const g = e.unit >= 0 && p.copies !== "each" ? `${set.units.group[e.unit]}|${piece}` : "";
    const into = byGroup.get(g);
    if (into && g) into.push(e);
    else merged.push(g ? (byGroup.set(g, [e]).get(g) as Entry[]) : [e]);
  }
  // Each entry of `positions`: a data ref, or `{ data, filter }`, whose covered rows are those the filter keeps.
  const positions = (Array.isArray(p.positions) ? p.positions : p.positions == null ? [] : [p.positions]).flatMap(
    (e: string | { data: string; filter: { field?: string | null; values: unknown[] } }) => {
      const [r, f] = typeof e === "string" ? [e, null] : [e.data, e.filter];
      const ann = kindOf(set, r) === KeyKind.annotation ? set.annotations[t.id] : undefined;
      if (!ann || ann.geometry === "rows") return [];
      const keeps = f && matcher(f.values);
      const field = f?.field ?? "category";
      return [{ ann, keep: keeps ? (row: number) => keeps(ann, set, row, field) : undefined }];
    },
  );
  const blocks: BlockRow[] = [];
  const leftOut: BlockRow[] = [];
  for (const group of merged) {
    const e = group[0];
    const aln = set.alignments[e.alignment];
    const uniq = [...new Set(group.map((x) => x.unit))];
    const members = uniq.length > 1 ? (Int32Array.from(uniq).sort() as SortedArray<UnitId>) : null;
    const mol = aln ? null : set.molecules[set.units.molecule[e.unit]];
    const ligand = mol?.kind === MoleculeKind.ligand && p.ligand_form !== "atoms" ? 1 : 0;
    const unitLabel = (u: number) => set.strings.get(set.units.label[u]);
    let label = e.params.label as string | null;
    if (label == null && members) {
      label = `${set.strings.get(mol?.label ?? -1)} ×${uniq.length} (${uniq.map(unitLabel).join(", ")})`;
    }
    const nodes = [...new Set(group.flatMap((x) => (x.node ? [x.node] : [])))];
    // The piece: `range` of the domain or all of it (a ligand's one position), laid out in `direction`. On a whole
    // cyclic unit a range may wrap the origin (start > stop), and `origin` starts the whole piece elsewhere.
    const base = aln ? 0 : set.units.rangeStart[e.unit];
    const length = aln ? aln.columnCount : ligand ? 1 : set.units.rangeStop[e.unit] - base;
    const cyclic = !aln && !ligand && cyclicUnit(set, e.unit);
    const range = ligand ? null : (e.params.range as [number, number] | null);
    const params = `${e.node ? doc.nodes.get(e.node)?.path : n.path}/params`;
    const wraps = range != null && range[1] < range[0] && cyclic && range[0] < length && range[1] >= 0;
    if (range && (range[0] < base || range[1] > base + length || (range[1] <= range[0] && !wraps))) {
      const why = range[1] < range[0] && !cyclic ? "; start > stop wraps only a cyclic sequence's origin" : "";
      say(
        "svs.position.range",
        `${params}/range`,
        `[${range}] is outside the unit's domain [${base}, ${base + length})${why}: left out`,
      );
      continue;
    }
    let origin = range ? 0 : ((e.params.origin as number | undefined) ?? 0);
    if (origin !== 0 && !(cyclic && origin > 0 && origin < length)) {
      const where = cyclic ? `within [0, ${length})` : "on a cyclic sequence's whole unit";
      say("svs.position.range", `${params}/origin`, `origin ${origin} is ${where} only: ignored`);
      origin = 0;
    }
    const [from, to] = origin ? [origin, origin] : range ? range : [base, base + length];
    const [start, stop, dir] = [from, to, dirOf(e)];
    // VL-28: the layout, the viewer's axis and the model's tables hold a block's positions as Int32 (a section of
    // blocks may hold more: a whole genome), so a block past `MaxBlockPositions` gets an issue here, not a wrap later.
    const count = start === stop ? length : stop > start ? stop - start : length - start + stop;
    if (count > MaxBlockPositions) {
      const f = (x: number) => x.toLocaleString("en-US");
      say(
        "svs.position.range",
        params,
        `the block holds ${f(count)} positions, past the ${f(MaxBlockPositions)} a block lays out: left out; show a \`range\` of it`,
      );
      continue;
    }
    const piece = layoutOf(aln ? columnsSystem(aln.id) : unitSystem(e.unit), start, stop, dir, length);
    const block: BlockRow = {
      nodes,
      unit: e.unit,
      alignment: aln ? aln.id : -1,
      members,
      base,
      start,
      stop,
      dir,
      first: piece.first,
      layout: piece.map,
      length,
      label: str(
        label ??
          (aln
            ? (doc.nodes.get(set.strings.get(aln.key))?.params.label ?? set.strings.get(aln.key))
            : unitLabel(e.unit)),
      ),
      ligand,
      positions: null,
      varying: members ? varyingOf(set, members) : null,
      numbering: null,
      numberingPath: -1,
      params: e.params,
    };
    let shown: OrderedSet | null = null;
    if (positions.length > 0) {
      const mask = new Uint8Array(aln ? aln.columnCount : (mol?.length ?? 0));
      for (const { ann, keep } of positions) {
        for (const end of endsOf(ann as never))
          reach(set, end.system, block, (_m, path) => cover(ann as never, end, path, mask, keep));
      }
      shown = marked(mask, start, stop);
      if (!shown) {
        leftOut.push(block);
        continue;
      }
    }
    const at = e.node ? (doc.nodes.get(e.node)?.path ?? n.path) : n.path;
    const num = numberingOf(set, block, e.params.numbering ?? p.numbering, paths, say, at);
    blocks.push({ ...block, positions: ligand ? null : shown, numbering: num.numbering, numberingPath: num.path });
  }

  // Groups and tracks in order; a `rows` track's lanes (a lazy alignment's all rows: a run, placed after the stored).
  const groups: GroupRow[] = [];
  const tracks: TrackRow[] = [];
  const layers: { ref: string; track: number; type: string }[] = [];
  const lazy: (LazyLanes & { track: number })[] = [];
  const lanes = { row: [] as number[], depth: [] as number[], lazy: lazy as LazyLanes[] };
  const columns = blocks.length === 1 ? set.alignments[blocks[0].alignment] : undefined;
  let rowsLayer: string | null = null;
  const walk = (parent: IndexedNode, group: number) => {
    for (const c of childrenOf(doc, parent)) {
      const cp = c.params;
      if (c.kind === "group") {
        groups.push({
          ref: str(c.ref),
          parent: group,
          label: str(cp.label),
          collapsed: flag(cp.collapsed),
          hidden: flag(cp.hidden),
          params: paramsOf(c),
        });
        walk(c, groups.length - 1);
      }
      if (c.kind !== "track") continue;
      // The section's first `rows` layer wins; a later one is left out.
      const ls = childrenOf(doc, c).filter((l) => {
        if (l.kind !== "layer" || l.params.type !== "rows") return l.kind === "layer";
        rowsLayer ??= l.ref;
        if (rowsLayer === l.ref) return true;
        say("svs.layer.rows", l.path, `a section has one rows layer; "${rowsLayer}" is drawn, this one left out`);
        return false;
      });
      const laneStart = lanes.row.length;
      const rowsNode = ls.find((l) => l.params.type === "rows");
      const rows = rowsNode?.params;
      let lazyRows = 0;
      if (rows && columns) {
        const subset = rows.row_subset ? new Set(rows.row_subset as number[]) : null;
        const given = (rows.filter ?? []) as KeyFilterValue[];
        const asked = rows.order as OrderKeyValue | OrderKeyValue[] | null;
        // A streamed alignment is ordered by tree or data only: its tiles are fetched by lane, so a permutation by a key
        // table would draw one row's letters under another's label. The key order and filter are left out, and said so.
        const keys = columns.lazy && asked != null && asked !== "tree" && asked !== "data";
        const tell = (where: string, what: string) =>
          say(
            "svs.provider.lazy",
            `${rowsNode?.path}/params/${where}`,
            `a streamed alignment is ordered by tree or data only: its ${what} over key tables is left out`,
          );
        if (keys) tell("order", "order");
        if (columns.lazy && given.length > 0) tell("filter", "filter");
        const order = keys ? null : asked;
        const filter = columns.lazy ? [] : given;
        // Rows with a parent are indented under it only where they stay under it: in tree or data order.
        const first = Array.isArray(order) ? order[0] : order;
        const indent = first == null || first === "tree" || first === "data";
        // A lazy alignment's rows are never listed: a run, in tree or data order.
        if (columns.lazy && !subset) {
          lazyRows = columns.rowCount;
          const byTree = (order ?? (columns.tree ? "tree" : "data")) === "tree" && columns.tree;
          lazy.push({
            track: tracks.length,
            start: -1,
            length: lazyRows,
            order: byTree ? columns.tree.leafOrder : null,
            parent: columns.parent,
          });
        } else {
          const members = alignmentRows(columns.id, columns.rowCount, columns.tree ? columns.tree.leafOrder : null);
          const k = keyed(set, reg, members, { order, filter, subset }, (where, message) =>
            say("svs.layer.data", `${rowsNode?.path}/params/${where}`, message),
          );
          for (const r of k.permutation) {
            if (!k.kept[r]) continue;
            lanes.row.push(r);
            lanes.depth.push(indent ? depthOf(columns.parent, r) : 0);
          }
        }
      }
      const first =
        ls[0] && (refsOf(ls[0].params.data)[0] ?? childrenOf(doc, ls[0]).find((x) => x.kind === "annotation")?.ref);
      const label =
        cp.label != null
          ? str(cp.label)
          : kindOf(set, first) === KeyKind.annotation
            ? (set.annotations[t.id]?.label ?? -1)
            : -1;
      const split = ls.findIndex((l) => l.params.split != null);
      tracks.push({
        ref: str(c.ref),
        group,
        label,
        height: typeof cp.height === "number" ? cp.height : Number.NaN,
        hidden: flag(cp.hidden),
        laneStart,
        laneCount: lazyRows || lanes.row.length - laneStart,
        split: split < 0 ? -1 : layers.length + split,
        expanded: flag(cp.expanded ?? ls[split]?.params.type === "variants"),
        params: paramsOf(c),
      });
      for (const l of ls) layers.push({ ref: l.ref, track: tracks.length - 1, type: String(l.params.type) });
    }
  };
  walk(n, -1);
  let next = lanes.row.length;
  for (let i = 0; i < lazy.length; i++) {
    lazy[i] = { ...lazy[i], start: next };
    tracks[lazy[i].track] = { ...tracks[lazy[i].track], laneStart: next };
    next += lazy[i].length;
  }

  const axis = blocks.some((b) => b.alignment >= 0) ? SectionAxis.columns : SectionAxis.blocks;
  const focusNode = kids.filter((c) => c.kind === "focus").pop();
  const style = resolveStyle(p.style);
  const row: SectionRow = {
    ref: str(n.ref),
    axis,
    arrangement: p.arrangement === "stacked" ? Arrangement.stacked : Arrangement.side,
    ligandForm: p.ligand_form === "atoms" ? LigandForm.atoms : LigandForm.position,
    height: p.height === "content" ? Height.content : Height.box,
    heightRows: typeof p.height_rows === "number" ? p.height_rows : -1,
    rowSizing: p.row_sizing === "natural" ? RowSizing.natural : RowSizing.uniform,
    hidden: flag(p.hidden),
    style: str(style.preset),
    render: style.render,
    scrollGroup: str(p.scroll_group),
    focus: focusNode ? focusOf(set, focusNode, blocks, say) : null,
    params: paramsOf(n),
  };
  const layout = view?.params.layout === "text" ? Layout.text : Layout.tracks;
  const axisKey = [
    axis,
    layout,
    ...blocks.map((b) => `${b.nodes}|${b.unit}|${b.alignment}|${b.members}|${b.length}|${b.ligand}`),
    lanes.row.join(","),
    // A run's order by its contents too, as the stored lanes' rows are: a changed permutation is another axis.
    lazy.map((z) => `${z.start}+${z.length}${z.order ? `t${orderHash(z.order)}` : "d"}`).join(","),
    // A layer reads its track only for lanes: which layers sit on a track with lanes, never the order of tracks.
    tracks
      .flatMap((x, i) => {
        const on = x.laneCount > 0 ? layers.filter((l) => l.track === i).map((l) => l.ref) : [];
        return on.length ? [`${x.laneStart}:${x.laneCount}:${on}`] : [];
      })
      .join(","),
  ].join(";");
  const out = { ref: n.ref, path: n.path, layout, row, blocks, leftOut, groups, tracks, lanes, layers, paths, axisKey };
  return { section: out, issues };
}

/** The range a section opens on, in a block's positions; null (reported) when outside. */
function focusOf(set: SeqSet, f: IndexedNode, blocks: BlockRow[], say: Say) {
  const fp = f.params;
  const fail = (why: string) => {
    say("svs.focus.range", f.path, `the focus ${why}`);
    return null;
  };
  if (fp.annotation != null) {
    if (kindOf(set, fp.annotation) !== KeyKind.annotation) return null;
    const ann = set.annotations[t.id];
    const items = fp.items as string[] | null;
    const keep = (r: number) => !items || items.includes(set.strings.get(ann?.itemKey[ann.item[r]] ?? -1));
    for (let block = 0; ann && block < blocks.length; block++) {
      const mask = new Uint8Array(
        blocks[block].alignment >= 0
          ? blocks[block].length
          : (set.molecules[set.units.molecule[blocks[block].unit]]?.length ?? 0),
      );
      for (const end of endsOf(ann))
        reach(set, end.system, blocks[block], (_m, path) => cover(ann, end, path, mask, keep));
      const s = marked(mask, blocks[block].start, blocks[block].stop);
      if (s && blocks[block].alignment < 0 && cyclicUnit(set, blocks[block].unit))
        return { block, ...arcOf(s, blocks[block].length) };
      if (s) return { block, start: OrderedSet.min(s), stop: OrderedSet.max(s) + 1 };
    }
    return fail(`annotation "${fp.annotation}" covers nothing in this section`);
  }
  const k = kindOf(set, fp.target);
  // An unavailable target is ignored quietly; a derived ref its node did not yield is reported.
  if (k < 0 && kindOf(set, head(String(fp.target))) >= 0)
    say("svs.ref.derived_unknown", `${f.path}/params/target`, `"${fp.target}" names nothing its node yielded`);
  if (k < 0) return null;
  const id = t.id;
  const row = t.row;
  const aln = k === KeyKind.alignment || k === KeyKind.row ? set.alignments[id] : null;
  const mol = k === KeyKind.unit ? set.units.molecule[id] : k === KeyKind.molecule ? id : -1;
  const shows = (b: BlockRow) =>
    aln
      ? b.alignment === id
      : b.unit >= 0 &&
        (b.unit === id || b.members?.includes(id) || (k === KeyKind.molecule && set.units.molecule[b.unit] === mol));
  if (!blocks.some(shows)) return fail(`target "${fp.target}" is not shown in this section`);
  // The target's domain: a unit's is its range.
  const [lo, hi] = aln
    ? [0, row >= 0 ? aln.rows.lengthOf(row as never) : aln.columnCount]
    : k === KeyKind.unit
      ? [set.units.rangeStart[id], set.units.rangeStop[id]]
      : [0, set.molecules[mol]?.length ?? 0];
  // A whole cyclic molecule's or unit's range may wrap its origin (start > stop).
  const cyclic = !aln && !!set.molecules[mol]?.cyclic && lo === 0 && hi === set.molecules[mol]?.length;
  let [start, stop] = [fp.start as number | null, fp.stop as number | null];
  if (fp.beg != null) {
    const name = set.strings.intern(String(fp.numbering));
    const num =
      k === KeyKind.row
        ? aln?.numberingsOf(row as never)?.byName(name)
        : k === KeyKind.unit
          ? SeqSets.numberingOf(set, id as UnitId, name)
          : set.molecules[mol]?.numberings.byName(name);
    const out = new Int32Array(2);
    const end = fp.end ?? fp.beg;
    const okBeg = !!num && toLabel(set, num, fp.beg, la);
    const ok = okBeg && toLabel(set, num, end, lb);
    if (!num || !ok || !Numbering.rangeOf(num, la, lb, out, cyclic)) {
      // Written backwards (CF-01); else the label that is not there.
      if (num && ok && Numbering.backwards(num, la, lb))
        return fail(
          cyclic
            ? `beg "${fp.beg}" comes right after end "${end}" in numbering "${fp.numbering}": a range across the origin covers less than a whole lap`
            : `beg "${fp.beg}" comes after end "${end}" in numbering "${fp.numbering}": write beg ≤ end`,
        );
      if (num && ok) Numbering.placeRange(num, la, lb, out);
      const missing = !okBeg || (ok && out[0] < 0) ? fp.beg : end;
      return fail(`label ${missing} is not in numbering "${fp.numbering}" of "${fp.target}"`);
    }
    [start, stop] = [out[0], out[1]];
  }
  start ??= lo;
  stop ??= hi;
  const wraps = cyclic && start > stop && start < hi && stop >= 0;
  if (!(start >= lo && (start < stop || wraps) && stop <= hi))
    return fail(`[${start}, ${stop}) is outside "${fp.target}" (${lo}..${hi})`);
  if (k === KeyKind.row && aln) {
    const map = aln.rows.posMap(row as never) as Int32Array;
    const cols = Array.from(map.subarray(start, stop)).filter((c) => c >= 0);
    if (!cols.length) return fail("covers only inserts");
    [start, stop] = [cols[0], (cols.at(-1) as number) + 1];
  }
  // The first block whose piece shows part of it (two pieces of one unit are two blocks).
  const [s0, s1] = [start, stop];
  const block = blocks.findIndex((b) => shows(b) && meets(b, s0, s1));
  return block < 0 ? fail(`[${start}, ${stop}) is outside the pieces shown of "${fp.target}"`) : { block, start, stop };
}
