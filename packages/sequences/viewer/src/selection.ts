/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Selection (seqspec:interaction#selection): one loci per kind, set by behaviors and hosts with an origin, carried to
// every new set (lineage.ts: rebound or remapped within a lineage, remapped by key into a new one). `close` crosses
// coordinate systems through the edge-rate projector (seqspec:interaction#channels); the verbs take `Coordinates` too,
// converted to loci through `ctx.systems` (what names no registered system is dropped).

import { type Closure, SelectionGranularity } from "@molstar/seq-core/format";
import {
  AlignedLocation,
  AlignedLoci,
  BoundaryLocation,
  type BoundaryLoci,
  ColumnLocation,
  type ColumnLoci,
  Conversion,
  ItemLocation,
  type ItemLoci,
  type Location,
  Loci,
  OrderedSet,
  RowLocation,
  RowLoci,
  SeqLocation,
  type SeqLoci,
  type SeqSet,
  SortedArray,
} from "@molstar/seq-core/model";
import { carry, type SeqLociAny } from "./lineage.ts";
import { type Coordinates, isCoordinates, type Projector } from "./projection.ts";
import { live, type Origin, type SeqState, state, type WritableLive, type WritableSeqState } from "./state.ts";

export interface Selection {
  readonly version: number;
  readonly origin: Origin;
  readonly seq: SeqLoci | null;
  readonly aligned: AlignedLoci | null;
  readonly columns: ColumnLoci | null;
  readonly rows: RowLoci | null;
  readonly items: ItemLoci | null;
  readonly boundaries: BoundaryLoci | null;
}
/** Reused; mutated per move by a drag, drawn by the overlay. */
export interface PendingSelection {
  active: boolean;
  /** The drag removes its range from the selection (a Shift-drag pressed on a selected position, `shiftPositions: "toggle"`). */
  remove: boolean;
  section: number;
  kind: number;
  block: number;
  startShown: number;
  endShown: number;
  rowStart: number;
  rowEnd: number;
  /**
   * A dot plot's rectangle (seqspec:interaction#selection): the plot's layer, and its y range, a block's shown range
   * beside the x range; -1: none.
   */
  plotLayer: number;
  yBlock: number;
  yStartShown: number;
  yEndShown: number;
}
export interface SelectionManager {
  readonly current: SeqState<Selection>;
  set(loci: SeqLociAny | Coordinates | null, origin: Origin): void;
  add(loci: SeqLociAny | Coordinates, origin: Origin): void;
  remove(loci: SeqLociAny | Coordinates, origin: Origin): void;
  toggle(loci: SeqLociAny | Coordinates, origin: Origin): void;
  clear(origin: Origin): void;
  /** Closure over the current view (seq loci; other kinds as given); edge rate, allocates. */
  close(loci: SeqLociAny, closure: Closure): SeqLociAny;
  /** `SelectionGranularity`. */
  widen(loci: SeqLociAny, granularity: number): SeqLociAny;
  /** The loci is the whole selection of its kind. */
  isWhole(loci: SeqLociAny): boolean;
  readonly pending: { readonly current: PendingSelection; readonly version: SeqState<number> };
}

type Field = "seq" | "aligned" | "columns" | "rows" | "items" | "boundaries";
const FIELD: Record<SeqLociAny["kind"], Field> = {
  "seq-loci": "seq",
  "aligned-loci": "aligned",
  "column-loci": "columns",
  "row-loci": "rows",
  "item-loci": "items",
  "boundary-loci": "boundaries",
};
const FIELDS = Object.values(FIELD);
type Parts = { -readonly [k in Field]: SeqLociAny | null };
const NONE: Parts = { seq: null, aligned: null, columns: null, rows: null, items: null, boundaries: null };
/** Per loci kind: the element's owner and set fields, and the location's value field. */
const SHAPE = {
  "seq-loci": ["unit", "indices", "index"],
  "boundary-loci": ["unit", "indices", "at"],
  "row-loci": ["row", "indices", "index"],
  "aligned-loci": ["row", "columns", "column"],
  "item-loci": ["annotation", "items", "item"],
} as const;
const LOCATION = {
  "seq-loci": SeqLocation,
  "boundary-loci": BoundaryLocation,
  "row-loci": RowLocation,
  "aligned-loci": AlignedLocation,
  "column-loci": ColumnLocation,
  "item-loci": ItemLocation,
};
const locationOf = (loci: SeqLociAny) => LOCATION[loci.kind].create(loci.set) as Location;
const ofList = (xs: number[]) => OrderedSet.ofSortedArray(SortedArray.ofUnsortedArray(Int32Array.from(xs)));
const empty = (l: SeqLociAny | null) => !l || Loci.isEmpty(l);
const sameOwner = (a: SeqLociAny, b: SeqLociAny) =>
  a.kind === b.kind &&
  a.set.lineage === b.set.lineage &&
  (!("alignment" in a) || a.alignment === (b as never as typeof a).alignment);

/** The locations of `a` that pass `keep`, as loci of its kind; edge-rate. */
function filter(a: SeqLociAny, keep: (l: Location) => boolean): SeqLociAny {
  const l = locationOf(a);
  if (a.kind === "column-loci") {
    const cols: number[] = [];
    Loci.forEachLocation(a, l, (x) => keep(x) && cols.push((x as ColumnLocation).column), null);
    return { ...a, columns: ofList(cols) } as SeqLociAny;
  }
  const [owner, values, value] = SHAPE[a.kind];
  const by = new Map<number, number[]>();
  Loci.forEachLocation(
    a,
    l,
    (x) => {
      if (!keep(x)) return;
      const r = x as unknown as Record<string, number>;
      const list = by.get(r[owner]) ?? [];
      by.set(r[owner], list);
      list.push(r[value]);
    },
    null,
  );
  const elements = [...by.keys()].sort((x, y) => x - y).map((k) => ({ [owner]: k, [values]: ofList(by.get(k) ?? []) }));
  return { ...a, elements } as unknown as SeqLociAny;
}
const contains = (a: SeqLociAny, b: SeqLociAny) => {
  let missing = 0;
  Loci.forEachLocation(b, locationOf(b), (x) => !Loci.has(a, x) && missing++, null);
  return missing === 0;
};

export class Selections implements SelectionManager {
  readonly current: WritableSeqState<Selection>;
  readonly pending: WritableLive<PendingSelection> = live({
    active: false,
    remove: false,
    section: -1,
    kind: 0,
    block: -1,
    startShown: -1,
    endShown: -1,
    rowStart: -1,
    rowEnd: -1,
    plotLayer: -1,
    yBlock: -1,
    yStartShown: -1,
    yEndShown: -1,
  });
  private version = 0;
  /** The edge-rate projector over the shown view (`Viewports.projector`, with `ctx.systems`). */
  private readonly projector: () => Projector | null;

  constructor(projector: () => Projector | null = () => null) {
    this.projector = projector;
    this.current = state<Selection>({ version: 0, origin: "api", ...NONE } as Selection);
  }

  set(target: SeqLociAny | Coordinates | null, origin: Origin) {
    const list = target ? this.lociOf(target) : [];
    if (!list.length) return this.clear(origin);
    const parts: Parts = { ...NONE };
    for (const l of list) parts[FIELD[l.kind]] = l;
    this.emit(parts, origin);
  }
  add(target: SeqLociAny | Coordinates, origin: Origin) {
    for (const loci of this.lociOf(target)) {
      const prev = this.part(loci);
      this.put(loci, prev && sameOwner(prev, loci) ? Loci.union(prev, loci as never) : loci, origin);
    }
  }
  remove(target: SeqLociAny | Coordinates, origin: Origin) {
    for (const loci of this.lociOf(target)) {
      const prev = this.part(loci);
      if (prev && sameOwner(prev, loci))
        this.put(
          loci,
          filter(prev, (x) => !Loci.has(loci, x)),
          origin,
        );
    }
  }
  toggle(target: SeqLociAny | Coordinates, origin: Origin) {
    for (const loci of this.lociOf(target)) {
      const prev = this.part(loci);
      if (prev && sameOwner(prev, loci) && contains(prev, loci)) this.remove(loci, origin);
      else this.add(loci, origin);
    }
  }
  close(loci: SeqLociAny, closure: Closure): SeqLociAny {
    return this.projector()?.closeLoci(loci, closure) ?? loci;
  }
  clear(origin: Origin) {
    if (FIELDS.some((f) => this.current.value[f] !== null)) this.emit({ ...NONE }, origin);
  }
  isWhole(loci: SeqLociAny): boolean {
    const prev = this.part(loci);
    return prev !== null && Loci.areEqual(prev, loci);
  }

  /** Widens to a granularity: an item's extents (range), whole units or rows, copy groups, columns. */
  widen(loci: SeqLociAny, granularity: number): SeqLociAny {
    const G = SelectionGranularity;
    const set = loci.set;
    const c = Conversion.create();
    const units = set.units;
    if (granularity === G.range && loci.kind === "item-loci") {
      const by = new Map<number, number[]>();
      for (const e of loci.elements) {
        const ann = set.annotations[e.annotation];
        for (let r = 0; ann && r < ann.count; r++) {
          if (!OrderedSet.has(e.items, ann.item[r]) || ann.target.alignment[r] >= 0) continue;
          for (let u = 0; u < units.count; u++) {
            const own =
              ann.target.unit[r] === u || (ann.target.unit[r] < 0 && units.molecule[u] === ann.target.molecule[r]);
            if (!own || !units.alive[u]) continue;
            const list = by.get(u) ?? [];
            by.set(u, list);
            for (let i = ann.start[r]; i < ann.stop[r]; i++) list.push(i);
          }
        }
      }
      const elements = [...by].sort((a, b) => a[0] - b[0]).map(([unit, xs]) => ({ unit, indices: ofList(xs) }));
      return { kind: "seq-loci", set, elements } as unknown as SeqLoci;
    }
    if (loci.kind === "seq-loci" && (granularity === G.unit || granularity === G.copyGroup)) {
      const by = new Map<number, OrderedSet>();
      for (const e of loci.elements) {
        const whole = OrderedSet.ofRange(units.rangeStart[e.unit], units.rangeStop[e.unit]);
        const g = units.group[e.unit];
        const copies = granularity === G.copyGroup && g >= 0 ? Array.from(set.groups.members[g]) : [e.unit];
        for (const u of copies)
          by.set(u, OrderedSet.union(by.get(u) ?? OrderedSet.Empty, granularity === G.unit ? whole : e.indices));
      }
      const elements = [...by].sort((a, b) => a[0] - b[0]).map(([unit, indices]) => ({ unit, indices }));
      return { ...loci, elements } as unknown as SeqLoci;
    }
    const aligned =
      loci.kind === "row-loci" ? RowLoci.toAlignedLoci(loci, c) : loci.kind === "aligned-loci" ? loci : null;
    if (aligned && granularity === G.column) return AlignedLoci.toColumnLoci(aligned, c);
    if (aligned && (granularity === G.row || granularity === G.unit)) {
      const a = set.alignments[aligned.alignment];
      const rows = new Set<number>();
      if (Array.isArray(aligned.elements)) for (const e of aligned.elements) rows.add(e.row);
      else OrderedSet.forEach((aligned.elements as { rows: OrderedSet }).rows, (r) => rows.add(r), null);
      const sorted = [...rows].sort((x, y) => x - y);
      const elements = sorted.map((row) => ({
        row,
        indices: OrderedSet.ofRange(0, a?.rows.lengthOf(row as never) ?? 0),
      }));
      return { kind: "row-loci", set, alignment: aligned.alignment, elements } as unknown as RowLoci;
    }
    return loci;
  }

  /**
   * On every assembled set: each part carried to it (`carry`), origin `graph`; null (a failed update) changes nothing,
   * so the next set is carried from the last one.
   */
  update(set: SeqSet | null) {
    if (!set) return;
    const cur = this.current.value;
    const next: Parts = { ...NONE };
    let changed = false;
    for (const f of FIELDS) {
      const l = cur[f] as SeqLociAny | null;
      if (!l) continue;
      next[f] = carry(l, set);
      changed ||= next[f] !== l;
    }
    if (changed) this.emit(next, "graph");
  }

  /** Loci as given; `Coordinates` through the projector's systems into loci per kind (nothing registered: none). */
  private lociOf(target: SeqLociAny | Coordinates): SeqLociAny[] {
    if (!isCoordinates(target)) return [target];
    return this.projector()?.lociOf(target) ?? [];
  }
  private part(loci: SeqLociAny): SeqLociAny | null {
    return this.current.value[FIELD[loci.kind]] as SeqLociAny | null;
  }
  private put(loci: SeqLociAny, value: SeqLociAny, origin: Origin) {
    this.emit({ ...(this.current.value as Parts), [FIELD[loci.kind]]: empty(value) ? null : value }, origin);
  }
  private emit(parts: Parts, origin: Origin) {
    this.current.set({ ...parts, version: ++this.version, origin } as Selection, origin);
  }

  dispose() {
    this.current.complete();
    this.pending.complete();
  }
}
