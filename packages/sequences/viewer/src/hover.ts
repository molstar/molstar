/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The pointer and the hover (seqspec:interaction#hover): the record the dispatcher writes, the hover the `hover`
// behavior, the labels and the channels (channels.ts) read, an external owner's target (`setLoci`: loci or
// `Coordinates`, carried to every new set), and `toLoci`, what is lit. Reused records; no move inside an element
// touches either. A header's hover carries what it stands for (`header`, header-hover.ts) as the `headerHover` behavior
// describes it, counted once when it lands; `toLoci` and `size` are getters, resolved when read, cached per version.

import type { Closure } from "@molstar/seq-core/format";
import {
  type AlignedLocation,
  type ColumnLocation,
  type CoordinateSystem,
  Extents,
  type ItemLocation,
  Loci,
  type RowLocation,
  type SeqLocation,
  type SeqSet,
} from "@molstar/seq-core/model";
import type { ViewModel } from "./graph.js";
import {
  describeHeader,
  type HeaderHoverParams,
  type HeaderSet,
  headerLoci,
  headerSet,
  type LociSize,
  lociSize,
} from "./header-hover.js";
import { Hit, HitKey, HitRegion, HoverKind } from "./hit.js";
import { carry, type SeqLociAny } from "./lineage.js";
import { type Coordinates, isCoordinates, type Projector } from "./projection.js";
import { live, type SeqState, type WritableLive } from "./state.js";

/** Read once: a member read in a hover-rate path allocates under the test runner's module wrapper. */
const HEADER = HoverKind.header;

/** `ctx.pointer.current`: written by the input dispatcher. `x` and `y` start as doubles (a Smi field turning Double). */
export interface PointerState {
  section: number;
  x: number;
  y: number;
  inside: boolean;
}

/** Reused. */
export interface HoverState {
  owner: string;
  section: number;
  /** `HitRegion`, chrome included. */
  region: number;
  /** `HoverKind`. */
  kind: number;
  /** The hit's block; -1: external. */
  block: number;
  /** The model index on the block's axis. */
  index: number;
  /** The element's coordinate system (seqspec:interaction#channels: the element). */
  system: CoordinateSystem;
  closure: Closure;
  seq: SeqLocation;
  row: RowLocation;
  aligned: AlignedLocation;
  column: ColumnLocation;
  item: ItemLocation;
  /** Resolved-view indices; -1: none. */
  track: number;
  layer: number;
  lane: number;
  treeNode: number;
  dataColumn: number;
  /** A matrix's state (band) under the pointer; -1: none. */
  state: number;
  /** A split track's row: its value's, -2 its label's entry; -1: none. */
  split: number;
  /** A group's header: the group; -1: none. */
  group: number;
  /** A header's description (the `headerHover` behavior's): its size, whether it lights, its runs; `on` false else. */
  readonly header: HeaderSet;
  /** A renderer's part of the element (a fingerprint's kind); -1: none. */
  part: number;
  /**
   * A dot plot's partner under y (seqspec:interaction#hover), lit with the element: its block (the hit's section's),
   * unit and index; -1: none.
   */
  pairBlock: number;
  pairUnit: number;
  pairIndex: number;
  /** A streamed item's [start, stop) on the hovered unit (the hit's `span`, ST4b); -1, -1: none. */
  readonly span: Int32Array;
  /** An external owner's target, as given. */
  external: SeqLociAny | Coordinates | null;
}
export interface HoverManager {
  readonly live: { readonly current: HoverState; readonly version: SeqState<number> };
  /** "pointer", "keyboard". `closure` defaults to "containing". */
  setFromHit(hit: Hit, owner: string, closure?: Closure): void;
  /** External ("3d", a Context row); `closure` defaults to "element"; null clears it if the owner holds the hover. */
  setLoci(target: SeqLociAny | Coordinates | null, owner: string, opts?: { closure?: Closure }): void;
  /** Only if the owner still holds the hover. */
  clear(owner: string): void;
  /**
   * What is lit: the element closed, a molecule residue in every alive unit; closed already. A header's: its shown items
   * as item loci when the behavior forwards them, else null. Resolved when read, cached per version.
   */
  toLoci(): SeqLociAny | null;
  /** The size of what `toLoci` gives (items, covered positions), known before it resolves; null: nothing. */
  size(): LociSize | null;
  /**
   * The `headerHover` behavior's hook: its params for a section (null: the section's headers light nothing). Without
   * it, a header's hover is the header alone.
   */
  setHeaderHover(params: ((section: number) => HeaderHoverParams | null) | null): void;
  /** The pointer's hover again, its version bumped: what the element is may have changed (a folded clade's label). */
  refresh(): void;
}

export class Hover implements HoverManager {
  readonly live: WritableLive<HoverState>;
  /** The hovered element's `HitKey`; -2: none (a leave, a clear). */
  private readonly key = new Int32Array(HitKey.SIZE).fill(-2);
  /** The edge-rate projector over the shown view (`Viewports.projector`), for `toLoci`. */
  private readonly projector: () => Projector | null;
  /** The view shown (a header's description reads it). */
  private readonly model: () => ViewModel | null;
  /** The `headerHover` behavior's params per section; null: no behavior. */
  private headerParams: ((section: number) => HeaderHoverParams | null) | null = null;
  /** The params the current header was described with. */
  private described: HeaderHoverParams | null = null;
  /** A header's track and group by ref (indices change with the view): read again on a new set. */
  private trackRef: string | null = null;
  private groupRef: string | null = null;
  /** `toLoci` and `size` of `version`. */
  /**
   * `toLoci` and `size` per hover version, and per what resolving reads besides the hover: the projector's view and set
   * (`epoch`), the host's systems (`systems`, their version: a `Coordinates` hover lights once its system registers).
   */
  private cache = {
    version: -1,
    loci: null as SeqLociAny | null,
    size: null as LociSize | null,
    sized: -1,
    epoch: -1,
    systems: -1,
  };
  private readonly systemsVersion: () => number;

  constructor(
    projector: () => Projector | null = () => null,
    model: () => ViewModel | null = () => null,
    systemsVersion: () => number = () => 0,
  ) {
    this.systemsVersion = systemsVersion;
    this.projector = projector;
    this.model = model;
    const { seq, row, aligned, column, itemLoc } = Hit.create();
    this.live = live<HoverState>({
      owner: "",
      section: -1,
      region: HitRegion.none,
      kind: HoverKind.none,
      block: -1,
      index: -1,
      system: { kind: -1, id: -1, row: -1 },
      closure: "containing",
      seq,
      row,
      aligned,
      column,
      item: itemLoc,
      track: -1,
      layer: -1,
      lane: -1,
      treeNode: -1,
      dataColumn: -1,
      state: -1,
      split: -1,
      group: -1,
      header: headerSet(),
      part: -1,
      pairBlock: -1,
      pairUnit: -1,
      pairIndex: -1,
      span: new Int32Array([-1, -1]),
      external: null,
    });
  }

  /** The same element: `HitKey`, the index, and the set its location is in (a new set is a new element). */
  private same(hit: Hit): boolean {
    const h = this.live.current;
    if (h.index !== hit.index || !HitKey.same(hit, this.key)) return false;
    switch (hit.kind) {
      case HoverKind.seq:
        return h.seq.set === hit.seq.set;
      case HoverKind.row:
        return h.row.set === hit.row.set;
      case HoverKind.aligned:
        return h.aligned.set === hit.aligned.set;
      case HoverKind.column:
        return h.column.set === hit.column.set;
      case HoverKind.item:
        return h.item.set === hit.itemLoc.set;
      default:
        return true;
    }
  }

  setFromHit(hit: Hit, owner: string, closure: Closure = "containing") {
    const h = this.live.current;
    if (h.owner === owner && h.external === null && h.closure === closure && this.same(hit)) return;
    h.owner = owner;
    h.external = null;
    h.closure = closure;
    h.section = hit.section;
    h.region = hit.region;
    h.kind = hit.kind;
    h.block = hit.block;
    h.index = hit.index;
    h.track = hit.track;
    h.layer = hit.layer;
    h.lane = hit.lane;
    h.treeNode = hit.treeNode;
    h.dataColumn = hit.dataColumn;
    h.state = hit.state;
    h.split = hit.split;
    h.group = hit.group;
    h.part = hit.part;
    h.pairBlock = hit.pairBlock;
    h.pairUnit = hit.pairUnit;
    h.pairIndex = hit.pairIndex;
    h.span[0] = hit.span[0];
    h.span[1] = hit.span[1];
    h.system.kind = hit.system.kind;
    h.system.id = hit.system.id;
    h.system.row = hit.system.row;
    // The kind's location, field by field (a table of names would be a keyed store per field).
    switch (hit.kind) {
      case HoverKind.seq:
        h.seq.set = hit.seq.set;
        h.seq.unit = hit.seq.unit;
        h.seq.index = hit.seq.index;
        // The item under a residue hover (`UI.OwnRowHover`), for its label and card; none: -1.
        h.item.set = hit.itemLoc.set;
        h.item.annotation = (hit.annotation >= 0 ? hit.itemLoc.annotation : -1) as never;
        h.item.item = (hit.annotation >= 0 ? hit.itemLoc.item : -1) as never;
        break;
      case HoverKind.row:
        h.row.set = hit.row.set;
        h.row.alignment = hit.row.alignment;
        h.row.row = hit.row.row;
        h.row.index = hit.row.index;
        break;
      case HoverKind.aligned:
        h.aligned.set = hit.aligned.set;
        h.aligned.alignment = hit.aligned.alignment;
        h.aligned.row = hit.aligned.row;
        h.aligned.column = hit.aligned.column;
        break;
      case HoverKind.column:
        h.column.set = hit.column.set;
        h.column.alignment = hit.column.alignment;
        h.column.column = hit.column.column;
        break;
      case HoverKind.item:
        h.item.set = hit.itemLoc.set;
        h.item.annotation = hit.itemLoc.annotation;
        h.item.item = hit.itemLoc.item;
        break;
    }
    // A header: what it stands for, once per header landed on (a move inside it returned above).
    this.describe(hit.kind === HEADER, true);
    HitKey.write(hit, this.key);
    this.live.bump();
  }

  setLoci(target: SeqLociAny | Coordinates | null, owner: string, opts: { closure?: Closure } = {}) {
    const h = this.live.current;
    const closure = opts.closure ?? "element";
    if (!target) return this.clear(owner);
    if (h.owner === owner && h.external === target && h.closure === closure) return;
    this.reset(owner, target);
    h.closure = closure;
    this.live.bump();
  }

  clear(owner: string) {
    if (this.live.current.owner !== owner) return;
    this.reset("", null);
    this.live.bump();
  }

  setHeaderHover(params: ((section: number) => HeaderHoverParams | null) | null) {
    this.headerParams = params;
    if (this.live.current.kind !== HEADER) return;
    this.describe(true, false);
    this.live.bump();
  }

  /**
   * The header's description under the behavior's params: on landing from the hit's indices (their refs kept), else
   * again from its refs in the view now shown (a filter changed, data arrived, a new document).
   */
  private describe(header: boolean, landing: boolean) {
    const h = this.live.current;
    const s = h.header;
    s.on = false;
    s.extents.count = 0;
    const params = header ? (this.headerParams?.(h.section) ?? null) : null;
    this.described = params;
    const model = this.model();
    if (!params || !model) return;
    const { resolved: rv, set } = model;
    if (landing) {
      s.track = h.track;
      s.split = h.split;
      s.lane = h.lane;
      s.group = h.group;
      this.trackRef = s.track >= 0 ? set.strings.get(rv.tracks.ref[s.track]) : null;
      this.groupRef = s.group >= 0 ? set.strings.get(rv.groups.ref[s.group]) : null;
    } else {
      const [track, group] = [this.trackRef, this.groupRef];
      s.track = track === null ? -1 : rv.tracks.ref.indexOf(set.strings.find(track));
      s.group = group === null ? -1 : rv.groups.ref.indexOf(set.strings.find(group));
      if ((track !== null && s.track < 0) || (group !== null && s.group < 0)) return;
    }
    describeHeader(model, params, s);
  }

  /** Cached per version: a header's counted when it landed (nothing resolved), any other's from its loci. */
  size(): LociSize | null {
    const h = this.live.current;
    const c = this.fresh();
    const v = this.live.version.value;
    if (c.sized !== v) {
      const s = h.header;
      c.size =
        h.kind !== HEADER || h.external
          ? lociSize(this.toLoci())
          : s.on && s.forward
            ? { items: s.items, length: s.length }
            : null;
      c.sized = v;
    }
    return c.size;
  }

  /** Resolved on the first read of a version, then kept (seqspec:interaction#hover: "once a version"). */
  toLoci(): SeqLociAny | null {
    const c = this.fresh();
    const v = this.live.version.value;
    if (c.version !== v) {
      c.version = v;
      c.loci = this.resolve();
    }
    return c.loci;
  }

  /** The cache, emptied when the projector's view or set or the host's systems changed since it was filled. */
  private fresh() {
    const c = this.cache;
    const epoch = this.projector()?.epoch ?? -1;
    const systems = this.systemsVersion();
    if (c.epoch !== epoch || c.systems !== systems) {
      c.epoch = epoch;
      c.systems = systems;
      c.version = c.sized = -1;
    }
    return c;
  }

  private resolve(): SeqLociAny | null {
    const h = this.live.current;
    if (h.kind === HEADER && !h.external)
      return h.header.on && h.header.forward && this.described
        ? headerLoci(this.model(), this.described, h.header)
        : null;
    const p = this.projector();
    const ext = h.external;
    if (ext && !isCoordinates(ext)) return p ? p.closeLoci(ext, h.closure) : ext;
    // Rows, cells, columns and items as their own loci: span-1 maps close nothing, an item's closure is none.
    const k = h.kind;
    const loc =
      k === HoverKind.row ? h.row : k === HoverKind.aligned ? h.aligned : k === HoverKind.column ? h.column : null;
    if (!ext && (loc || k === HoverKind.item)) return Loci.fromLocation(loc ?? h.item);
    const e = Extents.create(4);
    if (ext) p?.fromLoci(ext, e);
    else p?.fromHover(h, e);
    if (!p || e.count === 0) return null;
    if (h.closure === "containing") p.close(e);
    return p.toSeqLoci(e);
  }

  /**
   * On every assembled set (seqspec:interaction#selection): an external loci carried to it; nothing left of it clears
   * the hover.
   */
  update(set: SeqSet | null) {
    const h = this.live.current;
    // A header's description again in the view now shown, by its refs (what its track shows may have changed).
    if (set && h.kind === HEADER && !h.external) {
      this.describe(true, false);
      this.live.bump();
      return;
    }
    if (!set || !h.external || isCoordinates(h.external)) return;
    const next = carry(h.external, set);
    if (next === h.external) return;
    if (next) h.external = next;
    else this.reset("", null);
    this.live.bump();
  }

  private reset(owner: string, external: SeqLociAny | Coordinates | null) {
    const h = this.live.current;
    h.section = h.block = h.index = -1;
    h.track = h.layer = h.lane = h.treeNode = h.dataColumn = h.state = h.split = h.group = h.part = -1;
    h.pairBlock = h.pairUnit = h.pairIndex = -1;
    h.header.on = false;
    h.header.extents.count = 0;
    h.span[0] = h.span[1] = -1;
    h.region = HitRegion.none;
    h.kind = HoverKind.none;
    h.system.kind = h.system.id = h.system.row = -1;
    h.external = external;
    h.owner = owner;
    this.key.fill(-2);
  }

  refresh() {
    if (this.live.current.owner === "pointer") this.live.bump();
  }
}
