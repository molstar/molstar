/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Refs → loci (seqspec:interaction#focus): what commands, highlight channels and hosts name. An edge path: each call
// allocates. `items` names items by id; unknown and empty ids are dropped, and with none left the ref resolves to null.

import { type DocNode, indexDoc, type SeqViewSpecDoc, toLabel } from "@molstar/seq-core/format";
import {
  Annotations,
  CoordinateSystemKind,
  type ItemLoci,
  KeyKind,
  type KeyTarget,
  Loci,
  NO_LABEL,
  Numbering,
  OrderedSet,
  type PositionLabel,
  type SeqSet,
  Systems,
  type UnitId,
} from "@molstar/seq-core/model";
import type { SeqLociAny } from "./lineage.js";

/**
 * A ref narrowed as the format's range args do (seqspec:format#descriptions-and-invocations): indices, labels in a
 * numbering, or item ids.
 */
export interface RangeArgs {
  readonly start?: number;
  readonly stop?: number;
  readonly beg?: string;
  readonly end?: string;
  readonly numbering?: string;
  readonly items?: readonly string[];
}

export interface RefManager {
  /**
   * Null where it names nothing, or where range args fall outside its domain (index args, or labels not there or
   * resolving outside it; a complex's outside every member: never clipped, as a document's focus is refused): then
   * `why` hears why, as `[start, stop) is outside "<ref>" (lo..hi)`.
   */
  toLoci(ref: string, args?: RangeArgs, why?: (message: string) => void): SeqLociAny | null;
  /** The set's keys and the document's refs (a data view's included). */
  list(): readonly string[];
}
export interface RefHost {
  set(): SeqSet | null;
  doc(): SeqViewSpecDoc | null;
}

const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const la: PositionLabel = { num: NO_LABEL, suffix: -1 as never, name: -1 as never };
const lb: PositionLabel = { ...la };
const range = new Int32Array(2);
/** A layer's or block's refs to the data it shows. */
const DATA = ["unit", "sequence", "ligand", "complex", "alignment"];

export class Refs implements RefManager {
  private readonly host: RefHost;
  private nodes = { doc: null as SeqViewSpecDoc | null, map: new Map<string, DocNode>() };

  constructor(host: RefHost) {
    this.host = host;
  }

  /**
   * A unit or molecule to its positions, a complex to its units', a row ref (`R/<k>`) to its residues, an alignment to
   * its columns, an annotation or computed node to its items, a track or layer to its data's, a block to its unit;
   * `args` narrows it. Null when it names nothing in the set.
   */
  toLoci(ref: string, args: RangeArgs = {}, why?: (message: string) => void): SeqLociAny | null {
    const set = this.host.set();
    if (!set) return null;
    if (set.byName(ref, t)) {
      const outside = rangeOutside(set, t.kind, t.id, t.row, args);
      if (outside === null) {
        refused.message = "";
        const l = fromKey(set, t.kind, t.id, t.row, args);
        if (!l && refused.message) why?.(refused.message.replace(REF, ref));
        return l;
      }
      why?.(
        `[${args.start ?? outside[0]}, ${args.stop ?? outside[1]}) is outside "${ref}" (${outside[0]}..${outside[1]})`,
      );
      return null;
    }
    const n = this.node(ref);
    if (n?.kind === "block") {
      const k = DATA.find((p) => n.params?.[p] != null);
      return k ? this.toLoci(String(n.params?.[k]), args, why) : null;
    }
    if (n?.kind !== "track" && n?.kind !== "layer") return null;
    const refs: string[] = [];
    const walk = (x: DocNode) => {
      if (x.kind === "annotation" && x.ref) refs.push(x.ref);
      else if (x.kind === "layer") refs.push(...[x.params?.data ?? []].flat().map(String));
      for (const c of x.kind === "layer" && x.params?.data != null ? [] : (x.children ?? [])) walk(c);
    };
    walk(n);
    let out: SeqLociAny | null = null;
    for (const r of refs) {
      const l = this.toLoci(r, args, why);
      const same =
        out && l?.kind === out.kind && (!("alignment" in out) || out.alignment === (l as typeof out).alignment);
      if (l && (!out || same)) out = out ? (Loci.union(out, l as never) as SeqLociAny) : l;
    }
    return out;
  }

  list(): readonly string[] {
    const out = new Set<string>();
    const set = this.host.set();
    if (set) {
      for (const list of [set.molecules, set.complexes, set.alignments, set.mappings, set.annotations])
        for (const x of list) if (x) out.add(set.strings.get(x.key));
      for (let u = 0; u < set.units.count; u++) if (set.units.alive[u]) out.add(set.strings.get(set.units.key[u]));
    }
    this.node("");
    for (const r of this.nodes.map.keys()) out.add(r);
    return [...out];
  }

  private node(ref: string): DocNode | undefined {
    const doc = this.host.doc();
    if (this.nodes.doc !== doc) {
      // Core's index: every node by ref, the implied `root` and `data` included.
      const map = new Map<string, DocNode>();
      if (doc) for (const [ref, n] of indexDoc(doc).nodes) map.set(ref, n.node);
      this.nodes = { doc, map };
    }
    return this.nodes.map.get(ref);
  }
}

/** Why the last `narrow` refused its args (`REF` for the ref); "" when it did not. */
const refused = { message: "" };
const REF = "\u0000";

/**
 * Indices [lo, hi) narrowed by `start`/`stop`, or by `beg`/`end` labels in a numbering; null when outside. `cyclic` (a
 * unit over all of a cyclic molecule): `start > stop`, or `beg` after `end`, runs across the origin
 * (seqspec:format#positions). `strict` (the whole domain the ref names, as a unit's or a row's): a range past [lo, hi),
 * or labels not in the numbering, are refused (`refused`), never clipped, as a document's focus is; else (one piece of
 * a molecule's) the range is cut to the piece.
 */
function narrow(
  set: SeqSet,
  n: Numbering | null,
  lo: number,
  hi: number,
  a: RangeArgs,
  cyclic = false,
  strict = false,
): OrderedSet | null {
  let [s, e] = [a.start ?? lo, a.stop ?? hi];
  if (a.beg != null) {
    const labels = `${a.beg}..${a.end ?? a.beg}`;
    if (
      !n ||
      !toLabel(set, n, a.beg, la) ||
      !toLabel(set, n, a.end ?? a.beg, lb) ||
      !Numbering.rangeOf(n, la, lb, range, cyclic)
    ) {
      if (strict) refused.message = `the labels ${labels} are not positions of "${REF}"${n ? "" : " (no numbering)"}`;
      return null;
    }
    [s, e] = [range[0], range[1]];
  }
  const wraps = cyclic && s > e && s < hi && e >= lo;
  if (strict && !wraps && (s < lo || e > hi || e <= s)) {
    refused.message = `[${s}, ${e})${a.beg != null ? ` (labels ${a.beg}..${a.end ?? a.beg})` : ""} is outside "${REF}" (${lo}..${hi})`;
    return null;
  }
  if (wraps) return OrderedSet.union(OrderedSet.ofRange(s, hi), OrderedSet.ofRange(lo, e));
  [s, e] = [Math.max(s, lo), Math.min(e, hi)];
  return e > s ? OrderedSet.ofRange(s, e) : null;
}

/**
 * Index args (`start`, `stop`) outside the key's domain, as a document's focus checks them (a unit's range, a
 * molecule's or row's length, an alignment's columns; across a whole cyclic molecule's origin `start > stop`): the
 * domain [lo, hi); null when they fit (or are not index args).
 */
function rangeOutside(set: SeqSet, kind: number, id: number, row: number, a: RangeArgs): [number, number] | null {
  if ((a.start == null && a.stop == null) || a.beg != null) return null;
  const units = set.units;
  const mol = kind === KeyKind.unit ? units.molecule[id] : kind === KeyKind.molecule ? id : -1;
  const aln = kind === KeyKind.row || kind === KeyKind.alignment ? set.alignments[id] : undefined;
  if (mol < 0 && !aln) return null;
  const [lo, hi] = aln
    ? [0, kind === KeyKind.row ? aln.rows.lengthOf(row as never) : aln.columnCount]
    : kind === KeyKind.unit
      ? [units.rangeStart[id], units.rangeStop[id]]
      : [0, set.molecules[mol]?.length ?? 0];
  const cyclic = !aln && !!set.molecules[mol]?.cyclic && lo === 0 && hi === set.molecules[mol]?.length;
  const [s, e] = [a.start ?? lo, a.stop ?? hi];
  const wraps = cyclic && s > e && s < hi && e >= 0;
  return s >= lo && (s < e || wraps) && e <= hi ? null : [lo, hi];
}

function fromKey(set: SeqSet, kind: number, id: number, row: number, args: RangeArgs): SeqLociAny | null {
  const units = set.units;
  // A numbering by name (a namespace alone names the one under it), as a document's focus reads it.
  const numbering = (k: number, x: number, r: number) =>
    args.numbering ? (Systems.numberingOf(set, k, x, r, args.numbering) ?? null) : null;
  /** `strict`: each unit whole (a unit's, a complex's members: one that does not hold the range is left out). */
  const unitsOf = (us: readonly number[], strict: boolean): SeqLociAny | null => {
    const elements = us.flatMap((u) => {
      const n = numbering(CoordinateSystemKind.unit, u, -1);
      const mol = set.molecules[units.molecule[u]];
      const [lo, hi] = [units.rangeStart[u], units.rangeStop[u]];
      const whole = !!mol?.cyclic && lo === 0 && hi === mol.length;
      const indices = units.alive[u] ? narrow(set, n, lo, hi, args, whole, strict) : null;
      return indices ? [{ unit: u as UnitId, indices }] : [];
    });
    return elements.length ? ({ kind: "seq-loci", set, elements } as unknown as SeqLociAny) : null;
  };
  const ranged = args.start != null || args.stop != null || args.beg != null;
  switch (kind) {
    case KeyKind.unit:
      return unitsOf([id], true);
    case KeyKind.molecule:
      return unitsOf(
        Array.from({ length: units.count }, (_, u) => u).filter((u) => units.molecule[u] === id),
        false,
      );
    case KeyKind.complex: {
      // Range args on a complex: the members that hold the range; none: refused, as a unit's would be.
      const l = unitsOf(Array.from(set.complexes[id]?.units ?? []), ranged);
      if (!l && ranged) refused.message = `no member of "${REF}" holds the range`;
      return l;
    }
    case KeyKind.row: {
      const a = set.alignments[id];
      const n = a ? numbering(CoordinateSystemKind.row, id, row) : null;
      const indices = a ? narrow(set, n, 0, a.rows.lengthOf(row as never), args, false, true) : null;
      return (
        indices && ({ kind: "row-loci", set, alignment: id, elements: [{ row, indices }] } as unknown as SeqLociAny)
      );
    }
    case KeyKind.alignment: {
      const columns = narrow(set, null, 0, set.alignments[id]?.columnCount ?? 0, args, false, true);
      return columns && ({ kind: "column-loci", set, alignment: id, columns } as unknown as SeqLociAny);
    }
    case KeyKind.annotation: {
      const ann = set.annotations[id];
      if (!ann) return null;
      let items = OrderedSet.ofRange(0, ann.itemCount);
      if (args.items) {
        const picked = new Set<number>();
        for (const id of args.items) {
          const item = id ? Annotations.itemOf(ann, set.strings.find(id) as never) : -1;
          if (item >= 0) picked.add(item);
        }
        if (picked.size === 0) return null;
        items = OrderedSet.ofSortedArray(Int32Array.from(picked).sort());
      }
      return { kind: "item-loci", set, elements: [{ annotation: id as never, items }] } as ItemLoci;
    }
    default:
      return null;
  }
}
