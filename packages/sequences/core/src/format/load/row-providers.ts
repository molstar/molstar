/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The kernel providers `row_stats` and `redundancy` (seqspec:format#computed-nodes, seqspec:model#providers,
// seqspec:format#keyed-sets): key tables over an alignment's rows, `rows` annotations keyed by `row`. The statistics
// are `@molstar/seq-core/algorithms`' (pure, over the alignment's cells); this file reads the set, counts the items of an
// annotation on each row through the mapping paths, and builds the annotation. `redundancy` yields to the event loop
// and honors the stage's abort signal, and its result is cached by the rows' content and its options, as the node
// providers' are.

import {
  type Clusters,
  redundancy as cluster,
  RedundancySizeError,
  type RowCodes,
  rowStats,
} from "@molstar/seq-core/algorithms";
import {
  type Alignment,
  type Annotation,
  type AnnotationInput,
  Annotations,
  CoordinateSystemKind,
  Issue,
  KeyKind,
  type KeyTarget,
  Paths,
  type RowIndex,
  type SeqSet,
} from "@molstar/seq-core/model";
import type { StageContext } from "../contract.js";
import type { ProviderEntry } from "../registry.js";
import { alignmentRows, keyed, type OrderKeyValue } from "../resolve/keyed.js";
import { endsOf, pathTo, steps } from "../resolve/reach.js";
import type { UnionParamsSchema } from "../schema/params-schema.js";
import { SvsTreeSchema } from "../tables/svs-tree.js";
import { report } from "./issues.js";
import { cached, contentKey } from "./node-providers.js";

const cases = (SvsTreeSchema.nodes.computed.params as UnionParamsSchema).cases;
const kind: KeyTarget = { kind: -1, id: -1, row: -1 };
const latin1 = new TextDecoder("latin1");

/** The alignment `ref` names, and its cells as the algorithms read them; null with the reason reported. */
function cellsOf(ctx: StageContext, set: SeqSet, ref: string, what: string): { a: Alignment; cells: RowCodes } | null {
  const fail = (code: string, message: string) => {
    ctx.report(Issue.error(code, "/params/alignment", `${what}: ${message}`));
    return null;
  };
  const a = set.byName(ref, kind) && kind.kind === KeyKind.alignment ? set.alignments[kind.id] : undefined;
  if (!a) return fail("model.provider.input", `no alignment "${ref}".`);
  if (a.lazy) return fail("model.provider.lazy", "a progressive alignment's rows are not all held.");
  // The held rows are in row order; a different slot order is copied into it.
  const C = a.columnCount;
  let codes = a.rows.codes;
  if (
    a.rows.count !== a.rowCount ||
    Array.from({ length: a.rowCount }, (_, r) => a.rows.slotOf(r as RowIndex)).some((s, r) => s !== r)
  ) {
    codes = new Uint8Array(a.rowCount * C).fill(45);
    for (let r = 0; r < a.rowCount; r++) {
      const s = a.rows.slotOf(r as RowIndex);
      if (s >= 0) codes.set(a.rows.codes.subarray(s * C, (s + 1) * C), r * C);
    }
  }
  return { a, cells: { codes, rows: a.rowCount, columns: C, gap: 45 } };
}

interface ItemOptions {
  items: string;
  category: string | null;
  columns: [number, number] | null;
  score: string | null;
}

/**
 * Per row, the number of items of `ann` on it and the largest `score` among them. An item is on a row when one of its rows
 * covers a position of the row (through a mapping path that binds the row, or directly) whose column
 * lies in `columns`; on the alignment's columns, a column where the row has a letter.
 */
function itemsOn(set: SeqSet, a: Alignment, cells: RowCodes, ann: Annotation, o: ItemOptions, signal: AbortSignal) {
  const count = new Float64Array(a.rowCount);
  const score = new Float64Array(a.rowCount).fill(Number.NaN);
  const category = o.category === null ? -1 : (set.strings.find(o.category) ?? -2);
  const [lo, hi] = o.columns ?? [0, a.columnCount];
  const seen = new Int32Array(ann.itemCount).fill(-1);
  const ends = endsOf(ann).filter((e) => e.side === 0);
  for (let r = 0; r < a.rowCount; r++) {
    signal.throwIfAborted();
    const row = { kind: CoordinateSystemKind.row, id: a.id, row: r };
    for (const end of ends) {
      const own = end.system.kind === CoordinateSystemKind.row;
      if (own && (end.system.id !== a.id || end.system.row !== r)) continue;
      const wide = end.system.kind === CoordinateSystemKind.columns;
      if (wide && end.system.id !== a.id) continue;
      const path = own || wide ? null : pathTo(set, end.system, row);
      if (path === undefined || path === "ambiguous") continue;
      // Only the rows a mapping binds: a path through the alignment's columns to another row is not on that row.
      if (path && steps.some((st) => st.mapping < 0 && st.unit < 0)) continue;
      const through = path?.heldSpans ? path : null;
      for (const i of end.rows) {
        if (category !== -1 && ann.category?.[i] !== category) continue;
        const item = ann.item[i];
        for (let p = ann.start[i]; p < Math.max(ann.stop[i], ann.start[i] + 1); p++) {
          // The position's image: a column (wide), else the row's position through the path, then its column.
          let on = false;
          const test = (q: number) => {
            const c = wide ? q : q < 0 ? -1 : a.rows.colOf(r as RowIndex, q);
            if (c >= lo && c < hi && c < a.columnCount && (!wide || cells.codes[r * a.columnCount + c] !== cells.gap))
              on = true;
          };
          if (through) for (let k = 0; k < 3; k++) test(Paths.spanAt(through, p, k));
          else test(path ? Paths.at(path, p) : p);
          if (!on) continue;
          if (seen[item] !== r) {
            seen[item] = r;
            count[r]++;
          }
          if (o.score !== null) {
            const x = Annotations.rowValue(ann, set, i, o.score);
            if (typeof x === "number" && !(score[r] >= x)) score[r] = x;
          }
          break;
        }
      }
    }
  }
  return { count, score };
}

/** `row_stats`: per row of an alignment its identity to a reference row, coverage, length, and the items on it. */
const rowStatsProvider: ProviderEntry = {
  name: "row_stats",
  geometry: "rows",
  options: cases.row_stats,
  async compute(ctx, options) {
    const o = options as { alignment: string; reference_row: number | null } & ItemOptions;
    if (!ctx.set) return { key: "", geometry: "rows", columns: {} };
    const found = cellsOf(ctx, ctx.set, o.alignment, "row_stats");
    if (!found) return null;
    const { a, cells } = found;
    const ref = o.reference_row ?? a.referenceRow;
    if (!Number.isInteger(ref) || ref < 0 || ref >= a.rowCount) {
      ctx.report(Issue.error("model.provider.input", "/params/reference_row", `row_stats: no reference row ${ref}.`));
      return null;
    }
    const stats = rowStats(cells, ref);
    const length = Float64Array.from(stats.length, (n, r) =>
      a.rows.lengthOf(r as RowIndex) >= 0 ? a.rows.lengthOf(r as RowIndex) : n,
    );
    const columns: Record<string, ArrayLike<unknown>> = {
      row: Int32Array.from({ length: a.rowCount }, (_, r) => r),
      identity: stats.identity,
      coverage: stats.coverage,
      length,
    };
    if (o.items !== null) {
      const ok = ctx.set.byName(o.items, kind) && kind.kind === KeyKind.annotation;
      const ann = ok ? ctx.set.annotations[kind.id] : undefined;
      if (!ann) {
        ctx.report(Issue.error("model.provider.input", "/params/items", `row_stats: no annotation "${o.items}".`));
        return null;
      }
      const r = itemsOn(ctx.set, a, cells, ann, o, ctx.signal);
      columns.items = r.count;
      if (o.score !== null) columns.score = r.score;
    }
    const key = a.key;
    return { key: "", geometry: "rows", target: ctx.set.strings.get(key), columns } as AnnotationInput;
  },
};

/** `redundancy`: greedy clustering of an alignment's rows at an identity threshold, in the given order. */
const redundancyProvider: ProviderEntry = {
  name: "redundancy",
  geometry: "rows",
  options: cases.redundancy,
  async compute(ctx, options) {
    const o = options as { alignment: string; threshold: number; order: OrderKeyValue | OrderKeyValue[] | null };
    if (!ctx.set) return { key: "", geometry: "rows", columns: {} };
    const found = cellsOf(ctx, ctx.set, o.alignment, "redundancy");
    if (!found) return null;
    const { a, cells } = found;
    // The rows in the order asked for (a key table's fields, the tree's or data order): the same resolver a rows layer uses.
    const members = alignmentRows(a.id, a.rowCount, a.tree ? a.tree.leafOrder : null);
    const order = keyed(ctx.set, ctx.registries, members, { order: o.order ?? "data" }, (where, message) =>
      report(ctx, "svs.layer.data", `/params/${where}`, message),
    ).permutation;
    const rows: string[] = [];
    for (let r = 0; r < cells.rows; r++)
      rows.push(latin1.decode(cells.codes.subarray(r * cells.columns, (r + 1) * cells.columns)));
    const key = contentKey("redundancy", { threshold: o.threshold }, rows, [String(order.join(","))]);
    let result: Clusters | null;
    try {
      result = await cached<Clusters>(key, () => cluster(cells, o.threshold, order, ctx.signal));
    } catch (e) {
      if (!(e instanceof RedundancySizeError)) throw e;
      ctx.report(Issue.error("model.provider.size", "/params/alignment", `redundancy: ${e.message}`));
      return null;
    }
    if (!result) return null;
    const columns = {
      row: Int32Array.from({ length: a.rowCount }, (_, r) => r),
      cluster: result.cluster,
      representative: Array.from(result.representative, (x) => x === 1),
    };
    return { key: "", geometry: "rows", target: ctx.set.strings.get(a.key), columns } as AnnotationInput;
  },
};

export const RowProviders: readonly ProviderEntry[] = [rowStatsProvider, redundancyProvider];
