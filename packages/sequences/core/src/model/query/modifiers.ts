/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QR2a's modifiers and filters (seqspec:seqql#each): `expand`, `group`, `each` with `current`, `pick`, `first`, `ids`;
// and the item slot's properties, which `pick`, `items` and `partners` test. Modifiers that change runs drop the item
// tag; filters keep their groups as they are.

import { type Annotation, Annotations, ColumnRoles } from "../annotation.js";
import type { CompileContext, QueryFn } from "./compile.js";
import { F64, type QueryContext } from "./context.js";
import { addRuns, eachRun, SelectionBuilder, sameSystem } from "./kernels.js";
import { compileWhere } from "./properties.js";
import { itemExtent, rt, unionOf } from "./runtime.js";
import { SeqSelection } from "./selection.js";
import { circleOf, domainOf } from "./systems.js";

const dom = new Int32Array(2);
const query = (c: CompileContext, key = "0") => c.query(key) as QueryFn<SeqSelection>;

/**
 * Each run widened in its own system, clipped to its domain, or on a circle across its origin (CM-10); runs that meet
 * merge (the builder's).
 */
function expand(c: CompileContext): QueryFn<SeqSelection> {
  const q = query(c);
  const n = c.raw<number>("1") ?? 0;
  const [left, right] = [c.raw<number>("left") ?? n, c.raw<number>("right") ?? n];
  return (ctx) => {
    const [s, b, set] = [q(ctx), ctx.builder(), ctx.set];
    let circle = 0;
    for (let g = 0, last = -1; g < s.groupCount; g++) {
      b.begin();
      for (let k = s.groupStart[g]; k < s.groupStart[g + 1]; k++) {
        if (last < 0 || !sameSystem(s, last, k)) {
          domainOf(set, s.kind[k], s.id[k], s.row[k], dom);
          circle = circleOf(set, s.kind[k], s.id[k]);
        }
        last = k;
        let lo = s.start[k] - left;
        let hi = s.stop[k] + right;
        if (circle > 0) {
          if (hi - lo >= circle) {
            b.add(s.kind[k], s.id[k], s.row[k], 0, circle);
            continue;
          }
          if (lo < 0) {
            b.add(s.kind[k], s.id[k], s.row[k], lo + circle, circle);
            lo = 0;
          }
          if (hi > circle) {
            b.add(s.kind[k], s.id[k], s.row[k], 0, hi - circle);
            hi = circle;
          }
        }
        b.add(s.kind[k], s.id[k], s.row[k], lo < dom[0] ? dom[0] : lo, hi > dom[1] ? dom[1] : hi);
      }
      b.end();
    }
    return ctx.finish(b);
  };
}

/** Regrouped: none, one group (MolQL's union); run, each maximal run; system, one per system. Untagged. */
function group(c: CompileContext): QueryFn<SeqSelection> {
  const [q, by] = [query(c), c.raw<string>("by")];
  return (ctx) => {
    const u = unionOf(ctx, q(ctx));
    if (by === "run") return eachRun(u);
    if (by !== "system") return u;
    const b = ctx.builder();
    for (let k = 0; k < u.runCount; k++) {
      if (k > 0 && !sameSystem(u, k - 1, k)) b.end();
      if (k === 0 || !sameSystem(u, k - 1, k)) b.begin();
      b.add(u.kind[k], u.id[k], u.row[k], u.start[k], u.stop[k]);
    }
    if (u.runCount > 0) b.end();
    return ctx.finish(b);
  };
}

/** One group of `s` as a selection, tagged as it was: what `(current)` gives inside `each`. */
function groupOf(ctx: QueryContext, s: SeqSelection, g: number): SeqSelection {
  const b = ctx.builder();
  b.beginAs(s, g);
  addRuns(b, s, g);
  b.end();
  return ctx.finish(b);
}

/** The bind: E once per group of Q, `(current)` that group; the results' groups concatenated in order. */
function each(c: CompileContext): QueryFn<SeqSelection> {
  const [q, e] = [query(c), query(c, "1")];
  return (ctx) => {
    const [s, prev, out] = [q(ctx), ctx.current, ctx.builder()];
    for (let g = 0; g < s.groupCount; g++) {
      ctx.current = groupOf(ctx, s, g);
      const r = e(ctx);
      for (let h = 0; h < r.groupCount; h++) {
        out.beginAs(r, h);
        addRuns(out, r, h);
        out.end();
      }
    }
    ctx.current = prev;
    return ctx.finish(out);
  };
}

/**
 * The groups passing `:where`, which reads the group, its item (when tagged) and its system (when it has one); the
 * group tests' selections are evaluated once, before the groups are walked (seqspec:seqql#elements-and-groups).
 */
function pick(c: CompileContext): QueryFn<SeqSelection> {
  const [q, { where, prepare }] = [query(c), compileWhere(c)];
  return (ctx) => {
    const [s, b, slot, missed] = [q(ctx), ctx.builder(), ctx.slot, ctx.missed];
    for (let i = 0; s.groupCount > 0 && i < prepare.length; i++) prepare[i](ctx);
    for (let g = 0; g < s.groupCount; g++) {
      const k = s.groupStart[g];
      const one = sameSystem(s, k, s.groupStart[g + 1] - 1);
      slot.selection = s;
      slot.group = g;
      slot.annotation = s.groupAnnotation[g];
      slot.item = s.groupItem[g];
      slot.kind = one ? s.kind[k] : -1;
      slot.id = one ? s.id[k] : -1;
      slot.row = one ? s.row[k] : -1;
      if (where(ctx) !== true) continue;
      b.beginAs(s, g);
      addRuns(b, s, g);
      b.end();
    }
    slot.selection = null;
    slot.group = -1;
    if (b.groupCount === 0) ctx.tested(c.path, s.groupCount, "group", missed);
    return ctx.finish(b);
  };
}

/** The first n groups, as they are. */
function first(c: CompileContext): QueryFn<SeqSelection> {
  const [q, n] = [query(c), c.raw<number>("1") ?? 1];
  return (ctx) => {
    const [s, b] = [q(ctx), ctx.builder()];
    for (let g = 0; g < s.groupCount && g < n; g++) {
      b.beginAs(s, g);
      addRuns(b, s, g);
      b.end();
    }
    return ctx.finish(b);
  };
}

/** The item groups whose item's id is listed; an id no group has is a warning naming some that are there. */
function ids(c: CompileContext): QueryFn<SeqSelection> {
  const [q, list] = [query(c), (c.raw<(string | number)[]>("1") ?? []).map(String)];
  return (ctx) => {
    const [s, b, set] = [q(ctx), ctx.builder(), ctx.set];
    const keys = Int32Array.from(list, (x) => set.strings.find(x) ?? -1);
    const found = new Uint8Array(list.length);
    const seen = new Set<string>();
    for (let g = 0; g < s.groupCount; g++) {
      const a = set.annotations[s.groupAnnotation[g]];
      const key = a ? a.itemKey[s.groupItem[g]] : -1;
      if (key < 0) continue;
      const i = keys.indexOf(key);
      if (seen.size < 6) seen.add(set.strings.get(key));
      if (i < 0) continue;
      // Every copy of a listed id is found (CM-13: `indexOf` marked the first only, and the rest were said missing).
      for (let k = i; k < keys.length; k++) if (keys[k] === key) found[k] = 1;
      b.beginAs(s, g);
      addRuns(b, s, g);
      b.end();
    }
    const missing = list.filter((_, i) => !found[i]);
    if (missing.length && s.groupCount > 0) {
      const have = seen.size ? `its items' ids include ${[...seen].join(", ")}` : "its groups are no items with ids";
      const say = `no item has the id ${missing.map((x) => `"${x}"`).join(", ")}; ${have}; ids are the annotation's own: match labels with (pick Q :where (matches (label) …))`;
      ctx.warn("seqql.ids.unknown", `${c.path}/args/1`, say);
    }
    return ctx.finish(b);
  };
}

/** `(current)`: the group `each` runs on. */
const current =
  (c: CompileContext): QueryFn<SeqSelection> =>
  (ctx) =>
    ctx.current ?? ctx.fail("seqql.current", c.path, "(current) is only meaningful inside each");

// The item slot's properties: what an item's first row says, its annotation's, or its extent's size.

const itemOf = (ctx: QueryContext) => {
  const a = ctx.set.annotations[ctx.slot.annotation];
  return a && ctx.slot.item >= 0 ? a : null;
};
const firstRow = (a: Annotation, item: number) => {
  const rows = Annotations.itemRows(a);
  return rows.start[item] < rows.start[item + 1] ? rows.rows[rows.start[item]] : -1;
};
/** A column or field of the item's first row; undefined: missing. */
const column = (ctx: QueryContext, name: string) => {
  const a = itemOf(ctx);
  const r = a ? firstRow(a, ctx.slot.item) : -1;
  return a && r >= 0 ? Annotations.rowValue(a, ctx.set, r, name) : undefined;
};
const itemProp = (name: string, read: (ctx: QueryContext) => unknown) => rt(name, () => read);

/** Positions in the item's extent: its group's when the slot is that group, else built once here. */
function itemLength(ctx: QueryContext): number | undefined {
  // Per item: no destructuring, whose array is made per call where the call is not inlined (8.4 MB a run of 50,000).
  const a = itemOf(ctx);
  const slot = ctx.slot;
  if (!a) return undefined;
  const s = slot.selection;
  if (s && slot.group >= 0 && s.groupAnnotation[slot.group] === a.id && s.groupItem[slot.group] === slot.item)
    return SeqSelection.positionCount(s, slot.group);
  const b = scratch.reset(null);
  b.begin();
  itemExtent(b, ctx.set, a, slot.item);
  b.end();
  let n = 0;
  for (let k = 0; k < b.runCount; k++) n += b.stop[k] - b.start[k];
  return n;
}
const scratch = new SelectionBuilder();

/** `(field name)`: any column or field; a name the annotation has none of is an error that lists what it has. */
function itemField(c: CompileContext): QueryFn<unknown> {
  const name = String(c.raw("0"));
  return (ctx) => {
    const v = column(ctx, name);
    const a = itemOf(ctx);
    if (v !== undefined || !a || ColumnRoles[a.geometry][name] !== undefined) return v;
    if (a.fields.some((f) => ctx.set.strings.get(f.name) === name)) return v;
    const cols = [a.category && "category", a.rowLabel && "label", a.tooltip && "tooltip", a.sparse && "value"];
    const has = [...cols.filter((x) => !!x), ...a.fields.map((f) => ctx.set.strings.get(f.name))];
    const what = ctx.set.strings.get(a.key);
    ctx.fail("seqql.field.unknown", c.path, `${what} has no field "${name}"; its fields: ${has.join(", ") || "none"}`);
  };
}

export const MODIFIER_RUNTIMES = [
  rt("seq.generator.current", current),
  rt("seq.modifier.expand", expand),
  rt("seq.modifier.group", group),
  rt("seq.modifier.each", each),
  rt("seq.filter.pick", pick),
  rt("seq.filter.first", first),
  rt("seq.filter.ids", ids),
  itemProp("seq.item.id", (ctx) => {
    const a = itemOf(ctx);
    return a && a.itemKey[ctx.slot.item] >= 0 ? ctx.set.strings.get(a.itemKey[ctx.slot.item]) : undefined;
  }),
  itemProp("seq.item.label", (ctx) => column(ctx, "label")),
  itemProp("seq.item.category", (ctx) => column(ctx, "category")),
  {
    ...itemProp("seq.item.value", (ctx) => column(ctx, "value")),
    // seqspec:seqql#performance: a sparse value read inline for the relations (QR2b), never boxed
    number: () => (ctx: QueryContext) => {
      const a = itemOf(ctx);
      const r = a?.sparse ? firstRow(a, ctx.slot.item) : -1;
      const v = a?.sparse ? 0 : column(ctx, "value");
      F64[0] = a?.sparse && r >= 0 ? a.sparse[r] : typeof v === "number" && !a?.sparse ? v : Number.NaN;
    },
  },
  itemProp("seq.item.type", (ctx) => {
    const a = itemOf(ctx);
    return a && a.semanticType >= 0 ? ctx.set.strings.get(a.semanticType) : undefined;
  }),
  itemProp("seq.item.annotation", (ctx) => {
    const a = itemOf(ctx);
    return a ? ctx.set.strings.get(a.key) : undefined;
  }),
  itemProp("seq.item.length", itemLength),
  rt("seq.item.field", itemField),
];
