/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Shared scales (seqspec:format#themes, `domain_group`; EN1): the layers of a view that name one group share every
// normalized domain they resolve, the union per (group, field), so one number takes one color (or one height) in all of
// them.
// Each member's `resolveLayer` asks for its group's domains, so a viewer that resolves layers one at a time (the
// graph's layer nodes, which list the members' data as their inputs: `groupMembers`) has them in its themes; the
// union is memoized per set and group. `value` goes to the layer's `domain` (what a `values` y scale and a palette
// over the value read), a palette over another field has its `null` ends filled, and `composeView` collects the
// table (`ResolvedView.domains`). An explicit `value_domain` end wins in its layer; so does a layer's own `domain`:
// when any member gives one, the group's value domain is the union of those, for all of its members. Layers of a
// template track (instances are not nodes) are left out.

import type { Annotation, SeqSet } from "@molstar/seq-core/model";
import { KeyKind, type KeyTarget } from "@molstar/seq-core/model";
import type { ThemeRef } from "../doc.js";
import { type DocIndex, type StreamedData, StreamProp } from "../load/compile.js";
import { numberColumn } from "../palettes/fields.js";
import type { SvsRegistries } from "../registry.js";
import type { IndexedNode } from "../schema/tree-schema.js";
import { childrenOf, refsOf } from "./section.js";
import { forEachOf } from "./template.js";

/** One shared domain of a layer's group, for the resolved table: a field and its union. */
export interface SharedField {
  readonly field: string;
  readonly min: number;
  readonly max: number;
}

/** The range of the value over annotations (a streamed one's manifest domain too); lo > hi: none. */
export function dataDomain(anns: readonly Annotation[]): [number, number] {
  let [lo, hi] = [Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY];
  const see = (v: number) => {
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  };
  for (const a of anns) {
    (a.custom.props.get(StreamProp) as StreamedData | undefined)?.manifest.domain?.forEach(see);
    for (const d of a.dense ?? []) d.forEach(see);
    for (let r = 0; a.sparse && r < a.count; r++) see(a.sparse[r]);
    a.matrix?.values.forEach(see);
    for (const row of a.grid?.rows ?? []) row.forEach(see);
  }
  return [lo, hi];
}

const target: KeyTarget = { kind: -1, id: -1, row: -1 };
/** A layer node's data: its `data` refs, else its annotation children. */
export const dataRefs = (doc: DocIndex, n: IndexedNode): string[] =>
  n.params.data != null
    ? refsOf(n.params.data)
    : childrenOf(doc, n)
        .filter((c) => c.kind === "annotation")
        .map((c) => c.ref);
const annotationsOf = (doc: DocIndex, set: SeqSet, n: IndexedNode): Annotation[] =>
  dataRefs(doc, n).flatMap((r) =>
    set.byName(r, target) && target.kind === KeyKind.annotation ? [set.annotations[target.id] as Annotation] : [],
  );

const viewOf = (doc: DocIndex, n: IndexedNode) => {
  let x: IndexedNode | undefined = n;
  while (x && x.kind !== "view") x = x.parent === null ? undefined : doc.nodes.get(x.parent);
  return x;
};
const grouped = new WeakMap<DocIndex, Map<string, IndexedNode[]>>();
/** Every `domain_group`'s layers of a document, by view and group (template layers left out): built once per document. */
function groupsOf(doc: DocIndex) {
  let out = grouped.get(doc);
  if (out) return out;
  out = new Map();
  grouped.set(doc, out);
  for (const n of doc.order) {
    const g = n.kind === "layer" ? n.params.domain_group : null;
    const view = typeof g === "string" && g !== "" ? viewOf(doc, n) : undefined;
    if (!view) continue;
    let template = false;
    for (
      let x: IndexedNode | undefined = n;
      x && x !== view;
      x = x.parent === null ? undefined : doc.nodes.get(x.parent)
    )
      template ||= forEachOf(x) !== null;
    if (template) continue;
    const key = `${view.ref}\0${g}`;
    out.set(key, [...(out.get(key) ?? []), n]);
  }
  return out;
}
/** The layers of `n`'s view in its `domain_group`, itself included, in document order; none when it is alone in it. */
export function groupMembers(doc: DocIndex, n: IndexedNode): readonly IndexedNode[] {
  const g = n.params.domain_group;
  const view = typeof g === "string" && g !== "" ? viewOf(doc, n) : undefined;
  const list = view ? groupsOf(doc).get(`${view.ref}\0${g}`) : undefined;
  return list && list.length > 1 && list.includes(n) ? list : [];
}

type Union = [number, number];
interface Group {
  readonly value: Union;
  readonly explicit: boolean;
  readonly fields: Map<string, Union>;
}
const unions = new WeakMap<readonly IndexedNode[], WeakMap<SeqSet, Group>>();
function groupOf(doc: DocIndex, set: SeqSet, members: readonly IndexedNode[]): Group {
  const bySet = unions.get(members) ?? new WeakMap<SeqSet, Group>();
  unions.set(members, bySet);
  let g = bySet.get(set);
  if (g) return g;
  // The value: the explicit `domain`s when any member has one, else the data's; a layer without data adds nothing.
  const explicit = members.filter((m) => Array.isArray(m.params.domain));
  const ranges = (explicit.length > 0 ? explicit : members).map((m) =>
    Array.isArray(m.params.domain) ? (m.params.domain as Union) : dataDomain(annotationsOf(doc, set, m)),
  );
  const value: Union = [Math.min(...ranges.map((r) => r[0])), Math.max(...ranges.map((r) => r[1]))];
  g = { value, explicit: explicit.length > 0, fields: new Map() };
  bySet.set(set, g);
  return g;
}

/** What `resolveLayer` takes from its group: the value domain (lo > hi: none), whether it is explicit, the rest below. */
export interface GroupDomains {
  readonly lo: number;
  readonly hi: number;
  readonly explicit: boolean;
  /** The union of another field over the members' data (lo > hi: none), read through the registries' field providers. */
  field(name: string): readonly [number, number];
}
export function sharedDomains(
  doc: DocIndex,
  set: SeqSet,
  members: readonly IndexedNode[],
  reg: SvsRegistries,
): GroupDomains | null {
  if (members.length === 0) return null;
  const g = groupOf(doc, set, members);
  return {
    lo: g.value[0],
    hi: g.value[1],
    explicit: g.explicit,
    field(name) {
      let held = g.fields.get(name);
      if (!held) {
        const found: Union = [Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY];
        for (const m of members)
          for (const a of annotationsOf(doc, set, m)) {
            const column = numberColumn(a, set, name, reg);
            for (let r = 0; column && r < a.count; r++) {
              if (column[r] < found[0]) found[0] = column[r];
              if (column[r] > found[1]) found[1] = column[r];
            }
          }
        g.fields.set(name, found);
        held = found;
      }
      return held;
    },
  };
}

type PaletteParams = {
  palette?: { kind?: string; mode?: string; value_domain?: readonly (number | null)[] };
  field?: unknown;
  rules?: { theme?: ThemeRef | null }[];
  default?: ThemeRef | null;
  domain?: unknown;
  map?: unknown;
};
/**
 * A theme with the `null` ends of its normalized palettes over another field than the value filled from the group's
 * union of that field (a `field_palette`, also inside `rules`); each field it used is added to `used`.
 */
export function fillShared(r: ThemeRef | null, g: GroupDomains, used: SharedField[], depth = 0): ThemeRef | null {
  const q = r?.params as PaletteParams | undefined;
  if (!r || !q || depth > 4) return r;
  if (r.theme === "field_palette" && typeof q.field === "string" && q.field !== "value") {
    const given = q.palette?.kind === "categorical" ? null : q.palette?.value_domain;
    if (!given || q.palette?.mode === "absolute" || !given.includes(null)) return r;
    const [min, max] = g.field(q.field);
    if (!(min <= max)) return r;
    if (!used.some((u) => u.field === q.field)) used.push({ field: q.field, min, max });
    return { ...r, params: { ...q, palette: { ...q.palette, value_domain: [given[0] ?? min, given[1] ?? max] } } };
  }
  // A `field_scale` over another field with no domain of its own takes the group's.
  if (r.theme === "field_scale" && typeof q.field === "string" && q.field !== "value" && !Array.isArray(q.domain)) {
    const [min, max] = g.field(q.field);
    if (!(min <= max) || q.map != null) return r;
    if (!used.some((u) => u.field === q.field)) used.push({ field: q.field, min, max });
    return { ...r, params: { ...q, domain: [min, max] } };
  }
  if (r.theme !== "rules") return r;
  const rules = (q.rules ?? []).map((x) => (x.theme ? { ...x, theme: fillShared(x.theme, g, used, depth + 1) } : x));
  return { ...r, params: { ...q, rules, default: fillShared(q.default ?? null, g, used, depth + 1) } };
}
