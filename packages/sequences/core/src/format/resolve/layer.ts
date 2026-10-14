/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `resolveLayer` (seqspec:format#placement, seqspec:format#layer-types, seqspec:format#themes): the type drawn after
// fallbacks, a theme per channel after defaults and fallbacks with named rule tables expanded, the layer's data and the
// rows its filter keeps, its domain over all its rows, and its placements onto the section's blocks and lanes. It reads
// only the section's axis (`axisKey`).

import {
  type Alphabet,
  type Annotation,
  Annotations,
  type ComposedPath,
  type CoordinateSystem,
  CoordinateSystemKind,
  type Issue,
  KeyKind,
  type KeyTarget,
  type NameId,
  type SeqSet,
} from "@molstar/seq-core/model";
import type { ThemeRef } from "../doc.ts";
import type { DocIndex } from "../load/compile.ts";
import { issueOf } from "../load/issues.ts";
import { streamedSources } from "../load/nodes.ts";
import { paletteValue } from "../palettes/fields.ts";
import { type SchemeSubject, schemeFits, schemeOf, schemeSubject } from "../palettes/schemes.ts";
import type { ChannelKind, LayerTypeEntry, SvsRegistries } from "../registry.ts";
import { isNamespaced, isPlainObject, TEMPLATE_NODE } from "../schema/field-schema.ts";
import type { IndexedNode } from "../schema/tree-schema.ts";
import { numericDefaultColor } from "../tables/kernel.ts";
import { themeRef } from "../tables/param-types.ts";
import { isStreamedSource } from "../validate.ts";
import { Layout, SectionAxis, Sizing } from "./codes.ts";
import { dataDomain, fillShared, groupMembers, type SharedField, sharedDomains } from "./domains.ts";
import { type KeyFilterValue, keyed, type Members, type OrderKeyValue } from "./keyed.ts";
import { endsOf, matcher, reach } from "./reach.ts";
import { childrenOf, paramsOf, type ResolvedSection, refsOf } from "./section.ts";
import { instanceRef } from "./template.ts";

/** One (layer, block, annotation, slot, member, lane) with its path; indices local to the section and layer. */
export interface PlacementRow {
  readonly block: number;
  readonly annotation: number;
  readonly slot: number;
  readonly member: number;
  readonly lane: number;
  readonly path: number;
}
/** One layer's rows of the resolved tables. */
export interface ResolvedLayer {
  readonly ref: string;
  /** Its track when resolved; `composeView` takes the section's (`ResolvedSection.layers`), so reuse stays right. */
  readonly track: number;
  /** NameId of the type drawn; -1: omitted. */
  readonly type: number;
  readonly declared: number;
  readonly sizing: number;
  /** Per channel of its type, in order (`color` first): after defaults and fallbacks, expanded; null: the viewer's. */
  readonly themes: readonly (ThemeRef | null)[];
  readonly domainMin: number;
  readonly domainMax: number;
  /** The domain is given, by the layer's `domain` or its group's: a `values` y scale uses it as it is, not down to 0. */
  readonly domainGiven: boolean;
  /** Its `domain_group` when it shares its domains with other layers of the view (domains.ts), else null. */
  readonly group: string | null;
  /** The group's other fields' domains its themes used (the value's are `domainMin/Max`). */
  readonly shared: readonly SharedField[];
  /** AnnotationIds. */
  readonly data: readonly number[];
  /**
   * A `matrix`'s states drawn, top to bottom, as indices into its data's states: the kept ones in their order
   * (seqspec:format#keyed-sets); null: every state, in data order.
   */
  readonly members: Int32Array | null;
  readonly kept: readonly (Uint8Array | null)[];
  /** Sorted by block, lane. */
  readonly placements: readonly PlacementRow[];
  readonly paths: readonly ComposedPath[];
  readonly params: Readonly<Record<string, unknown>>;
  /**
   * Its data names only dynamic nodes not held yet (seqspec:sources#dynamic-nodes): its track is not shown until one
   * arrives.
   */
  readonly pending: boolean;
  /** A template instance's node (what `"$node"` names); null: a document layer. */
  readonly node: string | null;
}

/** A template layer's params for one node: `"$node"` in `data` and `of` replaced by its ref. */
function withNode(p: Readonly<Record<string, unknown>>, node: string): Readonly<Record<string, unknown>> {
  const sub = (v: unknown) =>
    v === TEMPLATE_NODE ? node : Array.isArray(v) ? v.map((x) => (x === TEMPLATE_NODE ? node : x)) : v;
  return { ...p, data: sub(p.data), of: sub(p.of) };
}

/** Every registered layer type and theme, and the `char` and `path` glyphs: a viewer that draws everything. */
export function supportedAll(reg: SvsRegistries): Set<string> {
  const kinds = ["layerTypes", "themes"] as const;
  return new Set([...kinds.flatMap((k) => reg.list(k).map((e) => e.name)), "char", "path"]);
}

const t: KeyTarget = { kind: -1, id: -1, row: -1 };

/**
 * A split layer's rows: its values, `order`'s then the others in data order, and per annotation each kept row's (-1:
 * not kept). Values are strings; a missing or empty one is `-`. With `empty: "hide"`, only the values
 * a kept row has; `by` may be a field provider's (`reg`). Both viewers read rows through this.
 */
export function splitRows(
  set: SeqSet,
  anns: readonly Annotation[],
  kept: readonly (Uint8Array | null)[],
  split: { by: string; order?: readonly string[] | null; empty?: "show" | "hide" },
  reg?: SvsRegistries,
): { values: string[]; rowOf: Int32Array[] } {
  const values = [...(split.order ?? [])];
  const at = new Map(values.map((v, i) => [v, i]));
  const rowOf = anns.map((a, d) =>
    Int32Array.from({ length: a.count }, (_, r) => {
      if (kept[d] && !kept[d]?.[r]) return -1;
      const x = paletteValue(a, set, r, split.by, reg);
      const v = x == null || x === "" ? "-" : String(x);
      if (!at.has(v)) at.set(v, values.push(v) - 1);
      return at.get(v) as number;
    }),
  );
  if (split.empty !== "hide") return { values, rowOf };
  // Renumber the values some row has, in their order.
  const to = new Int32Array(values.length).fill(-1);
  for (const of of rowOf) for (const v of of) if (v >= 0) to[v] = 0;
  const shown: string[] = [];
  for (let i = 0; i < values.length; i++) if (to[i] === 0) to[i] = shown.push(values[i]) - 1;
  for (const of of rowOf) for (let r = 0; r < of.length; r++) if (of[r] >= 0) of[r] = to[of[r]];
  return { values: shown, rowOf };
}

/**
 * The `letter_map` param a layer draws through: its own; an `emphasis`'s, that of the `letters` (without `of`) or
 * `rows` layer before it in its track, whose letters it redraws.
 */
export function letterMapParam(doc: DocIndex, n: IndexedNode): unknown {
  if (n.params.type !== "emphasis") return n.params.letter_map;
  const track = n.parent === null ? undefined : doc.nodes.get(n.parent);
  const siblings = track ? childrenOf(doc, track) : [];
  let m: unknown;
  for (const x of siblings.slice(0, siblings.indexOf(n)))
    if ((x.params.type === "letters" && x.params.of == null) || x.params.type === "rows") m = x.params.letter_map;
  return m;
}

/**
 * The `letter_map` as its table, so a viewer reads one thing: a named map as its dict; none for a name no registry
 * holds (static validation reports it) or a map whose `alphabets` miss one the layer's letters are of (`alphabets`;
 * a ligand's null skipped).
 */
function letterMapped(
  params: Readonly<Record<string, unknown>>,
  m: unknown,
  reg: SvsRegistries,
  alphabets: () => readonly (Alphabet | null)[],
) {
  if (m == null) return params;
  const e = typeof m === "string" ? reg.get("letterMaps", m) : undefined;
  const fits = !e?.alphabets || alphabets().every((a) => a === null || e.alphabets?.includes(a));
  return { ...params, letter_map: typeof m === "string" ? (fits && e ? e.map : null) : m };
}

/**
 * One layer of a resolved section. `supported`: the drawing viewer's types, themes and glyphs. `node`: a template
 * layer's instance, `"$node"` naming it; its ref is `instanceRef(layer, node)`, and it is never unplaced.
 */
export function resolveLayer(
  section: ResolvedSection,
  layer: string,
  doc: DocIndex,
  set: SeqSet,
  reg: SvsRegistries,
  supported: ReadonlySet<string>,
  node?: string,
): { layer: ResolvedLayer; issues: Issue[] } {
  const issues: Issue[] = [];
  const n = doc.nodes.get(layer) as IndexedNode;
  const p = node === undefined ? n.params : withNode(n.params, node);
  // A section with no blocks (its data unavailable) reports nothing more: nothing cascades.
  const say = (code: string, param: string, message: string) =>
    section.blocks.length > 0 && issues.push(issueOf(reg, code, `${n.path}/params/${param}`, message));
  const track = section.layers.find((l) => l.ref === layer)?.track ?? -1;

  // The type drawn: the declared one, else down its fallback chain (the layer's first step, then the types'), each
  // checked against the viewer, the section's axis and its layout.
  const axis = section.row.axis === SectionAxis.columns ? "columns" : "blocks";
  const layout = section.layout === Layout.text ? "text" : "tracks";
  let type: string | null = p.type as string;
  for (let k = 0; type !== null && type !== "omit"; k++) {
    const e: LayerTypeEntry | undefined = reg.get("layerTypes", type);
    if (e && supported.has(type) && e.axes.includes(axis) && e.layouts.includes(layout)) break;
    type = k > 8 ? null : ((k === 0 && typeof p.fallback === "string" ? p.fallback : e?.fallback) ?? null);
  }
  const entry = type === null ? undefined : reg.get("layerTypes", type);
  if (!entry) say("svs.render.omitted", "type", `"${p.type}" cannot be drawn here; its fallbacks end in omit`);
  else if (type !== p.type) say("svs.render.fallback", "type", `"${p.type}" cannot be drawn here; drawn as "${type}"`);

  // Themes after defaults and fallbacks; a named rule table expanded (its name kept as `table`), so a viewer implements
  // `rules` alone, the themes nested in a `rules` theme resolved alike; glyphs a viewer cannot draw as their fallback
  // shape.
  const b0 = section.blocks[0];
  const aln0 = b0 ? set.alignments[b0.alignment] : undefined;
  const alphabet = b0 && (aln0?.alphabet ?? set.molecules[set.units.molecule[b0.unit]]?.alphabet);
  // The color's subject, set before the themes resolve (below).
  let subject: SchemeSubject | null = null;
  const theme = (
    given: ThemeRef | null,
    kind: ChannelKind,
    param: string,
    depth = 0,
    byDefault: ThemeRef | null = null,
  ): ThemeRef | null => {
    let r = given;
    for (let k = 0; r && k < 5; k++) {
      const e = reg.get("themes", r.theme);
      // Clustal X reads a protein column's composition: not on nucleotides, nor a lazy alignment's.
      const nt = r.theme === "clustal_x" && (alphabet !== "protein" || aln0?.lazy === true);
      // A color scheme of another subject (Viridis on letters), or unknown here: the layer's default
      // (seqspec:format#color-schemes).
      const def = r.theme === "scheme" ? schemeOf(reg, r.params ?? {}) : undefined;
      const ns = r.theme === "scheme" && !(def && subject && schemeFits(def, subject));
      // Only namespaced themes are checked against `supported`: a kernel theme is registered, so it counts as drawn.
      if (e?.params[kind] && !nt && !ns && (!isNamespaced(r.theme) || supported.has(r.theme))) break;
      r = nt
        ? { theme: "letter" }
        : ns
          ? (r.fallback ?? (byDefault !== r ? byDefault : null))
          : (r.fallback ?? e?.fallback ?? null);
    }
    if (r !== given)
      say("svs.render.theme", param, `"${given?.theme}" is drawn as ${r ? `"${r.theme}"` : "the default"}`);
    const expand = r && depth < 4 ? reg.get("themes", r.theme)?.expand : undefined;
    if (r && expand) r = { ...expand(r.params ?? {}, kind), table: r.theme };
    // Normalized, so a viewer reads every param (a type's default is written short).
    r = r && (themeRef(kind).normalize?.(r, reg) as ThemeRef);
    if (r?.theme !== "rules" || depth >= 4) return r;
    const q = r.params as { rules: { theme: ThemeRef | null }[]; default: ThemeRef | null };
    const nested = (x: ThemeRef | null) => x && theme(x, kind, param, depth + 1);
    const rules = q.rules.map((x) => (x.theme ? { ...x, theme: nested(x.theme) } : x));
    return { ...r, params: { ...q, rules, default: nested(q.default) } };
  };
  // Each kind of glyph drawn as its fallback is reported once per layer.
  const fellBack = new Set<string>();
  const shape = (g: unknown) => {
    if (
      !isPlainObject(g) ||
      !((g.char != null && !supported.has("char")) || (g.path != null && !supported.has("path")))
    )
      return g;
    const kind = g.char != null ? "char" : "path";
    if (!fellBack.has(kind) && fellBack.add(kind))
      say("svs.render.glyph", "glyph", `a ${kind} glyph is drawn as its fallback`);
    return g.fallback ?? "circle";
  };
  const glyphs = (glyph: ThemeRef | null) => {
    const gp = glyph?.params;
    const one = gp && glyph?.theme === "uniform" ? shape(gp.glyph) : undefined;
    if (gp && glyph?.theme === "uniform" && one !== gp.glyph) return { ...glyph, params: { ...gp, glyph: one } };
    if (gp && glyph?.theme === "field_map" && isPlainObject(gp.map)) {
      const map = Object.fromEntries(Object.entries(gp.map).map(([k, g]) => [k, shape(g)]));
      const fallback = shape(gp.default);
      if (Object.keys(map).some((k) => map[k] !== (gp.map as Record<string, unknown>)[k]) || fallback !== gp.default)
        return { ...glyph, params: { ...gp, map, default: fallback } };
    }
    return glyph;
  };

  // Data (a matrix's first data's states win), the rows the filter keeps, the domain over every row.
  const own = childrenOf(doc, n).filter((c) => c.kind === "annotation");
  const refs = p.data != null ? refsOf(p.data) : own.map((c) => c.ref);
  const all = refs.flatMap((r) =>
    set.byName(r, t) && t.kind === KeyKind.annotation ? [set.annotations[t.id] as Annotation] : [],
  );
  const states = all.find((a) => a.matrix)?.matrix?.states;
  const other = (a: Annotation) =>
    !!states &&
    !!a.matrix &&
    (a.matrix.states.length !== states.length || a.matrix.states.some((x, i) => x !== states[i]));
  // Valid data the layer cannot draw beside the first matrix: a loss, a warning (seqspec:format#validation).
  if (all.some(other))
    say("svs.render.omitted", "data", "data whose states differ from the first matrix's are left out");
  const anns = all.filter((a) => !other(a));
  // A rows or matrix layer's `filter` is a list of key filters over its members (keyed.ts), not a filter of the data's rows.
  const members = p.type === "matrix" && states ? matrixStates(set, reg, states, p, say) : null;
  const filter =
    p.type === "rows" || p.type === "matrix" ? null : (p.filter as { field: string | null; values: unknown[] } | null);
  const keeps = filter && matcher(filter.values);
  const kept = anns.map((a) => {
    if (!filter || !keeps) return null;
    const out = new Uint8Array(a.count);
    for (let r = 0; r < a.count; r++) out[r] = +keeps(a, set, r, filter.field ?? "category");
    return out;
  });
  // A split by a field no row of the data has (streamed data is not read here).
  const by = (p.split as { by: string } | null)?.by;
  const lacks = (a: Annotation) => {
    for (let r = 0; r < a.count; r++) if (paletteValue(a, set, r, by ?? "", reg) != null) return false;
    return !a.streamed;
  };
  if (by && anns.length > 0 && anns.every(lacks))
    say("svs.layer.split", "split", `no row of its data has "${by}": split by a field it has (category, alt, …)`);
  // The domain: the layer's own `domain`, else its data's; a `domain_group` shares the union of its members' (domains.ts).
  const explicit = Array.isArray(p.domain);
  let [lo, hi] = explicit ? (p.domain as number[]) : dataDomain(anns);
  const shared = node === undefined ? sharedDomains(doc, set, groupMembers(doc, n), reg) : null;
  if (shared && shared.lo <= shared.hi) [lo, hi] = [shared.lo, shared.hi];

  // Placements: every block the data reaches, on its row's lane in a `rows` track.
  const placements: PlacementRow[] = [];
  const paths: ComposedPath[] = [];
  const index = new Map<ComposedPath, number>();
  const tr = section.tracks[track];
  // A lazy run's lane is its row's rank: in data order the row itself; in an order ranked below (one scan of it), a row it leaves out dropped.
  const run = tr && tr.laneCount > 0 ? section.lanes.lazy.find((z) => z.start === tr.laneStart) : undefined;
  const ranked: number[] = [];
  const lanes = new Map<number, number>();
  for (let i = tr?.laneStart ?? 0; tr && !run && i < tr.laneStart + tr.laneCount; i++)
    lanes.set(section.lanes.row[i], i);
  /** A row's lane; -1: not on a track with lanes; -2: not on one of its lanes (left out). Tree order: the row, ranked below. */
  const laneOf = (row: number) => {
    if (row < 0) return -1;
    // With an order the row stands in for its lane until ranked below, where a row the order leaves out is dropped.
    if (run) return run.order ? row : row >= run.length ? -2 : run.start + row;
    return lanes.size > 0 ? (lanes.get(row) ?? -2) : -1;
  };
  const shown = section.blocks.flatMap((b, i) =>
    !Array.isArray(p.blocks) || b.nodes.some((x) => (p.blocks as string[]).includes(x)) ? [i] : [],
  );
  const add = (
    block: number,
    annotation: number,
    slot: number,
    member: number,
    path: ComposedPath | null,
    row: number,
  ) => {
    const lane = laneOf(row);
    if (lane === -2) return;
    if (row >= 0 && run?.order) ranked.push(placements.length);
    if (path && !index.has(path)) index.set(path, paths.push(path) - 1);
    placements.push({ block, annotation, slot, member, lane, path: path ? (index.get(path) as number) : -1 });
  };
  let [ambiguous, inserts, merged] = [false, 0, false];
  // `letters` and `ruler` with `of`: other polymers' letters or numberings on every block their coordinate
  // systems reach, as data is; each placement's path's `from` names its polymer. `mapped`: every polymer a codon
  // mapping takes to the block's coordinate system, a row's included.
  const ofs = (type === "letters" || type === "ruler") && p.of != null ? ofSystems(set, p.of) : null;
  if (ofs) {
    for (const [from, to] of ofs)
      for (const i of shown) {
        if (!onTo(set, section.blocks[i], to)) continue;
        ambiguous ||= !reach(set, from, section.blocks[i], (member, path, row) => add(i, -1, -1, member, path, row));
      }
    if (placements.length === 0 && node === undefined)
      say("svs.layer.unplaced", "of", `"${p.of}" reaches no block of this section`);
  } else if (entry && refs.length === 0) {
    const where = type === "rows" ? shown.filter((i) => section.blocks[i].alignment >= 0) : shown;
    if (entry.accepts === "none" || isNamespaced(p.type)) for (const i of where) add(i, -1, -1, -1, null, -1);
  }
  for (const a of entry ? anns : []) {
    if (a.geometry === "rows") continue;
    for (const end of endsOf(a)) {
      const own = ownMolecule(set, section.blocks, end.system);
      for (const i of shown) {
        const block = section.blocks[i];
        if (own >= 0 && (block.unit < 0 || set.units.molecule[block.unit] !== own)) continue;
        // A streamed annotation reaches a block through a path too: its regions are asked back through it (ST7,
        // seqspec:sources#data-patches).
        ambiguous ||= !reach(set, end.system, block, (member, path, row) => {
          // Not the copies of a merged block each by its own path: their regions would land on one set of indices.
          if (a.streamed && path && member >= 0) {
            merged = true;
            return;
          }
          add(i, a.id, end.slot, member, path, row);
          if (!path || a.geometry !== "ranges") return;
          for (const r of end.rows) if (a.start[r] < a.stop[r] && insert(path, a.start[r], a.stop[r])) inserts++;
        });
      }
    }
  }
  if (merged)
    say(
      "svs.layer.streamed_merged",
      "data",
      'streamed data reaches the copies of a merged block by different paths: show the copies apart (copies: "each")',
    );
  if (ambiguous) say("svs.mapping.ambiguous", "data", "several shortest paths reach a block with different results");
  if (inserts > 0)
    say(
      "svs.layer.insert",
      "data",
      `${inserts} range(s) lie inside an insert (unpaired between paired positions): drawn at the boundary after it`,
    );
  // Data with no rows (an empty table) reaches nothing, and that is no fault of the layer's (VW2).
  const rowless = anns.length > 0 && anns.every((a) => a.count === 0 && a.dense === null && !a.streamed);
  // Data naming only a streamed source's dynamic nodes, none held yet, is pending until one arrives, noted once per
  // source with a count and the first refs (`svs.ref.pending`, info; seqspec:sources#dynamic-nodes).
  const waited = refs.map((r) => waitedFor(doc, set, r));
  const pending = refs.length > 0 && all.length === 0 && waited.every((w) => w !== null);
  if (pending) notePending(waited as Waited[], (m) => say("svs.ref.pending", "data", m));
  if (entry && refs.length > 0 && placements.length === 0 && !rowless && !pending && !merged && node === undefined)
    say("svs.layer.unplaced", "data", "its data reaches no block of this section");
  // The color's default: the providers' when all the data is one provider's that declares one; a numeric
  // one over the layer's domain for a matrix, for color over values and a dot plot over grids, sequential or diverging
  // at 0.
  const from = refs.map((r) => doc.nodes.get(r)).map((x) => (x?.kind === "computed" ? String(x.params.provider) : ""));
  const declared = from.length > 0 && from.every((x) => x === from[0]) ? reg.get("providers", from[0]) : undefined;
  const numeric =
    type === "matrix" ||
    (type === "color" && anns.length > 0 && anns.every((a) => a.geometry === "values")) ||
    (type === "dot_plot" && anns.length > 0 && anns.every((a) => a.geometry === "grid"));
  const byDefault = numeric ? numericDefaultColor(lo, hi, shared ? shared.explicit : explicit) : null;
  // The color's subject: letters of `of`'s polymer where the layer has one.
  const layerSubject = () => {
    const from = p.of != null && placements[0]?.path >= 0 ? paths[placements[0].path].from : null;
    const K = CoordinateSystemKind;
    const m = from?.kind === K.unit ? set.units.molecule[from.id] : from?.kind === K.molecule ? from.id : -1;
    return schemeSubject(set, section, entry, anns, from ? (set.molecules[m]?.alphabet ?? null) : undefined);
  };
  const used: SharedField[] = [];
  const themes = (entry?.channels ?? []).map((c) => {
    if (c.kind === "color") subject = layerSubject();
    const given = (p[c.name] ??
      (c.kind === "color" ? (declared?.defaultColor ?? byDefault) : null) ??
      c.default) as ThemeRef | null;
    const r = theme(
      given,
      c.kind,
      c.name,
      0,
      c.kind === "color" ? (declared?.defaultColor ?? byDefault ?? c.default) : null,
    );
    return c.kind === "glyph" ? glyphs(r) : shared ? fillShared(r, shared, used) : r;
  });
  if (run?.order && ranked.length > 0) {
    const rank = new Map<number, number>();
    for (const i of ranked) rank.set(placements[i].lane, -1);
    for (let k = 0; k < run.order.length; k++) if (rank.has(run.order[k])) rank.set(run.order[k], run.start + k);
    const out: PlacementRow[] = [];
    const rankedAt = new Set(ranked);
    placements.forEach((x, i) => {
      const lane = rankedAt.has(i) ? (rank.get(x.lane) as number) : x.lane;
      if (lane !== -1 || !rankedAt.has(i)) out.push(rankedAt.has(i) ? { ...x, lane } : x);
    });
    placements.length = 0;
    placements.push(...out);
  }
  placements.sort((x, y) => x.block - y.block || x.lane - y.lane);

  const str = (s: string) => set.strings.intern(s);
  const out: ResolvedLayer = {
    ref: node === undefined ? layer : instanceRef(layer, node),
    track,
    type: entry && type ? str(type) : -1,
    declared: str(p.type as string),
    sizing: entry?.sizing === "free" ? Sizing.free : Sizing.row,
    themes,
    domainMin: lo <= hi ? lo : Number.NaN,
    domainMax: lo <= hi ? hi : Number.NaN,
    domainGiven: shared ? shared.explicit : explicit,
    group: shared ? String(p.domain_group) : null,
    shared: used,
    data: anns.map((a) => a.id),
    members,
    kept,
    placements,
    paths,
    params: letterMapped(paramsOf(n), letterMapParam(doc, n), reg, () =>
      ofs
        ? ofs.map(([from]) => alphabetOf(set, from))
        : shown.map((i) => {
            const b = section.blocks[i];
            return b.alignment >= 0
              ? (set.alignments[b.alignment]?.alphabet ?? null)
              : (set.molecules[set.units.molecule[b.unit]]?.alphabet ?? null);
          }),
    ),
    pending,
    node: node ?? null,
  };
  return { layer: out, issues };
}

/** A dynamic node a layer's data waits for: its ref and the streamed source that would send it. */
interface Waited {
  readonly source: string;
  readonly ref: string;
}

/** The streamed source whose dynamic node `ref` names (`<source>/<key>`), else null. */
function sourceUnder(doc: DocIndex, set: SeqSet, ref: string): string | null {
  const sources = streamedSources(set);
  for (let k = ref.indexOf("/"); k > 0; k = ref.indexOf("/", k + 1)) {
    const head = ref.slice(0, k);
    if (sources.has(head) || isStreamedSource(doc.nodes.get(head))) return head;
  }
  return null;
}

/**
 * What a data ref waits for: a dynamic node it names, or the one a document annotation's target names (an item
 * coordinate system, ST7); null when it waits for nothing streamed.
 */
function waitedFor(doc: DocIndex, set: SeqSet, ref: string): Waited | null {
  const own = sourceUnder(doc, set, ref);
  if (own !== null) return { source: own, ref };
  const target = doc.nodes.get(ref)?.params.target;
  const source = typeof target === "string" ? sourceUnder(doc, set, target) : null;
  return source === null ? null : { source, ref: target as string };
}

/** Refs a pending note names before "and n more". */
const PENDING_SHOWN = 3;

/**
 * One note per source (GV-07's rule: one per cause and node, a count and the first examples). No nearest key: a
 * streamed source sends its items as the view reaches them, so a ref not sent yet is no typo until the source has sent
 * all it will, which no source says yet (todo.md).
 */
function notePending(waited: readonly Waited[], note: (message: string) => void): void {
  const by = new Map<string, string[]>();
  for (const w of waited) {
    const refs = by.get(w.source) ?? [];
    if (!refs.includes(w.ref)) refs.push(w.ref);
    by.set(w.source, refs);
  }
  for (const [source, refs] of by) {
    const shown = refs
      .slice(0, PENDING_SHOWN)
      .map((r) => `"${r}"`)
      .join(", ");
    const more = refs.length > PENDING_SHOWN ? ` and ${refs.length - PENDING_SHOWN} more` : "";
    const one = refs.length === 1;
    note(
      `${shown}${more} name${one ? "s" : ""} nothing "${source}" has sent yet: shown when ${one ? "it arrives" : "they arrive"}`,
    );
  }
}

/**
 * A matrix's states in the order and kept by the filter of its key tables (keyed by `state`: a row names the first
 * state of that name); null without either, every state in data order. A matrix has no tree: `tree` is reported and
 * ordered as data; a filter that keeps no state is noted (`svs.layer.empty`, info: the track is drawn empty).
 */
function matrixStates(
  set: SeqSet,
  reg: SvsRegistries,
  states: Int32Array,
  p: Readonly<Record<string, unknown>>,
  say: (code: string, param: string, message: string) => void,
): Int32Array | null {
  const asked = p.order as OrderKeyValue | OrderKeyValue[] | null | undefined;
  const keys = asked == null ? [] : Array.isArray(asked) ? asked : [asked];
  const tree = keys.indexOf("tree");
  if (tree >= 0)
    say(
      "svs.layer.data",
      Array.isArray(asked) ? `order/${tree}` : "order",
      'a matrix has no tree: its states are in data order; order by "data" or a key table\'s field',
    );
  const order = tree >= 0 ? keys.map((k) => (k === "tree" ? "data" : k)) : asked;
  const filter = (Array.isArray(p.filter) ? p.filter : []) as KeyFilterValue[];
  if (keys.every((k) => k === "data" || k === "tree") && filter.length === 0) return null;
  const at = new Map<string, number>();
  for (let i = 0; i < states.length; i++) {
    const name = set.strings.get(states[i] as NameId);
    if (!at.has(name)) at.set(name, i);
  }
  const members: Members = {
    count: states.length,
    tree: null,
    // A row with no state names none (GV-11: `String(undefined)` named a state called "undefined").
    of: (ann, i) => {
      const v = Annotations.rowValue(ann, set, i, "state");
      return v === undefined ? -1 : (at.get(String(v)) ?? -1);
    },
  };
  const k = keyed(set, reg, members, { order, filter }, (where, message) => say("svs.layer.data", where, message));
  const kept = Int32Array.from(k.permutation.filter((m) => k.kept[m] === 1));
  if (kept.length === 0 && states.length > 0)
    say("svs.layer.empty", "filter", "the filter keeps no state of the matrix: its track is drawn empty");
  return kept;
}

/** A coordinate system's letters' alphabet: a row's alignment's, else its molecule's (a ligand's: null). */
function alphabetOf(set: SeqSet, s: CoordinateSystem): Alphabet | null {
  const K = CoordinateSystemKind;
  if (s.kind === K.row) return set.alignments[s.id]?.alphabet ?? null;
  const m = s.kind === K.unit ? set.units.molecule[s.id] : s.kind === K.molecule ? s.id : -1;
  return set.molecules[m]?.alphabet ?? null;
}

/**
 * `of`'s polymers as (their coordinate system, the one a codon mapping takes it to; null: any): a ref, a list of refs,
 * or `mapped` (every polymer the from side of a `codon` mapping, placed only on blocks over its `to`).
 */
function ofSystems(set: SeqSet, of: unknown): [CoordinateSystem, CoordinateSystem | null][] {
  if (of === "mapped") {
    const out = new Map<string, [CoordinateSystem, CoordinateSystem]>();
    for (const m of set.mappings) {
      if (m?.method !== "codon") continue;
      out.set(`${m.from.kind}:${m.from.id}:${m.to.kind}:${m.to.id}:${m.to.row}`, [{ ...m.from }, { ...m.to }]);
    }
    return [...out.values()];
  }
  const kinds: Record<number, number> = {
    [KeyKind.molecule]: CoordinateSystemKind.molecule,
    [KeyKind.unit]: CoordinateSystemKind.unit,
    [KeyKind.row]: CoordinateSystemKind.row,
  };
  const out: [CoordinateSystem, null][] = [];
  for (const r of Array.isArray(of) ? of : [of]) {
    if (typeof r !== "string" || !set.byName(r, t) || kinds[t.kind] === undefined) continue;
    out.push([{ kind: kinds[t.kind], id: t.id, row: t.kind === KeyKind.row ? t.row : -1 }, null]);
  }
  return out;
}

/**
 * A unit or molecule end's molecule when a block of the section shows it, else -1 (direct placement wins): such an end
 * is placed only on its molecule's blocks, copies by the copy rule, never through an alignment or another path onto
 * another molecule's block.
 */
function ownMolecule(set: SeqSet, blocks: ResolvedSection["blocks"], s: CoordinateSystem): number {
  const mol =
    s.kind === CoordinateSystemKind.unit
      ? set.units.molecule[s.id]
      : s.kind === CoordinateSystemKind.molecule
        ? s.id
        : -1;
  for (let i = 0; mol >= 0 && i < blocks.length; i++)
    if (blocks[i].alignment < 0 && set.units.molecule[blocks[i].unit] === mol) return mol;
  return -1;
}

/** A block over `to` (null: any): its alignment's for a row or the columns, else its unit's molecule's. */
function onTo(set: SeqSet, b: { unit: number; alignment: number }, to: CoordinateSystem | null): boolean {
  if (!to) return true;
  if (to.kind === CoordinateSystemKind.row || to.kind === CoordinateSystemKind.columns) return b.alignment === to.id;
  const mol = to.kind === CoordinateSystemKind.unit ? set.units.molecule[to.id] : to.id;
  return b.unit >= 0 && set.units.molecule[b.unit] === mol;
}

/** [start, stop) unmapped, between positions the path pairs (an insert; beyond them, data draws nothing). */
function insert(path: ComposedPath, start: number, stop: number): boolean {
  const r = path.runs;
  return r.length > 0 && r[0] < start && r[r.length - 1] > stop && unmapped(path, start, stop);
}
/**
 * No element of [start, stop) maps through the path (a range inside an insert): its held elements read in place, from
 * the first run ending past `start` (a binary search).
 */
function unmapped(path: ComposedPath, start: number, stop: number): boolean {
  const { runs, slots, heldTo } = path;
  let lo = 0;
  let hi = slots.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (runs[2 * mid + 1] <= start) lo = mid + 1;
    else hi = mid;
  }
  for (let k = lo; k < slots.length && runs[2 * k] < stop; k++)
    for (let i = Math.max(start, runs[2 * k]); i < Math.min(stop, runs[2 * k + 1]); i++)
      if (heldTo[slots[k] + i - runs[2 * k]] >= 0) return false;
  return true;
}
