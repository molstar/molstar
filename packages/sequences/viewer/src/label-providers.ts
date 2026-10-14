/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The built-in label providers (seqspec:viewer#labels): structured lines for the hovered element, called by
// `labels.compute` once per hovered element (labels.ts), so their strings are the only allocation. `applies` checks
// numbers only. Numberings reached through mappings (another isoform's) are found once per unit and set (`lc.reached`),
// then read.

import { paletteValue, TILE } from "@molstar/seq-core/format";
import {
  type Alignment,
  type Annotation,
  Annotations,
  CoordinateSystemKind,
  type Numbering,
  OrderedSet,
  Paths,
  type SeqSet,
  SortedArray,
} from "@molstar/seq-core/model";
import { codonLabel } from "./codon-label.ts";
import type { ViewModel } from "./graph.ts";
import { copyAt, HoverKind, sourceAt } from "./hit.ts";
import type { HoverState } from "./hover.ts";
import type { LabelContext, LabelProvider, LabelSink } from "./labels.ts";
import type { PerPosition } from "./prepare.ts";
import type { MatrixState } from "./rows.ts";
import { RegionState, type StreamItems } from "./streams.ts";
import { TreeNode, type TreeRows } from "./tree.ts";

const K = HoverKind;
const name = (set: SeqSet, id: number) => (id >= 0 ? set.strings.get(id) : "");
/** A position's label in a numbering; "" outside it. */
const fmt = (lc: LabelContext, n: Numbering, i: number) => (i >= 0 && i < n.length ? lc.format(n, i as never) : "");

/** The hovered unit's block in the section: its own or a merged copy group's; -1: none. */
function blockOf(v: ViewModel, section: number, unit: number): number {
  const r = v.resolved;
  const b0 = r.sections.blockStart[section];
  for (let b = b0; b < b0 + r.sections.blockCount[section]; b++) {
    const m = r.blocks.members[b];
    if (r.blocks.unit[b] === unit || (m && SortedArray.indexOf(m, unit) >= 0)) return b;
  }
  return -1;
}
/** The hovered block (`local`, the section's) when it shows the unit: of two pieces of one unit, the one hovered. */
function blockAt(v: ViewModel, section: number, unit: number, local: number): number {
  const r = v.resolved;
  const b = r.sections.blockStart[section] + local;
  const m = local >= 0 && local < r.sections.blockCount[section] ? r.blocks.members[b] : null;
  const shows = local >= 0 && (r.blocks.unit[b] === unit || (m && SortedArray.indexOf(m, unit) >= 0));
  return shows ? b : blockOf(v, section, unit);
}

const one = new Uint8Array(1);
/**
 * A lazy molecule's letter at index i of resolved block b, from a `letters` layer of the hovered track; 0: none held
 * (the label is made again when a region arrives).
 */
function streamedLetter(h: HoverState, lc: LabelContext, b: number, i: number): number {
  const r = lc.view.resolved;
  const local = b - r.sections.blockStart[lc.section];
  for (let k = 0; b >= 0 && k < r.layers.count; k++) {
    const st = r.layers.track[k] === h.track && typeOf(lc, k) === "letters" ? lc.stream(k) : null;
    if (!st?.letters) continue;
    if (st.letters(local, i, i + 1, one) > 0) return one[0];
    lc.waiting();
  }
  return 0;
}

/** Layer k's type name ("" when omitted). */
const typeOf = (lc: LabelContext, k: number) => {
  const t = lc.view.resolved.layers.type[k];
  return t >= 0 ? lc.set.strings.get(t) : "";
};
type MapState = { readonly map?: Uint8Array | null; readonly of?: { readonly unit: Int32Array } | null } | undefined;
/**
 * The `letter_map` table layer k shows a residue of `unit` through (-1: a row's): null as stored; undefined when the
 * layer shows no letters of it (not `letters` or `rows`, or an `of` layer of other polymers).
 */
function mapOf(lc: LabelContext, k: number, unit: number): Uint8Array | null | undefined {
  const type = k >= 0 ? typeOf(lc, k) : "";
  const s = type === "letters" || type === "rows" ? (lc.view.prepared[k]?.state as MapState) : undefined;
  if (!s || (s.of && (unit < 0 || !s.of.unit.includes(unit)))) return undefined;
  return s.map ?? null;
}

/**
 * The letter shown for a stored residue code of `unit` (`letter_map`, seqspec:format#layer-types): the hovered layer's
 * map when it shows the residue, else the first layer of the section that does (a position under a features track reads
 * what the letters above it show); the code itself when that layer maps nothing. The model keeps the stored letter.
 */
function shownCode(h: HoverState, lc: LabelContext, code: number, unit: number): number {
  const r = lc.view.resolved;
  let map = mapOf(lc, h.layer, unit);
  for (let k = 0; map === undefined && k < r.layers.count; k++)
    if (r.tracks.section[r.layers.track[k]] === lc.section) map = mapOf(lc, k, unit);
  return map ? map[code & 127] : code;
}

/** The position label's working record (reused, cleared on return: it holds no set between labels). */
const P = {
  lc: null as LabelContext | null,
  out: null as LabelSink | null,
  bn: null as Numbering | null,
  number: "",
  named: false,
};
/** A line per other numbering whose label differs from the title's; `ordinal` only alone; never `label` (a structure's internal one). */
function numbered(n: Numbering, j: number) {
  const lc = P.lc as LabelContext;
  const key = name(lc.set, n.name);
  const v = key === "ordinal" || n === P.bn ? "" : fmt(lc, n, j);
  P.named ||= key !== "ordinal";
  if (key !== "label" && v !== "" && v !== P.number) (P.out as LabelSink).line(key, v);
}

/**
 * The position label: the title in the block's numbering, formatted once; every other numbering (the unit's, the
 * molecule's, those mappings reach) formatted once and written only when its string differs from the title's number;
 * never `label`; `ordinal` only when nothing named was seen; equal lines grouped by the sink.
 */
const position: LabelProvider = {
  name: "position",
  priority: 100,
  applies: (h) => h.kind === K.seq,
  label(h, lc, out) {
    const set = lc.set;
    const u = h.seq.unit;
    const i = h.seq.index;
    const mol = set.molecules[set.units.molecule[u]];
    if (!mol) return;
    const r = lc.view.resolved;
    const b = blockAt(lc.view, lc.section, u, h.block);
    const bn = b >= 0 ? r.blocks.numbering[b] : null;
    const np = b >= 0 ? r.blocks.numberingPath[b] : -1;
    const at = np >= 0 && r.paths[np] ? sourceAt(r.paths[np], i) : i;
    const label = name(set, b >= 0 && r.blocks.label[b] >= 0 ? r.blocks.label[b] : set.units.label[u]);
    // A lazy molecule's letter from the hovered track's letters regions (ST3a); none held: no letter.
    const stored = mol.element ? 0 : mol.lazy ? streamedLetter(h, lc, b, i) : mol.code[i];
    const shown = mol.element || stored === 0 ? 0 : shownCode(h, lc, stored, u);
    const residue = mol.element ? name(set, mol.chemCompId?.[0] ?? -1) : shown ? String.fromCharCode(shown) : "";
    const number = bn ? fmt(lc, bn, at) : `${i + 1}`;
    // A construct position whose source residue has several copies says which, past the first (the ruler's rule).
    const copy = np >= 0 && bn ? copyAt(r.paths[np], i) : 1;
    out.title(`${label} ${residue}${mol.element ? " " : ""}${number}${copy > 1 ? ` (copy ${copy})` : ""}`.trim());
    if (shown !== stored) out.line("stored", String.fromCharCode(stored));
    P.lc = lc;
    P.out = out;
    P.bn = bn;
    P.number = number;
    P.named = bn !== null && name(set, bn.name) !== "ordinal";
    for (const n of set.units.numberings[u].items) numbered(n, i);
    for (const n of mol.numberings.items) numbered(n, i);
    for (const x of lc.reached(u)) {
      const j = Paths.at(x.path, i);
      if (j >= 0) numbered(x.n, j);
    }
    if (!P.named) out.line("ordinal", number);
    P.lc = null;
    P.out = null;
    P.bn = null;
    const obs = set.units.observed[u];
    let observed = obs === null;
    for (let k = 0; obs && k < obs.length; k += 2) observed ||= i >= obs[k] && i < obs[k + 1];
    if (!observed) out.line("unobserved", "not in the model");
  },
};

const copies: LabelProvider = {
  name: "copies",
  priority: 95,
  applies: (h) => h.kind === K.seq,
  label(h, lc, out) {
    const r = lc.view.resolved;
    const b = blockOf(lc.view, lc.section, h.seq.unit);
    const m = b >= 0 ? r.blocks.members[b] : null;
    if (!m) return;
    const labels = Array.from(m, (u) => name(lc.set, lc.set.units.label[u]));
    out.line("copies", `× ${m.length} (${labels.join(", ")})`);
    const varying = r.blocks.varying[b];
    if (varying && OrderedSet.has(varying, h.seq.index)) out.line("copies", "varies");
  },
};

/**
 * A hovered tree node: the clade it leads into, named (its label, else "Clade"), its size in rows, and what a click
 * does. The tree is the section's rows alignment's (the header's hover names no alignment: its first block's).
 */
function clade(h: HoverState, lc: LabelContext, out: LabelSink) {
  const r = lc.view.resolved;
  const aln = lc.set.alignments[r.blocks.alignment[r.sections.blockStart[lc.section]] ?? -1];
  const t = aln?.tree;
  if (!t || h.treeNode < 0 || h.treeNode >= t.count) return;
  let n = 0;
  for (let k = h.treeNode; k < t.cladeEnd[h.treeNode]; k++) if (t.leafRow[k] >= 0) n++;
  let collapsed = false;
  for (let k = 0; k < r.layers.count && !collapsed; k++) {
    const g = r.layers.track[k] === h.track ? (lc.view.prepared[k]?.state as TreeRows | undefined)?.tree : null;
    collapsed = !!g && g.state[h.treeNode] === TreeNode.wedge;
  }
  out.title(name(lc.set, t.label[h.treeNode]) || "Clade");
  out.line("clade", `${n} ${n === 1 ? "row" : "rows"}`);
  out.line("click", collapsed ? "expand" : "collapse");
}

/** What a streamed row's cell is while its tile is not here: in flight, or the source could not serve it. */
const LOADING = -1;
const UNAVAILABLE = -2;
/**
 * A streamed (lazy) row's residue code at `column` from the finest tile the rows layer's region cache holds for it
 * (0: no residue there, a gap), `LOADING` while the tile is in flight (the label then waits on its arrival),
 * `UNAVAILABLE` when its source failed. The tile's row is the row's lane offset in its track, as the draw reads it.
 */
function streamedCell(h: HoverState, lc: LabelContext, row: number, column: number): number {
  const r = lc.view.resolved;
  let layer = -1;
  for (let k = 0; k < r.layers.count && layer < 0; k++)
    if (r.layers.track[k] === h.track && typeOf(lc, k) === "rows") layer = k;
  const stream = layer >= 0 ? lc.stream(layer) : null;
  if (!stream || column < 0) return LOADING;
  const k = (h.lane >= 0 ? h.lane : row) - (h.track >= 0 ? r.tracks.laneStart[h.track] : 0);
  const [ty, tx] = [k >> 8, column >> 8];
  const state = stream.state(0, 0, ty, tx);
  if (state === RegionState.empty) return 0;
  if (state === RegionState.failed) return UNAVAILABLE;
  const tile = state === RegionState.loaded ? stream.tile(0, 0, ty, tx) : null;
  if (tile) return tile.cls[(k & 255) * TILE + (column & 255)];
  lc.waiting();
  return LOADING;
}

const row: LabelProvider = {
  name: "row",
  priority: 90,
  applies: (h) => h.kind === K.aligned || h.kind === K.row || (h.kind === K.header && h.treeNode >= 0),
  label(h, lc, out) {
    const set = lc.set;
    if (h.kind === K.header) return clade(h, lc, out);
    const al = h.kind === K.aligned ? h.aligned.alignment : h.row.alignment;
    const aln = set.alignments[al];
    if (!aln) return;
    const r = h.kind === K.aligned ? h.aligned.row : h.row.row;
    if (aln.lazy) return lazyRow(h, lc, out, aln, r);
    out.title(name(set, aln.labelOf(r)) || `row ${r + 1}`);
    const index = h.kind === K.aligned ? aln.rows.colToPos(r, h.aligned.column) : h.row.index;
    const rowN = aln.numberingsOf(r)?.byName(set.strings.find("row"));
    const stored = index >= 0 ? aln.rows.residues[aln.rows.residueStart[aln.rows.slotOf(r)] + index] : 0;
    const shown = index >= 0 ? shownCode(h, lc, stored, -1) : 0;
    const residue = index >= 0 ? String.fromCharCode(shown) : "";
    out.line("residue", index < 0 ? "gap" : `${residue}${rowN ? fmt(lc, rowN, index) : index + 1}`);
    if (shown !== stored) out.line("stored", String.fromCharCode(stored));
    if (h.kind === K.aligned) out.line("column", `${h.aligned.column + 1}`);
  },
};

/**
 * A streamed row: its label once its source served it (the held one's, else the manifest's chunk), the residue of its
 * tile (the letter; its number in the row needs the whole row, so none), "loading" while the tile is in flight.
 */
function lazyRow(h: HoverState, lc: LabelContext, out: LabelSink, aln: Alignment, r: number) {
  const set = lc.set;
  const column = h.kind === K.aligned ? h.aligned.column : -1;
  const held = aln.labelOf(r as never);
  let label = held >= 0 ? name(set, held) : "";
  for (let k = 0; label === "" && k < lc.view.resolved.layers.count; k++)
    if (lc.view.resolved.layers.track[k] === h.track && typeOf(lc, k) === "rows") label = lc.stream(k)?.label(r) ?? "";
  out.title(label || `row ${r + 1}`);
  const stored = column >= 0 ? streamedCell(h, lc, r, column) : 0;
  if (stored === LOADING) out.line("residue", "loading");
  else if (stored === UNAVAILABLE) out.line("residue", "not available");
  else if (stored <= 0 || stored === 45 || stored === 46) out.line("residue", "gap");
  else {
    const shown = shownCode(h, lc, stored, -1);
    out.line("residue", String.fromCharCode(shown));
    if (shown !== stored) out.line("stored", String.fromCharCode(stored));
  }
  if (h.kind === K.aligned) out.line("column", `${column + 1}`);
}

/** The first row of an item. */
function firstRow(item: Int32Array, it: number) {
  for (let r = 0; r < item.length; r++) if (item[r] === it) return r;
  return -1;
}

/** A row's field as a label shows it: the data's own, else a field provider's (`ns:name`); undefined: none. */
function fieldOf(lc: LabelContext, a: Annotation, row: number, field: string): unknown {
  const v = paletteValue(a, lc.set, row, field, lc.registry.svs);
  return v === null || (typeof v === "number" && Number.isNaN(v)) ? undefined : v;
}

/** The fields a layer's channel themes read (`field` anywhere in their params: palettes, scales, rule tables). */
const encodedOf = new WeakMap<object, string[]>();
function encoded(params: Readonly<Record<string, unknown>>): string[] {
  let out = encodedOf.get(params);
  if (out) return out;
  out = [];
  const walk = (x: unknown, key: string) => {
    if (key === "field" && typeof x === "string" && !out?.includes(x)) out?.push(x);
    else if (x && typeof x === "object") for (const [k, y] of Object.entries(x)) walk(y, k);
  };
  for (const v of Object.values(params)) if (v && typeof v === "object" && "theme" in v) walk(v, "");
  encodedOf.set(params, out);
  return out;
}

/** `byName`'s target, reused; the item's fields shown after its title. */
const target = { kind: -1, id: -1, row: -1 };
const FIELDS = ["category", "tooltip", "alt"];
const item: LabelProvider = {
  name: "item",
  priority: 80,
  applies: (h) =>
    h.kind === K.item ||
    (h.kind === K.seq && h.item.annotation >= 0) ||
    (h.dataColumn >= 0 && (h.kind === K.row || h.kind === K.aligned)),
  label(h, lc, out) {
    const set = lc.set;
    if (h.dataColumn >= 0 && h.layer >= 0) {
      const c = (
        lc.view.resolved.layers.params[h.layer].data_columns as { data: string; field: string; label?: string }[]
      )?.[h.dataColumn];
      const a = c && set.byName(c.data, target) ? set.annotations[target.id] : null;
      const r = h.kind === K.aligned ? h.aligned.row : h.row.row;
      const v = a && a.count > r ? fieldOf(lc, a, r, c.field) : undefined;
      if (c && v !== undefined) out.line(c.label ?? c.field, String(v));
      if (h.kind !== K.item) return;
    }
    const a = set.annotations[h.item.annotation];
    const r = a ? firstRow(a.item, h.item.item) : -1;
    if (!a || r < 0) return;
    const field = (f: string) => fieldOf(lc, a, r, f);
    out.title(String(field("label") ?? field("id") ?? (name(set, a.label) || "item")));
    if (a.label >= 0) out.line("annotation", name(set, a.label));
    for (const f of FIELDS) if (field(f) !== undefined) out.line(f, String(field(f)));
    // The fields the hovered layer's themes encode (its color or size by a field: "score 64").
    const params = h.layer >= 0 ? lc.view.resolved.layers.params[h.layer] : null;
    for (const f of params ? encoded(params) : [])
      if (!FIELDS.includes(f) && field(f) !== undefined) out.line(f, String(field(f)));
    for (let k = 0; k < a.links.length; k++) {
      const url = Annotations.linkOf(a, set, k, r);
      if (url) out.line(name(set, a.links[k].label) || "link", url);
    }
  },
};

/** A value line, to 4 significant digits; none for NaN (no value there). */
const valueLine = (out: LabelSink, key: string, v: number) =>
  Number.isNaN(v) || out.line(key || "value", `${Number(v.toPrecision(4))}`);

/**
 * The value drawn under the pointer: for each `values` and `color` layer of the hovered track, its prepared value at
 * (block, `index − base`) through LY2's `PerPosition`, so data on the unit, its molecule, a mapped polymer or through
 * codons reads what the draw reads (no `slotOf`); a matrix's at the hovered state.
 */
const value: LabelProvider = {
  name: "value",
  priority: 70,
  applies: (h) => h.kind === K.seq && h.track >= 0 && h.block >= 0,
  label(h, lc, out) {
    const r = lc.view.resolved;
    const b = h.block;
    const j = h.index - r.blocks.base[r.sections.blockStart[lc.section] + b];
    for (let k = 0; k < r.layers.count; k++) {
      if (r.layers.track[k] !== h.track) continue;
      const type = name(lc.set, r.layers.type[k]);
      const st = lc.view.prepared[k]?.state as (PerPosition & MatrixState) | undefined;
      if ((type === "values" || type === "color") && st?.value) {
        const v = st.value[b];
        const at = st.annotation[b];
        const a = v && j >= 0 && j < v.length ? lc.set.annotations[at ? at[j] : st.only[b]] : undefined;
        if (v && a) valueLine(out, name(lc.set, a.label), v[j]);
      } else if (type === "matrix" && st?.rows && h.state >= 0 && h.state < st.stride[b]) {
        const rows = st.rows[b];
        const m = st.ann[b]?.matrix;
        const row = rows && j >= 0 && j < rows.length ? rows[j] : -1;
        if (m && row >= 0) valueLine(out, st.labels[h.state], m.values[row * st.stride[b] + h.state]);
      }
    }
  },
};

/** The partner endpoint of a pair row, in the partner unit's own first named numbering. */
const pair: LabelProvider = {
  name: "pair",
  priority: 70,
  applies: (h) => h.kind === K.item,
  label(h, lc, out) {
    const set = lc.set;
    const a = set.annotations[h.item.annotation];
    if (!a?.pairs) return;
    const end = (u: number, m: number, i: number) => {
      const mol = set.molecules[u >= 0 ? set.units.molecule[u] : m];
      if (!mol || i < 0) return "";
      const ns = u >= 0 ? set.units.numberings[u].items : [];
      const n = ns[0] ?? mol.numberings.items[0];
      const who = name(set, u >= 0 ? set.units.label[u] : mol.label);
      const res = mol.element ? name(set, mol.chemCompId?.[0] ?? -1) : String.fromCharCode(mol.code[i]);
      return `${who} ${res}${n ? fmt(lc, n, i) : i + 1}`;
    };
    for (let r = 0; r < a.count; r++) {
      if (a.item[r] !== h.item.item) continue;
      const s = a.slot[r];
      const b = a.pairs.bTarget;
      out.line(
        "pair",
        `${end(a.slots.unit[s], a.slots.molecule[s], a.start[r])} – ${end(b.unit[r], b.molecule[r], a.pairs.bStart[r])}`,
      );
    }
  },
};

/** Consensus and conservation: annotations of the view over the alignment's columns, at the hovered column. */
const column: LabelProvider = {
  name: "column",
  priority: 60,
  applies: (h) => h.kind === K.column || h.kind === K.aligned,
  label(h, lc, out) {
    const set = lc.set;
    const loc = h.kind === K.column ? h.column : h.aligned;
    const space = { kind: CoordinateSystemKind.columns, id: loc.alignment, row: -1 };
    const data = lc.view.resolved.data;
    if (h.kind === K.column) out.title(`column ${loc.column + 1}`);
    for (let d = 0; d < data.length; d++) {
      // Each annotation once (several layers may draw it), without a set per label.
      if (data.indexOf(data[d]) !== d) continue;
      const a = set.annotations[data[d]];
      const slot = a ? Annotations.slotOf(a, space) : -1;
      if (!a || slot < 0) continue;
      if (a.geometry === "values") {
        const v = Annotations.valueAt(a, slot, loc.column);
        if (!Number.isNaN(v)) out.line(name(set, a.label) || "value", `${Number(v.toPrecision(3))}`);
      } else
        Annotations.itemsOver(
          a,
          slot,
          loc.column,
          loc.column + 1,
          (r: number, o: LabelSink) => {
            const c = Annotations.rowValue(a, set, r, "category") ?? Annotations.rowValue(a, set, r, "label");
            if (c !== undefined) o.line(name(set, a.label) || "consensus", String(c));
          },
          out,
        );
    }
  },
};

/** A renderer type's own lines for the hovered element of its layer (`LayerRenderer.label`). */
const layerEntry: LabelProvider = {
  name: "layer",
  priority: 75,
  applies: (h) => h.layer >= 0,
  label(h, lc, out) {
    const p = lc.view.prepared[h.layer];
    lc.registry
      .get("renderers", p?.renderer ?? "")
      ?.label?.(h, lc, out, p.state, lc.view.resolved.layers.params[h.layer]);
  },
};

const at = new Int32Array(6);
/** A streamed item's columns not written as fields (its place, its key, its map). */
const PLACE = new Set(["start", "stop", "id", "label", "map"]);
const FIELD_LINES = 8;
/**
 * Streamed layers of the hovered track (ST3b): "loading…" over a pending region (made again when it arrives), "failed"
 * over a failed one (Retry in the track's "…" menu); over a streamed item (the hover's `part`, its number in the
 * layer), its label and fields from its region's table, as `item` writes an annotation's.
 */
const stream: LabelProvider = {
  name: "stream",
  priority: 50,
  applies: (h) => h.kind === K.seq && h.track >= 0 && h.block >= 0,
  label(h, lc, out) {
    const r = lc.view.resolved;
    const t = h.track;
    const key = name(lc.set, r.tracks.label[t] >= 0 ? r.tracks.label[t] : r.tracks.ref[t]) || "data";
    for (let k = 0; k < r.layers.count; k++) {
      const st = r.layers.track[k] === t ? lc.stream(k) : null;
      const s = st?.stateAt ? st.stateAt(h.block, h.index) : RegionState.empty;
      if (s === RegionState.pending) {
        out.line(key, "loading…");
        lc.waiting();
      } else if (s === RegionState.failed) out.line(key, "failed: Retry in the track's … menu");
      if (s !== RegionState.loaded || h.part < 0 || k !== h.layer || !st?.items) continue;
      // The data's regions over the index (`at[0]` -1: a data read), the first holding the item the hit found
      // (through a path a split codon's parts may lie in two).
      at[0] = -1;
      let x: StreamItems | null = null;
      let q = -1;
      for (let j = 0; q < 0; j++) {
        x = st.items(h.block, h.index, h.index + 1, j, at);
        if (x === null) break;
        q = x.features && x.item ? x.item.indexOf(h.part) : -1;
      }
      const f = x?.features;
      if (!f || !x?.table || q < 0) continue;
      const cell = (c: string) => {
        const col = x.table?.[c];
        const v = col && f.row[q] < col.length ? col[f.row[q]] : null;
        return v == null || v === "" ? null : String(v);
      };
      out.line(key, cell("label") ?? cell("id") ?? "item");
      let n = 0;
      for (const c in x.table)
        if (!PLACE.has(c) && cell(c) !== null && n++ < FIELD_LINES) out.line(c, cell(c) as string);
    }
  },
};

export const DefaultLabelProviders: readonly LabelProvider[] = [
  position,
  copies,
  row,
  item,
  layerEntry,
  value,
  pair,
  column,
  stream,
  codonLabel,
];
