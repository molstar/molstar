/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The prepare of the alignment and matrix renderers (seqspec:viewer#extension-points): `rows` (an alignment's rows over
// its track's lanes, the tree beside them: `TreeRows`, tree.ts), `logo`, `matrix` and `seqv:map`. A rows layer's data
// by region is its `stream` (streams.ts, LY4), its data columns' prepare RN3's, its draws RN2's and RN3's. Prepare-time
// code allocates.

import {
  followsTree,
  letterMapTable,
  paletteValue,
  type StreamedData,
  StreamProp,
  type SvsRegistries,
} from "@molstar/seq-core/format";
import { type Annotation, type Field, KeyKind, type KeyTarget, Paths, type SeqSet } from "@molstar/seq-core/model";
import type { AxisView } from "./axis.ts";
import type { Hit } from "./hit.ts";
import type { LayerInput, PrepareContext } from "./registry.ts";
import type { StreamLayer } from "./streams.ts";
import { hitTree, TREE_INSET, type TreeRows, treeOf, treeSpan } from "./tree.ts";

/** A rows layer prepared: its alignment, its lanes (offset k from the track's first lane), the held codes, the tree. */
export interface RowsState extends TreeRows {
  readonly lazy: boolean;
  readonly columns: number;
  readonly lanes: number;
  /** Lane offset → its slot in the held rows' codes (-1: not held), and its depth under `parent_rows`. */
  readonly slot: Int32Array | null;
  readonly depth: Int32Array | null;
  readonly codes: Uint8Array | null;
  /** `letter_map` (see `LettersState.map`): applied to the residue a cell and a letter show; null: as stored. */
  readonly map: Uint8Array | null;
  /** Its data by region: tiles fetched from its streamed source, or reduced from the held codes. */
  readonly source: StreamLayer;
  /** Its data columns: the renderer drawn (`text` for an unknown type), its prepared state, its head. */
  readonly dataColumns: readonly { readonly renderer: string; readonly state: unknown; readonly label: string }[];
}

/** Content identities: one number per held-codes object, so an upsert's new cells change the stream's key. */
const stamps = new WeakMap<object, number>();
let stamped = 0;
const stampOf = (o: object) => {
  let s = stamps.get(o);
  if (s === undefined) {
    s = ++stamped;
    stamps.set(o, s);
  }
  return s;
};

export function prepareRows(input: LayerInput, rc: PrepareContext): RowsState {
  const { set, section, layer } = input;
  const a = set.alignments[section.blocks[0]?.alignment ?? -1];
  const track = section.tracks[layer.track];
  const lanes = track ? track.laneCount : 0;
  const start = track ? track.laneStart : 0;
  const lazy = !!a?.lazy;
  const run = section.lanes.lazy.find((z) => z.start === start);
  let rowOf: Int32Array | null = null;
  let slot: Int32Array | null = null;
  let depth: Int32Array | null = null;
  if (!run) {
    rowOf = Int32Array.from({ length: lanes }, (_, k) => section.lanes.row[start + k]);
    depth = Int32Array.from({ length: lanes }, (_, k) => section.lanes.depth[start + k]);
    slot = Int32Array.from(rowOf, (r) => (a ? a.rows.slotOf(r as never) : -1));
  }
  const remote = a ? ((a.custom.props.get(StreamProp) as StreamedData | undefined) ?? null) : null;
  const columns = a?.columnCount ?? 0;
  const codes = a && !lazy ? a.rows.codes : null;
  // The rows in lane order, the shape, and the content: the held cells' identity or the remote source's address. The
  // order as two independent 32-bit mixes (a base-31 sum collided on small row indices: a reorder kept the old tiles).
  let h1 = 0x811c9dc5 | 0;
  let h2 = lanes | 0;
  for (let k = 0; rowOf && k < rowOf.length; k++) {
    h1 = Math.imul(h1 ^ rowOf[k], 0x01000193);
    h2 = Math.imul(h2 ^ rowOf[k], 0x5bd1e995);
    h2 ^= h2 >>> 13;
  }
  const order = `${(h1 >>> 0).toString(36)}.${(h2 >>> 0).toString(36)}`;
  const src = codes ? stampOf(codes) : remote ? (remote.source.uri ?? remote.source.stream ?? remote.source.ref) : "";
  const key = `${a?.id ?? -1}:${lanes}:${columns}:${order}:${src}`;
  return {
    set,
    alignment: a?.id ?? -1,
    lazy,
    columns,
    lanes,
    rowOf,
    order: run?.order ?? null,
    slot,
    depth,
    codes,
    map: letterMapTable(input.params.letter_map as never),
    source: {
      key,
      rows: lanes,
      columns,
      remote,
      local: codes && slot ? { codes, slot } : null,
      systems: null,
    },
    tree: a?.tree && input.params.show_tree !== false && followsTree(input.params.order, true) ? treeOf(a.tree) : null,
    dataColumns: (
      (input.params.data_columns ?? []) as { data: string; field: string; label?: string; type?: string }[]
    ).map((c) => {
      const renderer = rc.registry.get("columns", String(c.type)) ? String(c.type) : "text";
      const t = { kind: -1, id: -1, row: -1 };
      const ann = set.byName(String(c.data), t) && t.kind === KeyKind.annotation ? set.annotations[t.id] : null;
      const name = set.strings.find(c.field);
      const field = ann?.fields.find((f) => f.name === name) ?? providerField(ann, set, c.field, rc.registry.svs);
      const rows = rowsOfAnnotation(ann, a?.rowCount ?? 0);
      const state = rc.registry.get("columns", renderer)?.prepare?.(field, rows, ann, set) ?? null;
      return { renderer, state, label: c.label ?? c.field };
    }),
  };
}

/** Per alignment row, the `rows` annotation's row for it; -1: none. */
function rowsOfAnnotation(ann: Annotation | null, n: number): Int32Array {
  const out = new Int32Array(n).fill(-1);
  for (let r = 0; ann && r < ann.count; r++) {
    const row = ann.target.row[r];
    if (row >= 0 && row < n) out[row] = r;
  }
  return out;
}
/** A field provider's values (`ns:name`) as a column's field, per annotation row; null: the name gives none. */
function providerField(ann: Annotation | null, set: SeqSet, name: string, svs: SvsRegistries): Field | null {
  if (!ann || ann.count === 0) return null;
  const data = Array.from({ length: ann.count }, (_, r) => paletteValue(ann, set, r, name, svs));
  return data.some((x) => x !== null && x !== undefined) ? { name: -1 as never, kind: "json", data } : null;
}

/** The rows entry's data by region (`stream`). */
export const rowsStream = (s: unknown): StreamLayer => (s as RowsState).source;

/** The rows header's headless half: its inset, its hit, the rows a node stands for. */
export const RowsHeader = {
  inset: (s: unknown) => ((s as RowsState).tree ? TREE_INSET : 0),
  hit: (v: AxisView, s: unknown, _p: unknown, out: Hit) => hitTree(v, s as RowsState, out),
  span: (s: unknown, node: number, out: Int32Array) => treeSpan(s as RowsState, node, out),
};

/** Per block, by `i − base`, a matrix annotation's row (-1: none); the states' letters and labels. */
export interface MatrixState {
  readonly rows: readonly (Int32Array | null)[];
  readonly ann: readonly (Annotation | null)[];
  /**
   * The bands: the states drawn (seqspec:format#keyed-sets), else the most states any block's annotation has (a block
   * with fewer has empty bands below its own).
   */
  readonly states: number;
  /** Per band, top to bottom, its state (an index into the data's states: `labels`, `codes`, the values' stride). */
  readonly stateOf: Int32Array;
  /** Per state, its band; -1: not drawn (a filter left it out). */
  readonly bandOf: Int32Array;
  /** How a reduced device row picks across its bands (the layer's `density`): 1 the largest value, 0 the most common color. */
  readonly maximum: number;
  /** Per block, its annotation's own state count: the stride of its values (`values[row · stride + state]`). */
  readonly stride: Int32Array;
  readonly codes: Uint8Array;
  /** Per state, its label in the header (`stateLabel`). */
  readonly labels: readonly string[];
  readonly lo: number;
  readonly hi: number;
  /** Logo: per block by `i − base`, its letters' heights (0..1 of the band) and their order bottom to top. */
  readonly heights: readonly (Float32Array | null)[];
  readonly order: readonly (Uint8Array | null)[];
}
export function prepareMatrix(input: LayerInput, logo = false): MatrixState {
  const { set, section, layer } = input;
  const rows: (Int32Array | null)[] = section.blocks.map(() => null);
  const anns: (Annotation | null)[] = section.blocks.map(() => null);
  let states: Int32Array = new Int32Array(0);
  const stride = new Int32Array(section.blocks.length);
  let lo = Number.POSITIVE_INFINITY;
  let hi = Number.NEGATIVE_INFINITY;
  for (const pl of layer.placements) {
    const a = set.annotations[pl.annotation];
    if (!a?.matrix || pl.lane >= 0) continue;
    const path = pl.path >= 0 ? layer.paths[pl.path] : null;
    const { length: n, base } = section.blocks[pl.block];
    const out = rows[pl.block] ?? new Int32Array(n).fill(-1);
    rows[pl.block] = out;
    anns[pl.block] = a;
    stride[pl.block] = a.matrix.states.length;
    if (a.matrix.states.length > states.length) states = a.matrix.states;
    for (let r = 0; r < a.count; r++) {
      if (a.slot[r] !== pl.slot) continue;
      const i = a.start[r];
      const j = (path ? Paths.at(path, i) : i) - base;
      if (j >= 0 && j < n && i >= 0) out[j] = r;
    }
    for (const v of a.matrix.values)
      if (!Number.isNaN(v)) {
        lo = Math.min(lo, v);
        hi = Math.max(hi, v);
      }
  }
  const S = states.length;
  const codes = Uint8Array.from(states, (id) => set.strings.get(id).charCodeAt(0) & 127);
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  const labels = Array.from(states, (id) => stateLabel(set, set.strings.get(id), t));
  const bits = input.params.scale !== "probability";
  const heights = rows.map((rr) => (logo && rr ? new Float32Array(rr.length * S) : null));
  const order = rows.map((rr) => (logo && rr ? new Uint8Array(rr.length * S) : null));
  const hmax = Math.log2(Math.max(2, S >= 20 ? 20 : S));
  const idx = Array.from({ length: S }, (_, k) => k);
  rows.forEach((rr, b) => {
    const h = heights[b];
    const o = order[b];
    const m = anns[b]?.matrix;
    const W = stride[b];
    for (let i = 0; h && o && m && rr && i < rr.length; i++) {
      const r = rr[i];
      if (r < 0) continue;
      let ent = 0;
      for (let s = 0; s < W; s++) {
        const p = m.values[r * W + s];
        if (p > 0) ent -= p * Math.log2(p);
      }
      const scale = bits ? Math.max(0, hmax - ent) / hmax : 1;
      for (let s = 0; s < W; s++) h[i * S + s] = Math.max(0, m.values[r * W + s] || 0) * scale;
      idx.sort((x, y) => h[i * S + x] - h[i * S + y]);
      for (let s = 0; s < S; s++) o[i * S + s] = idx[s];
    }
  });
  const range = { lo: lo <= hi ? lo : 0, hi: lo < hi ? hi : lo + 1 };
  const stateOf = layer.members ?? Int32Array.from({ length: S }, (_, k) => k);
  const bandOf = new Int32Array(S).fill(-1);
  stateOf.forEach((q, band) => {
    if (q < S) bandOf[q] = band;
  });
  const maximum = input.params?.density === "maximum" ? 1 : 0;
  return {
    rows,
    ann: anns,
    states: stateOf.length,
    stateOf,
    bandOf,
    maximum,
    stride,
    codes,
    labels,
    ...range,
    heights,
    order,
  };
}
/**
 * A matrix state's label: where it names a unit (a fingerprint's ligands and poses), that unit's label, or its key
 * where the label is only its name in its complex ("1iep/C", not "C"); else the state itself.
 */
function stateLabel(set: SeqSet, state: string, t: KeyTarget): string {
  if (!set.byName(state, t) || t.kind !== KeyKind.unit) return state;
  const label = set.strings.get(set.units.label[t.id]);
  return !label || state.endsWith(`/${label}`) ? state : label;
}
/** A matrix's states as bands over its track (top to bottom), labeled in the header. */
export const MatrixBands = {
  count: (s: unknown) => (s as MatrixState).states,
  label: (s: unknown, band: number) => (s as MatrixState).labels[(s as MatrixState).stateOf[band]] ?? "",
};

/** `seqv:map`: the tracks it draws in miniature (refs; empty: every other track of its section). */
export const prepareMap = (input: LayerInput) => ({
  tracks: ((input.params.options as { tracks?: string[] } | undefined)?.tracks ?? []) as readonly string[],
});
