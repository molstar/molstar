/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `interactions:fingerprint` renderer (seqspec:viewer#extension-points; seqspec:libraries#interactions), a registry
// module like the built-ins, on the renderer API only. `prepare` reads the `Fingerprint` of each block's matrix
// annotation (its runtime property) and keeps what the draw, the hit and the label read as typed arrays by `i − base`:
// each block's column per position, its contacted positions, the kinds a column holds in any row. A block without the
// property is drawn by the `matrix` renderer's tiles and hit and labeled as its bands. The draw is draw.ts; the integer
// geometry the draw and the hit share is geometry.ts.

import { SortedArray } from "@molstar/seq-core/model";
import {
  type AxisView,
  type Hit,
  HoverKind,
  type HoverState,
  type LabelContext,
  type LabelSink,
  type LayerInput,
  type LayerRenderer,
  type LegendSink,
  type MatrixState,
  MatrixTiles,
  prepareMatrix,
  sourceAt,
  type Theme,
} from "@molstar/seq-viewer/draw";
import type { Fingerprint } from "../fingerprint.ts";
import { fingerprintOf } from "../registry.ts";
import { KindMask, Kinds } from "../vocabulary.ts";
import { draw, overlay, tiny } from "./draw.ts";
import { compactRows, Geo, geo, LEAD, rowFrame, rowTop, STRIPES, stripeTop } from "./geometry.ts";

export const FINGERPRINT = "interactions:fingerprint";

/** The engine's cache per prepared state and skin: the fallback blocks' tiles. */
export interface FingerprintCache {
  readonly tiles: MatrixTiles;
}

/**
 * The prepared layer. Per-block lists are typed arrays by `i − base` (never by shown position) and empty, not null, for
 * a block with no fingerprint, so each per-block load is one shape.
 */
export interface FingerprintState {
  /** The provider's matrix as `matrix` prepares it: rows per position, a band (state) per ligand. */
  readonly matrix: MatrixState;
  /** Per block: the runtime property of its matrix annotation; null: none (the block is drawn as a matrix). */
  readonly fp: readonly (Fingerprint | null)[];
  /** Per block: the rows it draws, `min(states, its own rows)`; 0 without a fingerprint. */
  readonly rowsOf: Int32Array;
  /** Per block: its fingerprint's columns; its cells (rows × columns); its kinds per column in any row it draws. */
  readonly width: Int32Array;
  readonly cells: readonly Uint16Array[];
  readonly union: readonly Uint16Array[];
  /** Per block, by `i − base`: the fingerprint's column there; -1: none. */
  readonly colOf: readonly Int32Array[];
  /** Per block: the positions (`i − base`, ascending) that have a column with a contact, and those columns. */
  readonly contacted: readonly Int32Array[];
  readonly contactedCol: readonly Int32Array[];
  /** The kinds present in the rows drawn, in stripe order: a stripe each, at the same place in every row. */
  readonly kinds: Uint8Array;
  /** Kind → its stripe; -1: none. 256 long, so a mask's lead (255: none) is in bounds. */
  readonly stripeOf: Int8Array;
  /** The blocks with matrix rows and no fingerprint, as the matrix state `drawMatrixTiles` takes; null: none. */
  readonly fallback: MatrixState | null;
}

const NONE32 = new Int32Array(0);
const NONE16 = new Uint16Array(0);

export function prepareFingerprint(input: LayerInput): FingerprintState {
  const matrix = prepareMatrix(input);
  const S = matrix.states;
  const nb = matrix.rows.length;
  const fp: (Fingerprint | null)[] = [];
  const rowsOf = new Int32Array(nb);
  const width = new Int32Array(nb);
  const cells: Uint16Array[] = [];
  const union: Uint16Array[] = [];
  const colOf: Int32Array[] = [];
  const contacted: Int32Array[] = [];
  const contactedCol: Int32Array[] = [];
  let present = 0;
  for (let b = 0; b < nb; b++) {
    const rows = matrix.rows[b];
    const a = matrix.ann[b];
    const f = rows && a ? (fingerprintOf(a) ?? null) : null;
    fp.push(f);
    if (!f || !rows || !a) {
      cells.push(NONE16);
      union.push(NONE16);
      colOf.push(NONE32);
      contacted.push(NONE32);
      contactedCol.push(NONE32);
      continue;
    }
    const C = f.columns.length;
    const R = Math.min(S, f.rows.length);
    rowsOf[b] = R;
    width[b] = C;
    cells.push(f.cells);
    const u = new Uint16Array(C);
    for (let q = 0; q < R; q++) for (let c = 0; c < C; c++) u[c] |= f.cells[q * C + c];
    union.push(u);
    // The matrix's rows are the system's positions (its `index`): each one's column in the fingerprint.
    const col = Int32Array.from(rows, (r) => (r < 0 ? -1 : SortedArray.indexOf(f.columns, a.start[r])));
    colOf.push(col);
    const at: number[] = [];
    const cs: number[] = [];
    for (let j = 0; j < col.length; j++) {
      if (col[j] < 0 || u[col[j]] === 0) continue;
      at.push(j);
      cs.push(col[j]);
      present |= u[col[j]];
    }
    contacted.push(Int32Array.from(at));
    contactedCol.push(Int32Array.from(cs));
  }
  const kinds = Uint8Array.from(STRIPES.filter((k) => present & (1 << k)));
  const stripeOf = new Int8Array(256).fill(-1);
  kinds.forEach((k, j) => {
    stripeOf[k] = j;
  });
  const lacking = matrix.rows.some((rows, b) => rows !== null && fp[b] === null);
  const fallback = lacking ? { ...matrix, rows: matrix.rows.map((rows, b) => (fp[b] === null ? rows : null)) } : null;
  return { matrix, fp, rowsOf, width, cells, union, colOf, contacted, contactedCol, kinds, stripeOf, fallback };
}

/**
 * The row under the point is the hit's `state` (the position, already the hit's, is the system's residue); the stripe
 * under it, where the cell holds its kind, is its `part` (that kind), a compact row's the cell's first kind. A block
 * without the property gives the band, as `matrix` does. The geometry is the draw's (geometry.ts): O(1) a hit.
 */
function hit(v: AxisView, s: FingerprintState, _p: unknown, out: Hit): boolean {
  const S = s.matrix.states;
  const b = out.block;
  if (b < 0 || b >= s.fp.length || S === 0 || v.axis === null) return false;
  rowFrame(v, S);
  const a = geo[Geo.a];
  const h = geo[Geo.height];
  const gap = geo[Geo.gap];
  const yd = out.point[1] * v.dpr;
  if (h < 1 || !(yd >= a && yd < a + h)) return false;
  if (s.fp[b] === null) {
    if (s.fallback === null || s.fallback.rows[b] === null) return false;
    out.state = (((yd - a) * S) / h) | 0;
    return true;
  }
  // The last row whose top is at or above the point (a row's gap below it is its own).
  let q = (((yd - a) * S) / (h + gap)) | 0;
  if (q >= S) q = S - 1;
  while (q + 1 < S && rowTop(a, h, gap, S, q + 1) <= yd) q++;
  while (q > 0 && rowTop(a, h, gap, S, q) > yd) q--;
  if (q >= s.rowsOf[b]) return false;
  out.state = q;
  const j = out.index - v.axis.base[b];
  const cols = s.colOf[b];
  const col = j >= 0 && j < cols.length ? cols[j] : -1;
  if (col < 0) return true;
  const m = s.cells[b][q * s.width[b] + col];
  if (m === 0) return true;
  // In the row's gap (a compact row's too): the row, no kind.
  const y0 = rowTop(a, h, gap, S, q);
  const y1 = rowTop(a, h, gap, S, q + 1) - gap;
  if (y1 <= y0 || yd >= y1) return true;
  const K = s.kinds.length;
  if (compactRows(K)) {
    out.part = LEAD[m];
    return true;
  }
  const hh = y1 - y0;
  let t = (((yd - y0) * K) / hh) | 0;
  if (t >= K) t = K - 1;
  while (t + 1 < K && stripeTop(y0, hh, K, t + 1) <= yd) t++;
  while (t > 0 && stripeTop(y0, hh, K, t) > yd) t--;
  if (m & (1 << s.kinds[t])) out.part = s.kinds[t];
  return true;
}

/** A residue in its block's numbering (`T315`), else its 1-based index; a column on an alignment's axis. */
function positionText(h: HoverState, lc: LabelContext, g: number): string {
  const index = h.index;
  if (h.kind === HoverKind.column) return `column ${index + 1}`;
  const set = lc.set;
  const rv = lc.view.resolved;
  const mol = set.molecules[set.units.molecule[h.seq.unit]];
  const n = rv.blocks.numbering[g];
  const np = rv.blocks.numberingPath[g];
  const at = np >= 0 && rv.paths[np] ? sourceAt(rv.paths[np], index) : index;
  const letter = mol && !mol.element ? String.fromCharCode(mol.code[index]) : "";
  return `${letter}${n && at >= 0 && at < n.length ? lc.format(n, at as never) : index + 1}`;
}

/**
 * Over a fingerprint's cell: the ligand, the residue, the stripe's kind or else the cell's (`… · T315 · hbond`); over a
 * block without the property, the matrix's value at the state (`value`'s line, which keys on the type `matrix`).
 */
function label(h: HoverState, lc: LabelContext, out: LabelSink, s: FingerprintState): void {
  const b = h.block;
  if (b < 0 || b >= s.fp.length || h.state < 0) return;
  if (h.kind !== HoverKind.seq && h.kind !== HoverKind.column) return;
  const rv = lc.view.resolved;
  const g = rv.sections.blockStart[lc.section] + b;
  const j = h.index - rv.blocks.base[g];
  const ligand = s.matrix.labels[h.state] ?? "";
  if (s.fp[b] !== null) {
    if (h.state >= s.rowsOf[b]) return;
    const cols = s.colOf[b];
    const col = j >= 0 && j < cols.length ? cols[j] : -1;
    const mask = col >= 0 ? s.cells[b][h.state * s.width[b] + col] : 0;
    const kinds = h.part >= 0 ? Kinds[h.part] : mask ? KindMask.names(mask).join(", ") : "no contact";
    out.replaceTitle(`${ligand} · ${positionText(h, lc, g)} · ${kinds}`);
    out.line("ligand", ligand);
    return;
  }
  const rows = s.fallback?.rows[b];
  const m = s.matrix.ann[b]?.matrix;
  const row = rows && j >= 0 && j < rows.length ? rows[j] : -1;
  const stride = s.matrix.stride[b];
  const v = m && row >= 0 && h.state < stride ? m.values[row * stride + h.state] : Number.NaN;
  if (!Number.isNaN(v)) out.line(ligand, `${Number(v.toPrecision(4))}`);
}

const COUNTS = "contacts (number of kinds)";
/** The kinds present (the `color` channel's legend), and, if any block has no fingerprint, the counts' ramp, captioned. */
function legend(_p: unknown, themes: readonly (Theme | null)[], out: LegendSink, s: FingerprintState): void {
  themes[0]?.legend?.(out);
  if (s.fallback === null) return;
  themes[1]?.legend?.({
    item: (label, color, glyph) => out.item(label, color, glyph),
    ramp: (stops, lo, hi) => out.ramp(stops, lo, hi, COUNTS),
    size: (label, value, unit, of) => out.size(label, value, unit, of),
  });
}

/** The fingerprint entry: a band per ligand or pose, a stripe per kind. */
export const FingerprintRenderer: LayerRenderer<FingerprintState> = {
  name: FINGERPRINT,
  geometries: ["linear"],
  content: false,
  stable: false,
  minHeight: 4,
  prepare: prepareFingerprint,
  // A band per ligand with its 2 px gap, at least 2 px a kind present: 20 poses with 8 kinds keep every stripe.
  height: (s, _p, props) =>
    Math.max(24, s.matrix.states * Math.max(Math.min(18, props.rowHeight), 2 * s.kinds.length + 2)),
  // A row is never under a device pixel at any DPR: a track `height` hint is raised to a CSS px a row (and the borders).
  minTrackHeight: (s) => s.matrix.states + 2,
  hit,
  label,
  legend,
  bands: { count: (s) => s.matrix.states, label: (s, q) => s.matrix.labels[q] ?? "" },
  cache: { create: () => ({ tiles: new MatrixTiles() }), dispose: (c) => (c as FingerprintCache).tiles.dispose() },
  draw,
  tiny,
  overlay,
};
