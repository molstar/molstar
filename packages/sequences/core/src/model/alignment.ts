/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Alignments and trees (seqspec:model#alignments-and-trees): a column space, and per row a map between the columns and
// the row's own residues. Gapped strings are input syntax; rows are struct-of-arrays by slot.

import {
  type AlignmentId,
  type ColumnIndex,
  CustomData,
  Issue,
  type NameId,
  type RowIndex,
  type SeqIndex,
  type StringTable,
} from "./common.ts";
import { type Alphabet, Molecules, type ResidueCheck } from "./molecule.ts";
import { NO_LABEL, Numbering, type NumberingSet, type OffsetNumbering } from "./numbering.ts";

/** The code of a gap in `codes` (`-`). */
export const GAP = 45;
/** A checkpoint every 2^6 = 64 columns. */
const SHIFT = 6;
/** How many rows' `posMap` arrays are kept, least-recently-used evicted. */
const POS_CACHE = 256;

/** Rooted, struct-of-arrays over nodes in preorder (the root is 0; a clade is the nodes [node, cladeEnd[node])). */
export interface Tree {
  readonly count: number;
  /** -1: the root. */
  readonly parent: Int32Array;
  /** NaN: none. */
  readonly branchLength: Float64Array;
  /** RowIndex for a leaf, -1 inside. */
  readonly leafRow: Int32Array;
  /** NameId, -1: none. */
  readonly label: Int32Array;
  /** Rows in leaf order. */
  readonly leafOrder: Int32Array;
  /**
   * The end of each node's clade in preorder (an addition to seqspec:model#indices-and-labels: `cladeRows` reads it).
   */
  readonly cladeEnd: Int32Array;
  readonly custom: CustomData;
}

/** Struct-of-arrays over the rows held, by slot. Every operation on a row not held answers -1 or false. */
export interface AlignmentRows {
  readonly count: number;
  /** -1: not held. */
  slotOf(row: RowIndex): number;
  /** Each held row's own residues, inserts included, concatenated (upper case). */
  readonly residues: Uint8Array;
  /** count + 1. */
  readonly residueStart: Int32Array;
  /** count + 1, into `inserts` (an index into the array, two per insert). */
  readonly insertStart: Int32Array;
  /** (after column, length) pairs per row; after column -1: before the first column. */
  readonly inserts: Int32Array;
  /** The `row` numbering's first label; NO_LABEL: none. */
  readonly rowFirst: Int32Array;
  /** count × columnCount, GAP (45) at a gap. */
  readonly codes: Uint8Array;
  /** Residues before every 64th column, inserts counted: count × ceil(columnCount / 64). */
  readonly checkpoints: Int32Array;
  /** -1 at a gap, or a row not held. */
  colToPos(row: RowIndex, column: ColumnIndex): SeqIndex;
  /** The row's whole map into `out` (at least its length), -1 at an insert; false: not held. */
  posToCol(row: RowIndex, out: Int32Array): boolean;
  /** The row's own length; -1: not held (an addition to the spec, to size `posToCol`'s out). */
  lengthOf(row: RowIndex): number;
  /** The row's map, built once and cached up to a bound (LRU); null: not held. Allocates only on a miss. */
  posMap(row: RowIndex): Int32Array | null;
  /** The column of the row's residue `pos`; -1 in an insert, past its end, or a row not held. Allocation-free. */
  colOf(row: RowIndex, pos: number): number;
}

export interface Alignment {
  readonly id: AlignmentId;
  readonly key: NameId;
  readonly alphabet: Alphabet;
  readonly rowCount: number;
  readonly columnCount: number;
  /** 1 for an insert column; null: none. */
  readonly insertColumn: Uint8Array | null;
  /** Rows served by a progressive source. */
  readonly lazy: boolean;
  /** Always present; lazy: only the rows held. */
  readonly rows: AlignmentRows;
  /** A held row's label; -1 for a row not held (lazy) or unlabeled. */
  labelOf(row: RowIndex): NameId;
  /** Any row's label: a held row's, else the source's once served (never interned); null until then. */
  labelText(row: RowIndex): string | null;
  /** A held row's `ordinal` and `row` (from `rowFirst`, default when present); null: not held. Edges. */
  numberingsOf(row: RowIndex): NumberingSet | null;
  /** Per row: the row it is a variant of, or -1. */
  readonly parent: Int32Array;
  readonly referenceRow: RowIndex;
  readonly tree: Tree | null;
  readonly custom: CustomData;
}

/** A tree's input: any node order with one root; `leafRow` -1 names no row, -2 several (from `Trees.matchLeaves`). */
export interface TreeInput {
  parent: Int32Array;
  branchLength: Float64Array;
  label: (string | null)[];
  leafRow: Int32Array;
}
/** A parsed tree before its leaves are matched to rows. */
export type ParsedTree = Omit<TreeInput, "leafRow">;

export interface LazyAlignmentInput {
  rowCount: number;
  columnCount: number;
  /** The rows held, by row index: the label and the gapped string (in the alignment's flavor). */
  rows: Map<number, { label: string; gapped: string }>;
  /** null until the source has served it. */
  labelOf(row: number): string | null;
}

/** `gapped` is ignored when `lazy` is given. */
export interface AlignmentInput {
  key: string;
  gapped: string[];
  flavor: "aligned" | "a3m";
  labels?: string[];
  alphabet?: Alphabet;
  referenceRow?: number;
  parentRows?: (number | null)[];
  insertColumns?: Uint8Array;
  /** Per row; absent: from `NAME/27-137` headers in `labels`. */
  rowNumberingFirst?: (number | null)[];
  tree?: TreeInput;
  lazy?: LazyAlignmentInput;
  custom?: object;
  /**
   * Runtime properties (seqspec:model#custom-data), attached to the built alignment's `custom.props` (a lazy one's
   * progressive source).
   */
  props?: Record<string, unknown>;
}

interface RowData {
  count: number;
  held: Int32Array | null;
  columns: number;
  residues: Uint8Array;
  residueStart: Int32Array;
  insertStart: Int32Array;
  inserts: Int32Array;
  rowFirst: Int32Array;
  codes: Uint8Array;
  checkpoints: Int32Array;
}

class Rows implements AlignmentRows {
  readonly count: number;
  readonly residues: Uint8Array;
  readonly residueStart: Int32Array;
  readonly insertStart: Int32Array;
  readonly inserts: Int32Array;
  readonly rowFirst: Int32Array;
  readonly codes: Uint8Array;
  readonly checkpoints: Int32Array;
  /** The rows held, sorted (slot = its index); null: every row, slot = row. */
  private readonly held: Int32Array | null;
  private readonly columns: number;
  private readonly blocks: number;
  private readonly cacheSlot = new Map<number, number>();
  private readonly cacheRow = new Int32Array(POS_CACHE).fill(-1);
  private readonly cacheUse = new Float64Array(POS_CACHE);
  private readonly cacheMap: Int32Array[] = [];
  private clock = 0;

  constructor(d: RowData) {
    this.count = d.count;
    this.held = d.held;
    this.columns = d.columns;
    this.blocks = (d.columns + 63) >> SHIFT;
    this.residues = d.residues;
    this.residueStart = d.residueStart;
    this.insertStart = d.insertStart;
    this.inserts = d.inserts;
    this.rowFirst = d.rowFirst;
    this.codes = d.codes;
    this.checkpoints = d.checkpoints;
  }

  slotOf(row: RowIndex): number {
    const h = this.held;
    if (!h) return row >= 0 && row < this.count ? row : -1;
    let lo = 0;
    let hi = h.length - 1;
    while (lo <= hi) {
      const mid = (lo + hi) >> 1;
      if (h[mid] === row) return mid;
      if (h[mid] < row) lo = mid + 1;
      else hi = mid - 1;
    }
    return -1;
  }

  lengthOf(row: RowIndex): number {
    const s = this.slotOf(row);
    return s < 0 ? -1 : this.residueStart[s + 1] - this.residueStart[s];
  }

  colToPos(row: RowIndex, column: ColumnIndex): SeqIndex {
    const s = this.slotOf(row);
    const C = this.columns;
    if (s < 0 || column < 0 || column >= C) return -1 as SeqIndex;
    const base = s * C;
    const codes = this.codes;
    if (codes[base + column] === GAP) return -1 as SeqIndex;
    const from = (column >> SHIFT) << SHIFT;
    let pos = this.checkpoints[s * this.blocks + (column >> SHIFT)];
    for (let c = from; c < column; c++) if (codes[base + c] !== GAP) pos++;
    // Inserts after a column in [from, column): the first by binary search over the row's pairs.
    const ins = this.inserts;
    let lo = this.insertStart[s] >> 1;
    let hi = this.insertStart[s + 1] >> 1;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (ins[2 * mid] < from) lo = mid + 1;
      else hi = mid;
    }
    const end = this.insertStart[s + 1];
    for (let k = 2 * lo; k < end && ins[k] < column; k += 2) pos += ins[k + 1];
    return pos as SeqIndex;
  }

  posToCol(row: RowIndex, out: Int32Array): boolean {
    const s = this.slotOf(row);
    if (s < 0) return false;
    const C = this.columns;
    const base = s * C;
    const ins = this.inserts;
    const end = this.insertStart[s + 1];
    let k = this.insertStart[s];
    let pos = 0;
    for (let c = 0; c <= C; c++) {
      for (; k < end && ins[k] < c; k += 2) for (let j = 0; j < ins[k + 1]; j++) out[pos++] = -1;
      if (c < C && this.codes[base + c] !== GAP) out[pos++] = c;
    }
    return true;
  }

  colOf(row: RowIndex, pos: number): number {
    const s = this.slotOf(row);
    if (s < 0 || pos < 0) return -1;
    // The last checkpoint at or below `pos`, then the block's columns and the inserts after them, as `colToPos` walks.
    const cp = this.checkpoints;
    const b0 = s * this.blocks;
    let lo = 0;
    let hi = this.blocks - 1;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (cp[b0 + mid] <= pos) lo = mid;
      else hi = mid - 1;
    }
    const from = lo << SHIFT;
    const ins = this.inserts;
    let a = this.insertStart[s] >> 1;
    let z = this.insertStart[s + 1] >> 1;
    while (a < z) {
      const mid = (a + z) >> 1;
      if (ins[2 * mid] < from) a = mid + 1;
      else z = mid;
    }
    const end = this.insertStart[s + 1];
    const base = s * this.columns;
    let p = cp[b0 + lo];
    for (let c = from, k = 2 * a; c < this.columns && p <= pos; c++) {
      if (this.codes[base + c] !== GAP && p++ === pos) return c;
      for (; k < end && ins[k] === c; k += 2) p += ins[k + 1];
    }
    return -1;
  }

  posMap(row: RowIndex): Int32Array | null {
    const at = this.cacheSlot.get(row);
    if (at !== undefined) {
      this.cacheUse[at] = ++this.clock;
      return this.cacheMap[at];
    }
    const n = this.lengthOf(row);
    if (n < 0) return null;
    let victim = 0;
    for (let i = 1; i < POS_CACHE; i++) if (this.cacheUse[i] < this.cacheUse[victim]) victim = i;
    if (this.cacheRow[victim] >= 0) this.cacheSlot.delete(this.cacheRow[victim]);
    const map = new Int32Array(n);
    this.posToCol(row, map);
    this.cacheRow[victim] = row;
    this.cacheUse[victim] = ++this.clock;
    this.cacheMap[victim] = map;
    this.cacheSlot.set(row, victim);
    return map;
  }
}

/** Columns of a gapped row: every character but whitespace, and in A3M but lower case and `.` (inserts). */
function columnsOf(g: string, a3m: boolean): number {
  let n = 0;
  for (let i = 0; i < g.length; i++) {
    const c = g.charCodeAt(i);
    if (c > 32 && !(a3m && (c === 46 || (c >= 97 && c <= 122)))) n++;
  }
  return n;
}

/** The index of the first character outside ASCII, or -1: rows are never truncated to a byte. */
function nonAscii(g: string): number {
  for (let i = 0; i < g.length; i++) if (g.charCodeAt(i) > 127) return i;
  return -1;
}

/** Codes, residues, inserts and checkpoints of rows in slot order; `lower`/`upper` flag columns by case. */
function fillRows(gapped: readonly string[], a3m: boolean, C: number, lower: Uint8Array, upper: Uint8Array) {
  const count = gapped.length;
  const blocks = (C + 63) >> SHIFT;
  let total = 0;
  // Index loops, here and for the labels: `for…of` and `Int32Array.from` with a map allocate per row until optimized.
  for (let s = 0; s < count; s++) total += gapped[s].length;
  const residues = new Uint8Array(total);
  const residueStart = new Int32Array(count + 1);
  const insertStart = new Int32Array(count + 1);
  const inserts: number[] = [];
  const codes = new Uint8Array(count * C);
  const checkpoints = new Int32Array(count * blocks);
  let r = 0;
  for (let s = 0; s < count; s++) {
    const g = gapped[s];
    const base = s * C;
    residueStart[s] = r;
    insertStart[s] = inserts.length;
    let col = 0;
    let pos = 0;
    let open = -1;
    for (let i = 0; i < g.length; i++) {
      const c = g.charCodeAt(i);
      if (c <= 32) continue;
      const low = c >= 97 && c <= 122;
      if (a3m && (low || c === 46)) {
        if (c === 46) continue;
        if (open < 0) {
          inserts.push(col - 1, 0);
          open = inserts.length - 1;
        }
        inserts[open]++;
        residues[r++] = c - 32;
        pos++;
        continue;
      }
      open = -1;
      if ((col & 63) === 0) checkpoints[s * blocks + (col >> SHIFT)] = pos;
      if (c === 45 || c === 46) codes[base + col] = GAP;
      else {
        const u = low ? c - 32 : c;
        codes[base + col] = u;
        residues[r++] = u;
        pos++;
        if (low) lower[col] = 1;
        else upper[col] = 1;
      }
      col++;
    }
  }
  residueStart[count] = r;
  insertStart[count] = inserts.length;
  return {
    residues: residues.slice(0, r),
    residueStart,
    insertStart,
    inserts: Int32Array.from(inserts),
    codes,
    checkpoints,
  };
}

const HEADER = /\/(-?[0-9]+)-(-?[0-9]+)$/;

/** An alignment from its input, or null with errors. */
function build(
  input: AlignmentInput,
  id: number,
  strings: StringTable,
  issues: Issue[],
  path: string,
): Alignment | null {
  const err = (code: string, at: string, m: string) => {
    issues.push(Issue.error(code, `${path}${at}`, `Alignment "${input.key}": ${m}`));
    return null;
  };
  const a3m = input.flavor === "a3m";
  const lazy = input.lazy ?? null;
  let held: Int32Array | null = null;
  let gapped: readonly string[] = input.gapped;
  let heldLabels: readonly (string | null)[] = input.labels ?? [];
  let rowCount = gapped.length;
  let C = rowCount ? columnsOf(gapped[0], a3m) : 0;
  if (lazy) {
    rowCount = lazy.rowCount;
    C = lazy.columnCount;
    if (!(Number.isInteger(rowCount) && Number.isInteger(C) && rowCount > 0 && C > 0)) {
      return err("model.alignment.shape", "/lazy", "a lazy alignment needs positive row and column counts.");
    }
    const keys = [...lazy.rows.keys()];
    if (keys.some((r) => !(Number.isInteger(r) && r >= 0 && r < rowCount))) {
      return err("model.alignment.shape", "/lazy/rows", `held rows must be within 0..${rowCount - 1}.`);
    }
    held = Int32Array.from(keys).sort();
    const rows = Array.from(held, (r) => lazy.rows.get(r) as { label: string; gapped: string });
    gapped = rows.map((x) => x.gapped);
    heldLabels = rows.map((x) => x.label);
  } else if (rowCount === 0 || C === 0) return err("model.alignment.shape", "/gapped", "it has no rows or no columns.");
  // A row's path is made on the error path only: a string per row is 20k strings at 20k rows.
  const rowAt = (s: number) => (held ? `/lazy/rows/${held[s]}` : `/gapped/${s}`);
  for (let s = 0; s < gapped.length; s++) {
    const n = columnsOf(gapped[s], a3m);
    const row = held ? held[s] : s;
    if (n !== C) return err("model.alignment.ragged", rowAt(s), `row ${row} has ${n} columns, not ${C}.`);
    const bad = nonAscii(gapped[s]);
    if (bad >= 0) return err("model.residues.invalid", rowAt(s), `row ${row} has a non-ASCII character at ${bad}.`);
  }
  const perRow = {
    labels: lazy ? null : input.labels,
    parentRows: input.parentRows,
    rowNumberingFirst: input.rowNumberingFirst,
  };
  for (const [field, v] of Object.entries(perRow)) {
    if (v && v.length !== rowCount)
      return err("model.alignment.shape", `/${field}`, `${v.length} ${field} for ${rowCount} rows.`);
  }
  if (input.insertColumns && input.insertColumns.length !== C) {
    return err("model.alignment.shape", "/insertColumns", `${input.insertColumns.length} flags for ${C} columns.`);
  }

  const lower = new Uint8Array(C);
  const upper = new Uint8Array(C);
  const d = fillRows(gapped, a3m, C, lower, upper);
  // One counting pass over the residues detects and validates; the bad letter's index is looked up only on failure.
  const check: ResidueCheck = { at: -1, invalid: -1 };
  const alphabet = Molecules.checkResidues(d.residues, input.alphabet, check);
  if (check.at >= 0) {
    const m = `"${String.fromCharCode(d.residues[check.at])}" is no nucleotide code, so the rows are read as protein.`;
    issues.push(Issue.warning("model.alphabet.detected", path, `Alignment "${input.key}": ${m}`));
  }
  if (check.invalid >= 0) {
    let s = 0;
    while (d.residueStart[s + 1] <= check.invalid) s++;
    const c = String.fromCharCode(d.residues[check.invalid]);
    return err("model.residues.invalid", rowAt(s), `row ${held ? held[s] : s}: "${c}" is not a ${alphabet} letter.`);
  }

  // Insert columns: given (Stockholm's RF), else in aligned input the columns whose residues are all lower case.
  // Input that is lower case throughout is lower-case letters, not inserts.
  let insertColumn = input.insertColumns ?? null;
  let inserts = 0;
  let anyLower = false;
  let withResidues = 0;
  if (!insertColumn && !a3m && !lazy) {
    const flags = new Uint8Array(C);
    for (let c = 0; c < C; c++) {
      if (lower[c] || upper[c]) withResidues++;
      if (lower[c] && !upper[c]) flags[c] = 1;
      inserts += flags[c];
      anyLower ||= lower[c] === 1;
    }
    if (inserts > 0 && inserts < withResidues) insertColumn = flags;
    else inserts = 0;
  } else if (insertColumn) for (let c = 0; c < C; c++) inserts += insertColumn[c] ? 1 : 0;
  const what = a3m
    ? d.inserts.length > 0 && `${d.inserts.length / 2} lower-case inserts are not columns (A3M)`
    : input.insertColumns
      ? `${inserts} insert columns are given (Stockholm's RF)`
      : insertColumn
        ? `${inserts} columns with only lower-case residues are insert columns`
        : anyLower && "lower case marks no insert columns";
  if (what) {
    issues.push(
      Issue.info("model.alignment.inserts", path, `Alignment "${input.key}": ${what}; letters are upper case.`),
    );
  }

  // The row numbering: given, else from a `NAME/27-137` header; a descending header gives none.
  const rowFirst = new Int32Array(gapped.length).fill(NO_LABEL);
  let descending = 0;
  for (let s = 0; s < gapped.length; s++) {
    const row = held ? held[s] : s;
    if (input.rowNumberingFirst) rowFirst[s] = input.rowNumberingFirst[row] ?? NO_LABEL;
    else {
      const m = HEADER.exec(heldLabels[s] ?? "");
      if (m && Number(m[1]) <= Number(m[2])) rowFirst[s] = Number(m[1]);
      else if (m) descending++;
    }
  }
  if (descending) {
    const m = `${descending} descending row headers (NAME/137-27) give no row numbering.`;
    issues.push(Issue.warning("model.numbering.descending", path, `Alignment "${input.key}": ${m}`));
  }

  // A bad parent, or the parent that closes a cycle, is a warning and that row gets none.
  const parent = new Int32Array(rowCount).fill(-1);
  const parentRows = input.parentRows ?? [];
  // The first ten said, then a count (GV-07: a parent list off by one gave one issue per row).
  let bad = 0;
  const noParent = (r: number, m: string) => {
    parent[r] = -1;
    const message = `Alignment "${input.key}": row ${r}'s parent ${m}; it has none.`;
    if (bad++ < 10) issues.push(Issue.warning("model.alignment.parent", `${path}/parentRows/${r}`, message));
  };
  for (let r = 0; r < parentRows.length; r++) {
    const p = parentRows[r];
    if (p == null) continue;
    if (!(Number.isInteger(p) && p >= 0 && p < rowCount && p !== r)) noParent(r, `${p} is not another row`);
    else parent[r] = p;
  }
  // One walk up from each row, each row visited once: `seen` holds the walk that reached it.
  const seen = new Int32Array(parentRows.length ? rowCount : 0);
  for (let r = 0; r < seen.length; r++) {
    let x = r;
    while (!seen[x]) {
      seen[x] = r + 1;
      if (parent[x] < 0) break;
      if (seen[parent[x]] === r + 1) noParent(x, `${parent[x]} closes a cycle`);
      else x = parent[x];
    }
  }
  if (bad > 10) {
    const message = `Alignment "${input.key}": ${bad - 10} more rows' parents left out.`;
    issues.push(Issue.warning("model.alignment.parent", `${path}/parentRows`, message));
  }
  const referenceRow = input.referenceRow ?? 0;
  if (!(Number.isInteger(referenceRow) && referenceRow >= 0 && referenceRow < rowCount)) {
    return err("model.alignment.shape", "/referenceRow", `reference row ${referenceRow} is not a row.`);
  }

  const rows = new Rows({ count: gapped.length, held, columns: C, rowFirst, ...d });
  const labelIds = new Int32Array(heldLabels.length);
  for (let s = 0; s < labelIds.length; s++) {
    const l = heldLabels[s];
    labelIds[s] = l == null ? -1 : strings.intern(l);
  }
  const labelOf = (row: RowIndex) => {
    const s = rows.slotOf(row);
    return (s >= 0 ? labelIds[s] : -1) as NameId;
  };
  const labelText = (row: RowIndex) => {
    const s = rows.slotOf(row);
    if (s >= 0) return labelIds[s] < 0 ? null : strings.get(labelIds[s]);
    return lazy && row >= 0 && row < rowCount ? lazy.labelOf(row) : null;
  };
  const rowNumberings = new Map<number, NumberingSet>();
  const offset = (name: string, first: number, length: number): OffsetNumbering => {
    const n = { name: strings.intern(name), isDefault: name === "row", custom: CustomData.create(), length };
    return { ...n, method: "offset", first, skipZero: false };
  };
  const numberingsOf = (row: RowIndex) => {
    const s = rows.slotOf(row);
    if (s < 0) return null;
    let set = rowNumberings.get(s);
    if (!set) {
      const n = rows.lengthOf(row);
      const own = rowFirst[s] === NO_LABEL ? [] : [offset("row", rowFirst[s], n)];
      set = Numbering.setOf([offset("ordinal", 1, n), ...own]);
      rowNumberings.set(s, set);
    }
    return set;
  };
  return {
    id: id as AlignmentId,
    key: strings.intern(input.key),
    numberingsOf,
    alphabet,
    rowCount,
    columnCount: C,
    insertColumn,
    lazy: !!lazy,
    rows,
    labelOf,
    labelText,
    parent,
    referenceRow: referenceRow as RowIndex,
    tree: input.tree ? buildTree(input.tree, rowCount, strings, issues, `${path}/tree`) : null,
    custom: CustomData.create(input.custom, input.props),
  };
}

/** A tree in preorder (children in input order), or null with errors: one root, every row exactly one leaf. */
function buildTree(input: TreeInput, rowCount: number, strings: StringTable, issues: Issue[], path: string) {
  const err = (code: string, m: string) => {
    issues.push(Issue.error(code, path, `Tree: ${m}`));
    return null;
  };
  const { parent: inParent } = input;
  const n = inParent.length;
  if (input.branchLength.length !== n || input.label.length !== n || input.leafRow.length !== n) {
    return err("model.tree.shape", `its columns must all have ${n} nodes.`);
  }
  const childStart = new Int32Array(n + 1);
  let root = -1;
  for (let i = 0; i < n; i++) {
    const p = inParent[i];
    if (p === -1 && root >= 0) return err("model.tree.shape", "it has more than one root.");
    if (p === -1) root = i;
    else if (p < 0 || p >= n || p === i) return err("model.tree.shape", `node ${i}'s parent ${p} is not a node.`);
    else childStart[p + 1]++;
  }
  if (root < 0) return err("model.tree.shape", "it has no root.");
  for (let i = 0; i < n; i++) childStart[i + 1] += childStart[i];
  const children = new Int32Array(n);
  const fill = childStart.slice();
  for (let i = 0; i < n; i++) if (inParent[i] >= 0) children[fill[inParent[i]]++] = i;

  const order = new Int32Array(n);
  const stack = [root];
  let k = 0;
  while (stack.length) {
    const v = stack.pop() as number;
    order[k++] = v;
    for (let j = childStart[v + 1] - 1; j >= childStart[v]; j--) stack.push(children[j]);
  }
  if (k !== n) return err("model.tree.shape", "some nodes do not descend from the root.");
  const newId = new Int32Array(n);
  for (let i = 0; i < n; i++) newId[order[i]] = i;

  const parent = Int32Array.from(order, (o) => (inParent[o] < 0 ? -1 : newId[inParent[o]]));
  const cladeEnd = Int32Array.from(order, (_, i) => i + 1);
  for (let i = n - 1; i > 0; i--) cladeEnd[parent[i]] = Math.max(cladeEnd[parent[i]], cladeEnd[i]);
  const leafRow = new Int32Array(n).fill(-1);
  const seen = new Uint8Array(rowCount);
  const leafOrder: number[] = [];
  const bad: string[] = [];
  for (let i = 0; i < n; i++) {
    const r = input.leafRow[order[i]];
    const name = `"${input.label[order[i]] ?? ""}"`;
    if (cladeEnd[i] !== i + 1) {
      if (r !== -1) bad.push(`the inner node ${name} names a row`);
    } else if (r === -2) bad.push(`the leaf ${name} names several rows`);
    else if (r < 0 || r >= rowCount) bad.push(`the leaf ${name} names no row`);
    else if (seen[r]) bad.push(`row ${r} has several leaves`);
    else {
      seen[r] = 1;
      leafRow[i] = r;
      leafOrder.push(r);
    }
  }
  if (leafOrder.length < rowCount) bad.push(`${rowCount - leafOrder.length} rows have no leaf`);
  if (bad.length) {
    const more = bad.length > 5 ? `, and ${bad.length - 5} more` : "";
    return err("model.tree.leaves", `every row must be exactly one leaf: ${bad.slice(0, 5).join("; ")}${more}.`);
  }
  return {
    count: n,
    parent,
    branchLength: Float64Array.from(order, (o) => input.branchLength[o]),
    leafRow,
    label: Int32Array.from(order, (o) => (input.label[o] == null ? -1 : strings.intern(input.label[o] as string))),
    leafOrder: Int32Array.from(leafOrder),
    cladeEnd,
    custom: CustomData.create(),
  } satisfies Tree;
}

/**
 * Leaf rows for a parsed tree: by `order`, leaf k (in node order) is row k; by `labels`, a leaf names the row
 * whose label equals it, else whose label before its last `/` does. -1: no row, -2: several (reported by build).
 */
function matchLeaves(tree: ParsedTree, rowLabels: readonly (string | null)[], by: "labels" | "order"): Int32Array {
  const n = tree.parent.length;
  const inner = new Uint8Array(n);
  for (let i = 0; i < n; i++) if (tree.parent[i] >= 0) inner[tree.parent[i]] = 1;
  const out = new Int32Array(n).fill(-1);
  const exact = new Map<string, number>();
  const prefix = new Map<string, number>();
  const add = (m: Map<string, number>, key: string, row: number) => m.set(key, m.has(key) ? -2 : row);
  if (by === "labels") {
    rowLabels.forEach((l, r) => {
      if (l == null) return;
      add(exact, l, r);
      const slash = l.lastIndexOf("/");
      if (slash > 0) add(prefix, l.slice(0, slash), r);
    });
  }
  let k = 0;
  for (let i = 0; i < n; i++) {
    if (inner[i]) continue;
    const name = tree.label[i] ?? "";
    out[i] = by === "order" ? k++ : (exact.get(name) ?? prefix.get(name) ?? -1);
  }
  return out;
}

/** Visits a clade's rows (in leaf order) without allocating. */
function cladeRows<C>(tree: Tree, node: number, cb: (row: RowIndex, ctx: C) => void, ctx: C): void {
  const end = tree.cladeEnd[node];
  for (let i = node; i < end; i++) if (tree.leafRow[i] >= 0) cb(tree.leafRow[i] as RowIndex, ctx);
}

export const Alignments = { build };
export const Trees = { build: buildTree, matchLeaves, cladeRows };
