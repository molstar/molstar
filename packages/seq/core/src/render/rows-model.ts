/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { TreeRead } from '@molstar/seqx-core/data/alignment-readers';
import type { MeasureEntry, ThemeDatum } from '@molstar/seqx-core/theme/index';

// The rows' model (_spec/README.md "Render", The rows layer): what the rows layer and its readers (the viewer's
// commands, a navigator) compute from an alignment without drawing: its clades (a tree's nodes over its rows), the
// display rows (an order over the rows, folded clades as one row, pinned rows in a band at the top, hidden rows out),
// the column measures registered for themes, a clade's summary and the profile. Pure; arrays sized by the row count are
// made when the order, the folds or the data change, never a frame.

/** A row's height when its band gives none (px: a letter row). Fixed: the letters' step needs 20 px. */
export const RowHeight = 20;

/** How near a tree's branch a point hits its node (px; decided). */
export const TreeHitPx = 6;

/** A gap byte: `-` or `.`. */
export function isGap(code: number): boolean {
  return code === 45 || code === 46 || code === 0;
}

/** A letter byte's uppercase (a lowercase insert reads as its letter). */
export function upper(code: number): number {
  return code >= 97 && code <= 122 ? code - 32 : code;
}

/** An alignment held whole: its rows × columns as bytes, row-major, its labels and first residue numbers. */
export interface HeldCells {
  readonly rows: number;
  readonly columns: number;
  readonly codes: Uint8Array;
  readonly labels: readonly string[];
  readonly first: Int32Array;
  /** Moves with the data it was joined from: what the measures, the summaries and the pyramid are kept by. */
  readonly version: number;
}

/** A tree's nodes over the rows: each node's leaves a run of the tree's leaf order (`lo`, `hi`). */
export interface Clades {
  readonly nodes: number;
  readonly parent: Int32Array;
  /** Each node's run of leaves `[lo, hi)` in the tree's leaf order. */
  readonly lo: Int32Array;
  readonly hi: Int32Array;
  readonly depth: Int32Array;
  /** Each node's distance from the root by branch length, and its height (the farthest leaf below it). */
  readonly distance: Float64Array;
  readonly height: Float64Array;
  /**
   * Each node's name: its label where unique in the tree and not a number, else its first and last leaves' labels as
   * `first..last` (`(clade :of [a b])`).
   */
  readonly name: readonly string[];
  /** The row each leaf (in leaf order) names, -1 for none; and the leaf of each row, -1 for none. */
  readonly leafRow: Int32Array;
  readonly rowLeaf: Int32Array;
  /** Rows in tree order: the leaves' rows, then the rows no leaf names, in data order. */
  readonly order: Int32Array;
  /** Children of each node, as the tree reads them. */
  readonly children: readonly (readonly number[])[];
  /** Each name's node (the first of a name): a fold's lookup. */
  readonly byName: ReadonlyMap<string, number>;
}

/** The clades of a tree over `rows` labels (leaves matched by `matchLeaves`); undefined for an empty tree. */
export function buildClades(tree: TreeRead, leafRow: Int32Array, rows: number): Clades | undefined {
  const nodes = tree.parent.length;
  if (!nodes) return undefined;
  const lo = new Int32Array(nodes).fill(0x7fffffff);
  const hi = new Int32Array(nodes);
  const depth = new Int32Array(nodes);
  const distance = new Float64Array(nodes);
  const height = new Float64Array(nodes);
  const leafAt = new Int32Array(nodes).fill(-1);
  tree.leaves.forEach((n, k) => (leafAt[n] = k));
  // pre-order: a parent before its children
  for (let n = 1; n < nodes; n++) {
    const p = tree.parent[n];
    depth[n] = depth[p] + 1;
    distance[n] = distance[p] + tree.length[n];
  }
  for (let n = nodes - 1; n >= 0; n--) {
    if (leafAt[n] >= 0) {
      lo[n] = leafAt[n];
      hi[n] = leafAt[n] + 1;
    }
    const p = tree.parent[n];
    if (p >= 0) {
      if (lo[n] < lo[p]) lo[p] = lo[n];
      if (hi[n] > hi[p]) hi[p] = hi[n];
      const h = height[n] + tree.length[n];
      if (h > height[p]) height[p] = h;
    }
  }
  const rowLeaf = new Int32Array(rows).fill(-1);
  for (let k = 0; k < leafRow.length; k++) if (leafRow[k] >= 0) rowLeaf[leafRow[k]] = k;
  const order = new Int32Array(rows);
  let at = 0;
  for (let k = 0; k < leafRow.length; k++) if (leafRow[k] >= 0) order[at++] = leafRow[k];
  for (let r = 0; r < rows; r++) if (rowLeaf[r] < 0) order[at++] = r;
  const leafLabel = (k: number) => tree.label[tree.leaves[k]] ?? '';
  // a node's id names its clade where it is unique in the tree and not a number (a support value, `((A,B)95,(C,D)95)`,
  // names no clade): else its first and last leaves (_spec/README.md "Render", The rows layer, Clades)
  const uses = new Map<string, number>();
  for (let n = 0; n < nodes; n++) {
    const l = tree.label[n];
    if (l) uses.set(l, (uses.get(l) ?? 0) + 1);
  }
  const name: string[] = [];
  const byName = new Map<string, number>();
  for (let n = 0; n < nodes; n++) {
    const l = tree.label[n];
    const own = l && leafAt[n] < 0 && uses.get(l) === 1 && !isNumeric(l);
    const s = own ? l : `${leafLabel(lo[n])}..${leafLabel(hi[n] - 1)}`;
    name.push(s);
    if (!byName.has(s)) byName.set(s, n);
  }
  return {
    nodes,
    parent: tree.parent,
    lo,
    hi,
    depth,
    distance,
    height,
    name,
    leafRow,
    rowLeaf,
    order,
    children: tree.children,
    byName,
  };
}

/**
 * A label that reads as support values (a bootstrap count, IQ-TREE's `87.3/95`: every `/`- or `|`-separated part a
 * number): never a clade's name, since a reinferred tree changes it.
 */
function isNumeric(label: string): boolean {
  return label.split(/[/|]/).every((part) => part.trim() !== '' && Number.isFinite(Number(part)));
}

/** A clade by its name (`cladeNamed`), -1 for none. */
export function cladeNamed(clades: Clades, name: string): number {
  return clades.byName.get(name) ?? -1;
}

/**
 * The display rows: what each display row shows (a data row, or `-1 − clade` for a folded clade's summary row), and
 * where each data row is shown (its own display row, its folded clade's, or -1 when hidden). Pinned rows first, in a
 * band of their own (`pinned`), and kept in their place too (`inPlace`). Data order with nothing pinned, hidden or
 * folded is arithmetic (`identity`): no arrays, so a million lazy rows cost what ten do.
 */
export class DisplayRows {
  /** Display rows; the first `pinned` of them are the pinned band. */
  count = 0;
  pinned = 0;
  identity = true;
  /** Rows a hidden-rows chip counts. */
  hidden = 0;
  entry: Int32Array = new Int32Array(0);
  displayOf: Int32Array = new Int32Array(0);
  /** Each data row's place in the order (a pinned row's in-place row, else its display row): the tree's, the emphasis'. */
  placeOf: Int32Array = new Int32Array(0);
  /** The rows in the order the display was built in (data rows; undefined: data order): a hidden row's neighbours. */
  ordered: Int32Array | undefined;
  /** Data rows in the order. */
  rows = 0;
  /** Each folded clade shown (by display row of its summary): its shown rows, in display order. */
  readonly members = new Map<number, Int32Array>();
  /** Moves with every rebuild: what a layer's kept marks (summaries, the pyramid in display order) are keyed by. */
  version = 0;

  /** The data row a display row shows, or `-1 − clade` for a summary row; -1 − nothing past the end is never asked. */
  rowAt(i: number): number {
    return this.identity ? i : this.entry[i];
  }

  /** The display row a data row is shown at (its folded clade's summary row), -1 when hidden or filtered. */
  displayRow(row: number): number {
    return this.identity ? row : (this.displayOf[row] ?? -1);
  }

  /**
   * Where a data row is in the order: a pinned row kept in place at its in-place row (its band's copy is `displayRow`),
   * else as `displayRow`. What the tree, the hover and the selection read.
   */
  placeRow(row: number): number {
    return this.identity ? row : (this.placeOf[row] ?? -1);
  }

  /**
   * The display row a data row stands at, clamped (a kept window's top row): its own place, its folded clade's summary
   * row, else (hidden) the nearest row shown in the order (the one after it on a tie); -1 when nothing is shown.
   */
  nearestRow(row: number): number {
    if (this.identity) return this.count ? Math.max(0, Math.min(this.count - 1, row)) : -1;
    if (row < 0 || row >= this.rows) return -1;
    const own = this.placeOf[row];
    if (own >= 0) return own;
    const order = this.ordered;
    let at = row;
    if (order) for (let i = 0; i < this.rows; i++) if (order[i] === row) at = i;
    const placeAt = (i: number) => this.placeOf[order ? order[i] : i];
    for (let d = 1; d < this.rows; d++) {
      const a = at + d < this.rows ? placeAt(at + d) : -1;
      if (a >= 0) return a;
      const b = at - d >= 0 ? placeAt(at - d) : -1;
      if (b >= 0) return b;
    }
    return -1;
  }

  /** Data order over `n` rows: no arrays. */
  setIdentity(n: number): void {
    this.identity = true;
    this.ordered = undefined;
    this.rows = n;
    this.count = n;
    this.pinned = 0;
    this.hidden = 0;
    this.members.clear();
    this.version++;
  }

  /**
   * Rows in `order` (data rows in display order), clades of `clades` folded (`folds`, outermost wins; only in tree
   * order), `pinned` rows (data rows) in a band at the top, `hidden` rows left out (a clade's "+n" counts shown rows;
   * a clade with none shown goes, one with one folds to that row).
   */
  set(
    n: number,
    order: Int32Array | undefined,
    clades: Clades | undefined,
    folds: readonly number[],
    pinned: readonly number[],
    hidden: ReadonlySet<number>,
    inPlace: boolean,
  ): void {
    this.identity = false;
    this.members.clear();
    this.version++;
    if (this.entry.length < n + pinned.length) this.entry = new Int32Array(n + pinned.length);
    if (this.displayOf.length < n) {
      this.displayOf = new Int32Array(n);
      this.placeOf = new Int32Array(n);
    }
    this.displayOf.fill(-1, 0, n);
    this.placeOf.fill(-1, 0, n);
    this.rows = n;
    const place = this.placeOf;
    const entry = this.entry;
    let at = 0;
    const pin = new Set<number>();
    for (const r of pinned) if (r >= 0 && r < n && !hidden.has(r) && !pin.has(r)) pin.add(r);
    for (const r of pin) {
      entry[at] = r;
      this.displayOf[r] = at;
      if (!inPlace) place[r] = at;
      at++;
    }
    this.pinned = at;
    // the outermost folded clade of each leaf (pre-order: an ancestor is seen first)
    let foldOf: Int32Array | undefined;
    if (clades && folds.length) {
      foldOf = new Int32Array(clades.leafRow.length).fill(-1);
      const sorted = [...folds].sort((a, b) => a - b);
      for (const c of sorted) {
        if (c < 0 || c >= clades.nodes) continue;
        for (let k = clades.lo[c]; k < clades.hi[c]; k++) if (foldOf[k] < 0) foldOf[k] = c;
      }
    }
    const rows = order ?? clades?.order;
    this.ordered = rows;
    let hiddenCount = 0;
    for (let i = 0; i < n; i++) {
      const r = rows ? rows[i] : i;
      if (hidden.has(r)) {
        hiddenCount++;
        continue;
      }
      const leaf = foldOf && clades ? clades.rowLeaf[r] : -1;
      const c = leaf >= 0 ? foldOf![leaf] : -1;
      if (c < 0) {
        if (pin.has(r) && !inPlace) continue;
        entry[at] = r;
        if (!pin.has(r)) this.displayOf[r] = at;
        place[r] = at;
        at++;
        continue;
      }
      // a folded clade: its shown rows (pinned ones count too) as one display row, once
      if (at > 0 && entry[at - 1] === -1 - c) {
        if (!pin.has(r)) this.displayOf[r] = at - 1;
        place[r] = at - 1;
        continue;
      }
      const shown: number[] = [];
      for (let k = clades!.lo[c]; k < clades!.hi[c]; k++) {
        const rr = clades!.leafRow[k];
        if (rr >= 0 && !hidden.has(rr)) shown.push(rr);
      }
      if (shown.length === 1) {
        entry[at] = shown[0];
        if (!pin.has(shown[0])) this.displayOf[shown[0]] = at;
        place[shown[0]] = at;
        at++;
        continue;
      }
      if (!shown.length) continue;
      entry[at] = -1 - c;
      this.members.set(at, Int32Array.from(shown));
      for (const rr of shown) {
        if (!pin.has(rr)) this.displayOf[rr] = at;
        place[rr] = at;
      }
      at++;
    }
    this.count = at;
    this.hidden = hiddenCount;
  }
}

/** The clades `collapseAll` folds: every clade but the root (the root's children fold, each to its row). */
export function allClades(clades: Clades): number[] {
  const out: number[] = [];
  for (let n = 1; n < clades.nodes; n++) if (clades.hi[n] - clades.lo[n] > 1) out.push(n);
  return out;
}

/** The clades `collapseTo { depth }` folds: those at that depth with two leaves or more. */
export function cladesAtDepth(clades: Clades, depth: number): number[] {
  const out: number[] = [];
  for (let n = 0; n < clades.nodes; n++) if (clades.depth[n] === depth && clades.hi[n] - clades.lo[n] > 1) out.push(n);
  return out;
}

/**
 * The clades `collapseTo { height }` folds: every largest clade whose leaves lie within `height` of its node by branch
 * length (its parent's do not), of two leaves or more.
 */
export function cladesWithin(clades: Clades, height: number): number[] {
  const out: number[] = [];
  for (let n = 0; n < clades.nodes; n++) {
    if (clades.hi[n] - clades.lo[n] < 2 || clades.height[n] > height) continue;
    const p = clades.parent[n];
    if (p >= 0 && clades.height[p] <= height) continue;
    out.push(n);
  }
  return out;
}

/** The clades holding a row (its ancestors in the tree), innermost last. */
export function cladesHolding(clades: Clades, row: number): number[] {
  const leaf = clades.rowLeaf[row];
  if (leaf < 0) return [];
  const out: number[] = [];
  for (let n = 0; n < clades.nodes; n++)
    if (clades.children[n].length && clades.lo[n] <= leaf && leaf < clades.hi[n]) out.push(n);
  return out;
}

/**
 * A clade's per-column consensus over `rows` (exact: one pass over its cells): the majority letter byte a column
 * (gaps counted as a letter: a column mostly gaps reads as a gap) and its share 0–255, and the non-gap share. Into
 * `out` (3 bytes a column: letter, share, non-gap share).
 */
export function consensus(held: HeldCells, rows: ArrayLike<number>, out: Uint8Array, counts?: Uint32Array) {
  const C = held.columns;
  const n = rows.length;
  // rows outer, columns inner (the cells row-major): a count a letter a column
  const k = counts && counts.length >= C * 128 ? counts : new Uint32Array(C * 128);
  k.fill(0, 0, C * 128);
  const filled = new Uint32Array(C);
  const codes = held.codes;
  for (let i = 0; i < n; i++) {
    const row = rows[i] * C;
    for (let j = 0; j < C; j++) {
      const c = upper(codes[row + j]) & 127;
      k[j * 128 + c]++;
      if (!isGap(c)) filled[j]++;
    }
  }
  for (let j = 0; j < C; j++) {
    let best = 45;
    let most = 0;
    for (let c = 0; c < 128; c++) {
      const v = k[j * 128 + c];
      if (v > most) {
        most = v;
        best = c;
      }
    }
    out[3 * j] = best;
    out[3 * j + 1] = n ? Math.round((255 * most) / n) : 0;
    out[3 * j + 2] = n ? Math.round((255 * filled[j]) / n) : 0;
  }
  return out;
}

/** The column measures of a held alignment (_spec/README.md "Render", Column measures), computed once, kept by version. */
export interface ColumnMeasures {
  /** The most frequent letter's share a column (gaps left out of both counts with `gaps: ignore`). */
  readonly share: Float32Array;
  /** The share of rows equal to the reference row's letter a column (NaN without a reference, or a gap there). */
  readonly identity: Float32Array;
  /** One value a row: identity to the reference row, else its non-gap share (what a navigator draws, `rowsValue`). */
  readonly rowsValue: Float32Array;
  /** The reference row, -1 for none. */
  readonly reference: number;
}

/** The measures of `held` against its reference row (by label; -1 none). */
export function columnMeasures(
  held: HeldCells,
  reference: number,
  gaps: 'ignore' | 'count' = 'ignore',
): ColumnMeasures {
  const { rows: R, columns: C, codes } = held;
  const share = new Float32Array(C);
  const identity = new Float32Array(C).fill(Number.NaN);
  const counts = new Uint32Array(256);
  for (let j = 0; j < C; j++) {
    let most = 0;
    let n = 0;
    for (let i = 0; i < R; i++) {
      const c = upper(codes[i * C + j]);
      if (gaps === 'ignore' && isGap(c)) continue;
      n++;
      const k = ++counts[c];
      if (k > most) most = k;
    }
    for (let i = 0; i < R; i++) counts[upper(codes[i * C + j])] = 0;
    share[j] = n ? most / n : 0;
    if (reference >= 0) {
      const ref = upper(codes[reference * C + j]);
      if (isGap(ref)) continue;
      let same = 0;
      let m = 0;
      for (let i = 0; i < R; i++) {
        const c = upper(codes[i * C + j]);
        if (isGap(c)) continue;
        m++;
        if (c === ref) same++;
      }
      identity[j] = m ? same / m : 0;
    }
  }
  const rowsValue = new Float32Array(R);
  for (let i = 0; i < R; i++) {
    let same = 0;
    let filled = 0;
    let both = 0;
    for (let j = 0; j < C; j++) {
      const c = upper(codes[i * C + j]);
      if (isGap(c)) continue;
      filled++;
      if (reference < 0) continue;
      const ref = upper(codes[reference * C + j]);
      if (isGap(ref)) continue;
      both++;
      if (ref === c) same++;
    }
    rowsValue[i] = reference >= 0 ? (both ? same / both : 0) : C ? filled / C : 0;
  }
  return { share, identity, rowsValue, reference };
}

/**
 * The measures as theme registry entries: `share` and `identity` read a datum's column (`start`), `matches-reference`
 * its letter against the reference row's at that column (1 or 0; NaN without a reference).
 */
export function measureEntries(held: HeldCells, m: ColumnMeasures): ReadonlyMap<string, MeasureEntry> {
  const C = held.columns;
  const column = (d: ThemeDatum) => (d.start >= 0 && d.start < C ? d.start : -1);
  return new Map<string, MeasureEntry>([
    ['share', { name: 'share', read: (d) => (column(d) < 0 ? Number.NaN : m.share[column(d)]) }],
    ['identity', { name: 'identity', read: (d) => (column(d) < 0 ? Number.NaN : m.identity[column(d)]) }],
    [
      'matches-reference',
      {
        name: 'matches-reference',
        read: (d) => {
          const j = column(d);
          if (j < 0 || m.reference < 0 || !d.letter) return Number.NaN;
          return upper(held.codes[m.reference * C + j]) === upper(d.letter.charCodeAt(0)) ? 1 : 0;
        },
      },
    ],
  ]);
}

/**
 * The profile (the `matrix` shape on the columns): each column's letter frequencies over `rows` (every row when
 * undefined), gaps left out (`ignore`) or counted. `letters` the letter bytes seen (uppercase, sorted), `freq` a
 * column-major matrix (`columns × letters.length`), and each column's counted rows.
 */
export interface Profile {
  readonly letters: Uint8Array;
  readonly freq: Float32Array;
  readonly counted: Float32Array;
}

export function profileOf(held: HeldCells, gaps: 'ignore' | 'count' = 'ignore', rows?: ArrayLike<number>): Profile {
  const { columns: C, codes } = held;
  const R = rows ? rows.length : held.rows;
  const seen = new Uint8Array(256);
  for (let i = 0; i < R; i++) {
    const r = rows ? rows[i] : i;
    for (let j = 0; j < C; j++) seen[upper(codes[r * C + j])] = 1;
  }
  const letters: number[] = [];
  for (let c = 0; c < 256; c++) if (seen[c] && (gaps === 'count' || !isGap(c))) letters.push(c);
  const index = new Int16Array(256).fill(-1);
  letters.forEach((c, k) => (index[c] = k));
  const K = letters.length;
  const freq = new Float32Array(C * K);
  const counted = new Float32Array(C);
  for (let i = 0; i < R; i++) {
    const r = rows ? rows[i] : i;
    for (let j = 0; j < C; j++) {
      const k = index[upper(codes[r * C + j])];
      if (k < 0) continue;
      freq[j * K + k]++;
      counted[j]++;
    }
  }
  for (let j = 0; j < C; j++) if (counted[j]) for (let k = 0; k < K; k++) freq[j * K + k] /= counted[j];
  return { letters: Uint8Array.from(letters), freq, counted };
}

/**
 * A column's letter heights in a logo: `bits` (log2(k) − H, less the small-sample correction `(k − 1) / (2 ln 2 n)`
 * with `small-sample`), each letter its frequency times that; or `probability` (its frequency). `k` is 20 for protein,
 * 4 for nucleotides. Into `out` (one height a letter of the profile).
 */
export function logoHeights(
  p: Profile,
  column: number,
  scale: 'bits' | 'probability',
  correction: 'none' | 'small-sample',
  alphabetSize: number,
  out: Float32Array,
): Float32Array {
  const K = p.letters.length;
  if (scale === 'probability') {
    for (let k = 0; k < K; k++) out[k] = p.freq[column * K + k];
    return out;
  }
  let h = 0;
  for (let k = 0; k < K; k++) {
    const f = p.freq[column * K + k];
    if (f > 0) h -= f * Math.log2(f);
  }
  const n = p.counted[column];
  const e = correction === 'small-sample' && n > 0 ? (alphabetSize - 1) / (2 * Math.LN2 * n) : 0;
  const info = Math.max(0, Math.log2(alphabetSize) - h - e);
  for (let k = 0; k < K; k++) out[k] = p.freq[column * K + k] * info;
  return out;
}
