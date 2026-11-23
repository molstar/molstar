/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The shared 2D reduction (_spec/README.md "Render", The shared reduction): one pyramid over a typed grid of samples,
// its levels by a factor and a pairing the caller gives (square levels for a grid; a row level and a column level
// apart for rows, so a row stays a row zoomed out on columns), and a reducer (mean, max, min, or the majority with its
// share). Levels are built on demand, each from the finer level beside it, into buffers kept across bases of any size;
// at once (`level`) or a budget of cells at a time across calls (`levelWithin`: a frame builds a level across frames).

/** How a sample of a coarser level reads its finer samples. */
export type Reducer = 'mean' | 'max' | 'min' | 'majority';

export interface PyramidOptions {
  /** Samples a coarser sample takes along an axis (2 for a grid, 4 for rows). */
  readonly factor: number;
  /** `square`: one level for both axes; `independent`: a row level and a column level apart. */
  readonly pairing: 'square' | 'independent';
  readonly reducer: Reducer;
}

/**
 * A level of the pyramid: `rows` × `columns` samples, row-major. Numbers (`values`, NaN missing) for `mean`, `max`
 * and `min`; for `majority` a class a sample (`classes`, 255 none) and its share (`shares`, 0 to 255 for 0 to 1).
 */
export interface PyramidLevel {
  readonly rowLevel: number;
  readonly columnLevel: number;
  readonly rows: number;
  readonly columns: number;
  readonly values: Float32Array;
  readonly classes: Uint8Array;
  readonly shares: Uint8Array;
}

/** A sample with no class (a majority level's empty cell). */
export const NoClass = 255;

/**
 * A majority level whose sample covers fewer base rows than this counts its cells exactly, from the base; a larger one
 * is the majority of the majorities below it (_spec/README.md "Render", the rows layer). Why: a majority of majorities
 * misses the letter in small blocks (a lead up to 0.41 over seeds 1 to 12 of an unbiased generator), and a few rows a
 * sample are cheap to count.
 */
export const ExactMajorityRows = 64;

interface MutableLevel {
  rowLevel: number;
  columnLevel: number;
  rows: number;
  columns: number;
  values: Float32Array;
  classes: Uint8Array;
  shares: Uint8Array;
  /** The base's version it was built for; -1 never. */
  built: number;
  /** The base's version it is being built for, and its rows built so far (a build across calls). */
  building: number;
  done: number;
}

const Empty32 = new Float32Array(0);
const Empty8 = new Uint8Array(0);

/** The level whose sample is at least a device pixel: `ceil(log_f(units a pixel))`, 0 at a unit a pixel or more. */
export function pyramidLevelFor(unitsPerPixel: number, factor: number): number {
  if (!(unitsPerPixel > 1)) return 0;
  return Math.max(0, Math.ceil(Math.log(unitsPerPixel) / Math.log(factor) - 1e-9));
}

/**
 * The pyramid: `setBase` gives the finest samples (by reference: the caller keeps them unchanged while it reads the
 * pyramid), `level` a coarser level, built when first asked since the base changed. Buffers are kept by level and
 * grown, never shrunk: a base after a larger one builds into what was there.
 */
export class Pyramid {
  readonly options: PyramidOptions;
  private _rows = 0;
  private _columns = 0;
  private _version = 0;
  private readonly _levels = new Map<number, MutableLevel>();
  private readonly _base: MutableLevel = {
    rowLevel: 0,
    columnLevel: 0,
    rows: 0,
    columns: 0,
    values: Empty32,
    classes: Empty8,
    shares: Empty8,
    built: 0,
    building: 0,
    done: 0,
  };
  /** Votes a class while a majority sample is reduced, the sample each was last counted in (the last: the stamp). */
  private readonly _votes = new Float64Array(256);
  private readonly _stamps = new Uint32Array(257);
  /** The cells each finer sample covers, by row and by column (kept). */
  private readonly _weights = { rows: new Float64Array(0), columns: new Float64Array(0) };

  constructor(options: PyramidOptions) {
    if (!(options.factor >= 2) || !Number.isInteger(options.factor)) throw new Error('A pyramid’s factor is 2 or more');
    this.options = options;
  }

  get rows(): number {
    return this._rows;
  }

  get columns(): number {
    return this._columns;
  }

  /** Moves with every base: what a cache of drawn levels is keyed by. */
  get version(): number {
    return this._version;
  }

  /** Numbers a sample (`mean`, `max`, `min`), NaN missing. */
  setValues(rows: number, columns: number, values: Float32Array): void {
    if (this.options.reducer === 'majority') throw new Error('A majority pyramid reads classes');
    this._set(rows, columns);
    this._base.values = values;
  }

  /** Classes a sample (`majority`), each with its share (255 when absent: a cell is wholly its class). */
  setClasses(rows: number, columns: number, classes: Uint8Array, shares?: Uint8Array): void {
    if (this.options.reducer !== 'majority') throw new Error('Only a majority pyramid reads classes');
    this._set(rows, columns);
    this._base.classes = classes;
    this._base.shares = shares ?? Empty8;
  }

  /** How many levels an axis has: the last is one sample along it. */
  levels(axis: 'rows' | 'columns'): number {
    let n = axis === 'rows' ? this._rows : this._columns;
    let k = 0;
    while (n > 1) {
      n = Math.ceil(n / this.options.factor);
      k++;
    }
    return k + 1;
  }

  /**
   * A level: `level(k)` for a square pyramid, `level(rowLevel, columnLevel)` for an independent one; clamped to the
   * levels there are. The base itself at (0, 0), its shares 255 when it was given none. Built whole now when it is not.
   */
  level(rowLevel: number, columnLevel: number = rowLevel): PyramidLevel {
    this._clamp(rowLevel, columnLevel);
    const r = this._r;
    const c = this._c;
    this._advance(r, c, Number.POSITIVE_INFINITY);
    return this._slot(r, c);
  }

  /**
   * A level as `level` gives it when it is built; else its build goes on by about `cells` cells read (its finer levels
   * first; at least a row of samples a call) and undefined until it is whole: what a frame asks while the view moves.
   */
  levelWithin(rowLevel: number, columnLevel: number, cells: number): PyramidLevel | undefined {
    this._clamp(rowLevel, columnLevel);
    const r = this._r;
    const c = this._c;
    return this._advance(r, c, cells) < 0 ? undefined : this._slot(r, c);
  }

  /** Whether a level (clamped) is built for the current base. */
  isBuilt(rowLevel: number, columnLevel: number = rowLevel): boolean {
    this._clamp(rowLevel, columnLevel);
    return (this._r === 0 && this._c === 0) || this._slot(this._r, this._c).built === this._version;
  }

  // ---- internals ----

  /** The clamped level pair `_clamp` writes (fields: a pair returned would be an allocation a call). */
  private _r = 0;
  private _c = 0;

  private _clamp(rowLevel: number, columnLevel: number): void {
    const clamp = (k: number, n: number) => Math.max(0, Math.min(Math.floor(k), n - 1));
    if (this.options.pairing === 'square')
      this._r = this._c = clamp(rowLevel, Math.max(this.levels('rows'), this.levels('columns')));
    else {
      this._r = clamp(rowLevel, this.levels('rows'));
      this._c = clamp(columnLevel, this.levels('columns'));
    }
  }

  /** A level's kept buffers (the base at (0, 0)). */
  private _slot(r: number, c: number): MutableLevel {
    if (r === 0 && c === 0) return this._base;
    const key = r * 65536 + c;
    let l = this._levels.get(key);
    if (!l) {
      l = {
        rowLevel: r,
        columnLevel: c,
        rows: 0,
        columns: 0,
        values: Empty32,
        classes: Empty8,
        shares: Empty8,
        built: -1,
        building: -1,
        done: 0,
      };
      this._levels.set(key, l);
    }
    return l;
  }

  /**
   * Builds level (r, c) on by at most about `cells` cells read: from the base when its block is small (a majority's
   * exact count), else from the finer level beside it (a square level from the one below on both axes, else rows
   * first), which is built first. The budget left, or −1 while the level is not whole.
   */
  private _advance(r: number, c: number, cells: number): number {
    if (r === 0 && c === 0) return cells;
    const l = this._slot(r, c);
    if (l.built === this._version) return cells;
    const f = this.options.factor;
    let from: PyramidLevel;
    let fr = 1;
    let fc = 1;
    if (this.options.reducer === 'majority' && f ** r < ExactMajorityRows) {
      // a small block counts its cells exactly, from the base (_spec/README.md "Render", the rows layer)
      from = this._base;
      fr = f ** r;
      fc = f ** c;
    } else {
      const square = this.options.pairing === 'square';
      const fr0 = square || r > 0 ? r - 1 : r;
      const fc0 = square || r === 0 ? c - 1 : c;
      cells = this._advance(fr0, fc0, cells);
      if (cells < 0) return -1;
      from = this._slot(fr0, fc0);
      fr = fr0 === r ? 1 : f;
      fc = fc0 === c ? 1 : f;
    }
    if (l.building !== this._version) {
      l.building = this._version;
      l.done = 0;
      this._size(from, l, fr, fc);
    }
    // a row of samples reads `fr` rows of the finer level, each its columns
    const perRow = Math.max(1, fr * from.columns);
    const rows = cells === Number.POSITIVE_INFINITY ? l.rows : Math.max(1, Math.floor(cells / perRow));
    const i1 = Math.min(l.rows, l.done + rows);
    this._reduce(from, l, fr, fc, l.done, i1);
    cells -= (i1 - l.done) * perRow;
    l.done = i1;
    if (l.done < l.rows) return -1;
    l.built = this._version;
    return Math.max(0, cells);
  }

  private _set(rows: number, columns: number): void {
    if (!(rows >= 0 && columns >= 0) || !Number.isInteger(rows) || !Number.isInteger(columns))
      throw new Error('A pyramid’s base is whole rows and columns');
    this._rows = this._base.rows = rows;
    this._columns = this._base.columns = columns;
    this._version++;
    this._base.built = this._version;
  }

  /** `to`'s size from `from`, `fr` × `fc` samples a sample, its kept buffers grown to hold it. */
  private _size(from: PyramidLevel, to: MutableLevel, fr: number, fc: number): void {
    to.rows = Math.ceil(from.rows / fr);
    to.columns = Math.ceil(from.columns / fc);
    const n = to.rows * to.columns;
    if (this.options.reducer === 'majority') {
      if (to.classes.length < n) {
        to.classes = new Uint8Array(n);
        to.shares = new Uint8Array(n);
      }
    } else if (to.values.length < n) to.values = new Float32Array(n);
  }

  /** `to`'s rows `i0` to `i1` from `from`, `fr` × `fc` samples a sample. */
  private _reduce(from: PyramidLevel, to: MutableLevel, fr: number, fc: number, i0: number, i1: number): void {
    if (this.options.reducer === 'majority') {
      const f = this.options.factor;
      reduceMajority(
        from,
        to,
        fr,
        fc,
        this._votes,
        f ** from.rowLevel,
        f ** from.columnLevel,
        this._rows,
        this._columns,
        this._stamps,
        this._weights,
        i0,
        i1,
      );
    } else reduceNumbers(from, to, fr, fc, this.options.reducer, i0, i1);
  }
}

function reduceNumbers(
  from: PyramidLevel,
  to: MutableLevel,
  fr: number,
  fc: number,
  reducer: Reducer,
  i0: number,
  i1: number,
): void {
  const src = from.values;
  const out = to.values;
  const W = from.columns;
  const H = from.rows;
  for (let i = i0; i < i1; i++) {
    const r0 = i * fr;
    const r1 = Math.min(H, r0 + fr);
    for (let j = 0; j < to.columns; j++) {
      const c0 = j * fc;
      const c1 = Math.min(W, c0 + fc);
      let acc = reducer === 'max' ? -Infinity : reducer === 'min' ? Infinity : 0;
      let n = 0;
      for (let r = r0; r < r1; r++)
        for (let c = c0; c < c1; c++) {
          const v = src[r * W + c];
          if (Number.isNaN(v)) continue; // missing, not counted
          n++;
          if (reducer === 'mean') acc += v;
          else if (reducer === 'max') {
            if (v > acc) acc = v;
          } else if (v < acc) acc = v;
        }
      out[i * to.columns + j] = n === 0 ? Number.NaN : reducer === 'mean' ? acc / n : acc;
    }
  }
}

/**
 * The majority of majorities: each finer sample votes its class by its share times the cells it covers (a sample at an
 * edge covers fewer); the winner is the class with the most votes (the lower class on a tie), its share the votes over
 * the cells (an approximation of its share of the cells: a majority is not mergeable).
 */
function reduceMajority(
  from: PyramidLevel,
  to: MutableLevel,
  fr: number,
  fc: number,
  votes: Float64Array,
  unitRows: number,
  unitColumns: number,
  baseRows: number,
  baseColumns: number,
  stamps: Uint32Array,
  weights: { rows: Float64Array; columns: Float64Array },
  i0: number,
  i1: number,
): void {
  const classes = from.classes;
  const shares = from.shares;
  const whole = shares.length === 0;
  const W = from.columns;
  const H = from.rows;
  // the cells each finer sample covers: whole units but at the far edges
  if (weights.rows.length < H) weights.rows = new Float64Array(H);
  if (weights.columns.length < W) weights.columns = new Float64Array(W);
  const rw = weights.rows;
  const cw = weights.columns;
  for (let r = 0; r < H; r++) rw[r] = Math.min(baseRows, (r + 1) * unitRows) - r * unitRows;
  for (let c = 0; c < W; c++) cw[c] = Math.min(baseColumns, (c + 1) * unitColumns) - c * unitColumns;
  let stamp = stamps[256];
  for (let i = i0; i < i1; i++) {
    const r0 = i * fr;
    const r1 = Math.min(H, r0 + fr);
    for (let j = 0; j < to.columns; j++) {
      const c0 = j * fc;
      const c1 = Math.min(W, c0 + fc);
      // a new stamp a sample: a class's votes count from 0 when first seen in it (no pass to clear them)
      stamp = (stamp + 1) >>> 0;
      let n = 0;
      let best = NoClass;
      let most = 0;
      for (let r = r0; r < r1; r++) {
        const h = rw[r];
        const row = r * W;
        for (let c = c0; c < c1; c++) {
          const k = row + c;
          const cells = h * cw[c];
          n += cells;
          const cl = classes[k];
          if (cl === NoClass) continue;
          if (stamps[cl] !== stamp) {
            stamps[cl] = stamp;
            votes[cl] = 0;
          }
          const v = (votes[cl] += (whole ? 255 : shares[k]) * cells);
          // votes only grow: the argmax followed as they do, the lower class on a tie
          if (v > most || (v === most && cl < best)) {
            most = v;
            best = cl;
          }
        }
      }
      const o = i * to.columns + j;
      to.classes[o] = best;
      to.shares[o] = n === 0 || best === NoClass ? 0 : Math.round(most / n);
    }
  }
  stamps[256] = stamp;
}
