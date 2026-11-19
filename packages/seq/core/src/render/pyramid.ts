/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The shared 2D reduction (_spec/README.md "Render", The shared reduction): one pyramid over a typed grid of samples,
// its levels by a factor and a pairing the caller gives (square levels for a grid; a row level and a column level
// apart for rows, so a row stays a row zoomed out on columns), and a reducer (mean, max, min, or the majority with its
// share). Levels are built on demand, each from the finer level beside it, into buffers kept across bases of any size.

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
  };
  /** Votes a class while a majority sample is reduced. */
  private readonly _votes = new Float64Array(256);

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
   * levels there are. The base itself at (0, 0), its shares 255 when it was given none.
   */
  level(rowLevel: number, columnLevel: number = rowLevel): PyramidLevel {
    const square = this.options.pairing === 'square';
    const clamp = (k: number, n: number) => Math.max(0, Math.min(Math.floor(k), n - 1));
    let r: number;
    let c: number;
    if (square) r = c = clamp(rowLevel, Math.max(this.levels('rows'), this.levels('columns')));
    else {
      r = clamp(rowLevel, this.levels('rows'));
      c = clamp(columnLevel, this.levels('columns'));
    }
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
      };
      this._levels.set(key, l);
    }
    if (l.built === this._version) return l;
    // from the finer level beside it: a square level from the one below on both axes, else rows first
    let from: PyramidLevel;
    let fr = 1;
    let fc = 1;
    if (square) {
      from = this.level(r - 1, c - 1);
      fr = fc = this.options.factor;
    } else if (r > 0) {
      from = this.level(r - 1, c);
      fr = this.options.factor;
    } else {
      from = this.level(r, c - 1);
      fc = this.options.factor;
    }
    this._reduce(from, l, fr, fc);
    l.built = this._version;
    return l;
  }

  // ---- internals ----

  private _set(rows: number, columns: number): void {
    if (!(rows >= 0 && columns >= 0) || !Number.isInteger(rows) || !Number.isInteger(columns))
      throw new Error('A pyramid’s base is whole rows and columns');
    this._rows = this._base.rows = rows;
    this._columns = this._base.columns = columns;
    this._version++;
    this._base.built = this._version;
  }

  /** `to` from `from`, `fr` × `fc` samples a sample, into its kept buffers. */
  private _reduce(from: PyramidLevel, to: MutableLevel, fr: number, fc: number): void {
    const rows = Math.ceil(from.rows / fr);
    const columns = Math.ceil(from.columns / fc);
    const n = rows * columns;
    to.rows = rows;
    to.columns = columns;
    if (this.options.reducer === 'majority') {
      if (to.classes.length < n) {
        to.classes = new Uint8Array(n);
        to.shares = new Uint8Array(n);
      }
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
      );
    } else {
      if (to.values.length < n) to.values = new Float32Array(n);
      reduceNumbers(from, to, fr, fc, this.options.reducer);
    }
  }
}

function reduceNumbers(from: PyramidLevel, to: MutableLevel, fr: number, fc: number, reducer: Reducer): void {
  const src = from.values;
  const out = to.values;
  const W = from.columns;
  const H = from.rows;
  for (let i = 0; i < to.rows; i++) {
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
): void {
  const classes = from.classes;
  const shares = from.shares;
  const whole = shares.length === 0;
  const W = from.columns;
  const H = from.rows;
  for (let i = 0; i < to.rows; i++) {
    const r0 = i * fr;
    const r1 = Math.min(H, r0 + fr);
    for (let j = 0; j < to.columns; j++) {
      const c0 = j * fc;
      const c1 = Math.min(W, c0 + fc);
      let n = 0;
      for (let r = r0; r < r1; r++) {
        const h = Math.min(baseRows, (r + 1) * unitRows) - r * unitRows;
        for (let c = c0; c < c1; c++) {
          const k = r * W + c;
          const cl = classes[k];
          const cells = h * (Math.min(baseColumns, (c + 1) * unitColumns) - c * unitColumns);
          n += cells;
          if (cl !== NoClass) votes[cl] += (whole ? 255 : shares[k]) * cells;
        }
      }
      let best = NoClass;
      let most = 0;
      for (let r = r0; r < r1; r++)
        for (let c = c0; c < c1; c++) {
          const cl = classes[r * W + c];
          if (cl === NoClass) continue;
          const v = votes[cl];
          if (v > most || (v === most && cl < best)) {
            most = v;
            best = cl;
          }
        }
      // clear what was voted (only the classes seen: no pass over 256 a sample)
      for (let r = r0; r < r1; r++) for (let c = c0; c < c1; c++) votes[classes[r * W + c]] = 0;
      const o = i * to.columns + j;
      to.classes[o] = best;
      to.shares[o] = n === 0 || best === NoClass ? 0 : Math.round(most / n);
    }
  }
}
