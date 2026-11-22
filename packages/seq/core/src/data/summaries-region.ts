/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, categoryFromCodes, issueOf } from '@molstar/core/repr/tree';
import { StoreIssues } from './issues.js';
import type {
  Coverage,
  GridDense,
  OpenedSource,
  PairsBins,
  PairsList,
  RegionAnswer,
  RegionRequest,
  ShapeInfo,
} from './protocol.js';
import { binPairsRange, coarsenPairBins, pairsBinner } from './summaries-pairs.js';
import { cutGrid, type GridPlane, gridPlane, halveGrid, type WholeBlock, wholePlane } from './summaries-grid.js';
import { checkAnswer } from './tables.js';
import { regionBytes, sliceRegion } from './tables-region.js';

// The store's pyramids for one-level sources' pairs and grids (_spec/README.md "Data", Summaries): read whole once
// within `maxBytes`, built in time slices of `sliceMs`, owned by the attached source and counted against the store's
// bound, gone with it.

/** What the region pyramids are bounded by and how they are built. */
export interface RegionSummarySettings {
  /** Each pairs level's bin width over the one below. */
  readonly factor: number;
  /** A whole read and its pyramid, bytes: 16 MB. */
  readonly maxBytes: number;
  /** The time a slice of a pyramid's build runs before it yields, ms: 4. */
  readonly sliceMs: number;
}

/** How a build yields between slices (a frame by default): tests run it at once. */
export interface SliceClock {
  now(): number;
  yield(): Promise<void>;
}

export const FrameSlices: SliceClock = {
  now: () => performance.now(),
  yield: () => new Promise((r) => setTimeout(r, 0)),
};

/**
 * A build started once and shared: `start` makes it the first time, `building` is it until it settles (the store is not
 * settled while a source's summaries build, EC-125: a row drawn from them never comes after "settled").
 */
export function sharedBuild<T>(make: () => Promise<T>) {
  let built: Promise<T> | undefined;
  let running: Promise<T> | undefined;
  return {
    start(): Promise<T> {
      if (built) return built;
      const p = (built = running = make());
      const end = () => {
        if (running === p) running = undefined;
      };
      p.then(end, end);
      return p;
    },
    /** Forgets a failed build: the store's retry asks again. */
    forget(p: Promise<T>) {
      if (built === p) built = undefined;
    },
    get building(): Promise<T> | undefined {
      return running;
    },
  };
}

/** A failure answer of a summaries' code. */
function refused(code: 'seq.store.pairs.levels' | 'seq.store.grid.levels', message: string): RegionAnswer {
  return { ok: false, kind: 'too-large', issue: issueOf(StoreIssues, code, message) };
}

/** A pyramid's answer to one request: a table, or why not. */
type Built<T> = T | RegionAnswer;

/** The levels a pyramid of `factor` adds over `[1]`, up to `longest` (the first whose bin covers it, excluded). */
export function pyramidLevels(longest: number, factor: number): number[] {
  const out: number[] = [];
  for (let l = factor; l < longest; l *= factor) out.push(l);
  return out;
}

/** The largest of `own` dividing `level`, no wider: the level a request is cut from (then coarsened). */
function sourceLevel(own: readonly number[], level: number): number | undefined {
  let at: number | undefined;
  for (const l of own) if (l <= level && level % l === 0) at = l;
  return at;
}

const NoLevels: readonly number[] = [];

/** The levels a set leaves out. */
function without(levels: readonly number[], kept: ReadonlySet<number>): number[] {
  return levels.filter((l) => !kept.has(l));
}

/** A read of a whole system's pairs (its side's, for a sided source). */
export interface PairsRead {
  readonly system: string;
  readonly side?: 'a' | 'b';
  readonly list: PairsList;
}

/**
 * Lists read on each system merged, each pair once: a pair is in the read of each system an end is on, and kept from
 * the read of its `a` end's system (of side `a`), else, that system not read, from its `b` end's. No identity is
 * compared: a pair's place decides.
 */
export function mergePairs(reads: readonly PairsRead[], field?: string): PairsList {
  const steps = mergeSteps(reads, field);
  for (;;) {
    const r = steps.next();
    if (r.done) return r.value;
  }
}

/** Rows a merge step takes before it yields (a slice's grain). */
const MergeStep = 4096;

/**
 * `mergePairs` a step at a time (each at most `MergeStep` rows), so a pyramid's merge is built in time slices; `field`
 * a numeric field carried (summed by the bins: Decided 2).
 */
export function* mergeSteps(reads: readonly PairsRead[], field?: string): Generator<void, PairsList> {
  const readA = new Set(reads.filter((r) => r.side !== 'b').map((r) => r.system));
  const dict: string[] = [];
  const index = new Map<string, number>();
  const code = (name: string) => {
    let c = index.get(name);
    if (c === undefined) {
      c = dict.length;
      dict.push(name);
      index.set(name, c);
    }
    return c;
  };
  const typeDict: string[] = [];
  const typeIndex = new Map<string, number>();
  let n = 0;
  const keep: { read: PairsRead; rows: Int32Array; aCodes: Int32Array; bCodes: Int32Array; tCodes?: Int32Array }[] = [];
  for (const read of reads) {
    const t = read.list;
    // each list's dictionaries mapped once to the merged ones
    const own = code(read.system);
    const aMap = t.aSystem ? t.aSystem.dict.map(code) : undefined;
    const bMap = t.bSystem ? t.bSystem.dict.map(code) : undefined;
    const tMap = t.type?.dict.map((v) => {
      let c = typeIndex.get(v);
      if (c === undefined) typeIndex.set(v, (c = typeDict.push(v) - 1));
      return c;
    });
    const len = t.aStart.length;
    const rows = new Int32Array(len);
    const aCodes = new Int32Array(len);
    const bCodes = new Int32Array(len);
    const ac = t.aSystem?.codes;
    const bc = t.bSystem?.codes;
    // whether each merged system is read on side `a` (by code, filled as codes come)
    const aRead = dict.map((d) => readA.has(d));
    const takeA = read.side !== 'b';
    const takeB = read.side !== 'a';
    let m = 0;
    for (let i = 0; i < len; i++) {
      if ((i & (MergeStep - 1)) === MergeStep - 1) yield;
      const as = aMap ? aMap[ac![i]] : own;
      const bs = bMap ? bMap[bc![i]] : own;
      aCodes[i] = as;
      bCodes[i] = bs;
      const fromA = takeA && as === own;
      if (fromA || (takeB && bs === own && !(aRead[as] ?? (aRead[as] = readA.has(dict[as]))))) rows[m++] = i;
    }
    const tCodes = tMap && Int32Array.from(t.type!.codes, (c) => (c < 0 ? -1 : tMap[c]));
    keep.push({ read, rows: rows.subarray(0, m), aCodes, bCodes, tCodes });
    n += m;
  }
  const strands = reads.some((r) => r.list.aStrand || r.list.bStrand);
  const types = reads.some((r) => r.list.type);
  const out = {
    aSys: new Int32Array(n),
    bSys: new Int32Array(n),
    aStart: new Int32Array(n),
    aEnd: new Int32Array(n),
    bStart: new Int32Array(n),
    bEnd: new Int32Array(n),
    aStrand: new Int8Array(n),
    bStrand: new Int8Array(n),
    type: new Int32Array(n).fill(-1),
    field: new Float32Array(field === undefined ? 0 : n),
  };
  const carried = field !== undefined && reads.some((r) => fieldOf(r.list, field));
  let k = 0;
  for (const { read, rows, aCodes, bCodes, tCodes } of keep) {
    const t = read.list;
    const { aStart, aEnd, bStart, bEnd, aStrand, bStrand } = t;
    const f = field === undefined ? undefined : fieldOf(t, field);
    for (let j = 0; j < rows.length; j++, k++) {
      if ((j & (MergeStep - 1)) === MergeStep - 1) yield;
      const i = rows[j];
      if (carried) out.field[k] = f ? f[i] : Number.NaN;
      out.aSys[k] = aCodes[i];
      out.bSys[k] = bCodes[i];
      out.aStart[k] = aStart[i];
      out.aEnd[k] = aEnd[i];
      out.bStart[k] = bStart[i];
      out.bEnd[k] = bEnd[i];
      if (aStrand) out.aStrand[k] = aStrand[i];
      if (bStrand) out.bStrand[k] = bStrand[i];
      if (tCodes) out.type[k] = tCodes[i];
    }
  }
  return {
    form: 'list',
    aSystem: categoryFromCodes(dict, out.aSys),
    bSystem: categoryFromCodes(dict, out.bSys),
    aStart: out.aStart,
    aEnd: out.aEnd,
    bStart: out.bStart,
    bEnd: out.bEnd,
    ...(strands && { aStrand: out.aStrand, bStrand: out.bStrand }),
    ...(types && { type: categoryFromCodes(typeDict, out.type) }),
    ...(carried && { fields: [{ name: field!, column: out.field }] }),
  };
}

/** A list's numeric field by name. */
function fieldOf(t: PairsList, name: string): ArrayLike<number> | undefined {
  const c = t.fields?.find((x) => x.name === name)?.column;
  return c && !(c instanceof CategoryColumn) && typeof c !== 'string' ? (c as ArrayLike<number>) : undefined;
}

/** A pyramid's bins at a level cut to a request's region (an end's bin meeting it on the request's system and side). */
function cutPairs(bins: PairsBins, own: number, request: RegionRequest): PairsBins {
  const at = request.level === own ? bins : coarsenPairBins(bins, request.level / own);
  const { system, start, end, side } = request;
  return sliceRegion(at, request.level, { start, end, system, ...(side && { side }) });
}

/**
 * The pairs pyramid of a one-level source: every pairs system read whole at level 1 (shared by every caller), merged
 * once a pair (merged and binned in time slices, the shape's numeric field summed), binned at the first level and
 * coarsened for each further one, only the levels that halve the bins kept: a dropped level leaves the coverage once
 * built (`dropped`). Past `maxBytes`, a coarse tile is its items read alone when its `readCost` says they fit, else
 * refused before reading.
 */
export function pairsPyramid(
  inner: OpenedSource,
  coverage: Coverage,
  levels: readonly number[],
  settings: RegionSummarySettings,
  maxRows: number,
  clock: SliceClock,
) {
  const shapes = coverage.shapes.filter((s) => s.shape === 'pairs');
  const category = shapes.some((s) => s.binBy === 'type');
  // the shape's numeric field the bins sum (Decided 2): the first the coverage lists
  const field = shapes.flatMap((s) => s.fields ?? []).find((f) => f.type === 'float32' || f.type === 'int32')?.name;
  /** The levels built and kept (those that halve the rows below), once the pyramid is built. */
  let kept: ReadonlySet<number> | undefined;
  let dropped: readonly number[] | undefined;
  let bytes = 0;
  const reading = new AbortController();

  const tooMany = () =>
    refused('seq.store.pairs.levels', `The pairs hold more than ${settings.maxBytes} bytes: no levels`);

  const build = async (): Promise<Map<number, PairsBins> | RegionAnswer> => {
    const reads: PairsRead[] = [];
    let held = 0;
    for (const sys of coverage.systems) {
      if (!shapes.some((s) => s.system === sys.name) || sys.length === undefined) continue;
      const request: RegionRequest = {
        system: sys.name,
        start: 0,
        end: sys.length,
        level: 1,
        shapes: ['pairs'],
        ...(sys.side && { side: sys.side }),
      };
      const a = await inner.read(request, reading.signal);
      const bad = checkAnswer(a, request, maxRows);
      if (bad) return refused('seq.store.pairs.levels', `The pairs read whole are malformed: ${bad}`);
      if (!a.ok) return a;
      const list = a.tables.pairs;
      if (list?.form !== 'list') continue;
      // the reads hold a pair once an end: twice the bound, before they are merged
      held += regionBytes(list);
      if (held > 2 * settings.maxBytes) return tooMany();
      reads.push({ system: sys.name, list, ...(sys.side && { side: sys.side }) });
      // a system's read checked: one task a system
      await clock.yield();
    }
    // merged in slices too: 600,000 pairs are a few dozen of them, never one long task
    const steps = mergeSteps(reads, field);
    let merged = steps.next();
    for (;;) {
      const until = clock.now() + settings.sliceMs;
      while (!merged.done && clock.now() < until) merged = steps.next();
      if (merged.done) break;
      await clock.yield();
    }
    const items = merged.value;
    if (regionBytes(items) > settings.maxBytes) return tooMany();
    const made = new Map<number, PairsBins>();
    let below: PairsBins | undefined;
    let belowLevel = 1;
    let belowRows = items.aStart.length;
    for (const level of levels) {
      let bins: PairsBins;
      if (!below) {
        const run = binPairsRange(
          items,
          level,
          { category, systems: true, ...(field !== undefined && { field }) },
          0,
          items.aStart.length,
          pairsBinner(),
        );
        // binned in slices: a 600,000-pair read is a few dozen slices, never one long task
        let at = 0;
        for (;;) {
          const until = clock.now() + settings.sliceMs;
          while (at < items.aStart.length && clock.now() < until) run.step((at += 4096));
          if (at >= items.aStart.length) break;
          await clock.yield();
        }
        bins = run.finish();
      } else bins = coarsenPairBins(below, level / belowLevel);
      // each level its own task
      await clock.yield();
      // only the levels that halve the rows below are kept: sparse pairs barely compress at fine levels
      if (bins.count.length * 2 <= belowRows) {
        made.set(level, bins);
        bytes += regionBytes(bins);
      }
      below = bins;
      belowLevel = level;
      belowRows = bins.count.length;
    }
    kept = new Set(made.keys());
    return made;
  };

  /**
   * A coarse tile from its own items: past `maxBytes` (`bounded`), read only when the source tells the read's cost
   * (`readCost`) within `maxBytes` before reading, else refused at once (the layer reads coarser, then says zoom in),
   * never read to be refused after; a level asked while the pyramid was building that it then dropped, read (the
   * pairs read whole fit `maxBytes`).
   */
  const fromItems = async (
    request: RegionRequest,
    signal: AbortSignal,
    bounded: boolean,
  ): Promise<Built<PairsBins | PairsList>> => {
    const one: RegionRequest = { ...request, level: 1, shapes: ['pairs'] };
    const cost = bounded ? inner.readCost?.(one) : 0;
    if (cost === undefined || cost > settings.maxBytes)
      return refused('seq.store.pairs.levels', 'A coarse pairs tile past the summaries’ bytes: zoom in');
    const a = await inner.read(one, signal);
    if (!a.ok) return a;
    return (
      a.tables.pairs ?? {
        form: 'list',
        aStart: new Int32Array(0),
        aEnd: new Int32Array(0),
        bStart: new Int32Array(0),
        bEnd: new Int32Array(0),
      }
    );
  };

  const shared = sharedBuild(build);
  return {
    get bytes() {
      return bytes;
    },
    /** The build while it runs (`sharedBuild`). */
    get building() {
      return shared.building;
    },
    /** The pairs levels the pyramid dropped (they did not halve the rows below), once built: left out of the coverage. */
    /** Whether the pyramid is built: its `dropped` is then final. */
    get built(): boolean {
      return kept !== undefined;
    },
    get dropped(): readonly number[] {
      // made once the pyramid is built: the coverage reads it every frame (a store's `coverage`), allocating nothing
      if (dropped) return dropped;
      if (kept) dropped = without(levels, kept);
      return dropped ?? NoLevels;
    },
    async read(request: RegionRequest, signal: AbortSignal): Promise<Built<PairsBins | PairsList>> {
      const p = shared.start();
      // a failure is not kept: the store's retry asks again
      p.then((m) => {
        if (!(m instanceof Map)) shared.forget(p);
      });
      const made = await untilAborted(p, signal);
      if (!(made instanceof Map)) {
        if (made.ok === false && made.issue.code === 'seq.store.pairs.levels') return fromItems(request, signal, true);
        return made;
      }
      // a dropped level asked while the pyramid was building: from a kept finer level when there is one
      const own = sourceLevel([...made.keys()], request.level);
      // none: its items (the pairs read whole fit `maxBytes`, so a tile's do)
      if (own === undefined) return fromItems(request, signal, false);
      return cutPairs(made.get(own)!, own, request);
    },
    dispose() {
      reading.abort();
      bytes = 0;
    },
  };
}

/** Rejects with an `AbortError` as soon as `signal` aborts, whatever `p` does. */
export function untilAborted<T>(p: Promise<T>, signal: AbortSignal): Promise<T> {
  if (signal.aborted) return Promise.reject(abortErrorOf());
  return new Promise<T>((resolve, reject) => {
    const onAbort = () => reject(abortErrorOf());
    signal.addEventListener('abort', onAbort, { once: true });
    p.then(
      (v) => {
        signal.removeEventListener('abort', onAbort);
        resolve(v);
      },
      (e) => {
        signal.removeEventListener('abort', onAbort);
        reject(e);
      },
    );
  });
}

function abortErrorOf(): Error {
  const e = new Error('The read was aborted');
  e.name = 'AbortError';
  return e;
}

/** A grid's planes by level for one pair of systems. */
type Planes = Map<number, GridPlane>;

/**
 * The grid pyramid of a one-level source: each listed pair of systems (a symmetric one's served triangle) read whole
 * at level 1 within `maxBytes` together, summarized at factor 2 a side by its `aggregate`; a tile cut from its level.
 * A grid with a `whole` system gets its whole-assembly plane at a level from its blocks' planes. Past `maxBytes`,
 * coarse tiles are `too-large` (`seq.store.grid.levels`: a format with its own levels).
 */
export function gridPyramid(
  inner: OpenedSource,
  coverage: Coverage,
  levels: readonly number[],
  settings: RegionSummarySettings,
  maxRows: number,
) {
  const grids = coverage.shapes.filter((s) => s.shape === 'grid');
  const lengthOf = (name: string) => coverage.systems.find((s) => s.name === name)?.length;
  let bytes = 0;
  const wholes = new Map<string, GridPlane>();
  const reading = new AbortController();
  const order = (name: string) => coverage.systems.findIndex((s) => s.name === name);

  /** The pairs of systems a grid shape lists (a symmetric one's at or above the diagonal). */
  // a whole-assembly system is never read: its plane is made from its blocks'
  const wholeSystems = new Set(grids.flatMap((g) => (g.whole ? [g.whole.system] : [])));
  const pairsOf = (g: ShapeInfo): [string, string][] => {
    if (wholeSystems.has(g.system)) return [];
    const ys = g.y === '*' ? coverage.systems.map((s) => s.name) : (g.y ?? [g.system]);
    return ys
      .filter((y) => !g.symmetric || order(g.system) <= order(y))
      .filter((y) => !wholeSystems.has(y))
      .map((y) => [g.system, y]);
  };

  const build = async (): Promise<Map<string, Planes> | RegionAnswer> => {
    const all = grids.flatMap(pairsOf);
    let cells = 0;
    for (const [x, y] of all) cells += (lengthOf(x) ?? Infinity) * (lengthOf(y) ?? Infinity);
    if (4 * cells > settings.maxBytes)
      return refused('seq.store.grid.levels', `The grid's cells are past ${settings.maxBytes} bytes: no levels`);
    const out = new Map<string, Planes>();
    for (const [x, y] of all) {
      const nx = lengthOf(x)!;
      const ny = lengthOf(y)!;
      const request: RegionRequest = {
        system: x,
        start: 0,
        end: nx,
        level: 1,
        shapes: ['grid'],
        y: { system: y, start: 0, end: ny },
      };
      const a = await inner.read(request, reading.signal);
      const bad = checkAnswer(a, request, maxRows);
      if (bad) return refused('seq.store.grid.levels', `The grid read whole is malformed: ${bad}`);
      if (!a.ok) return a;
      if (!a.tables.grid) continue;
      const aggregate = grids.find((g) => g.system === x)?.aggregate ?? 'mean';
      const planes: Planes = new Map();
      let below = gridPlane(a.tables.grid, nx, ny);
      let belowLevel = 1;
      for (const level of levels) {
        for (let l = belowLevel * 2; l <= level; l *= 2) below = halveGrid(below, aggregate);
        belowLevel = level;
        planes.set(level, below);
        bytes += below.value.byteLength;
      }
      out.set(`${x}\u0000${y}`, planes);
    }
    return out;
  };

  /** The whole-assembly plane at `level`, from the blocks' planes at the finest pyramid level dividing it. */
  const wholeAt = (g: ShapeInfo, made: Map<string, Planes>, level: number): GridPlane | undefined => {
    const w = g.whole!;
    const key = `${w.system}\u0000${level}`;
    const at = wholes.get(key);
    if (at) return at;
    const offsets = new Map<string, number>();
    let length = 0;
    for (const s of w.sequences) {
      offsets.set(s, length);
      length += lengthOf(s) ?? 0;
    }
    const blocks: WholeBlock[] = [];
    for (const [k, planes] of made) {
      const [x, y] = k.split('\u0000');
      const own = sourceLevel([...planes.keys()], level);
      if (own === undefined || !offsets.has(x) || !offsets.has(y)) continue;
      const plane = planes.get(own)!;
      blocks.push({ x: offsets.get(x)!, y: offsets.get(y)!, plane, level: own });
      // a symmetric grid's other triangle, from the served one
      if (g.symmetric && x !== y)
        blocks.push({ x: offsets.get(y)!, y: offsets.get(x)!, plane: transposed(plane), level: own });
    }
    const plane = wholePlane(length, level, blocks, g.aggregate ?? 'mean');
    if (bytes + plane.value.byteLength > settings.maxBytes) return undefined;
    bytes += plane.value.byteLength;
    wholes.set(key, plane);
    return plane;
  };

  const shared = sharedBuild(build);
  return {
    get bytes() {
      return bytes;
    },
    /** The build while it runs (`sharedBuild`). */
    get building() {
      return shared.building;
    },
    async read(request: RegionRequest, signal: AbortSignal): Promise<Built<GridDense>> {
      const p = shared.start();
      p.then((m) => {
        if (!(m instanceof Map)) shared.forget(p);
      });
      const made = await untilAborted(p, signal);
      if (!(made instanceof Map)) return made;
      const y = request.y!;
      const L = request.level;
      const [x0, x1] = [Math.floor(request.start / L), Math.ceil(request.end / L)];
      const [y0, y1] = [Math.floor(y.start / L), Math.ceil(y.end / L)];
      const whole = grids.find((g) => g.whole?.system === request.system);
      if (whole) {
        const plane = wholeAt(whole, made, L);
        if (!plane) return refused('seq.store.grid.levels', 'The whole-assembly level is past the summaries’ bytes');
        return cutGrid(plane, x0, x1, y0, y1);
      }
      const planes = made.get(`${request.system}\u0000${y.system}`);
      const own = planes && sourceLevel([...planes.keys()], L);
      if (!planes || own === undefined) return refused('seq.store.grid.levels', 'No level of the grid serves the tile');
      return cutGrid(planes.get(own)!, x0, x1, y0, y1);
    },
    dispose() {
      reading.abort();
      wholes.clear();
      bytes = 0;
    },
  };
}

function transposed(p: GridPlane): GridPlane {
  const value = new Float32Array(p.value.length);
  for (let i = 0; i < p.nx; i++) for (let j = 0; j < p.ny; j++) value[j * p.nx + i] = p.value[i * p.ny + j];
  return { nx: p.ny, ny: p.nx, value };
}
