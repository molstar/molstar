/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import { SourceIssues } from './issues.js';
import {
  abortError,
  type Coverage,
  type DynamicNode,
  type OpenedSource,
  type RegionAnswer,
  type RegionRequest,
  type Tables,
  type ValuesBins,
  type ValuesPositions,
} from './protocol.js';
import { checkAnswer, DefaultMaxRows } from './tables.js';

// Value summaries (_spec/README.md "Data", Summaries): a source that serves level 1 only gets coarser levels for its
// values, so a values layer zoomed out reads a coarse level and stays within the store's tile bound. A system of known
// length up to `maxLength` gets a pyramid: read whole once at level 1 (checked as any answer), every level built from
// the one below, then the level-1 values let go; its bytes count against the store's bound. A longer system's coarse
// tile is binned on the fly from a level-1 read of that tile alone, at most `maxRows` positions.

/** How the store summarizes a level-1 source's values: a level `factor` times the one below, up to `maxLength`. */
export interface SummarySettings {
  /** Each level's bin width over the one below: 16. */
  readonly factor: number;
  /** The longest system summarized (positions); longer ones are read at level 1 only: 2^24. */
  readonly maxLength: number;
}

export const DefaultSummarySettings: SummarySettings = { factor: 16, maxLength: 1 << 24 };

function emptyBins(start: number, n: number): ValuesBins {
  const min = new Float32Array(n).fill(Number.NaN);
  const max = new Float32Array(n).fill(Number.NaN);
  const mean = new Float32Array(n).fill(Number.NaN);
  return { form: 'bins', start, min, max, mean, count: new Float32Array(n) };
}

/**
 * Per-position values summarized into bins of `width` over `[from, length)` (`from` aligned to `width`, 0 by default);
 * NaN is not counted.
 */
export function binValues(t: ValuesPositions, width: number, length: number, from = 0): ValuesBins {
  const out = emptyBins(from, Math.max(0, Math.ceil((length - from) / width)));
  const sum = new Float64Array(out.count.length);
  for (let i = 0; i < t.value.length; i++) {
    const v = t.value[i];
    const p = t.start + i;
    if (Number.isNaN(v) || p < from || p >= length) continue;
    const b = Math.floor((p - from) / width);
    if (out.count[b] === 0 || v < out.min[b]) out.min[b] = v;
    if (out.count[b] === 0 || v > out.max[b]) out.max[b] = v;
    out.count[b]++;
    sum[b] += v;
  }
  for (let b = 0; b < sum.length; b++) if (out.count[b] > 0) out.mean[b] = sum[b] / out.count[b];
  return out;
}

/** Bins `factor` times wider from finer ones starting at 0: min of mins, max of maxes, means weighted by count. */
export function coarsenBins(t: ValuesBins, factor: number): ValuesBins {
  const out = emptyBins(0, Math.ceil(t.count.length / factor));
  const sum = new Float64Array(out.count.length);
  for (let i = 0; i < t.count.length; i++) {
    const c = t.count[i];
    if (c === 0) continue;
    const b = Math.floor(i / factor);
    if (out.count[b] === 0 || t.min[i] < out.min[b]) out.min[b] = t.min[i];
    if (out.count[b] === 0 || t.max[i] > out.max[b]) out.max[b] = t.max[i];
    out.count[b] += c;
    sum[b] += t.mean[i] * c;
  }
  for (let b = 0; b < sum.length; b++) if (out.count[b] > 0) out.mean[b] = sum[b] / out.count[b];
  return out;
}

/** The bins of `[start, end)` (aligned to `width`) cut from a whole system's. */
function cutBins(t: ValuesBins, width: number, start: number, end: number): ValuesBins {
  const a = Math.min(t.count.length, Math.floor(start / width));
  const b = Math.min(t.count.length, Math.max(a, Math.ceil(end / width)));
  return {
    form: 'bins',
    start: a * width,
    min: t.min.subarray(a, b),
    max: t.max.subarray(a, b),
    mean: t.mean.subarray(a, b),
    count: t.count.subarray(a, b),
  };
}

/** The coarse levels over a source's: `[1, f, f², …]` up to the longest values system of known length. */
function summaryLevels(coverage: Coverage, settings: SummarySettings): number[] | undefined {
  if (coverage.levels.length !== 1 || coverage.levels[0] !== 1 || settings.factor < 2) return undefined;
  let longest = 0;
  for (const s of coverage.shapes) {
    if (s.shape !== 'values') continue;
    const length = coverage.systems.find((x) => x.name === s.system)?.length;
    if (length !== undefined) longest = Math.max(longest, length);
  }
  const levels = [1];
  for (let l = settings.factor; l < longest; l *= settings.factor) levels.push(l);
  return levels.length > 1 ? levels : undefined;
}

/** A source with summaries: what its pyramids hold, which the store counts against its byte bound. */
export interface SummarizedSource extends OpenedSource {
  readonly summaryBytes: number;
}

/** What a source holds besides the store's entries: its pyramids' bytes (0 for one without). */
export function summaryBytes(source: OpenedSource): number {
  return (source as Partial<SummarizedSource>).summaryBytes ?? 0;
}

/** Rejects with an `AbortError` as soon as `signal` aborts, whatever `p` does. */
function untilAborted<T>(p: Promise<T>, signal: AbortSignal): Promise<T> {
  if (signal.aborted) return Promise.reject(abortError());
  return new Promise<T>((resolve, reject) => {
    const onAbort = () => reject(abortError());
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

type Pyramid = Map<number, ValuesBins>;

/**
 * A level-1 source with coarse levels for its values: its coverage lists them; a read at level 1 is the source's own;
 * at a coarser level the values come from the system's pyramid (up to `maxLength`) or are binned from the tile's level-1
 * values (longer), and the other shapes are the source's at level 1 over the same region. The source itself when there
 * is nothing to summarize.
 */
export function withSummaries(
  inner: OpenedSource,
  settings: SummarySettings = DefaultSummarySettings,
  maxRows = DefaultMaxRows,
): OpenedSource {
  const levels = summaryLevels(inner.coverage, settings);
  if (!levels) return inner;
  const coverage: Coverage = { ...inner.coverage, levels };
  /** Pyramids being built (the shared whole-system read), and built: the level-1 values are not kept. */
  const building = new Map<string, Promise<Pyramid | RegionAnswer>>();
  const pyramids = new Map<string, Pyramid>();
  let bytes = 0;
  const reading = new AbortController();
  const infoOf = (system: string) => coverage.systems.find((s) => s.name === system || s.aliases?.includes(system));

  const failure = (why: string): RegionAnswer => ({
    ok: false,
    kind: 'malformed',
    issue: issueOf(SourceIssues, 'seq.source.malformed', `The source answered a region malformed: ${why}`),
  });

  /** Reads the whole system once (by none of the callers' signals: it is shared) and builds every level from it. */
  const pyramid = (system: string, length: number): Promise<Pyramid | RegionAnswer> => {
    const built = pyramids.get(system);
    if (built) return Promise.resolve(built);
    let p = building.get(system);
    if (p) return p;
    const request: RegionRequest = { system, start: 0, end: length, level: 1, shapes: ['values'] };
    p = inner.read(request, reading.signal).then((a) => {
      const bad = checkAnswer(a, request, maxRows);
      if (bad) return failure(bad);
      if (!a.ok) return a;
      const v = a.tables.values;
      if (v?.form !== 'positions') return failure('no values at level 1');
      const made: Pyramid = new Map();
      let below: ValuesBins | undefined;
      for (const level of levels.slice(1)) {
        below = below ? coarsenBins(below, settings.factor) : binValues(v, level, length);
        made.set(level, below);
        bytes += 4 * 4 * below.count.length;
      }
      pyramids.set(system, made);
      return made;
    });
    // kept while it is being built; a failure is not kept, so the store's retry asks again
    const done = () => building.delete(system);
    p.then(done, done);
    building.set(system, p);
    return p;
  };

  const values = async (request: RegionRequest, signal: AbortSignal): Promise<ValuesBins | RegionAnswer> => {
    const info = infoOf(request.system);
    const { level, start } = request;
    if (info?.length !== undefined && info.length <= settings.maxLength) {
      const made = await untilAborted(pyramid(info.name, info.length), signal);
      if (!(made instanceof Map)) return made;
      return cutBins(made.get(level)!, level, start, request.end);
    }
    // longer (or of unknown length): this tile's own level-1 values, binned on the fly
    const end = info?.length === undefined ? request.end : Math.min(request.end, info.length);
    if (end - start > maxRows) return failure(`a tile of ${end - start} positions to summarize, more than ${maxRows}`);
    const one: RegionRequest = {
      system: request.system,
      start,
      end: Math.max(start, end),
      level: 1,
      shapes: ['values'],
    };
    const a = await inner.read(one, signal);
    const bad = checkAnswer(a, one, maxRows);
    if (bad) return failure(bad);
    if (!a.ok) return a;
    const v = a.tables.values;
    return v?.form === 'positions' ? binValues(v, level, end, start) : failure('no values at level 1');
  };

  const read = async (request: RegionRequest, signal: AbortSignal): Promise<RegionAnswer> => {
    if (request.level === 1 || !levels.includes(request.level)) return inner.read(request, signal);
    const others = request.shapes.filter((s) => s !== 'values');
    const tables: { -readonly [K in keyof Tables]: Tables[K] } = {};
    let nodes: readonly DynamicNode[] | undefined;
    if (others.length) {
      const a = await inner.read({ ...request, level: 1, shapes: others }, signal);
      if (!a.ok) return a;
      Object.assign(tables, a.tables);
      if (a.nodes?.length) nodes = a.nodes;
    }
    if (request.shapes.includes('values')) {
      const v = await values(request, signal);
      if ('ok' in v) return v;
      tables.values = v;
    }
    return nodes ? { ok: true, tables, nodes } : { ok: true, tables };
  };

  const out: SummarizedSource = {
    coverage,
    get summaryBytes() {
      return bytes;
    },
    read(request, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      return read(request, signal);
    },
    dispose() {
      reading.abort();
      pyramids.clear();
      bytes = 0;
      inner.dispose();
    },
  };
  return out;
}
