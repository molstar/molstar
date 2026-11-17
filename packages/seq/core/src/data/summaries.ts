/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

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

// Value summaries (_spec/README.md "Data", Summaries): a source that serves level 1 only gets a pyramid of coarser
// levels for its values, built per level (each from the level below, never from level 1 again) the first time a level
// is read, so a values layer zoomed out reads a coarse level and stays within the store's tile bound. Only for systems
// of known length up to `maxLength`: the pyramid is built over a whole system, read once at level 1.

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

/** Per-position values summarized into bins of `width` from position 0 over `length` positions; NaN is not counted. */
export function binValues(t: ValuesPositions, width: number, length: number): ValuesBins {
  const out = emptyBins(0, Math.ceil(length / width));
  const sum = new Float64Array(out.count.length);
  for (let i = 0; i < t.value.length; i++) {
    const v = t.value[i];
    const p = t.start + i;
    if (Number.isNaN(v) || p < 0 || p >= length) continue;
    const b = Math.floor(p / width);
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

/** The pyramid's levels over a source's: `[1, f, f², …]` up to the longest values system it summarizes. */
function pyramidLevels(coverage: Coverage, settings: SummarySettings): number[] | undefined {
  if (coverage.levels.length !== 1 || coverage.levels[0] !== 1 || settings.factor < 2) return undefined;
  let longest = 0;
  for (const s of coverage.shapes) {
    if (s.shape !== 'values') continue;
    const length = coverage.systems.find((x) => x.name === s.system)?.length;
    if (length !== undefined && length <= settings.maxLength) longest = Math.max(longest, length);
  }
  const levels = [1];
  for (let l = settings.factor; l < longest; l *= settings.factor) levels.push(l);
  return levels.length > 1 ? levels : undefined;
}

/**
 * A level-1 source with a value pyramid: its coverage lists the pyramid's levels; a read at level 1 is the source's
 * own; at a coarser level, values come from the pyramid (the whole system read once at level 1, each level built from
 * the one below on first use) and the other shapes from the source at level 1 over the same region. The source itself
 * when there is nothing to summarize.
 */
export function withSummaries(inner: OpenedSource, settings: SummarySettings = DefaultSummarySettings): OpenedSource {
  const levels = pyramidLevels(inner.coverage, settings);
  if (!levels) return inner;
  const coverage: Coverage = { ...inner.coverage, levels };
  const whole = new Map<string, Promise<RegionAnswer>>();
  const pyramids = new Map<string, Map<number, ValuesBins>>();
  const reading = new AbortController();
  const lengthOf = (system: string) => coverage.systems.find((s) => s.name === system || s.aliases?.includes(system));

  /** The whole system's level-1 values, read once (by none of the callers' signals: it is shared). */
  const level1 = (system: string, length: number) => {
    let p = whole.get(system);
    if (!p) {
      const request: RegionRequest = { system, start: 0, end: length, level: 1, shapes: ['values'] };
      p = inner.read(request, reading.signal);
      // a failed read is not kept: the store's retry asks again
      p.then(
        (a) => !a.ok && whole.delete(system),
        () => whole.delete(system),
      );
      whole.set(system, p);
    }
    return p;
  };

  const bins = (system: string, values: ValuesPositions, length: number, level: number): ValuesBins => {
    let byLevel = pyramids.get(system);
    if (!byLevel) pyramids.set(system, (byLevel = new Map()));
    const hit = byLevel.get(level);
    if (hit) return hit;
    const below = level / settings.factor;
    const made =
      below === 1
        ? binValues(values, level, length)
        : coarsenBins(bins(system, values, length, below), settings.factor);
    byLevel.set(level, made);
    return made;
  };

  const read = async (request: RegionRequest, signal: AbortSignal): Promise<RegionAnswer> => {
    if (request.level === 1 || !levels.includes(request.level)) return inner.read(request, signal);
    const info = lengthOf(request.system);
    const others = request.shapes.filter((s) => s !== 'values');
    const tables: { -readonly [K in keyof Tables]: Tables[K] } = {};
    let nodes: readonly DynamicNode[] | undefined;
    if (others.length) {
      const a = await inner.read({ ...request, level: 1, shapes: others }, signal);
      if (!a.ok) return a;
      Object.assign(tables, a.tables);
      if (a.nodes?.length) nodes = a.nodes;
    }
    if (request.shapes.includes('values') && info?.length !== undefined && info.length <= settings.maxLength) {
      const a = await level1(info.name, info.length);
      if (signal.aborted) throw abortError();
      if (!a.ok) return a;
      const v = a.tables.values;
      if (v?.form === 'positions')
        tables.values = cutBins(
          bins(info.name, v, info.length, request.level),
          request.level,
          request.start,
          request.end,
        );
    }
    return nodes ? { ok: true, tables, nodes } : { ok: true, tables };
  };

  return {
    coverage,
    read(request, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      return read(request, signal);
    },
    dispose() {
      reading.abort();
      inner.dispose();
    },
  };
}
