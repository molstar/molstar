/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Stage measurement for the pipeline benches: wall time, GC pauses, JS-heap allocation, retained memory after a full
// GC, peak RSS, and an optional CPU or allocation profile per stage. The same method as bench.ts (its header has the
// column definitions); many.ts uses this module.

import { writeFileSync } from 'node:fs';
import { Session } from 'node:inspector/promises';
import { PerformanceObserver } from 'node:perf_hooks';
import v8 from 'node:v8';

export type StageResult = {
  ms: number;
  allocMB: number;
  retainedHeapMB: number;
  retainedAbMB: number;
  rssMB: number;
  gcMs: number;
  /** what one stage is divided by in the per-item column (structures, hovers, ...) */
  count?: number;
  info?: string;
  /** the stage threw: its time and memory columns measure the failure, not the work (perf-audit-3 12/B4) */
  error?: string;
};

const v8gc = (globalThis as any).gc as (options?: object) => void;
if (!v8gc) throw new Error('run with --expose-gc');
/**
 * A full collection whose freed array-buffer backing stores are also accounted: they are swept concurrently, and
 * `memoryUsage().arrayBuffers` drops only after the sweep, so without the wait the next stage gets the credit
 */
export async function gc() {
  v8gc({ type: 'major', execution: 'sync', flavor: 'last-resort' });
  await new Promise((r) => setTimeout(r, 20));
  v8gc({ type: 'major', execution: 'sync', flavor: 'last-resort' });
}

// cumulative JS-heap allocation: `total_allocated_bytes` is in Node 25+; without it the Alloc column is NaN
const hasAllocCounter = typeof (v8.getHeapStatistics() as any).total_allocated_bytes === 'number';
if (!hasAllocCounter) console.error('warning: Node 25+ is needed for the Alloc column (v8 total_allocated_bytes)');
function mem() {
  const m = process.memoryUsage();
  return {
    heap: m.heapUsed,
    ab: m.arrayBuffers,
    alloc: hasAllocCounter ? ((v8.getHeapStatistics() as any).total_allocated_bytes as number) : NaN,
  };
}

/** GC pauses (start, duration); a stage sums those that start inside it */
const gcPauses: [number, number][] = [];
new PerformanceObserver((list) => {
  for (const e of list.getEntries()) gcPauses.push([e.startTime, e.duration]);
}).observe({ entryTypes: ['gc'] });
const gcTime = (t0: number, t1: number) => gcPauses.reduce((a, [t, d]) => (t >= t0 && t <= t1 ? a + d : a), 0);

const matches = (name: string, pattern: string) =>
  pattern.endsWith('*') || pattern.endsWith(':') ? name.startsWith(pattern.replace(/\*$/, '')) : name === pattern;

export type StageRunner = ReturnType<typeof createStageRunner>;

/**
 * `stage(name, results, fn, { required, count })`: a `required` stage runs unmeasured when --only/--skip leave it out,
 * so that later stages do not pay for it; `count` divides the time in the per-item column.
 */
export function createStageRunner(options: { only: string[]; skip: string[]; profileDir?: string; alloc: boolean }) {
  const { only, skip, profileDir, alloc } = options;
  let session: Session | undefined;
  const profiles: Record<string, any> = {};
  let failures = 0;
  const enabled = (name: string) =>
    !skip.some((p) => matches(name, p)) && (!only.length || only.some((p) => matches(name, p)));

  async function stage(
    name: string,
    results: Record<string, StageResult>,
    fn: () => Promise<string | void> | string | void,
    opts: { required?: boolean; count?: number } = {},
  ) {
    if (!enabled(name)) {
      if (opts.required) await fn();
      return;
    }
    await gc();
    const before = mem();
    if (session && alloc) {
      await session.post('HeapProfiler.startSampling', {
        samplingInterval: 64 * 1024,
        includeObjectsCollectedByMajorGC: true,
        includeObjectsCollectedByMinorGC: true,
      });
    } else if (session) await session.post('Profiler.start');
    const t0 = performance.now();
    let info: string | void = undefined;
    let error: string | undefined;
    try {
      info = await fn();
    } catch (e) {
      error = (e as Error)?.message ?? String(e);
      failures++;
      console.error(e);
    }
    const t1 = performance.now();
    await new Promise((r) => setTimeout(r, 0)); // gc entries are delivered asynchronously
    const stageGcMs = gcTime(t0, t1);
    gcPauses.length = 0;
    if (session) {
      const { profile } = alloc ? await session.post('HeapProfiler.stopSampling') : await session.post('Profiler.stop');
      profiles[name] = profile;
    }
    const mid = mem();
    await gc();
    const after = mem();
    const r: StageResult = {
      ms: t1 - t0,
      allocMB: (mid.alloc - before.alloc) / 2 ** 20,
      retainedHeapMB: (after.heap - before.heap) / 2 ** 20,
      retainedAbMB: (after.ab - before.ab) / 2 ** 20,
      rssMB: process.resourceUsage().maxRSS / 1024,
      gcMs: stageGcMs,
      count: opts.count,
      info: info || undefined,
      error,
    };
    results[name] = r;
    if (!session) {
      if (r.error) {
        console.error(`${name.padEnd(36)} FAILED: ${r.error}`);
        return;
      }
      const per = r.count ? `${(r.ms / r.count).toFixed(2).padStart(8)} ms/${r.count}` : ''.padStart(14);
      console.error(
        `${name.padEnd(36)} ${r.ms.toFixed(1).padStart(9)} ms ${per} ${r.allocMB.toFixed(0).padStart(6)} MB alloc ${(r.retainedHeapMB + r.retainedAbMB).toFixed(0).padStart(6)} MB ret ${r.rssMB.toFixed(0).padStart(6)} MB rss ${r.gcMs.toFixed(0).padStart(6)} ms gc  ${r.info ?? ''}`,
      );
    }
  }

  async function startProfiling() {
    if (!profileDir) return;
    session = new Session();
    session.connect();
    if (alloc) await session.post('HeapProfiler.enable');
    else {
      await session.post('Profiler.enable');
      await session.post('Profiler.setSamplingInterval', { interval: 100 });
    }
  }

  function writeProfiles(results: Record<string, StageResult>) {
    if (!profileDir) return;
    for (const [name, p] of Object.entries(profiles))
      writeFileSync(
        `${profileDir}/${name.replace(/[:/]/g, '_')}.${alloc ? 'heapprofile' : 'cpuprofile'}`,
        JSON.stringify(p),
      );
    writeFileSync(`${profileDir}/stages.json`, JSON.stringify(results, null, 1));
    console.log('profiles written to', profileDir);
  }

  /** stages that threw so far; a run with any exits non-zero (`exitOnFailures`) */
  const failureCount = () => failures;

  return { stage, enabled, startProfiling, writeProfiles, failureCount };
}

const median = (xs: number[]) => {
  const a = [...xs].sort((x, y) => x - y);
  return a.length % 2 ? a[(a.length - 1) / 2] : (a[a.length / 2 - 1] + a[a.length / 2]) / 2;
};

/**
 * Median of runs per stage, the JSON row format of bench.ts. A stage that failed in any run has no numbers (null), so
 * that a comparison does not read the failure as a speedup, and its `error`.
 */
export function summarize(all: Record<string, StageResult[]>) {
  return Object.entries(all).map(([k, v]) => {
    const error = v.find((x) => x.error)?.error;
    if (error) return { stage: k, error, coldMs: null, medMs: null, allocMB: null, info: v[0].info ?? '' };
    return summarizeRow(k, v);
  });
}

function summarizeRow(k: string, v: StageResult[]) {
  return {
    stage: k,
    coldMs: +v[0].ms.toFixed(1),
    medMs: +median(v.map((x) => x.ms)).toFixed(1),
    perItemMs: v[0].count ? +(median(v.map((x) => x.ms)) / v[0].count).toFixed(3) : undefined,
    allocMB: +median(v.map((x) => x.allocMB)).toFixed(1),
    retHeapMB: +median(v.map((x) => x.retainedHeapMB)).toFixed(1),
    retAbMB: +median(v.map((x) => x.retainedAbMB)).toFixed(1),
    rssMB: +Math.max(...v.map((x) => x.rssMB)).toFixed(0),
    gcMs: +median(v.map((x) => x.gcMs)).toFixed(0),
    info: v[0].info ?? '',
  };
}

/** Exits non-zero after a run in which a stage threw: the table is still written, but the run is not a result */
export function exitOnFailures(failures: number) {
  if (failures === 0) return;
  console.error(`${failures} stage${failures === 1 ? '' : 's'} failed`);
  process.exitCode = 1;
}
