/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Canvas3DContext } from '@molstar/graphics/canvas3d/canvas3d';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { TimerResult } from '@molstar/graphics/gl/timer';
import { isWebGLContext } from '@molstar/graphics/gl/webgl/context';

/** CPU, or the GPU path of a backend */
export type BenchmarkPath = 'cpu' | 'webgl' | 'webgpu';

/** One run of a path: wall times of its stages (GPU work completed before each is taken) and its output */
export interface PathRun {
  readonly stages: { readonly label: string; readonly ms: number }[];
  readonly triangles: number;
  /** byte count of the result, where known */
  readonly bytes?: number;
  /** draw calls of a frame, of render cases */
  readonly drawCalls?: number;
  /** GPU time per frame of the timed pass groups (`--gpu-timing`), by path of group labels, ms */
  readonly gpuPasses?: { readonly [path: string]: number };
  /** what a case measured besides the times, printed with them */
  readonly notes?: readonly string[];
}

/** The context of the page, for cases that render with a `Canvas3D`; `gpuTiming` with timing mode on */
export const BenchmarkEnv: {
  context?: Canvas3DContext;
  gpuTiming?: boolean;
  /** every render case with draw lists or without (`RendererProps.drawLists`, `--draw-lists`, `--no-draw-lists`) */
  drawLists?: boolean;
  /** the canvas after the last run of each render case, as PNG data URLs by case id (`--images`) */
  images?: { [id: string]: string };
} = {};

/**
 * Waits until the GPU completed the submitted work: on WebGL with a 1-pixel readback (the fence of
 * `waitForGpuCommandsComplete` is polled and adds about 90 ms), on WebGPU until the submitted work is done
 */
export async function waitForGpu(graphics: GraphicsContext) {
  if (isWebGLContext(graphics)) graphics.waitForGpuCommandsCompleteSync();
  else await graphics.waitForGpuCommandsComplete();
}

/** Times the stages of a run, waiting for the GPU before taking each time (see `waitForGpu`) */
export class StageTimer {
  readonly stages: { label: string; ms: number }[] = [];
  private last = performance.now();

  constructor(private readonly graphics?: GraphicsContext) {}

  async stage(label: string) {
    if (this.graphics) await waitForGpu(this.graphics);
    const now = performance.now();
    this.stages.push({ label, ms: now - this.last });
    this.last = now;
  }
}

/**
 * With `--gpu-timing`, times the GPU passes of `fn` in a timer group: the GPU time per pass label summed over the run
 * (repeated passes, e.g. one per unit, add up, their count in the label), ms; waits for the timestamp queries
 */
export async function withGpuPasses<T>(
  graphics: GraphicsContext,
  fn: () => Promise<T>,
): Promise<{ result: T; gpuPasses?: { [label: string]: number } }> {
  if (!BenchmarkEnv.gpuTiming) return { result: await fn() };
  const label = 'benchmark run';
  graphics.timer.group(label);
  let result: T;
  try {
    result = await fn();
  } finally {
    graphics.timer.groupEnd(label);
  }
  const sums: { [label: string]: { ms: number; count: number } } = {};
  const add = (results: readonly TimerResult[]) => {
    for (const r of results) {
      if (r.children.length === 0 && !isNaN(r.gpuElapsed)) {
        const s = (sums[r.label] ??= { ms: 0, count: 0 });
        s.ms += r.gpuElapsed / 1e6;
        s.count += 1;
      }
      add(r.children);
    }
  };
  for (let i = 0; i < 100; ++i) {
    await waitForGpu(graphics);
    const results = graphics.timer.resolve().filter((r) => r.label === label);
    if (results.length) {
      add(results);
      break;
    }
    await new Promise((resolve) => setTimeout(resolve, 10));
  }
  const gpuPasses: { [label: string]: number } = {};
  let total = 0;
  for (const [k, { ms, count }] of Object.entries(sums)) {
    gpuPasses[count > 1 ? `${k} (${count}x)` : k] = ms;
    total += ms;
  }
  gpuPasses.total = total;
  return { result: result!, gpuPasses };
}

export interface BenchmarkCase<I = unknown> {
  /** `group/name`, e.g. `mc/volume/blobs-128` */
  readonly id: string;
  readonly description: string;
  /** The input of all paths, created once (not timed) */
  readonly setup: () => Promise<I> | I;
  /** Paths of the case; a path returns a skip reason when it can't handle the input */
  readonly paths: {
    readonly [P in BenchmarkPath]?: (graphics: GraphicsContext, input: I) => Promise<PathRun | string>;
  };
}

export interface PathResult {
  readonly id: string;
  readonly path: BenchmarkPath;
  readonly status: 'ok' | 'skipped' | 'error';
  readonly reason?: string;
  /** total wall time of the first run, ms */
  readonly first?: number;
  readonly median?: number;
  readonly min?: number;
  /** median per stage, ms */
  readonly stages?: { readonly [label: string]: number };
  readonly triangles?: number;
  readonly bytes?: number;
  readonly drawCalls?: number;
  readonly gpuPasses?: { readonly [path: string]: number };
  readonly notes?: readonly string[];
}

function median(values: number[]) {
  const s = [...values].sort((a, b) => a - b);
  const m = Math.floor(s.length / 2);
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}

const total = (run: PathRun) => run.stages.reduce((sum, s) => sum + s.ms, 0);

/** A run of `fn` after a full garbage collection (with `--expose-gc`), so runs don't pay for earlier garbage */
function collected<T>(fn: () => Promise<T>) {
  (globalThis as { gc?: () => void }).gc?.();
  return fn();
}

/** The first run, `warmup` runs, then `runs` measured runs of a path */
export async function runPath(
  c: BenchmarkCase<any>,
  path: BenchmarkPath,
  graphics: GraphicsContext,
  input: unknown,
  options: { runs: number; warmup: number },
): Promise<PathResult> {
  const fn = c.paths[path]!;
  try {
    const first = await collected(() => fn(graphics, input));
    if (typeof first === 'string') return { id: c.id, path, status: 'skipped', reason: first };
    for (let i = 0; i < options.warmup; ++i) await collected(() => fn(graphics, input));
    const runs: PathRun[] = [];
    for (let i = 0; i < options.runs; ++i) runs.push((await collected(() => fn(graphics, input))) as PathRun);
    const stages: { [label: string]: number } = {};
    for (const { label } of first.stages)
      stages[label] = median(runs.map((r) => r.stages.find((s) => s.label === label)!.ms));
    return {
      id: c.id,
      path,
      status: 'ok',
      first: total(first),
      median: median(runs.map(total)),
      min: Math.min(...runs.map(total)),
      stages,
      triangles: first.triangles,
      bytes: first.bytes,
      // of the measured runs: the first may draw otherwise while what it needs compiles
      drawCalls: runs[runs.length - 1]?.drawCalls ?? first.drawCalls,
      gpuPasses: runs[runs.length - 1]?.gpuPasses,
      notes: runs[runs.length - 1]?.notes,
    };
  } catch (e) {
    return { id: c.id, path, status: 'error', reason: e instanceof Error ? e.message : String(e) };
  }
}
