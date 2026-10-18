/**
 * Copyright (c) 2022-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { isDebugMode } from '@molstar/core/util/debug';
import { now } from '@molstar/core/util/now';
import type { RenderPassTarget } from './context.js';
import type { GraphicsStats } from './stats.js';

export type TimerCalls = {
  drawInstanced: number;
  counts: number;
};

export type TimerResult = {
  readonly label: string;
  /** nanoseconds, `NaN` if no timed pass began in the group */
  readonly gpuElapsed: number;
  readonly gpuAvg: number;
  readonly cpuElapsed: number;
  readonly cpuAvg: number;
  /** nested groups, then the group's own passes (labelled `pass <label>`) */
  readonly children: TimerResult[];
  readonly calls?: TimerCalls;
  readonly note?: string;
};

export type TimerOptions = { captureStats?: boolean; note?: string };

/**
 * GPU times of passes, summed in groups (`group` to `groupEnd`). A group times the passes that begin inside it; it
 * begins and ends outside of passes (checked in debug mode).
 */
export type GraphicsTimer = {
  /** Check with GPU for finished timers. */
  resolve: () => TimerResult[];
  group: (label: string, options?: TimerOptions) => void;
  groupEnd: (label: string) => void;
  stats: () => { gpu: Record<string, number>; cpu: Record<string, number> };
  formatedStats: () => Record<string, string>;
  clear: () => void;
  destroy: () => void;
};

/** GPU time measurement of a backend: one segment per pass, begun and ended at pass boundaries */
export interface TimerSource {
  /** Begins timing a pass and returns the segment id, or -1 if it can't be timed */
  begin: () => number;
  /** Ends the segment begun last */
  end: () => void;
  /** Segments finished since the last call: id and elapsed nanoseconds, `undefined` if the time is invalid */
  poll: () => Iterable<readonly [number, number | undefined]>;
  /** Forgets all segments */
  clear: () => void;
}

/** The timer of a backend: the `GraphicsTimer` plus the pass boundaries reported by the backend's context */
export type BackendTimer = GraphicsTimer & {
  /** A render pass begins on `target` (ending the open one) */
  passBegin: (target: RenderPassTarget | 'compute', label?: string) => void;
  /** The open render pass ends */
  passEnd: () => void;
};

function movingAverage(avg: number, sample: number, count: number) {
  avg -= avg / count;
  avg += sample / count;
  return avg;
}

class MovingAverage {
  private readonly avgs = new Map<string, number>();

  add(label: string, sample: number) {
    let avg = this.avgs.get(label) || sample;
    avg = movingAverage(avg, sample, this.count);
    this.avgs.set(label, avg);
    return avg;
  }

  get(label: string) {
    return this.avgs.get(label);
  }

  stats() {
    return Object.fromEntries(this.avgs.entries());
  }

  clear() {
    this.avgs.clear();
  }

  constructor(private count: number) {}
}

function getDrawCalls(stats: GraphicsStats) {
  const { calls } = stats;
  return calls.drawInstanced + calls.drawInstancedBase + calls.multiDrawInstancedBase;
}

function clearStats(stats: GraphicsStats) {
  stats.calls.drawInstanced = 0;
  stats.calls.drawInstancedBase = 0;
  stats.calls.multiDrawInstancedBase = 0;
  stats.calls.counts = 0;

  stats.culled.lod = 0;
  stats.culled.frustum = 0;
  stats.culled.occlusion = 0;

  stats.cacheHits.cull = 0;
}

function getPassLabel(target: RenderPassTarget | 'compute', label?: string) {
  if (label) return label;
  if (target === 'drawing-buffer' || target === 'compute') return target;
  return 'framebuffer' in target ? `render-target ${target.id}` : `framebuffer ${target.id}`;
}

type Segment = { elapsed?: number; done: boolean; invalid: boolean };

type Group = {
  label: string;
  parent: Group | undefined;
  children: Group[];
  /** passes that began with this group innermost */
  passes: Pass[];
  /** segments of all passes that began inside the group, including nested groups */
  segments: number[];
  cpu: { start: number; end: number };
  captureStats: boolean;
  calls?: TimerCalls;
  note?: string;
  /** root only: passes begun inside it whose timing has not ended */
  openPasses: number;
};

type Pass = {
  label: string;
  /** innermost group when the pass began */
  group: Group;
  root: Group;
  /** groups open when the pass began */
  groups: Group[];
  segment: number;
  cpu: { start: number; end: number };
};

function sameGroups(a: readonly Group[], b: readonly Group[]) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; ++i) if (a[i] !== b[i]) return false;
  return true;
}

/**
 * Timer core shared by the backends: the passes timed by `source` are summed in the groups open when they began.
 * A pass's timing ends with the pass, or when its innermost group ends (it has no draws after that).
 *
 * In debug mode, draws are checked against the placement rule: a pass draws only while exactly the groups open at its
 * begin are open (`stats.calls` changes are noted at group and pass events).
 *
 * Without `source` (no GPU timing) every method does nothing.
 */
export function createGraphicsTimer(
  source: TimerSource | undefined,
  stats: GraphicsStats,
  options?: { avgCount?: number },
): BackendTimer {
  const avgCount = options?.avgCount ?? 30;

  const segments = new Map<number, Segment>();
  const pending = new Map<string, Group>();
  const stack: Group[] = [];
  const gpuAvgs = new MovingAverage(avgCount);
  const cpuAvgs = new MovingAverage(avgCount);

  /** ended root groups waiting for their GPU times */
  let roots: Group[] = [];
  let capturingStats = false;
  /** the timed pass */
  let currentPass: Pass | undefined = undefined;
  /** label of the open pass on the GPU (timed or not), for the placement check */
  let openPassLabel: string | undefined = undefined;

  // draw calls are monotonic here, `clearStats` resets the stats
  let drawOffset = 0;
  let lastDraws = 0;
  const warned = new Set<string>();

  const warnOnce = (message: string) => {
    if (warned.has(message)) return;
    warned.add(message);
    console.warn(message);
  };

  /** Placement check of the draws since the last event, under the groups before the event changes them */
  const checkDraws = () => {
    if (!isDebugMode) return;
    const draws = drawOffset + getDrawCalls(stats);
    if (draws === lastDraws) return;
    lastDraws = draws;
    if (!stack.length || openPassLabel === undefined) return;
    if (!currentPass || !sameGroups(currentPass.groups, stack)) {
      const group = stack[stack.length - 1].label;
      warnOnce(`Timer group '${group}' begins or ends inside pass '${openPassLabel}'`);
    }
  };

  const endPassTiming = () => {
    const pass = currentPass;
    if (!pass) return;
    if (pass.segment !== -1) source!.end();
    pass.cpu.end = now();
    pass.root.openPasses -= 1;
    currentPass = undefined;
  };

  const clear = () => {
    pending.clear();
    stack.length = 0;
    gpuAvgs.clear();
    cpuAvgs.clear();
    segments.clear();
    roots = [];
    capturingStats = false;
    currentPass = undefined;
    openPassLabel = undefined;
    source?.clear();
  };

  const getGpuElapsed = (ids: readonly number[]) => {
    if (ids.length === 0) return NaN;
    let elapsed = 0;
    for (const id of ids) elapsed += segments.get(id)!.elapsed!;
    return elapsed;
  };

  const getResult = (g: Group, path: string): TimerResult => {
    const gpuElapsed = getGpuElapsed(g.segments);
    const cpuElapsed = g.cpu.end - g.cpu.start;
    const key = path ? `${path}/${g.label}` : g.label;
    const children = g.children.map((c) => getResult(c, key));
    for (const p of g.passes) {
      const passElapsed = getGpuElapsed(p.segment === -1 ? [] : [p.segment]);
      const passCpu = p.cpu.end - p.cpu.start;
      const passKey = `${key}/pass ${p.label}`;
      children.push({
        label: `pass ${p.label}`,
        gpuElapsed: passElapsed,
        gpuAvg: isNaN(passElapsed) ? NaN : gpuAvgs.add(passKey, passElapsed),
        cpuElapsed: passCpu,
        cpuAvg: cpuAvgs.add(passKey, passCpu),
        children: [],
      });
    }
    return {
      label: g.label,
      gpuElapsed,
      gpuAvg: isNaN(gpuElapsed) ? NaN : gpuAvgs.add(g.label, gpuElapsed),
      cpuElapsed,
      cpuAvg: cpuAvgs.add(g.label, cpuElapsed),
      children,
      calls: g.calls,
      note: g.note,
    };
  };

  return {
    resolve: () => {
      const results: TimerResult[] = [];
      if (!source || !roots.length || capturingStats) return results;

      for (const [id, elapsed] of source.poll()) {
        const segment = segments.get(id);
        if (!segment) continue;
        segment.done = true;
        if (elapsed === undefined) segment.invalid = true;
        else segment.elapsed = elapsed;
      }

      const unresolved: Group[] = [];
      for (const root of roots) {
        // a root's segments include those of all its passes
        if (root.openPasses > 0 || !root.segments.every((id) => segments.get(id)!.done)) {
          unresolved.push(root);
          continue;
        }
        // results with invalid (disjoint) times are dropped
        if (!root.segments.some((id) => segments.get(id)!.invalid)) results.push(getResult(root, ''));
        for (const id of root.segments) segments.delete(id);
      }
      roots = unresolved;
      return results;
    },
    group: (label: string, options?: TimerOptions) => {
      if (!source) return;

      if (pending.has(label)) {
        throw new Error(`Timer group '${label}' already exists`);
      }

      checkDraws();
      const captureStats = options?.captureStats ?? false;
      const parent = stack.length ? stack[stack.length - 1] : undefined;
      const group: Group = {
        label,
        parent,
        children: [],
        passes: [],
        segments: [],
        cpu: { start: now(), end: -1 },
        captureStats,
        openPasses: 0,
      };
      if (options?.note) group.note = options.note;
      pending.set(label, group);
      parent?.children.push(group);
      stack.push(group);

      if (captureStats) {
        if (capturingStats) {
          throw new Error('Already capturing stats');
        }
        drawOffset += getDrawCalls(stats);
        clearStats(stats);
        capturingStats = true;
      }
    },
    groupEnd: (label: string) => {
      if (!source) return;

      const group = pending.get(label);
      if (!group) {
        throw new Error(`Timer group '${label}' does not exist`);
      }

      if (stack[stack.length - 1]?.label !== label) {
        throw new Error(`Timer group '${label}' has a pending nested group`);
      }

      checkDraws();
      // the pass has no draws after its innermost group (placement rule)
      if (currentPass?.group === group) endPassTiming();
      stack.pop();
      pending.delete(label);

      group.cpu.end = now();
      if (group.captureStats) {
        group.calls = { ...stats.calls };
        capturingStats = false;
      }

      if (!group.parent) roots.push(group);
    },
    passBegin: (target: RenderPassTarget | 'compute', label?: string) => {
      if (!source) return;
      checkDraws();
      endPassTiming();
      if (!stack.length) {
        openPassLabel = isDebugMode ? getPassLabel(target, label) : '';
        return;
      }
      const passLabel = getPassLabel(target, label);
      openPassLabel = passLabel;
      const group = stack[stack.length - 1];
      const root = stack[0];
      const segment = source.begin();
      currentPass = {
        label: passLabel,
        group,
        root,
        groups: stack.slice(),
        segment,
        cpu: { start: now(), end: -1 },
      };
      root.openPasses += 1;
      group.passes.push(currentPass);
      if (segment !== -1) {
        segments.set(segment, { done: false, invalid: false });
        for (const g of stack) g.segments.push(segment);
      }
    },
    passEnd: () => {
      if (!source) return;
      checkDraws();
      endPassTiming();
      openPassLabel = undefined;
    },
    stats: () => {
      return {
        gpu: gpuAvgs.stats(),
        cpu: cpuAvgs.stats(),
      };
    },
    formatedStats: () => {
      const stats: Record<string, string> = {};
      const gpu = gpuAvgs.stats();
      const cpu = cpuAvgs.stats();
      for (const l of Object.keys(cpu)) {
        const g = gpu[l] === undefined ? 'n/a' : `${(gpu[l] / 1000 / 1000).toFixed(2)}`;
        const c = `${cpu[l].toFixed(2)}`;
        stats[l] = `${g} ms | CPU: ${c} ms`;
      }
      return stats;
    },
    clear,
    destroy: () => {
      clear();
    },
  };
}

function formatTimerResult(result: TimerResult) {
  const ms = (ns: number) => (isNaN(ns) ? 'n/a' : `${(ns / 1000 / 1000).toFixed(2)}`);
  const gpu = ms(result.gpuElapsed);
  const gpuAvg = ms(result.gpuAvg);
  const cpu = `${result.cpuElapsed.toFixed(2)}`;
  const cpuAvg = `${result.cpuAvg.toFixed(2)}`;
  return `${result.label} ${gpu} ms (avg. ${gpuAvg} ms) | CPU: ${cpu} ms (avg. ${cpuAvg} ms)`;
}

export function printTimerResults(results: TimerResult[]) {
  results.map((r) => {
    const f = formatTimerResult(r);
    if (r.children.length || r.calls || r.note) {
      console.groupCollapsed(f);
      if (r.calls) console.log(r.calls);
      if (r.note) console.log(r.note);
      printTimerResults(r.children);
      console.groupEnd();
    } else {
      console.log(f);
    }
  });
}
