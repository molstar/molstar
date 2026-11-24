/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { TimerResult } from '../timer.js';
import { createTimerTrace, TimerTraceTracks } from '../timer-trace.js';

type Entry = { name: string; start: number; end: number; track: string };

/** a trace writing into `entries` */
function setup() {
  const entries: Entry[] = [];
  const trace = createTimerTrace((name, { start, end, detail }) => {
    entries.push({ name, start, end, track: (detail as { devtools: { track: string } }).devtools.track });
  });
  return { trace, entries };
}

/** a result (GPU times in ns on the GPU's clock, CPU in ms) */
function result(
  label: string,
  gpu: [number, number] | undefined,
  cpu: [number, number],
  children: TimerResult[] = [],
): TimerResult {
  return {
    label,
    gpuElapsed: gpu ? gpu[1] - gpu[0] : NaN,
    gpuAvg: NaN,
    ...(gpu ? { gpuStart: gpu[0], gpuEnd: gpu[1] } : undefined),
    cpuStart: cpu[0],
    cpuElapsed: cpu[1] - cpu[0],
    cpuAvg: NaN,
    children,
  };
}

describe('timer trace', () => {
  it("places the GPU's passes on the CPU's clock by the least lag, the CPU's spans as they are", () => {
    const { trace, entries } = setup();
    // the GPU's clock 5000 ms ahead: pass X ran 2 ms after its encoding began, Y as it began (the least lag)
    const frame = result(
      'frame',
      [5_102e6, 5_116e6],
      [100, 104],
      [result('pass X', [5_102e6, 5_110e6], [100, 101]), result('pass Y', [5_112e6, 5_116e6], [112, 113])],
    );
    trace.report([frame]);
    expect(trace.offset).toBe(5000);
    const gpu = entries.filter((e) => e.track === TimerTraceTracks.gpu);
    expect(gpu.map((e) => [e.name, e.start, e.end])).toEqual([
      ['frame', 102, 116],
      ['pass X', 102, 110],
      ['pass Y', 112, 116],
    ]);
    const cpu = entries.filter((e) => e.track === TimerTraceTracks.cpu);
    expect(cpu.map((e) => [e.name, e.start, e.end])).toEqual([
      ['frame', 100, 104],
      ['pass X', 100, 101],
      ['pass Y', 112, 113],
    ]);

    // a later frame lagging less: the offset tightened, its passes placed by it
    entries.length = 0;
    trace.report([result('frame', [5_199e6, 5_201e6], [200, 201])]);
    expect(trace.offset).toBe(4999);
    expect(entries.find((e) => e.track === TimerTraceTracks.gpu)).toMatchObject({ start: 200, end: 202 });
  });

  it("puts a result without GPU times on the CPU's track alone", () => {
    const { trace, entries } = setup();
    trace.report([result('frame', undefined, [10, 12], [result('pass X', undefined, [10, 11])])]);
    expect(trace.offset).toBe(Infinity);
    expect(entries.map((e) => [e.name, e.track])).toEqual([
      ['frame', TimerTraceTracks.cpu],
      ['pass X', TimerTraceTracks.cpu],
    ]);
  });
});
