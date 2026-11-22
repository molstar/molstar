/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { TimerResult } from './timer.js';

/** The tracks of the profile, in a group of their own */
export const TimerTraceTracks = { group: 'Mol* graphics', gpu: 'GPU passes', cpu: 'CPU, encoding' } as const;

type Measure = (name: string, options: { start: number; end: number; detail: unknown }) => void;

/**
 * The graphics timer's results (`GraphicsTimer.resolve`) as tracks of a browser's performance profile (Chrome
 * DevTools' extensibility API: `performance.measure` with `detail.devtools`): where the passes ran on the GPU, nested
 * in their groups (`TimerTraceTracks.gpu`), and when the CPU encoded them (`TimerTraceTracks.cpu`). Why: the profile's
 * own GPU track shows a WebGPU queue's submissions, not how long its compute and render passes run (a frame of 33 ms
 * of passes looked idle there).
 *
 * The GPU's clock is placed on `performance.now()`'s by the least lag seen: no pass starts on the GPU before the CPU
 * began encoding it, so the clock's offset is at most a pass's GPU start less its CPU start, the least of them over
 * the passes so far (a pass that ran as soon as it was encoded: in a frame the GPU waits for, about its submit). A
 * result without GPU times (another backend, an invalid time) only on the CPU's track. Chrome quantizes WebGPU
 * timestamps to 100 µs unless started with `--enable-webgpu-developer-features`.
 */
export function createTimerTrace(measure: Measure = defaultMeasure) {
  /** the GPU clock's milliseconds less `performance.now()`'s (`Infinity` none seen yet) */
  let offset = Infinity;

  const lag = (r: TimerResult) => {
    if (r.gpuStart !== undefined) offset = Math.min(offset, r.gpuStart / 1e6 - r.cpuStart);
    for (const c of r.children) lag(c);
  };

  const ms = (v: number) => (Number.isFinite(v) ? `${v.toFixed(3)} ms` : 'n/a');

  const emit = (r: TimerResult) => {
    const gpu = r.gpuElapsed / 1e6;
    const properties = [
      ['GPU', ms(gpu)],
      ['CPU', ms(r.cpuElapsed)],
    ];
    const pass = r.label.startsWith('pass ');
    if (r.gpuStart !== undefined && r.gpuEnd !== undefined && Number.isFinite(offset)) {
      measure(r.label, {
        start: r.gpuStart / 1e6 - offset,
        end: r.gpuEnd / 1e6 - offset,
        detail: {
          devtools: {
            dataType: 'track-entry',
            track: TimerTraceTracks.gpu,
            trackGroup: TimerTraceTracks.group,
            // a group's span holds gaps where the GPU ran other work: its own time in its properties
            color: pass ? 'tertiary' : 'primary-light',
            tooltipText: `${r.label}: ${ms(gpu)} on the GPU`,
            properties,
          },
        },
      });
    }
    measure(r.label, {
      start: r.cpuStart,
      end: r.cpuStart + r.cpuElapsed,
      detail: {
        devtools: {
          dataType: 'track-entry',
          track: TimerTraceTracks.cpu,
          trackGroup: TimerTraceTracks.group,
          color: pass ? 'secondary' : 'secondary-light',
          tooltipText: `${r.label}: ${ms(r.cpuElapsed)} encoding`,
          properties,
        },
      },
    });
    for (const c of r.children) emit(c);
  };

  return {
    /** the clock's offset placed so far (milliseconds, `Infinity` none), see `createTimerTrace` */
    get offset() {
      return offset;
    },
    /** `results` onto the tracks (each root once, as `resolve` returns it) */
    report(results: readonly TimerResult[]) {
      for (const r of results) lag(r);
      for (const r of results) emit(r);
    },
  };
}

function defaultMeasure(name: string, options: { start: number; end: number; detail: unknown }) {
  if (typeof performance === 'undefined' || !performance.measure) return;
  try {
    performance.measure(name, options);
  } catch {
    // an end before its start (a clock's jump): left out
  }
}
