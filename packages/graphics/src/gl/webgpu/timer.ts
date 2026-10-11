/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { TimerSource } from '../timer.js';
import { BufferUsage, MapMode } from './constants.js';

/** Timestamp queries per submission (two per pass encoder) */
const QueryCount = 4096;

/** A pass segment: the summed durations of its pass encoders */
type Segment = { sum: number; pending: number; ended: boolean; invalid: boolean };

/** `TimerSource` from pass timestamps, plus the hooks of the context that encodes the passes */
export interface WebGPUTimerSource extends TimerSource {
  /** Timestamp writes for a pass encoder about to begin, `undefined` outside of timed passes */
  getTimestampWrites: () => GPURenderPassTimestampWrites | undefined;
  /** Encodes the resolve of the timestamps written since the last submission; call the result after submitting */
  beforeSubmit: (encoder: GPUCommandEncoder) => (() => void) | undefined;
  destroy: () => void;
}

/**
 * Times passes with `timestampWrites` on their encoders: a pass can have several encoders (copies end an encoder and
 * the next draw continues the pass), its time is their sum. Needs the `timestamp-query` feature.
 */
export function createWebGPUTimerSource(device: GPUDevice): WebGPUTimerSource | undefined {
  if (!device.features.has('timestamp-query')) return undefined;

  // created on the first timed pass, nothing is allocated while timing is off
  let querySet: GPUQuerySet | undefined;
  let resolveBuffer: GPUBuffer | undefined;
  /** readback buffers not in use */
  const readbacks: GPUBuffer[] = [];

  const segments = new Map<number, Segment>();
  const finished: [number, number | undefined][] = [];
  /** encoders of the current submission: their segment and first query index */
  let encoders: { segment: number; index: number }[] = [];
  let used = 0;
  let nextId = 0;
  let current = -1;

  const finish = (id: number) => {
    const s = segments.get(id);
    if (!s || !s.ended || s.pending > 0) return;
    finished.push([id, s.invalid ? undefined : s.sum]);
    segments.delete(id);
  };

  return {
    begin: () => {
      querySet ??= device.createQuerySet({ type: 'timestamp', count: QueryCount });
      resolveBuffer ??= device.createBuffer({
        size: QueryCount * 8,
        usage: BufferUsage.QUERY_RESOLVE | BufferUsage.COPY_SRC,
      });
      const id = nextId++;
      segments.set(id, { sum: 0, pending: 0, ended: false, invalid: false });
      current = id;
      return id;
    },
    end: () => {
      const s = segments.get(current);
      if (s) s.ended = true;
      const id = current;
      current = -1;
      finish(id);
    },
    poll: () => finished.splice(0),
    clear: () => {
      // readbacks in flight ignore the forgotten segments
      segments.clear();
      finished.length = 0;
      current = -1;
    },
    getTimestampWrites: () => {
      const s = segments.get(current);
      if (!s) return undefined;
      if (used + 2 > QueryCount) {
        s.invalid = true;
        return undefined;
      }
      const index = used;
      used += 2;
      s.pending += 1;
      encoders.push({ segment: current, index });
      return { querySet: querySet!, beginningOfPassWriteIndex: index, endOfPassWriteIndex: index + 1 };
    },
    beforeSubmit: (encoder: GPUCommandEncoder) => {
      if (used === 0 || !querySet || !resolveBuffer) return undefined;
      const count = used;
      const records = encoders;
      encoders = [];
      used = 0;

      const size = count * 8;
      encoder.resolveQuerySet(querySet, 0, count, resolveBuffer, 0);
      const readback =
        readbacks.pop() ?? device.createBuffer({ size: QueryCount * 8, usage: BufferUsage.MAP_READ | BufferUsage.COPY_DST });
      encoder.copyBufferToBuffer(resolveBuffer, 0, readback, 0, size);

      const done = (times: BigUint64Array | undefined) => {
        for (const { segment, index } of records) {
          const s = segments.get(segment);
          if (!s) continue;
          const elapsed = times ? times[index + 1] - times[index] : -1n;
          // timestamps can be invalid (e.g. after power state changes)
          if (elapsed < 0n) s.invalid = true;
          else s.sum += Number(elapsed);
          s.pending -= 1;
          finish(segment);
        }
      };
      return () => {
        readback.mapAsync(MapMode.READ, 0, size).then(
          () => {
            done(new BigUint64Array(readback.getMappedRange(0, size)));
            readback.unmap();
            readbacks.push(readback);
          },
          () => {
            done(undefined);
            readback.destroy();
          },
        );
      };
    },
    destroy: () => {
      querySet?.destroy();
      resolveBuffer?.destroy();
      for (const b of readbacks) b.destroy();
      readbacks.length = 0;
    },
  };
}
