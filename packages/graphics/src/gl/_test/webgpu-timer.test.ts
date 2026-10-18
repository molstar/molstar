/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createWebGPUTimerSource } from '../webgpu/timer.js';

/**
 * A device whose timestamp resolves copy `timestamps` (by query index) into the readback buffers; `mapAsync` resolves
 * on the next microtask
 */
function createFakeDevice(timestamps: Map<number, bigint>, features = ['timestamp-query']) {
  const resolved: number[] = [];
  const device = {
    features: new Set(features),
    createQuerySet: () => ({ destroy: () => {} }),
    createBuffer: ({ size }: { size: number }) => {
      const data = new ArrayBuffer(size);
      return {
        data,
        mapAsync: () => Promise.resolve(),
        getMappedRange: (offset: number, length: number) => data.slice(offset, offset + length),
        unmap: () => {},
        destroy: () => {},
      };
    },
  };
  const encoder = {
    resolveQuerySet: (_set: unknown, first: number, count: number) => {
      for (let i = first; i < first + count; ++i) resolved.push(i);
    },
    copyBufferToBuffer: (
      _src: unknown,
      _srcOffset: number,
      dst: { data: ArrayBuffer },
      _dstOffset: number,
      size: number,
    ) => {
      const view = new BigUint64Array(dst.data, 0, size / 8);
      for (let i = 0; i < view.length; ++i) view[i] = timestamps.get(i) ?? 0n;
    },
  };
  return { device: device as unknown as GPUDevice, encoder: encoder as unknown as GPUCommandEncoder, resolved };
}

const flush = () => new Promise((resolve) => setTimeout(resolve, 0));

describe('webgpu timer source', () => {
  it('needs the timestamp-query feature', () => {
    const { device } = createFakeDevice(new Map(), []);
    expect(createWebGPUTimerSource(device)).toBeUndefined();
  });

  it('sums the encoders of a pass, across submissions', async () => {
    // pass 0: encoders at queries 0-1 (10 ns) and, after a submission, 0-1 again (5 ns); pass 1: queries 2-3 (7 ns)
    const timestamps = new Map<number, bigint>([
      [0, 100n],
      [1, 110n],
      [2, 200n],
      [3, 207n],
    ]);
    const { device, encoder } = createFakeDevice(timestamps);
    const source = createWebGPUTimerSource(device)!;

    const a = source.begin();
    expect(source.getTimestampWrites()).toMatchObject({ beginningOfPassWriteIndex: 0, endOfPassWriteIndex: 1 });
    source.beforeSubmit(encoder)!();
    await flush();
    // the pass is still open
    expect([...source.poll()]).toEqual([]);

    timestamps.set(1, 105n);
    expect(source.getTimestampWrites()).toMatchObject({ beginningOfPassWriteIndex: 0 });
    source.end();
    const b = source.begin();
    expect(source.getTimestampWrites()).toMatchObject({ beginningOfPassWriteIndex: 2, endOfPassWriteIndex: 3 });
    source.end();
    source.beforeSubmit(encoder)!();
    await flush();

    expect([...source.poll()]).toEqual([
      [a, 15],
      [b, 7],
    ]);
  });

  it('finishes passes without encoders right away', () => {
    const { device } = createFakeDevice(new Map());
    const source = createWebGPUTimerSource(device)!;
    const a = source.begin();
    source.end();
    expect([...source.poll()]).toEqual([[a, 0]]);
  });

  it('reports invalid timestamps', async () => {
    const { device, encoder } = createFakeDevice(
      new Map([
        [0, 100n],
        [1, 90n],
      ]),
    );
    const source = createWebGPUTimerSource(device)!;
    const a = source.begin();
    source.getTimestampWrites();
    source.end();
    source.beforeSubmit(encoder)!();
    await flush();
    expect([...source.poll()]).toEqual([[a, undefined]]);
  });

  it('times nothing outside of passes', () => {
    const { device, encoder } = createFakeDevice(new Map());
    const source = createWebGPUTimerSource(device)!;
    expect(source.getTimestampWrites()).toBeUndefined();
    expect(source.beforeSubmit(encoder)).toBeUndefined();
  });
});
