/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ComputeContext, StorageBuffer } from '../compute.js';

/** How long a scratch buffer is kept after its last use */
const IdleMs = 3000;

type Scratch = { buffer: StorageBuffer; timer: ReturnType<typeof setTimeout> | undefined };

const scratches = new WeakMap<ComputeContext, Map<string, Scratch>>();

/**
 * A storage buffer of at least `byteLength` bytes for the intermediate results of one computation, kept per context and
 * `key` between calls and grown when too small, released when unused for a few seconds. Creating large buffers costs
 * more than the computations using them (WebGPU: resource creation in the GPU process). Its contents are left over from
 * the last use: clear what the computation reads before writing. Use it only in commands recorded without awaiting in
between: another computation may get the same buffer meanwhile.
 */
export function getScratchBuffer(compute: ComputeContext, key: string, byteLength: number): StorageBuffer {
  let perContext = scratches.get(compute);
  if (!perContext) {
    perContext = new Map();
    scratches.set(compute, perContext);
  }
  let scratch = perContext.get(key);
  if (!scratch || scratch.buffer.byteLength < byteLength) {
    if (scratch) {
      clearTimeout(scratch.timer);
      scratch.buffer.destroy();
    }
    scratch = { buffer: compute.createBuffer(byteLength, key), timer: undefined };
    perContext.set(key, scratch);
  }
  const current = scratch;
  clearTimeout(current.timer);
  current.timer = setTimeout(() => {
    // destroying waits for the submission of commands recorded with it
    if (perContext.get(key) === current) perContext.delete(key);
    current.buffer.destroy();
  }, IdleMs);
  // not keeping Node (headless rendering, tests) alive: its timers are objects, browsers' numbers
  if (typeof current.timer === 'object') current.timer.unref();
  return current.buffer;
}
