/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RuntimeContext } from '../execution/runtime-context.js';

type UniformlyChunkedFn<S> = (chunkSize: number, state: S) => number;

/**
 * After an update, a chunk is a quarter of what was processed since the previous one, so that a chunk takes about a
 * quarter of the update interval at the current rate and a gap between yields stays within ~1.25× of the interval (a
 * shift in the rate costs one longer gap). Sizing a chunk to the whole interval overshoots it by up to one chunk, and
 * the summed times since the last update and the rate averaged over the whole task made the gaps reach several times
 * the interval (perf-audit-2 03/B8).
 */
const ChunksPerUpdate = 4;

async function chunkedSubtask<S>(
  ctx: RuntimeContext,
  initialChunk: number,
  state: S,
  f: UniformlyChunkedFn<S>,
  update: (ctx: RuntimeContext, state: S, processed: number) => Promise<void> | void,
): Promise<S> {
  let chunkSize = Math.max(initialChunk, 0);
  let globalProcessed = 0;

  if (ctx.isSynchronous) {
    f(Number.MAX_SAFE_INTEGER, state);
    return state;
  }

  let lastSize = 0,
    processedSinceUpdate = 0;

  while ((lastSize = f(chunkSize, state)) > 0) {
    globalProcessed += lastSize;
    processedSinceUpdate += lastSize;

    if (ctx.shouldUpdate) {
      await update(ctx, state, globalProcessed);

      // `update` must clear `shouldUpdate` (`ctx.update` does): with it left set, every chunk would be a quarter of
      // the one before
      chunkSize = Math.max(1, Math.round(processedSinceUpdate / ChunksPerUpdate));
      processedSinceUpdate = 0;
    }
  }
  if (ctx.shouldUpdate) {
    await update(ctx, state, globalProcessed);
  }
  return state;
}

export { chunkedSubtask };
