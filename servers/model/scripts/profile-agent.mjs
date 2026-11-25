/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 *
 * Preloaded into the server by `profile-memory.mjs` (`node --import`): answers memory, GC and heap-snapshot
 * requests over the IPC channel, so the server needs no profiling code of its own.
 */

import * as v8 from 'node:v8';

async function fullGc() {
  if (typeof globalThis.gc !== 'function') return;
  // Two rounds with a pause between them let weak callbacks and finalizers (external buffers) run.
  for (let i = 0; i < 2; i++) {
    globalThis.gc();
    await new Promise((resolve) => setTimeout(resolve, 50));
  }
  globalThis.gc();
}

if (process.send) {
  process.on('message', async (msg) => {
    if (!msg || typeof msg !== 'object' || msg.__profile !== true) return;
    try {
      if (msg.cmd === 'mem') {
        if (msg.gc) await fullGc();
        const heap = v8.getHeapStatistics();
        process.send({
          __profile: true,
          id: msg.id,
          mem: process.memoryUsage(),
          heapLimit: heap.heap_size_limit,
        });
      } else if (msg.cmd === 'snapshot') {
        await fullGc();
        process.send({ __profile: true, id: msg.id, file: v8.writeHeapSnapshot(msg.path) });
      }
    } catch (e) {
      process.send({ __profile: true, id: msg.id, error: '' + e });
    }
  });
  // The IPC channel must not keep the server alive on its own.
  process.channel?.unref();
}
