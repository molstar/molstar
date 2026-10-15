import assert from 'node:assert/strict';
import { Task } from '@molstar/core/task';
import { encodeVideoAnimation } from '@molstar/video-export-extension/encoder';
import { VideoHeadlessPluginContext } from '@molstar/video-export-extension/headless';
import { H264WasmBackend } from '@molstar/video-export-extension/backends/h264-wasm';
import { VideoExport } from '@molstar/video-export-extension';

const noop = async () => {};
const plugin = {
  animationLoop: { isAnimating: false, stop() {}, resetTime() {}, tick: noop },
  managers: { animation: { play: noop, stop: noop } },
};
const movie = await Task.create('Node MP4 smoke', (ctx) =>
  encodeVideoAnimation(plugin, ctx, {
    animation: { definition: {}, params: {}, customDurationMs: 100 },
    width: 16,
    height: 16,
    viewport: { x: 0, y: 0, width: 16, height: 16 },
    fps: 2,
    backend: H264WasmBackend,
    pass: { updateBackground: noop, getImageData: async () => ({ data: new Uint8Array(16 * 16 * 4).fill(255) }) },
  }),
).run();
assert(movie.data.length > 0);
assert.equal(new TextDecoder().decode(movie.data.subarray(4, 8)), 'ftyp', 'Encoder did not produce an MP4 container');

// Exercise the headless registration guard without requiring a native GL context.
const fakeHeadless = {
  state: { hasBehavior: (behavior) => behavior === VideoExport },
  runTask: async () => movie,
};
assert.equal(await VideoHeadlessPluginContext.prototype.getAnimation.call(fakeHeadless), movie);
await assert.rejects(
  VideoHeadlessPluginContext.prototype.getAnimation.call({
    ...fakeHeadless,
    state: { hasBehavior: () => false },
  }),
  /extension registered/,
);
console.log('Node ESM MP4 encoding and headless registration guard passed');
