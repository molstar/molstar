import assert from 'node:assert/strict';
import { Task } from '@molstar/core/task';
import { encodeMp4Animation } from '@molstar/mp4-export-extension/encoder';
import { Mp4HeadlessPluginContext } from '@molstar/mp4-export-extension/headless';
import { Mp4Export } from '@molstar/mp4-export-extension';

const noop = async () => {};
const plugin = {
  animationLoop: { isAnimating: false, stop() {}, resetTime() {}, tick: noop },
  managers: { animation: { play: noop, stop: noop } },
};
const movie = await Task.create('Node MP4 smoke', ctx => encodeMp4Animation(plugin, ctx, {
  animation: { definition: {}, params: {}, customDurationMs: 100 },
  width: 16, height: 16,
  viewport: { x: 0, y: 0, width: 16, height: 16 },
  fps: 2,
  pass: { updateBackground: noop, getImageData: async () => ({ data: new Uint8Array(16 * 16 * 4).fill(255) }) },
})).run();
assert(movie.length > 0);
assert.equal(new TextDecoder().decode(movie.subarray(4, 8)), 'ftyp', 'Encoder did not produce an MP4 container');

// Exercise the headless registration guard without requiring a native GL context.
const fakeHeadless = {
  state: { hasBehavior: behavior => behavior === Mp4Export },
  runTask: async () => movie,
};
assert.equal(await Mp4HeadlessPluginContext.prototype.getAnimation.call(fakeHeadless), movie);
await assert.rejects(Mp4HeadlessPluginContext.prototype.getAnimation.call({
  ...fakeHeadless, state: { hasBehavior: () => false },
}), /extension registered/);
console.log('Node ESM MP4 encoding and headless registration guard passed');
