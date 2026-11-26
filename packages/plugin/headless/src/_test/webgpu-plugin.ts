/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createWebGPUContext, requestWebGPUDevice, type WebGPUContext } from '@molstar/graphics/gl/contexts/webgpu';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import type { HeadlessScreenshotHelperOptions } from '@molstar/plugin-headless/screenshot';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { PluginSpec } from '@molstar/plugin/spec';
import { WebGPU } from './webgpu-module.js';

/**
 * One WebGPU context a worker for the tests' plugins (`createWebGPUPlugin`): the files of a worker share its module
 * registry (`isolate: false`), so its device is requested once and a pipeline compiled once, not again by every file
 * (a device about 0.35 s, the first frames' pipelines most of a test). Made again where a test lost it
 */
let shared: Promise<WebGPUContext> | undefined;

async function getSharedContext() {
  const current = shared && (await shared);
  if (current && !current.isContextLost) return current;
  shared = (async () => {
    const gpu = WebGPU.module!.create([]);
    const device = await requestWebGPUDevice('default', gpu);
    if (!device) throw new Error('no WebGPU adapter available');
    return createWebGPUContext(device, {
      width: 64,
      height: 64,
      // also keeps `gpu` alive with the context (collecting it while its device is in use crashes Dawn's bindings)
      requestDevice: () => requestWebGPUDevice('default', gpu),
      ownsDevice: true,
    });
  })();
  return shared;
}

export interface WebGPUPluginOptions extends Omit<HeadlessScreenshotHelperOptions, 'backend' | 'graphics'> {
  /**
   * a context of its own, not the worker's: a test that loses, replaces or releases its device, that reads what a new
   * context compiles, or whose context takes options (`webgpu`, which give one of its own as well)
   */
  readonly own?: boolean;
}

/**
 * A headless plugin on WebGPU, its canvas drawing with the worker's context (`getSharedContext`) unless `own`: its
 * disposal leaves the context. One plugin on the shared context at a time (one drawing buffer, sized to the plugin's)
 */
export async function createWebGPUPlugin(
  spec: PluginSpec = DefaultPluginSpec(),
  size: { width: number; height: number } = { width: 64, height: 64 },
  options: WebGPUPluginOptions = {},
) {
  const { own, ...rest } = options;
  const graphics = own || rest.webgpu ? undefined : await getSharedContext();
  if (graphics) {
    if (live) {
      throw new Error(
        'a plugin on the shared WebGPU context is not disposed yet: one at a time (its drawing buffer), else `own`',
      );
    }
    // what an earlier plugin left: its work done, its errors not this one's
    await graphics.resources.pending.whenReady();
    graphics.takeErrors();
  }
  const plugin = await HeadlessPluginContext.create({ webgpu: WebGPU.module! }, spec, size, {
    ...rest,
    backend: 'webgpu',
    ...(graphics ? { graphics } : {}),
  });
  if (graphics) {
    live = plugin;
    const dispose = plugin.dispose.bind(plugin);
    const root = { dispose: () => plugin.dispose() };
    plugin.dispose = (...args) => {
      if (live === plugin) live = undefined;
      roots()?.delete(root);
      dispose(...args);
    };
    // a file that leaves it: disposed after the file by the shared setup (`vitest.setup.ts`, its roots), not left to
    // fail the worker's next file
    roots()?.add(root);
  }
  return plugin;
}

/** the setup's registry of what a test file leaves (`vitest.setup.ts`): disposed after the file */
const roots = () =>
  (
    globalThis as {
      [k: symbol]: { add(r: { dispose(): void }): void; delete(r: object): void } | undefined;
    }
  )[Symbol.for('molstar.test.rootScopes')];

/** the plugin on the shared context, till it is disposed (`createWebGPUPlugin`) */
let live: HeadlessPluginContext | undefined;
