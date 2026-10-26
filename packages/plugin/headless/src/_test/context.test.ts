/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { describe, expect, it, vi } from 'vitest';
import { DrawPass } from '@molstar/graphics/canvas3d/passes/draw';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import type { ExternalModules } from '@molstar/plugin-headless/screenshot';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PluginStateObject } from '@molstar/plugin/state/objects';

const require = createRequire(import.meta.url);
let hasGl = true;
let glModule: typeof import('gl') | undefined;
try {
  glModule = require('gl') as typeof import('gl');
  const context = glModule(1, 1);
  if (!context) throw new Error('Optional package "gl" loaded but could not create a WebGL context.');
  context.getExtension('STACKGL_destroy_context')?.destroy();
} catch (e) {
  const error = e as NodeJS.ErrnoException;
  const message = error.message ?? String(e);
  const packageMissing = error.code === 'MODULE_NOT_FOUND' && /Cannot find module ['"]gl['"]/.test(message);
  const bindingMissing =
    /Could not locate the bindings file/i.test(message) ||
    (error.code === 'ERR_DLOPEN_FAILED' &&
      /NODE_MODULE_VERSION|compiled against|different Node\.js version/i.test(message)) ||
    (error.code === 'MODULE_NOT_FOUND' &&
      /\.node['"]/.test(message) &&
      (error.requireStack ?? []).some((path) => /[/\\]gl[/\\]/.test(path)));
  if (!packageMissing && !bindingMissing) throw e;
  hasGl = false;
  console.warn(`Skipping HeadlessPluginContext test: optional native package "gl" is unavailable (${message}).`);
}

describe('HeadlessPluginContext', () => {
  (hasGl ? it : it.skip)('HeadlessPluginContext', { tags: ['medium'] }, async () => {
    const externalModules: ExternalModules = { gl: glModule! };
    const spec = DefaultPluginSpec();
    const plugin = new HeadlessPluginContext(externalModules, spec, { height: 100, width: 100 }, {});
    expect(plugin).toBeTruthy();

    await plugin.init();

    plugin.dispose();
  });

  // the warm-up scheduled by a commit of a changed scene (a timer in Node) fired after `dispose` and planned a frame on
  // the destroyed context: an uncaught "Invalid GL context", which ends a plain Node script
  (hasGl ? it : it.skip)('plans no frame after dispose', { tags: ['medium'] }, async () => {
    const externalModules: ExternalModules = { gl: glModule! };
    const plugin = new HeadlessPluginContext(externalModules, DefaultPluginSpec(), { height: 64, width: 64 }, {});
    await plugin.init();
    const cif = readFileSync(new URL('../../../../../data/examples/1cbs_updated.cif', import.meta.url), 'utf8');
    const data = await plugin.builders.data.rawData({ data: cif, label: '1cbs' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
    await plugin.getImageRaw();

    const canvas = plugin.canvas3d!;
    const [cell] = plugin.state.data.selectQ((q) => q.ofType(PluginStateObject.Molecule.Structure.Representation3D));
    const repr = cell.obj!.data.repr;
    canvas.remove(repr);
    canvas.add(repr);
    canvas.commit(true);
    const plan = vi.spyOn(DrawPass.prototype, 'plan');
    plugin.dispose();
    // past the warm-up (a timer of 0 ms in Node)
    await new Promise((resolve) => setTimeout(resolve, 50));
    expect(plan).not.toHaveBeenCalled();
  });
});
