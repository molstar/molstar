/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it, vi } from 'vitest';
import { DrawPass } from '@molstar/graphics/canvas3d/passes/draw';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { WebGPU } from './webgpu-module.js';

const readExample = (name: string) =>
  readFileSync(new URL(`../../../../../data/examples/${name}`, import.meta.url), 'utf8');

// the warm-up scheduled by the scene's commit (a timer in Node) fired after `dispose` and planned a frame on the
// destroyed context (with native GL an uncaught "Invalid GL context", which ends a plain Node script)
describe.skipIf(!WebGPU.module)('timers after dispose on WebGPU', () => {
  it('plan no frame after dispose', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 64, height: 64 },
      { backend: 'webgpu' },
    );
    await plugin.init();
    const data = await plugin.builders.data.rawData({ data: readExample('1cbs_updated.cif'), label: '1cbs' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
    await plugin.getImageRaw();

    // a commit of a changed scene schedules a warm-up (as the render's does with native GL, which does not yield before
    // dispose): pending at dispose, it would plan a frame
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
  }, 60_000);
});
