/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = `ATOM      1  N   ALA A   1      -1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       0.551   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1      -0.200   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       0.540  -0.760  -1.210  1.00  0.00           C
ATOM      6  N   GLY A   2       1.870   1.550   0.000  1.00  0.00           N
ATOM      7  CA  GLY A   2       2.520   2.860   0.000  1.00  0.00           C
ATOM      8  C   GLY A   2       4.030   2.720   0.000  1.00  0.00           C
ATOM      9  O   GLY A   2       4.570   1.610   0.000  1.00  0.00           O
END
`;

// What a recorded pass (render bundle) relies on: the global bind group's dynamic offsets repeat in frames of the same
// passes; the global blocks have a ring of their own, so other blocks and compute kernels writing their ring between
// the passes of a frame don't move them (not provoked here)
describe.skipIf(!WebGPU.module)('global bindings on WebGPU', () => {
  it('bind the same offsets in frames of the same passes, the camera and the data changing', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 64, height: 64 },
      { backend: 'webgpu' },
    );
    // the dynamic offsets of the global bind group (0), per frame
    const encoder = (WebGPU.module!.globals as { GPURenderPassEncoder: { prototype: GPURenderPassEncoder } })
      .GPURenderPassEncoder.prototype;
    const setBindGroup = encoder.setBindGroup;
    let frame: number[] | undefined;
    encoder.setBindGroup = function (this: GPURenderPassEncoder, ...args: Parameters<typeof setBindGroup>) {
      const [index, , offsets] = args;
      if (frame && index === 0 && Array.isArray(offsets) && offsets.length) frame.push(offsets[0]);
      return (setBindGroup as (...a: unknown[]) => void).apply(this, args);
    } as typeof setBindGroup;
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      const { representation } = plugin.builders.structure;
      await representation.addRepresentation(structure, { type: 'ball-and-stick' });
      const cartoon = await representation.addRepresentation(structure, { type: 'cartoon' });
      const canvas = plugin.canvas3d!;
      const image = () => plugin.getImageRaw({ width: 64, height: 64 });
      // pipelines compile in the background: frames drawn before they are ready skip draws
      await image();
      await canvas.graphics.resources.pending.whenReady();
      await image();

      const frames: number[][] = [];
      for (let i = 0; i < 4; ++i) {
        // the camera moves: the global block changes in every frame
        const position = Vec3.add(
          Vec3(),
          canvas.camera.state.target,
          Vec3.create(20 * Math.cos(i), 5, 20 * Math.sin(i)),
        );
        canvas.camera.setState({ position });
        if (i === 2) {
          // data changes between frames (colors written to the items' textures and blocks)
          await plugin
            .build()
            .to(cartoon)
            .update({ colorTheme: { name: 'uniform', params: {} } })
            .commit();
          await canvas.graphics.resources.pending.whenReady();
          await image();
        }
        frame = [];
        await image();
        frames.push(frame);
        frame = undefined;
      }
      expect(frames[0].length).toBeGreaterThan(1);
      for (const f of frames) expect(f).toEqual(frames[0]);
    } finally {
      encoder.setBindGroup = setBindGroup;
      plugin.dispose();
    }
  });
});
