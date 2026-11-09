/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { now } from '@molstar/core/util/now';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

// Passes binding their inputs in turn every frame (Hi-Z's levels, bloom's blur) find the bind groups they made for
// those inputs before: a moving view's frames make none (36 a frame before, about 0.35 ms of CPU on Dawn).

const Pdb = `${Array.from({ length: 100 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const residue = String(i + 1).padStart(4);
  const [x, y, z] = [(i % 10) * 2.5, Math.floor(i / 10) * 2.5, (i % 7) - 3];
  return `ATOM  ${serial}  CA  ALA A${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

describe.skipIf(!WebGPU.module)('bind groups of passes on WebGPU', () => {
  it('are made once for the inputs a pass binds in turn: none in a moving view’s frames (Hi-Z, bloom)', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 128, height: 128 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      // emissive: bloom draws
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
        typeParams: { emissive: 0.5 } as never,
      });
      const canvas = plugin.canvas3d!;
      canvas.setProps({
        hiZ: { ...canvas.props.hiZ, enabled: true },
        multiSample: { ...canvas.props.multiSample, mode: 'off' },
      });
      const device = (canvas.graphics as WebGPUContext).device;
      let created = 0;
      const create = device.createBindGroup.bind(device);
      (device as { createBindGroup: GPUDevice['createBindGroup'] }).createBindGroup = (d) => {
        created += 1;
        return create(d);
      };
      let time = now();
      const snapshot = canvas.camera.getSnapshot();
      const radius = Vec3.distance(snapshot.position, snapshot.target);
      const frame = (i: number) => {
        const a = i * 0.05;
        const t = snapshot.target;
        canvas.camera.setState(
          { position: Vec3.create(t[0] + radius * Math.sin(a), t[1], t[2] + radius * Math.cos(a)) },
          0,
        );
        canvas.requestDraw();
        canvas.tick((time += 16) as now.Timestamp, { isSynchronous: true });
      };
      for (let i = 0; i < 10; ++i) frame(i);
      await canvas.graphics.resources.pending.whenReady();
      for (let i = 10; i < 20; ++i) frame(i);
      created = 0;
      for (let i = 20; i < 30; ++i) frame(i);
      expect(created).toBe(0);
    } finally {
      plugin.dispose();
    }
  });
});
