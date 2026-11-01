/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ColorNames } from '@molstar/core/util/color/names';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { createRenderObject } from '@molstar/graphics/gl/render-object';
import { Scene } from '@molstar/graphics/gl/scene';
import { KernelCode, StorageBufferSpec } from '@molstar/graphics/gl/compute';
import { UniformSpec } from '@molstar/graphics/gl/renderable/schema';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const CopySchema = {
  uOffset: UniformSpec('i'),
  uCount: UniformSpec('i'),
  page: StorageBufferSpec('u32', 'read'),
  out: StorageBufferSpec('u32', 'read-write'),
};

/** copies `uCount` words of `page` from `uOffset` on into `out` */
const CopyKernel = KernelCode('test-copy-data-layer', [64], {
  wgsl: /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) gid: vec3u) {
    if (gid.x >= u32(uParams.uCount)) { return; }
    out[gid.x] = page[u32(uParams.uOffset) + gid.x];
}
`,
});

// kernels read a render item's data layers where its programs do: the item's page, from the layer's word offset
describe.skipIf(!WebGPU.module)('data layers as kernel bindings on WebGPU', () => {
  it('read the positions and sizes the item draws with', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const compute = graphics.compute!;
      const centers = new Float32Array([0, 0, 0, 1.5, -2, 0.25, -3, 4, 5]);
      const spheres = Spheres.create(centers, new Float32Array([0, 1, 2]), 3);
      const props = PD.getDefaultValues(Spheres.Params);
      const simple = BaseGeometry.createSimple(ColorNames.orange as never, 1);
      const values = Spheres.Utils.createValues(
        spheres,
        simple.transform,
        simple.locationIterator,
        simple.theme,
        props,
      );
      const scene = Scene.create(graphics);
      scene.add(createRenderObject('spheres', values, Spheres.Utils.createRenderableState(props), -1));
      scene.commit();
      const item = scene.renderables[0].renderItem;

      const read = async (name: string, words: number) => {
        const layer = item.getDataLayer!(name)!;
        expect(layer).toBeDefined();
        const out = compute.createBuffer(words * 4);
        const kernel = compute.createKernel(CopyKernel, CopySchema, {
          uOffset: ValueCell.create(layer.offset),
          uCount: ValueCell.create(words),
          page: ValueCell.create(layer.buffer),
          out: ValueCell.create(out),
        });
        compute.beginComputePass({ label: 'copy' });
        kernel.dispatch(Math.ceil(words / 64));
        compute.endComputePass();
        const data = await out.readAsync();
        kernel.destroy();
        out.destroy();
        return data;
      };

      // a vec4 per sphere: the center and the group
      const positions = new Float32Array(await read('tPositionGroup', 3 * 4));
      expect(Array.from(positions)).toEqual(Array.from(values.tPositionGroup.ref.value.array.subarray(0, 12)));
      // a word per element: the packed size of the instance (one)
      const sizes = new Uint8Array(await read('tSize', 1));
      expect(Array.from(sizes)).toEqual(Array.from(values.tSize.ref.value.array.subarray(0, 4)));

      expect(item.getDataLayer!('tCylinder')).toBeUndefined();
      // the renderer's: not the kernels' to write
      expect(() => item.getDataLayer!('tPositionGroup')!.buffer.write(new Uint32Array(1))).toThrow();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  });
});
