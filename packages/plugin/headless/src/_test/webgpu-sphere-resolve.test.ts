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

const Size = 64;
/** What the view holds top to bottom with a `distance` (Å), about the structure's size */
const ViewHeight = 16;

/**
 * The structure as spacefill colored by element index, its spheres drawn by the raster (`on`) or as impostors; with
 * `distance` seen from that far, its size on screen kept by the field of view
 */
async function render(sphereRaster: 'on' | 'off' | 'auto' | 'only', distance?: number) {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: Size, height: Size },
    { backend: 'webgpu' },
  );
  try {
    await plugin.init();
    if (!plugin.canvas3d!.graphics.capabilities.atomicMinMax64) return undefined;
    plugin.canvas3d!.setProps({
      renderer: { sphereRaster, sphereRasterRadius: 32 },
      // the raster draws perspective views
      camera: { mode: 'perspective' } as never,
      cameraFog: { name: 'off', params: {} },
      postprocessing: {
        occlusion: { name: 'off', params: {} },
        outline: { name: 'off', params: {} },
      },
    });
    const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model);
    await plugin.builders.structure.representation.addRepresentation(structure, {
      type: 'spacefill',
      color: 'element-index',
    });
    if (distance) {
      const camera = plugin.canvas3d!.camera;
      const { target } = camera.getSnapshot();
      camera.setState(
        {
          target,
          position: Vec3.create(target[0], target[1], target[2] + distance),
          up: Vec3.create(0, 1, 0),
          fov: 2 * Math.atan(ViewHeight / 2 / distance),
          radius: ViewHeight,
          radiusMax: ViewHeight,
        },
        0,
      );
    }
    const image = () => plugin.getImageRaw({ width: Size, height: Size });
    // pipelines and kernels compile in the background: frames drawn before they are ready draw impostors
    for (let i = 0; i < 3; ++i) {
      await image();
      await plugin.canvas3d!.graphics.resources.pending.whenReady();
    }
    const pixels = Array.from((await image()).data);
    expect(plugin.canvas3d!.graphics.takeErrors?.() ?? []).toEqual([]);
    return pixels;
  } finally {
    plugin.dispose();
  }
}

// the raster's spheres resolved in the opaque pass draw the image the impostors draw
describe.skipIf(!WebGPU.module)('sphere raster resolve on WebGPU', () => {
  it('draws the image of the impostors', async () => {
    const device = (WebGPU.module!.globals as { GPUDevice: { prototype: GPUDevice } }).GPUDevice.prototype;
    const createShaderModule = device.createShaderModule;
    const shaders: string[] = [];
    device.createShaderModule = function (this: GPUDevice, descriptor: GPUShaderModuleDescriptor) {
      shaders.push(descriptor.code);
      return createShaderModule.call(this, descriptor);
    };
    try {
      const impostors = await render('off');
      if (!impostors) return;
      expect(shaders.some((c) => c.includes('fn resolveSphere'))).toBe(false);
      expect(new Set(impostors).size).toBeGreaterThan(2);
      // `auto` leaves objects not culled on the GPU (this one, unmerged) to the impostors
      expect(await render('auto')).toEqual(impostors);
      expect(shaders.some((c) => c.includes('fn resolveSphere'))).toBe(false);
      // the raster's spheres with the larger ones as impostors, and all by the raster (`only`, no impostor drawn)
      for (const raster of [(await render('on'))!, (await render('only'))!]) {
        // the resolve's program and the raster's kernels were made
        expect(shaders.some((c) => c.includes('fn resolveSphere'))).toBe(true);
        expect(shaders.some((c) => c.includes('atomicStoreMax'))).toBe(true);
        let differ = 0;
        for (let i = 0; i < impostors.length; i += 4) {
          const d = Math.max(
            Math.abs(impostors[i] - raster[i]),
            Math.abs(impostors[i + 1] - raster[i + 1]),
            Math.abs(impostors[i + 2] - raster[i + 2]),
          );
          if (d > 8) differ += 1;
        }
        // but where a pixel's center lies on a sphere's edge (the raster's arithmetic is not the fragment stage's)
        expect(differ).toBeLessThanOrEqual(Size * Size * 0.01);
      }
    } finally {
      device.createShaderModule = createShaderModule;
    }
  });

  // the distance squared is rounded by more than the radius squared far from the camera: the hits may not follow it
  it.each(['off', 'only'] as const)('shades alike from near and far (sphereRaster %s)', async (sphereRaster) => {
    const near = await render(sphereRaster, 400);
    if (!near) return;
    const far = (await render(sphereRaster, 30000))!;
    expect(new Set(near).size).toBeGreaterThan(2);
    let differ = 0;
    for (let i = 0; i < near.length; i += 4) {
      const d = Math.max(
        Math.abs(near[i] - far[i]),
        Math.abs(near[i + 1] - far[i + 1]),
        Math.abs(near[i + 2] - far[i + 2]),
      );
      if (d > 16) differ += 1;
    }
    expect(differ).toBeLessThanOrEqual(Size * Size * 0.05);
  });
});
