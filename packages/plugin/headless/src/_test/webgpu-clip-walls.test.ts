/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterAll, describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Clip } from '@molstar/core/util/clip';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { DefaultCanvas3DParams } from '@molstar/graphics/canvas3d/canvas3d';
import { clipObjectFrame, clipWallsMesh } from '@molstar/graphics/geo/clip-volumes';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { Box } from '@molstar/graphics/geo/primitive/box';
import * as Shape from '@molstar/graphics/geo/shape/shape';
import { ShapeRepresentation } from '@molstar/graphics/repr/shape/representation';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { createWebGPUPlugin } from './webgpu-plugin.js';
import { WebGPU } from './webgpu-module.js';

/** Where the walls enter the box the two offsets differ by a line of pixels: up to 40 at 256 px */
const MaxDiffering = (size: number) => (40 * size) / 256;

const Plane: Clip.Props['objects'][number] = {
  type: 'plane',
  invert: false,
  position: Vec3(),
  rotation: { axis: Vec3.create(0, 0, 1), angle: 30 },
  scale: Vec3.create(1, 1, 1),
  transform: Mat4.identity(),
};

function shape(name: string, mesh: Mesh, color: number) {
  return ShapeRepresentation(
    () =>
      Shape.create(
        name,
        {},
        mesh,
        () => Color(color),
        () => 1,
        () => name,
      ),
    Mesh.Utils,
  );
}

/** The pixels further than 10 levels apart (any channel) */
function differing(a: { data: Uint8ClampedArray }, b: { data: Uint8ClampedArray }) {
  let count = 0;
  for (let i = 0; i < a.data.length; i += 4) {
    if (Math.max(...[0, 1, 2].map((c) => Math.abs(a.data[i + c] - b.data[i + c]))) > 10) count += 1;
  }
  return count;
}

// a clip object's walls lie on the solid-interior cap it cuts: by the depth test alone (in the scene) or the overlay's
// occlusion they would flip per pixel
describe.skipIf(!WebGPU.module)('clip walls on a solid-interior cap', () => {
  // one context for both layers, each case's representations removed after it
  let shared: Promise<HeadlessPluginContext> | undefined;
  const sharedPlugin = () => {
    shared ??= createWebGPUPlugin(DefaultPluginSpec(), { width: 256, height: 256 }).then(async (plugin) => {
      await plugin.init();
      return plugin;
    });
    return shared;
  };
  afterAll(async () => {
    (await shared)?.dispose();
  });

  /** the differing pixels of the walls drawn at the cap, against the walls well in front of it, at `size` */
  async function differingAtCap(layer: 'scene' | 'overlay', size: number) {
    const plugin = await sharedPlugin();
    const canvas = plugin.canvas3d!;
    const reprs: ReturnType<typeof shape>[] = [];
    try {
      plugin.renderer.imagePass.setProps({
        multiSample: { ...DefaultCanvas3DParams.multiSample, mode: 'off' },
        overlay: layer === 'overlay',
      });
      const state = MeshBuilder.createState(64, 64);
      state.currentGroup = 0;
      MeshBuilder.addPrimitive(state, Mat4.fromUniformScaling(Mat4(), 6), Box());
      const box = shape('box', MeshBuilder.getMesh(state), 0x999999);
      const walls = shape('walls', clipWallsMesh(Plane, 6), 0x377eb8);
      await box
        .createOrUpdate(
          { ...PD.getDefaultValues(Mesh.Params), solidInterior: true, clip: { variant: 'pixel', objects: [Plane] } },
          {},
        )
        .run();
      const wallProps = {
        ...PD.getDefaultValues(Mesh.Params),
        alpha: 0.4,
        ignoreLight: true,
      };
      await walls.createOrUpdate(wallProps, {}).run();
      walls.setState({
        transform: clipObjectFrame(Plane, { center: Vec3(), radius: 6 } as Parameters<typeof clipObjectFrame>[1]),
        layer,
        ...(layer === 'overlay' ? { occludedAlpha: 0.25 } : {}),
      });
      canvas.add(box);
      canvas.add(walls);
      reprs.push(box, walls);
      canvas.commit(true);
      // from above the cap, obliquely: grazing enough that the two depths round differently
      canvas.camera.setState({ position: Vec3.create(3.5, 5.5, 3), target: Vec3(), up: Vec3.unitY }, 0);
      const image = async (depthOffset: number) => {
        await walls.createOrUpdate({ ...wallProps, depthOffset }, {}).run();
        canvas.update(walls);
        canvas.commit(true);
        const { data, width, height } = await plugin.getImageRaw({ width: size, height: size });
        return { data: new Uint8ClampedArray(data), width, height };
      };
      // well in front of the cap there is nothing to fight over: the walls are even, what the image has to match
      const even = await image(0.5);
      const drawn = await image(layer === 'overlay' ? 0 : 0.02);
      return differing(drawn, even);
    } finally {
      for (const r of reprs) {
        canvas.remove(r);
        r.destroy();
      }
      canvas.commit(true);
    }
  }

  // quick: the scene layer at 64 px (the walls at the cap without their offset: some 150 pixels flip)
  it('draws them evenly (scene, 64 px)', async () => {
    expect(await differingAtCap('scene', 64)).toBeLessThan(MaxDiffering(64));
  });

  // Medium: a 256² image of solid interiors, clip objects and the overlay, its pipelines compiled first (a second)
  for (const layer of ['scene', 'overlay'] as const) {
    it(`draws them evenly (${layer})`, { tags: ['medium'] }, async () => {
      expect(await differingAtCap(layer, 256)).toBeLessThan(MaxDiffering(256));
    });
  }
});
