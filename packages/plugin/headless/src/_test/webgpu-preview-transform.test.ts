/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { DefaultCanvas3DParams } from '@molstar/graphics/canvas3d/canvas3d';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { SpheresBuilder } from '@molstar/graphics/geo/geometry/spheres/spheres-builder';
import { Box } from '@molstar/graphics/geo/primitive/box';
import * as Shape from '@molstar/graphics/geo/shape/shape';
import { ShapeRepresentation } from '@molstar/graphics/repr/shape/representation';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 128;

/** The pixels further than 10 levels apart (any channel) */
function differing(a: Uint8ClampedArray, b: Uint8ClampedArray) {
  let count = 0;
  for (let i = 0; i < a.length; i += 4) {
    if (Math.max(...[0, 1, 2].map((c) => Math.abs(a[i + c] - b[i + c]))) > 10) count += 1;
  }
  return count;
}

// a transient transform (`previewTransform`, a drag before its release) is applied on the GPU: drawn as the same
// transform as a placement (`transform`, baked into the instances), culled and bounded where it is drawn
describe.skipIf(!WebGPU.module)('a transient transform on WebGPU', () => {
  it('draws a mesh and impostor spheres as the same transform as a placement', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      plugin.renderer.imagePass.setProps({ multiSample: { ...DefaultCanvas3DParams.multiSample, mode: 'off' } });
      const canvas = plugin.canvas3d!;

      const state = MeshBuilder.createState(64, 64);
      state.currentGroup = 0;
      MeshBuilder.addPrimitive(state, Mat4.fromUniformScaling(Mat4(), 4), Box());
      const box = ShapeRepresentation(
        () =>
          Shape.create(
            'box',
            {},
            MeshBuilder.getMesh(state),
            () => Color(0x999999),
            () => 1,
            () => 'box',
          ),
        Mesh.Utils,
      );
      const builder = SpheresBuilder.create(8, 8);
      for (let i = 0; i < 4; ++i) builder.add(i * 3 - 4.5, 4, 0, i);
      const spheres = ShapeRepresentation(
        () =>
          Shape.create(
            'spheres',
            {},
            builder.getSpheres(),
            () => Color(0x377eb8),
            () => 1,
            () => 'spheres',
          ),
        Spheres.Utils,
      );
      await box.createOrUpdate(PD.getDefaultValues(Mesh.Params), {}).run();
      await spheres.createOrUpdate(PD.getDefaultValues(Spheres.Params), {}).run();
      canvas.add(box);
      canvas.add(spheres);
      canvas.commit(true);
      canvas.camera.setState({ position: Vec3.create(0, 2, 30), target: Vec3(), up: Vec3.unitY }, 0);

      const image = async () => {
        canvas.update(box);
        canvas.update(spheres);
        canvas.commit(true);
        const { data } = await plugin.getImageRaw({ width: Size, height: Size });
        return new Uint8ClampedArray(data);
      };
      // turned and moved to the side, still in view
      const moved = Mat4.mul(
        Mat4(),
        Mat4.fromTranslation(Mat4(), Vec3.create(4, -2, 0)),
        Mat4.fromRotation(Mat4(), Math.PI / 5, Vec3.create(1, 1, 0)),
      );

      const still = await image();
      for (const r of [box, spheres]) r.setState({ transform: moved });
      const placed = await image();
      for (const r of [box, spheres]) r.setState({ transform: Mat4.identity(), previewTransform: moved });
      const previewed = await image();
      const transforms = box.renderObjects[0].values.aTransform.ref.version;
      for (const r of [box, spheres]) r.setState({ previewTransform: Mat4.identity() });
      const back = await image();

      // the move shows, as the placement draws it
      expect(differing(placed, still)).toBeGreaterThan(200);
      expect(differing(previewed, placed)).toBeLessThan(8);
      // put back: as before, no instance touched
      expect(differing(back, still)).toBe(0);
      expect(box.renderObjects[0].values.aTransform.ref.version).toBe(transforms);
    } finally {
      plugin.dispose();
    }
  });
});
