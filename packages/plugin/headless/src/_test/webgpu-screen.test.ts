/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { ButtonsType, ModifiersKeys } from '@molstar/core/util/input/input-observer';
import { now } from '@molstar/core/util/now';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { DefaultCanvas3DParams } from '@molstar/graphics/canvas3d/canvas3d';
import {
  CameraHelperAxis,
  CameraHelperParams,
  type CameraHelperProps,
  isCameraAxesLoci,
} from '@molstar/graphics/canvas3d/helper/camera-helper';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = 'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n';
const Size = 256;

/**
 * The corner axes on, bottom-left, thick enough to pick at the pick's scale: axes 60 CSS pixels long, 16 wide, the
 * origin inset by the bounding radius (75)
 */
const AxesOn: CameraHelperProps = {
  axes: {
    name: 'on',
    params: {
      ...(PD.getDefaultValues((CameraHelperParams.axes.map('on') as PD.Group<PD.Params>).params) as any),
      location: 'bottom-left',
      scale: 0.6,
      radiusScale: 0.15,
    },
  },
};
const AxesOff: CameraHelperProps = { axes: { name: 'off', params: {} } };
/** The axes' origin, in CSS pixels from the top-left */
const Origin = [75, Size - 75];

/** A plugin with an atom's spacefill, looking down -z with y up, and what the tests ask of its canvas */
async function setup() {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: Size, height: Size },
    { backend: 'webgpu' },
  );
  await plugin.init();
  const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  const spacefill = (await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' }))
    .data!.repr;
  const canvas = plugin.canvas3d!;
  const pixelRatio = canvas.graphics.pixelRatio;
  canvas.setProps({ cameraResetDurationMs: 0, camera: { ...canvas.props.camera, helper: AxesOn } });

  let time = now();
  const frame = async () => {
    canvas.commit(true);
    canvas.requestDraw();
    canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
    await canvas.graphics.resources.pending.whenReady();
  };
  /** the loci and representation picked at (x, y) in CSS pixels, after a drawn frame */
  const pickAt = async (x: number, y: number) => {
    await frame();
    const pick = await canvas.asyncIdentify(Vec2.create(x, y))?.wait();
    return pick ? canvas.getLoci(pick.id) : undefined;
  };
  // the camera reset to the scene first, then turned to look down -z
  await frame();
  canvas.camera.setState({ position: Vec3.create(0, 0, 50), target: Vec3(), up: Vec3.unitY }, 0);

  /** the color at (x, y) in CSS pixels of an image with the axes as `axes` */
  const colorAt = async (axes: CameraHelperProps, x: number, y: number) => {
    plugin.renderer.imagePass.setProps({ cameraHelper: axes });
    const { data, width } = await plugin.getImageRaw({ width: Size, height: Size });
    const i = (Math.round(y * pixelRatio) * width + Math.round(x * pixelRatio)) * 4;
    return [data[i], data[i + 1], data[i + 2]];
  };
  return { plugin, canvas, spacefill, frame, pickAt, colorAt, pixelRatio };
}

// the corner axes are in the screen layer: drawn and picked in their corner above the scene, out of its bounds; a
// click aligns the camera (the plugin's `CameraAxisHelper`); in images as their `cameraHelper` props say
describe.skipIf(!WebGPU.module)('the corner axes in the screen layer on WebGPU', () => {
  it('picks an axis or plane in the corner, and a click aligns the camera to it', async () => {
    const { plugin, canvas, spacefill, frame, pickAt, pixelRatio } = await setup();
    try {
      // the corner fits the axes (the test's sizes are in CSS pixels)
      expect(Size / pixelRatio).toBeGreaterThan(2 * 75 + 60);
      const sceneSphere = Sphere3D.clone(canvas.boundingSphere);

      // the x axis to the right of the origin, the xy plane up and right of it
      const x = await pickAt(Origin[0] + 35, Origin[1]);
      expect(x && isCameraAxesLoci(x.loci) && x.loci.elements[0].groupId).toBe(CameraHelperAxis.X);
      const xy = await pickAt(Origin[0] + 20, Origin[1] - 20);
      expect(xy && isCameraAxesLoci(xy.loci) && xy.loci.elements[0].groupId).toBe(CameraHelperAxis.XY);
      // the scene where the axes are not
      expect((await pickAt(Size / 2, Size / 2))?.repr).toBe(spacefill);
      // out of the scene's bounds
      expect(Sphere3D.equals(canvas.boundingSphere, sceneSphere)).toBe(true);

      const click = (loci: NonNullable<typeof x>) =>
        plugin.behaviors.interaction.click.next({
          current: loci,
          buttons: ButtonsType.Flag.Primary,
          button: ButtonsType.Flag.Primary,
          modifiers: ModifiersKeys.None,
        });
      // the x axis: up along it, the view direction kept
      click(x!);
      await frame();
      expect(Vec3.distance(canvas.camera.state.up, Vec3.unitX)).toBeLessThan(1e-3);
      const dir = Vec3.normalize(Vec3(), Vec3.sub(Vec3(), canvas.camera.state.target, canvas.camera.state.position));
      expect(Vec3.distance(dir, Vec3.negUnitZ)).toBeLessThan(1e-3);
      // the xy plane: looking down -z, y up
      const xy2 = await pickAt(Origin[0] - 20, Origin[1] - 20);
      // turned with the camera (x up now): the xy plane is up and left of the origin
      expect(xy2 && isCameraAxesLoci(xy2.loci) && xy2.loci.elements[0].groupId).toBe(CameraHelperAxis.XY);
      click(xy2!);
      await frame();
      expect(Vec3.distance(canvas.camera.state.up, Vec3.unitY)).toBeLessThan(1e-3);

      // off: the scene's pixels again
      canvas.setProps({ camera: { ...canvas.props.camera, helper: AxesOff } });
      const off = await pickAt(Origin[0] + 35, Origin[1]);
      expect(off === undefined || !isCameraAxesLoci(off.loci)).toBe(true);
    } finally {
      plugin.dispose();
    }
  });

  it('draws the axes into images as their props say', async () => {
    const { plugin, colorAt } = await setup();
    try {
      plugin.renderer.imagePass.setProps({ multiSample: { ...DefaultCanvas3DParams.multiSample, mode: 'off' } });
      const at = [Origin[0] + 35, Origin[1]] as const;
      const without = await colorAt(AxesOff, ...at);
      const withAxes = await colorAt(AxesOn, ...at);
      // the red x axis, translucent (alpha 0.51) over the background
      expect(withAxes[0] - withAxes[1]).toBeGreaterThan(60);
      expect(Math.abs(without[0] - without[1])).toBeLessThan(8);
    } finally {
      plugin.dispose();
    }
  });
});
