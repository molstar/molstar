/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { Color } from '@molstar/core/util/color';
import { now } from '@molstar/core/util/now';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { DefaultCanvas3DParams } from '@molstar/graphics/canvas3d/canvas3d';
import { DofParams } from '@molstar/graphics/canvas3d/passes/dof';
import type { PostprocessingProps } from '@molstar/graphics/canvas3d/passes/postprocessing';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { STYLIZED_POSTPROCESSING } from '@molstar/plugin-headless/screenshot';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

/** one atom: its spacefill hides its ball-and-stick */
const Pdb = 'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n';
const Size = 64;

function isRed([r, g, b]: number[]) {
  return r > 200 && g < 80 && b < 80;
}

/**
 * A plugin with the atom's spacefill and a red ball-and-stick inside it (`ballParams` added to its type params), on a
 * canvas of `size`, and what the tests ask of its canvas
 */
async function setup(size = Size, ballParams: { [k: string]: unknown } = {}) {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: size, height: size },
    { backend: 'webgpu' },
  );
  await plugin.init();
  const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  const { representation } = plugin.builders.structure;
  const spacefill = (await representation.addRepresentation(structure, { type: 'spacefill' })).data!.repr;
  const ball = (
    await representation.addRepresentation(structure, {
      type: 'ball-and-stick',
      // inside the spacefill (size factor 1), several pixels wide
      typeParams: { sizeFactor: 0.8, ...ballParams },
      color: 'uniform',
      colorParams: { value: Color(0xff0000) },
    })
  ).data!.repr;
  const canvas = plugin.canvas3d!;
  const pixelRatio = canvas.graphics.pixelRatio;

  let time = now();
  /** the last pick's position */
  const pickPosition = Vec3();
  /** the representation picked at (x, y) in pixels of the canvas, after a drawn frame (the pick renders again after one) */
  const pickAt = async (x: number, y: number) => {
    canvas.commit(true);
    canvas.requestDraw();
    canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
    await canvas.graphics.resources.pending.whenReady();
    const pick = await canvas.asyncIdentify(Vec2.create(x / pixelRatio, y / pixelRatio))?.wait();
    Vec3.copy(pickPosition, pick ? pick.position : Vec3.create(NaN, NaN, NaN));
    return pick && canvas.getLoci(pick.id).repr;
  };
  const pickCenter = () => pickAt(size / 2, size / 2);
  // the pick's pipelines may compile after the frame's
  for (let i = 0; i < 20 && !(await pickCenter()); i++);
  /** an image with the overlay or without it */
  const image = async (overlay: boolean, postprocessing?: Partial<PostprocessingProps>) => {
    plugin.renderer.imagePass.setProps({ overlay });
    return plugin.getImageRaw({ width: size, height: size }, postprocessing);
  };
  /** the color at the center of an image */
  const centerColor = async (overlay: boolean, postprocessing?: Partial<PostprocessingProps>) => {
    const { data } = await image(overlay, postprocessing);
    const i = ((size / 2) * size + size / 2) * 4;
    return [data[i], data[i + 1], data[i + 2]];
  };
  return { plugin, canvas, spacefill, ball, pickAt, pickCenter, pickPosition, image, centerColor };
}

/** The pixels of the center row of an image that are red */
function redWidth({ data, width, height }: { data: Uint8ClampedArray; width: number; height: number }) {
  let count = 0;
  const row = Math.floor(height / 2) * width * 4;
  for (let x = 0; x < width; ++x) {
    const i = row + x * 4;
    if (isRed([data[i], data[i + 1], data[i + 2]])) count += 1;
  }
  return count;
}

/** The pixels of the center row of an image that are part red: an edge's width, a few pixels sharp, more blurred */
function softWidth({ data, width, height }: { data: Uint8ClampedArray; width: number; height: number }) {
  let count = 0;
  const row = Math.floor(height / 2) * width * 4;
  for (let x = 0; x < width; ++x) {
    const i = row + x * 4;
    const redness = data[i] - (data[i + 1] + data[i + 2]) / 2;
    if (redness > 30 && !isRed([data[i], data[i + 1], data[i + 2]])) count += 1;
  }
  return count;
}

// a representation in the overlay layer is drawn and picked above the scene, is out of the scene's bounds, and is in
// images only when they ask for it
describe.skipIf(!WebGPU.module)('the overlay layer on WebGPU', () => {
  it('draws and picks above the scene, out of its bounds', async () => {
    const { plugin, canvas, spacefill, ball, pickCenter, pickPosition, centerColor } = await setup();
    try {
      // in the scene, inside the spacefill: hidden, not picked
      expect(await pickCenter()).toBe(spacefill);
      expect(isRed(await centerColor(true))).toBe(false);
      const sceneSphere = Sphere3D.clone(canvas.boundingSphere);

      ball.setState({ layer: 'overlay' });
      canvas.update(ball);
      expect(await pickCenter()).toBe(ball);

      // behind the scene, beyond the camera's far plane: drawn and picked with the overlay's own planes
      const { position, target } = canvas.camera.state;
      const behind = Vec3.scale(Vec3(), Vec3.normalize(Vec3(), Vec3.sub(Vec3(), target, position)), 30);
      ball.setState({ transform: Mat4.fromTranslation(Mat4(), behind) });
      canvas.update(ball);
      expect(await pickCenter()).toBe(ball);
      expect(Vec3.distance(behind, target)).toBeGreaterThan(canvas.camera.far - Vec3.distance(position, target));
      // on the ball's surface, its radius 0.8 * 1.7
      expect(Vec3.distance(pickPosition, behind)).toBeGreaterThan(1.1);
      expect(Vec3.distance(pickPosition, behind)).toBeLessThan(1.45);
      // the scene's bounds stay
      expect(Sphere3D.equals(canvas.boundingSphere, sceneSphere)).toBe(true);

      // on top in images that ask for it: multi-sampled (the default), single sample, post-processed
      expect(isRed(await centerColor(true))).toBe(true);
      plugin.renderer.imagePass.setProps({ multiSample: { ...DefaultCanvas3DParams.multiSample, mode: 'off' } });
      expect(isRed(await centerColor(true))).toBe(true);
      expect(isRed(await centerColor(true, STYLIZED_POSTPROCESSING))).toBe(true);
      // screenshots leave it out by default
      expect(isRed(await centerColor(false))).toBe(false);

      // back in the scene, its bounds grow
      ball.setState({ layer: 'scene' });
      canvas.update(ball);
      canvas.commit(true);
      expect(canvas.boundingSphere.radius).toBeGreaterThan(sceneSphere.radius + 10);
    } finally {
      plugin.dispose();
    }
  });

  // the scene's depth applies per object by its occluded alpha: on top (1), depth tested (0), x-ray (between); the
  // pick follows the alpha a click needs
  it('applies the scene depth by occluded alpha', async () => {
    const { plugin, canvas, spacefill, ball, pickCenter, centerColor } = await setup();
    try {
      plugin.renderer.imagePass.setProps({ multiSample: { ...DefaultCanvas3DParams.multiSample, mode: 'off' } });
      const sceneColor = await centerColor(true);
      ball.setState({ layer: 'overlay' });
      canvas.update(ball);
      const occlude = async (occludedAlpha: number) => {
        await ball.createOrUpdate({ occludedAlpha }).run();
        return { picked: await pickCenter(), color: await centerColor(true) };
      };

      let { picked, color } = await occlude(1);
      expect(picked).toBe(ball);
      expect(isRed(color)).toBe(true);

      // depth tested: the scene's
      ({ picked, color } = await occlude(0));
      expect(picked).toBe(spacefill);
      for (let i = 0; i < 3; ++i) expect(Math.abs(color[i] - sceneColor[i])).toBeLessThan(8);

      // x-ray: reddened, below the alpha a click needs
      ({ picked, color } = await occlude(0.25));
      expect(picked).toBe(spacefill);
      expect(color[0]).toBeGreaterThan(sceneColor[0] + 20);
      expect(isRed(color)).toBe(false);

      // in front of the scene: not occluded
      const { position, target } = canvas.camera.state;
      const front = Vec3.scale(Vec3(), Vec3.normalize(Vec3(), Vec3.sub(Vec3(), position, target)), 3);
      ball.setState({ transform: Mat4.fromTranslation(Mat4(), front) });
      canvas.update(ball);
      ({ picked, color } = await occlude(0));
      expect(picked).toBe(ball);
      expect(isRed(color)).toBe(true);
    } finally {
      plugin.dispose();
    }
  });

  // a transparent surface in front occludes as a click takes it: opaque enough to be picked, not fainter than that
  it('applies the scene depth of transparent objects opaque enough to be picked', async () => {
    const { plugin, canvas, spacefill, ball, pickCenter, centerColor } = await setup();
    try {
      plugin.renderer.imagePass.setProps({ multiSample: { ...DefaultCanvas3DParams.multiSample, mode: 'off' } });
      ball.setState({ layer: 'overlay' });
      canvas.update(ball);
      await ball.createOrUpdate({ occludedAlpha: 0 }).run();
      const seeThrough = async (alpha: number) => {
        await spacefill.createOrUpdate({ alpha }).run();
        canvas.update(spacefill);
        // a frame first: one held for what the update needs leaves the pick of the frame before
        await pickCenter();
        return { picked: await pickCenter(), color: await centerColor(true) };
      };

      // faint: not an occluder, so the ball is drawn and picked through it
      let { picked, color } = await seeThrough(0.2);
      expect(picked).toBe(ball);
      expect(isRed(color)).toBe(true);

      // picked itself: the ball behind it hidden, in the image and the pick
      ({ picked, color } = await seeThrough(0.8));
      expect(picked).toBe(spacefill);
      expect(isRed(color)).toBe(false);
    } finally {
      plugin.dispose();
    }
  });

  // a pick-only object (a hit area) is picked, with its position, but not drawn
  it('picks a pick-only object without drawing it', async () => {
    const { plugin, canvas, ball, pickCenter, pickPosition, centerColor } = await setup();
    try {
      plugin.renderer.imagePass.setProps({ multiSample: { ...DefaultCanvas3DParams.multiSample, mode: 'off' } });
      ball.setState({ layer: 'overlay' });
      canvas.update(ball);
      expect(isRed(await centerColor(true))).toBe(true);
      const drawnAt = Vec3.clone(pickPosition);
      expect(await pickCenter()).toBe(ball);
      Vec3.copy(drawnAt, pickPosition);

      // on its render objects: a structure representation does not take `State.pickOnly` (shapes do)
      for (const o of ball.renderObjects) o.state.pickOnly = true;
      canvas.requestDraw();
      expect(await pickCenter()).toBe(ball);
      // on its surface: the pick's depth has it too
      expect(Vec3.distance(pickPosition, drawnAt)).toBeLessThan(0.05);
      expect(isRed(await centerColor(true))).toBe(false);
      // the canvas does not draw it either
      expect(isRed(await centerColor(false))).toBe(false);
    } finally {
      plugin.dispose();
    }
  });

  // with depth of field the scene is blurred and the overlay is drawn after it: its edges stay sharp
  it('keeps the overlay out of the depth of field blur', async () => {
    const { plugin, canvas, ball, image } = await setup(128);
    try {
      plugin.renderer.imagePass.setProps({ multiSample: { ...DefaultCanvas3DParams.multiSample, mode: 'off' } });
      // in front of the atom, out of the focus at its center
      const { position, target } = canvas.camera.state;
      const front = Vec3.scale(Vec3(), Vec3.normalize(Vec3(), Vec3.sub(Vec3(), position, target)), 4);
      ball.setState({ transform: Mat4.fromTranslation(Mat4(), front) });
      canvas.update(ball);
      const dof = { dof: { name: 'on', params: { ...PD.getDefaultValues(DofParams), blurSize: 16, PPM: 1 } } } as const;
      const soft = async () => softWidth(await image(true, dof as Partial<PostprocessingProps>));

      const blurred = await soft();
      ball.setState({ layer: 'overlay' });
      canvas.update(ball);
      const sharp = await soft();
      // a few pixels of shaded and antialiased rim either side (7 measured), against a blur many pixels wide (19)
      expect(sharp).toBeLessThanOrEqual(9);
      expect(blurred).toBeGreaterThanOrEqual(sharp + 8);
    } finally {
      plugin.dispose();
    }
  });

  // a screen-sized mesh is as large on screen at any zoom: in images, and in the pick at its own scale (a quarter of
  // the canvas's pixels by default)
  it('keeps a screen-sized mesh at its size on screen', async () => {
    const size = 128;
    // 15 CSS pixels per unit
    const { plugin, canvas, ball, pickAt, image } = await setup(size, {
      tryUseImpostor: false,
      quality: 'custom',
      detail: 3,
      screenSize: 15,
    });
    try {
      plugin.renderer.imagePass.setProps({ multiSample: { ...DefaultCanvas3DParams.multiSample, mode: 'off' } });
      canvas.setProps({ camera: { ...canvas.props.camera, mode: 'perspective' } });
      ball.setState({ layer: 'overlay' });
      canvas.update(ball);
      const pixelRatio = canvas.graphics.pixelRatio;
      // the camera reset to the scene first
      await pickAt(0, 0);
      const position = Vec3.clone(canvas.camera.state.position);
      const target = Vec3.clone(canvas.camera.state.target);

      const widths: number[] = [];
      for (const factor of [1, 3]) {
        const p = Vec3.scaleAndAdd(Vec3(), target, Vec3.sub(Vec3(), position, target), factor);
        canvas.camera.setState({ position: p }, 0);
        // an image first: frames wait while the props (the camera mode) are prepared
        const width = redWidth(await image(true));
        widths.push(width);
        // inside the ball and beside it (beyond the pick padding of a few pixels), in CSS pixels
        expect(await pickAt(size / 2 + 0.3 * width * pixelRatio, size / 2)).toBe(ball);
        expect(await pickAt(size / 2 + 0.9 * width * pixelRatio, size / 2)).not.toBe(ball);
        expect(Vec3.distance(canvas.camera.state.position, target)).toBeCloseTo(
          factor * Vec3.distance(position, target),
        );
      }
      // 15 pixels per unit: the radius of the ball's mesh (about 1.1 units) is 16 pixels
      expect(widths[0]).toBeGreaterThan(25 * pixelRatio);
      expect(Math.abs(widths[1] - widths[0])).toBeLessThanOrEqual(2);
    } finally {
      plugin.dispose();
    }
  });
});
