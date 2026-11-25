/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterAll, describe, expect, it } from 'vitest';
import { Quat, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { OrderedSet } from '@molstar/core/data/int';
import { ShapeGroup } from '@molstar/model/model/shape/shape';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { GaussianSplats } from '@molstar/graphics/geo/geometry/gaussian-splats/gaussian-splats';
import { encodeSplatSh, setSplat, type SplatSh } from '@molstar/graphics/geo/geometry/gaussian-splats/splat-data';
import { int32ToFloatBits } from '@molstar/graphics/geo/geometry/center-group';
import * as Shape from '@molstar/graphics/geo/shape/shape';
import { ShapeRepresentation } from '@molstar/graphics/repr/shape/representation';
import type { PostprocessingProps } from '@molstar/graphics/canvas3d/passes/postprocessing';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { SplatSortKernels } from '@molstar/graphics/gl/compute/splat-sort';
import { WebGPU } from './webgpu-module.js';

const Size = 128;
/** The camera's distance from the splat: at the screenshot's 90° field of view about 2.1 px a unit */
const Distance = 30;

const NoPostprocessing: Partial<PostprocessingProps> = {
  occlusion: { name: 'off', params: {} },
  outline: { name: 'off', params: {} },
  antialiasing: { name: 'off', params: {} },
  shadow: { name: 'off', params: {} },
};

interface Splat {
  center?: Vec3;
  scale: Vec3;
  q?: Quat;
  color?: number;
  opacity?: number;
}

/** Splats as a representation, a group each (white, opacity 1 and at the origin by default) */
function splatRepr(list: Splat[], sh?: SplatSh) {
  const { positionGroup, splatData } = GaussianSplats.createArrays(list.length);
  list.forEach(({ center = Vec3(), scale, q = Quat.identity(), opacity = 1 }, i) => {
    positionGroup.set([center[0], center[1], center[2], int32ToFloatBits(i)], i * 4);
    setSplat(splatData, i, Math.log(scale[0]), Math.log(scale[1]), Math.log(scale[2]), q[0], q[1], q[2], q[3], opacity);
  });
  const splats = GaussianSplats.create(positionGroup, splatData, list.length, undefined, sh);
  return ShapeRepresentation(
    () =>
      Shape.create(
        'splats',
        {},
        splats,
        (g) => Color(list[g].color ?? 0xffffff),
        () => 1,
        () => 'splat',
      ),
    GaussianSplats.Utils,
  );
}

/** The center pixel's red, green and blue, 0 to 1 */
function centerColor(image: { data: Uint8ClampedArray | Uint8Array }) {
  const o = ((Size / 2) * Size + Size / 2) * 4;
  return [image.data[o] / 255, image.data[o + 1] / 255, image.data[o + 2] / 255];
}

/** The red channel of the image along its middle row (`axis` 'x') or column ('y'), a value per pixel, 0 to 1 */
function profile(image: { data: Uint8ClampedArray | Uint8Array }, axis: 'x' | 'y') {
  const values: number[] = [];
  const c = Size / 2;
  for (let i = 0; i < Size; ++i) {
    const [x, y] = axis === 'x' ? [i, c] : [c, i];
    values.push(image.data[(y * Size + x) * 4] / 255);
  }
  return values;
}

/**
 * `values` scaled to the falloff at their middle pixel (`expected`): the screenshot's jittered samples average the
 * peak down a little (4%), the shape is the splat's
 */
function normalized(values: number[], sigma: number) {
  const peak = values[Size / 2];
  expect(peak).toBeGreaterThan(0.85);
  const scale = expected(Size / 2, sigma) / peak;
  return values.map((v) => v * scale);
}

/** The Gaussian's falloff at pixel `i` of a profile: its center between the two middle pixels */
function expected(i: number, sigma: number) {
  const d = i + 0.5 - Size / 2;
  return Math.exp((-0.5 * d * d) / sigma / sigma);
}

describe.skipIf(!WebGPU.module)('gaussian splats on WebGPU', () => {
  // one context for the tests that draw a representation of their own, removed after
  let shared: Promise<HeadlessPluginContext> | undefined;
  const sharedPlugin = () => {
    shared ??= HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    ).then(async (plugin) => {
      await plugin.init();
      return plugin;
    });
    return shared;
  };
  afterAll(async () => {
    (await shared)?.dispose();
  });

  async function render(list: Splat[], eye = Vec3.create(0, 0, Distance), sh?: SplatSh) {
    const plugin = await sharedPlugin();
    const canvas = plugin.canvas3d!;
    canvas.setProps({ renderer: { backgroundColor: Color(0x000000) }, cameraFog: { name: 'off', params: {} } });
    const repr = splatRepr(list, sh);
    await repr.createOrUpdate(PD.getDefaultValues(GaussianSplats.Params), {}).run();
    canvas.add(repr);
    canvas.commit(true);
    try {
      canvas.camera.setState({ position: eye, target: Vec3(), up: Vec3.unitY }, 0);
      // the pixels per unit at the splat's distance
      const focal = Size / 2 / Math.tan(canvas.camera.state.fov / 2);
      const pxPerUnit = focal / Distance;

      // its pipelines compile in the background: rendered again once they are ready
      const { pending } = canvas.graphics.resources;
      // the order sorted on the GPU once its kernels compiled (the CPU's before): drawn until it is, then again
      const values = repr.renderObjects[0].values as { dSplatGpuOrder: { ref: { value: boolean } } };
      let image = await plugin.getImageRaw({ width: Size, height: Size }, NoPostprocessing);
      const sorted = () => list.length < 2 || values.dSplatGpuOrder.ref.value;
      for (let i = 0; i < 20 && (Math.max(...centerColor(image)) === 0 || !sorted()); i++) {
        await pending.whenReady();
        image = await plugin.getImageRaw({ width: Size, height: Size }, NoPostprocessing);
      }
      if (list.length > 1) expect(values.dSplatGpuOrder.ref.value).toBe(true);
      image = await plugin.getImageRaw({ width: Size, height: Size }, NoPostprocessing);
      expect(canvas.graphics.takeErrors?.() ?? []).toEqual([]);
      return { image, pxPerUnit };
    } finally {
      canvas.remove(repr);
      canvas.commit(true);
      repr.destroy();
    }
  }

  // the 2D Gaussian projected from the 3D one (EWA), white over black: a pixel's value is its coverage
  it('draws a splat as the projected Gaussian, to three standard deviations', async () => {
    const { image, pxPerUnit } = await render([{ scale: Vec3.create(2, 2, 2) }]);
    // the 0.3 px² dilation against aliasing widens it a little
    const sigma = Math.sqrt((2 * pxPerUnit) ** 2 + 0.3);
    expect(sigma).toBeGreaterThan(4);
    for (const axis of ['x', 'y'] as const) {
      const values = normalized(profile(image, axis), sigma);
      for (let i = 0; i < Size; ++i) {
        const d = Math.abs(i + 0.5 - Size / 2);
        // within the quad (three standard deviations); beyond nothing drawn
        if (d < 2.8 * sigma) expect(values[i]).toBeCloseTo(expected(i, sigma), 1);
        else if (d > 3.2 * sigma) expect(values[i]).toBe(0);
      }
    }
  }, 60_000);

  it('turns and stretches a splat by its rotation and scales', async () => {
    // twice as long along x, turned 90° about z: long along y on screen
    const q = Quat.setAxisAngle(Quat(), Vec3.create(0, 0, 1), Math.PI / 2);
    const { image, pxPerUnit } = await render([{ scale: Vec3.create(3, 1.5, 1.5), q }]);
    const long = Math.sqrt((3 * pxPerUnit) ** 2 + 0.3);
    const short = Math.sqrt((1.5 * pxPerUnit) ** 2 + 0.3);
    const x = normalized(profile(image, 'x'), short);
    const y = normalized(profile(image, 'y'), long);
    for (let i = 0; i < Size; ++i) {
      const d = Math.abs(i + 0.5 - Size / 2);
      if (d < 2.8 * short) expect(x[i]).toBeCloseTo(expected(i, short), 1);
      if (d < 2.8 * long) expect(y[i]).toBeCloseTo(expected(i, long), 1);
    }
    // a standard deviation of the long axis out: lit along y, dark along x
    const o = Math.round(Size / 2 + long);
    expect(y[o]).toBeGreaterThan(0.5);
    expect(x[o]).toBeLessThan(0.2);
  }, 60_000);

  // stored red first (in front of a camera on +z), blue behind it: drawn as stored, blue would cover red
  it('blends splats back to front for the view, from either side', async () => {
    const pair: Splat[] = [
      { center: Vec3.create(0, 0, 2), scale: Vec3.create(3, 3, 0.5), color: 0xff0000, opacity: 0.9 },
      { center: Vec3.create(0, 0, -2), scale: Vec3.create(3, 3, 0.5), color: 0x0000ff, opacity: 0.9 },
    ];
    const front = centerColor((await render(pair, Vec3.create(0, 0, Distance))).image);
    expect(front[0]).toBeGreaterThan(0.6);
    expect(front[2]).toBeLessThan(0.3);
    const back = centerColor((await render(pair, Vec3.create(0, 0, -Distance))).image);
    expect(back[2]).toBeGreaterThan(0.6);
    expect(back[0]).toBeLessThan(0.3);
  }, 60_000);

  // a splat's group picked where its coverage reaches the picking threshold (0.5), not in its faint tail
  it('picks a splat where it is at least half opaque, not by its tail', async () => {
    const plugin = await sharedPlugin();
    const canvas = plugin.canvas3d!;
    const sigma = 3;
    const repr = splatRepr([
      { center: Vec3.create(-8, 0, 0), scale: Vec3.create(sigma, sigma, sigma) },
      { center: Vec3.create(8, 0, 0), scale: Vec3.create(sigma, sigma, sigma), opacity: 0.4 },
    ]);
    await repr.createOrUpdate(PD.getDefaultValues(GaussianSplats.Params), {}).run();
    canvas.add(repr);
    canvas.commit(true);
    try {
      canvas.camera.setState({ position: Vec3.create(0, 0, Distance), target: Vec3(), up: Vec3.unitY }, 0);
      const pxPerUnit = Size / 2 / Math.tan(canvas.camera.state.fov / 2) / Distance;

      const pixelRatio = canvas.graphics.pixelRatio;
      const { pending } = canvas.graphics.resources;
      let time = 0;
      const pickAt = async (x: number) => {
        canvas.requestDraw();
        canvas.tick((time += 1000) as never, { isSynchronous: true });
        await pending.whenReady();
        const px = Size / 2 + x * pxPerUnit;
        const pick = await canvas.asyncIdentify(Vec2.create(px / pixelRatio, Size / 2 / pixelRatio))?.wait();
        if (!pick) return undefined;
        const loci = canvas.getLoci(pick.id).loci;
        return ShapeGroup.isLoci(loci) ? OrderedSet.getAt(loci.groups[0].ids, 0) : undefined;
      };
      for (let i = 0; i < 20 && (await pickAt(-8)) === undefined; i++) await pending.whenReady();

      expect(await pickAt(-8)).toBe(0);
      // a standard deviation out: coverage exp(-0.5) = 0.61, picked; two out: 0.14, not
      expect(await pickAt(-8 + sigma)).toBe(0);
      expect(await pickAt(-8 + 2 * sigma)).toBeUndefined();
      // opacity 0.4: below the threshold anywhere
      expect(await pickAt(8)).toBeUndefined();
      expect(canvas.graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      canvas.remove(repr);
      canvas.commit(true);
      repr.destroy();
    }
  }, 60_000);

  // degree 1, red on the z basis function (3D Gaussian splatting's `+C1 z`): z of the direction from the camera, -1 from
  // +z (red gone), +1 from -z (red full); green and blue the base color's either way
  it('adds the view-dependent color for the direction from the camera', async () => {
    const sh = encodeSplatSh(1, 1, (_, k, c) => (k === 1 && c === 0 ? 1 : 0));
    const splat: Splat[] = [{ scale: Vec3.create(2, 2, 2), color: 0x808080 }];
    const front = centerColor((await render(splat, Vec3.create(0, 0, Distance), sh)).image);
    const back = centerColor((await render(splat, Vec3.create(0, 0, -Distance), sh)).image);
    expect(front[0]).toBeLessThan(0.1);
    expect(back[0]).toBeGreaterThan(0.8);
    for (const c of [front, back]) {
      expect(c[1]).toBeGreaterThan(0.4);
      expect(c[1]).toBeLessThan(0.55);
    }
  }, 60_000);

  // the splats' programs blended whatever the scene's transparency (their own pass), the GPU sort's kernels with them
  it('a precompile list of splats prepares what their frames draw and sort them with, before they exist', {
    tags: ['medium'],
    timeout: 60_000,
  }, async () => {
    const create = async () => {
      const plugin = await HeadlessPluginContext.create(
        { webgpu: WebGPU.module! },
        DefaultPluginSpec(),
        { width: 32, height: 32 },
        { backend: 'webgpu' },
      );
      await plugin.init();
      return plugin;
    };
    const list: Splat[] = [
      { scale: Vec3.create(2, 2, 2) },
      { center: Vec3.create(1, 0, -1), scale: Vec3.create(1, 1, 1) },
    ];
    const add = async (plugin: HeadlessPluginContext) => {
      const repr = splatRepr(list);
      await repr.createOrUpdate(PD.getDefaultValues(GaussianSplats.Params), {}).run();
      plugin.canvas3d!.add(repr);
      plugin.canvas3d!.commit(true);
    };
    const a = await create();
    const b = await create();
    try {
      a.canvas3d!.resetPrecompileList();
      await add(a);
      const recorded = JSON.parse(JSON.stringify(a.canvas3d!.getPrecompileList()));
      expect(recorded.objects.some((o: { type: string }) => o.type === 'gaussian-splats')).toBe(true);

      const canvas = b.canvas3d!;
      // the scene's transparency other than the splats' own
      canvas.setProps({ transparency: 'wboit' });
      const { pending } = canvas.graphics.resources;
      await canvas.prepare({ interaction: true });
      await pending.whenReady();
      await canvas.precompile(recorded);
      expect(canvas.graphics.compute!.hasKernels(SplatSortKernels)).toBe(true);

      // added: their frames and the first interaction compile nothing more, sorted on the GPU at once
      const before = pending.trackedCount;
      await add(b);
      await canvas.prepare({ interaction: true });
      await pending.whenReady();
      expect(pending.trackedSince(before)).toEqual([]);
      await b.getImageRaw({ width: 32, height: 32 }, NoPostprocessing);
      const values = canvas.getRenderObjects()[0].values as { dSplatGpuOrder: { ref: { value: boolean } } };
      expect(values.dSplatGpuOrder.ref.value).toBe(true);
      expect(canvas.graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      a.dispose();
      b.dispose();
    }
  });
});
