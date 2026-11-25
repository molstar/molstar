/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterAll, describe, expect, it, vi } from 'vitest';
import { now } from '@molstar/core/util/now';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Clip } from '@molstar/core/util/clip';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ValueCell } from '@molstar/core/util/value-cell';

import { HiZPass } from '@molstar/graphics/canvas3d/passes/hi-z';
import { createPrepareBuilder } from '@molstar/graphics/canvas3d/passes/pass-builder';
import { Renderer } from '@molstar/graphics/gl/renderer';
import type { Camera } from '@molstar/graphics/canvas3d/camera';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 128;

/** a compact grid of atoms in the middle of the view, background around it */
const Pdb = `${Array.from({ length: 48 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const residue = String(i + 1).padStart(4);
  const [x, y, z] = [(i % 4) * 2.5, (Math.floor(i / 4) % 4) * 2.5, Math.floor(i / 16) * 2.5];
  return `ATOM  ${serial}  CA  ALA A${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

describe.skipIf(!WebGPU.module)('Hi-Z on WebGPU', () => {
  // the grid drawn in one context for the tests that build passes of their own on it, its camera as it framed the grid
  let shared: Promise<{ plugin: HeadlessPluginContext; camera: Camera.Snapshot }> | undefined;
  const sharedPlugin = async () => {
    shared ??= (async () => {
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
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      const canvas = plugin.canvas3d!;
      canvas.setProps({ multiSample: { ...canvas.props.multiSample, mode: 'off' } });
      // a first frame (which focuses the camera on the grid)
      canvas.commit(true);
      canvas.requestDraw();
      canvas.tick(now(), { isSynchronous: true });
      await canvas.graphics.resources.pending.whenReady();
      return { plugin, camera: canvas.camera.getSnapshot() };
    })();
    const { plugin, camera } = await shared;
    plugin.canvas3d!.camera.setState(camera, 0);
    return plugin;
  };
  afterAll(async () => {
    (await shared)?.plugin.dispose();
  });

  it('primes the first level from an occluder grid, what its clip objects remove hiding nothing', async () => {
    const plugin = await sharedPlugin();
    const canvas = plugin.canvas3d!;
    const { graphics } = canvas;
    const passes = plugin.renderer.canvas3dContext!.passes;
    const hiZ = new HiZPass(graphics, passes.draw, undefined, { enabled: true, minLevel: 1 });
    hiZ.setViewport(0, 0, Size, Size);
    // a dense grid all over the view (the atoms' box, 0 to 7.5 Å, and far around it): every ray opaque in it
    const n = 8;
    const grid = {
      min: Vec3.create(-200, -200, -50),
      max: Vec3.create(200, 200, 50),
      dimensions: [n, n, n] as [number, number, number],
      density: new Float32Array(n * n * n).fill(1),
    };
    /** the clip values of an object (`ClipStorageValues`, as a renderable's) of `objects` */
    const clipOf = (objects: Partial<PD.Values<typeof Clip.Params>['objects'][number]>[]) => {
      const defaults = PD.getDefaultValues(Clip.Params.objects.element);
      const clip = Clip.getClip({ variant: 'pixel', objects: objects.map((o) => ({ ...defaults, ...o })) }).objects;
      return {
        dClipObjects: ValueCell.create(clip.count > 0),
        dClipVariant: ValueCell.create('pixel'),
        dClipping: ValueCell.create(false),
        uClipObjectCount: ValueCell.create(clip.count),
        uClipObjectFrame: ValueCell.create(clip.frame),
        uClipObjectShape: ValueCell.create(clip.shape),
      };
    };
    let time = now();
    const prepare = createPrepareBuilder(graphics, Renderer.create(graphics));
    /** the first level's texels with a depth (not the background's 0), the left and right halves' */
    const primed = async (occluder: typeof grid | undefined, clip?: ReturnType<typeof clipOf>) => {
      hiZ.setOccluder(occluder, clip);
      for (let i = 0; i < 2; ++i) {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await graphics.resources.pending.whenReady();
        hiZ.setOccluder(occluder, clip);
        hiZ.plan(prepare, canvas.camera);
        hiZ.render(canvas.camera);
      }
      const pyramid = hiZ.gpuPyramid!;
      const { texture } = pyramid;
      const width = texture.getWidth();
      const height = texture.getHeight();
      const levels = new Float32Array(width * height);
      const reader = graphics.createTextureReader({ format: 'r32f' });
      reader.start([texture], 0, 0, width, height);
      for (let t = 0; t < 500 && reader.check([levels]) === 'pending'; ++t) {
        await new Promise((resolve) => setTimeout(resolve, 5));
      }
      reader.dispose();
      const ts = 1 << (pyramid.levelCount - 1 - pyramid.minLevel);
      let left = 0;
      let right = 0;
      for (let ty = 0; ty < ts; ++ty) {
        for (let tx = 0; tx < ts; ++tx) {
          if (!(levels[ty * width + tx] > 0)) continue;
          if (tx < ts / 2) left += 1;
          else right += 1;
        }
      }
      return { left, right, all: ts * ts };
    };
    const drawn = await primed(undefined);
    // the grid: every texel primed (its rays opaque within the grid)
    const whole = await primed(grid);
    expect(whole.left + whole.right).toBe(whole.all);
    // a plane through the camera's target clipping what is left of it (turned about z by 90°, its frame's y the
    // world's -x): the left half as drawn (a ray from the eye left of the target there all along)
    const at = Vec3.create(canvas.camera.state.target[0], 0, 0);
    const half = clipOf([{ type: 'plane', position: at, rotation: { axis: Vec3.create(0, 0, 1), angle: 90 } }]);
    const cut = await primed(grid, half);
    expect(cut.right).toBe(whole.right);
    expect(cut.left).toBeLessThan(whole.left / 2);
    expect(cut.left).toBeGreaterThanOrEqual(drawn.left);
    // a cube about all of the grid: all of it clipped away, the level as drawn
    const box = clipOf([{ type: 'cube', scale: Vec3.create(1000, 1000, 1000) }]);
    const gone = await primed(grid, box);
    expect(gone.left + gone.right).toBe(drawn.left + drawn.right);
    // inverted, a small sphere far out of view kept alone: nothing the rays pass kept, the level as drawn
    const far = clipOf([
      { type: 'sphere', invert: true, position: Vec3.create(5000, 0, 0), scale: Vec3.create(2, 2, 2) },
    ]);
    const none = await primed(grid, far);
    expect(none.left + none.right).toBe(drawn.left + drawn.right);
    // the cut moved (its uniforms alone, the grid as it was): the right half clipped instead
    const other = clipOf([{ type: 'plane', position: at, rotation: { axis: Vec3.create(0, 0, 1), angle: -90 } }]);
    const moved = await primed(grid, other);
    expect(moved.left).toBe(whole.left);
    expect(moved.right).toBeLessThan(whole.right / 2);

    // a grid updated in place (a view's step: its array written, a new version naming the cells written): its cells
    // alone uploaded, the level as a copy of it uploaded whole primes it
    const density = new Float32Array(n * n * n).fill(1);
    const base = { ...grid, density, version: 1 };
    expect(await primed(base)).toEqual(whole);
    const texture = (hiZ as unknown as { occluderTexture: { loadRows: (...a: unknown[]) => void } }).occluderTexture;
    const rows = vi.spyOn(texture, 'loadRows');
    for (let z = 0; z < n; ++z) {
      for (let y = 0; y < n; ++y) density.fill(0, (z * n + y) * n, (z * n + y) * n + n / 2);
    }
    const updated = { ...base, version: 2, update: { from: 1, boxes: Int32Array.of(0, 0, 0, n / 2, n, n) } };
    const partial = await primed(updated);
    expect(rows).toHaveBeenCalledTimes(1);
    rows.mockRestore();
    const copy = await primed({ ...grid, density: density.slice() });
    expect(partial).toEqual(copy);
    expect(partial.left).toBeLessThan(whole.left);
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
    hiZ.dispose();
  }, 60_000);

  it('levels whose depth is the farthest of the pixels they cover, frame after frame', async () => {
    const plugin = await sharedPlugin();
    const canvas = plugin.canvas3d!;
    const { graphics } = canvas;
    const passes = plugin.renderer.canvas3dContext!.passes;
    // a pass of its own on the canvas's draw pass: its levels from 1 (4 pixels a texel at this size)
    const hiZ = new HiZPass(graphics, passes.draw, undefined, { enabled: true, minLevel: 1 });
    hiZ.setViewport(0, 0, Size, Size);

    let time = now();
    // built a few times, prepared (planned without drawing) in between as Canvas3D does: a build after a prepare
    // set the first level's size into the second level's uniforms (`uInvSize` shared), wrong from there on
    const prepare = createPrepareBuilder(graphics, Renderer.create(graphics));
    for (let i = 0; i < 4; ++i) {
      canvas.commit(true);
      canvas.requestDraw();
      canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
      await graphics.resources.pending.whenReady();
      hiZ.plan(prepare, canvas.camera);
      hiZ.render(canvas.camera);
    }
    const image = await graphics.readRenderTarget(passes.draw.drawTarget, 0, 0, Size, Size);
    const pyramid = hiZ.gpuPyramid!;
    expect(pyramid).toBeDefined();
    const { texture } = pyramid;
    const width = texture.getWidth();
    const height = texture.getHeight();
    const levels = new Float32Array(width * height);
    const reader = graphics.createTextureReader({ format: 'r32f' });
    reader.start([texture], 0, 0, width, height);
    for (let t = 0; t < 500 && reader.check([levels]) === 'pending'; ++t) {
      await new Promise((resolve) => setTimeout(resolve, 5));
    }
    reader.dispose();

    // level `minLevel` is the first in the texture: a texel with a depth (not the background's 0) has geometry in all
    // the pixels it covers, the farthest depth of them
    const ts = 1 << (pyramid.levelCount - 1 - pyramid.minLevel);
    const px = Size / ts;
    const background = [image[0], image[1], image[2]];
    const isBackground = (x: number, y: number) => {
      const i = (y * Size + x) * 4;
      return image[i] === background[0] && image[i + 1] === background[1] && image[i + 2] === background[2];
    };
    let filled = 0;
    const wrong: number[][] = [];
    for (let ty = 0; ty < ts; ++ty) {
      for (let tx = 0; tx < ts; ++tx) {
        if (!(levels[ty * width + tx] > 0)) continue;
        filled += 1;
        let covered = true;
        for (let y = 0; y < px && covered; ++y) {
          for (let x = 0; x < px && covered; ++x) if (isBackground(tx * px + x, ty * px + y)) covered = false;
        }
        if (!covered) wrong.push([tx, ty]);
      }
    }
    expect(filled).toBeGreaterThan(4);
    expect(wrong).toEqual([]);
    hiZ.dispose();
  }, 60_000);

  // a frame drawn straight to the drawing buffer leaves the opaque depth texture as an earlier frame drew it: levels
  // built from it hid what is visible now (the GPU culling of an image pass, which has no second phase, lost instances)
  it('levels only of a frame that drew the opaque depth texture', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      const canvas = plugin.canvas3d!;
      const context = plugin.renderer.canvas3dContext!;
      // blended: the order-independent modes always draw the texture
      context.setProps({ transparency: 'blended' });
      const { graphics } = canvas;
      // the levels' passes (the canvas's Hi-Z pass is its own, not reachable from here)
      let levelPasses = 0;
      const begin = graphics.beginRenderPass.bind(graphics);
      vi.spyOn(graphics, 'beginRenderPass').mockImplementation((target, options) => {
        if (options?.label === 'hi-z') levelPasses += 1;
        begin(target, options);
      });
      let time = now();
      const draw = async (offscreen: boolean) => {
        canvas.setProps({
          multiSample: { ...canvas.props.multiSample, mode: 'off' },
          illumination: { ...canvas.props.illumination, enabled: false },
          hiZ: { ...canvas.props.hiZ, enabled: true },
          // post-processing (its antialiasing) draws offscreen, none of it straight to the drawing buffer
          postprocessing: { ...canvas.props.postprocessing, enabled: offscreen },
        });
        // the last of a few frames: the first ones skip draws while their pipelines compile (no levels of those)
        for (let i = 0; i < 3; ++i) {
          await graphics.resources.pending.whenReady();
          levelPasses = 0;
          canvas.requestDraw();
          canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        }
      };

      await draw(true);
      expect(context.passes.draw.opaqueDepthDrawn).toBe(true);
      expect(levelPasses).toBeGreaterThan(0);

      await draw(false);
      expect(context.passes.draw.opaqueDepthDrawn).toBe(false);
      expect(levelPasses).toBe(0);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it("asks for a frame again once the levels of a stopped camera's view arrive, once a camera pose", async () => {
    const plugin = await sharedPlugin();
    const canvas = plugin.canvas3d!;
    const { graphics, camera } = canvas;
    const passes = plugin.renderer.canvas3dContext!.passes;
    // a pass of its own on the canvas's draw pass, driven as the canvas's frames drive theirs
    const hiZ = new HiZPass(graphics, passes.draw, undefined, { enabled: true, minLevel: 1 });
    hiZ.setViewport(0, 0, Size, Size);
    const prepare = createPrepareBuilder(graphics, Renderer.create(graphics));
    let time = now();
    const { target } = camera.state;
    const at = (z: number) => {
      camera.setState({ target, position: Vec3.create(target[0], target[1], target[2] + z) }, 0);
      camera.update();
    };
    /** a frame: culled with the levels so far (noted), its levels built from a camera at `levelsZ` (the frame's) */
    const frame = async (levelsZ?: number) => {
      canvas.commit(true);
      canvas.requestDraw();
      canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
      await graphics.resources.pending.whenReady();
      hiZ.noteCulled(camera);
      const z = camera.state.position[2] - target[2];
      if (levelsZ !== undefined) at(levelsZ);
      hiZ.plan(prepare, camera);
      hiZ.render(camera);
      if (levelsZ !== undefined) at(z);
    };
    /** the levels read back (`tick` polled, as the canvas's ticks do) */
    const arrive = async () => {
      for (let i = 0; i < 200 && hiZ.isPending; ++i) {
        await new Promise((resolve) => setTimeout(resolve, 5));
        hiZ.tick();
      }
      expect(hiZ.isPending).toBe(false);
    };

    // the levels of a view, a frame of it culled with them: nothing to draw again
    at(40);
    await frame();
    await arrive();
    await frame();
    await arrive();
    expect(hiZ.takeCatchUp()).toBe(false);
    // the camera moved and stopped: its frame culled with the levels of the view before, those of its own arriving
    at(25);
    await frame();
    await arrive();
    expect(hiZ.takeCatchUp()).toBe(true);
    expect(hiZ.takeCatchUp()).toBe(false);
    // that frame drawn again, with the levels of its view: no more
    await frame();
    await arrive();
    expect(hiZ.takeCatchUp()).toBe(false);
    // levels never of the frame's view (built of another): drawn again once, not again until the camera moves
    at(30);
    await frame(31);
    await arrive();
    expect(hiZ.takeCatchUp()).toBe(true);
    await frame(31);
    await arrive();
    expect(hiZ.takeCatchUp()).toBe(false);
    await frame(31);
    await arrive();
    expect(hiZ.takeCatchUp()).toBe(false);
    at(35);
    await frame();
    await arrive();
    expect(hiZ.takeCatchUp()).toBe(true);
  }, 60_000);
});
