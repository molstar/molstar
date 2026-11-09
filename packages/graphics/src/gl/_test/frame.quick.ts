/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The quick allocation tier (`test:perf`) of a frame at rest (the view unchanged): the renderer's update, and its opaque
// pass (the pass walk over a scene, its render items' draws and uniform writes; on WebGPU the draw list's replay). On
// WebGL the GL shim (WebGL1: every uniform a classic one), on WebGPU a stub device (stub-gpu.ts); no GPU. Each with its
// positive control.

import { it } from 'vitest';
import { quickCase } from '@molstar/test-harness/quick';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Camera } from '../../canvas3d/camera.js';
import type { GraphicsContext } from '../context.js';
import { createFrame } from '../renderable.js';
import { Renderer } from '../renderer.js';
import { Scene } from '../scene.js';
import { createContext } from '../webgl/context.js';
import { createWebGPUContext } from '../webgpu/context.js';
import { createGl } from './gl.shim.js';
import { createMeshObject, createSpheresObject } from './frame-cases.js';
import { createStubDevice } from './stub-gpu.js';

const Size = 64;

/**
 * A renderer and a scene of four meshes and four sphere sets (two programs, alternating), its programs or pipelines
 * ready (WebGPU compiles them in the background: frames until none skips a draw)
 */
async function createFrameCase(backend: 'webgl' | 'webgpu') {
  const ctx: GraphicsContext =
    backend === 'webgl'
      ? createContext(createGl(Size, Size, { preserveDrawingBuffer: true }))
      : createWebGPUContext(createStubDevice(), {
          width: Size,
          height: Size,
          // the device probes these by drawing and reading back, which a stub can't: as a GPU has them
          fragDepthInBundles: true,
          viewportOutsideTarget: true,
        });
  const camera = new Camera({ position: Vec3.create(0, 0, 50) });
  camera.viewport = { x: 0, y: 0, width: Size, height: Size };
  camera.update();
  const renderer = Renderer.create(ctx);
  renderer.setViewport(0, 0, Size, Size);
  const scene = Scene.create(ctx);
  for (let i = 0; i < 4; ++i) {
    scene.add(createMeshObject(30, i * 2));
    scene.add(createSpheresObject(8, i * 2));
  }
  scene.commit();
  if (backend === 'webgl') ctx.resources.finalizePrograms(undefined, true);
  const frame = createFrame();
  const update = () => renderer.update(camera, scene, frame);
  // a frame's task ends with its commands submitted (WebGPU: in a microtask, which a measured loop never reaches)
  const opaque = () => {
    ctx.beginRenderPass('drawing-buffer');
    renderer.renderOpaque(scene.primitives, camera);
    ctx.endRenderPass();
    ctx.flush();
  };
  update();
  for (let k = 0; k < 20; ++k) {
    const skipped = ctx.stats.skippedDraws;
    opaque();
    await ctx.resources.pending.whenReady();
    if (k > 2 && ctx.stats.skippedDraws === skipped) break;
  }
  return { ctx, update, opaque };
}

for (const backend of ['webgl', 'webgpu'] as const) {
  it(`update at rest (${backend})`, async ({ task }) => {
    const { ctx, update } = await createFrameCase(backend);
    quickCase(task, `frame: update at rest (${backend})`, () => update());
    ctx.destroy();
  });

  it(`opaque pass at rest (${backend})`, async ({ task }) => {
    const { ctx, opaque } = await createFrameCase(backend);
    quickCase(task, `frame: opaque pass at rest (${backend})`, () => opaque());
    ctx.destroy();
  });
}
