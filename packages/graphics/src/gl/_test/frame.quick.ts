/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The quick allocation tier (`test:perf`) of a frame on WebGL (the GL shim, no GPU: WebGL1, every uniform a classic
// one), at rest (the view unchanged): the opaque pass (the pass walk over a scene, its render items' draws and uniform
// writes), with its positive control. Not yet the renderer's update: it boxes its camera uniforms anew every frame
// (`ValueCell.update`), and its bytes a step drift with the JIT's state.

import { it } from 'vitest';
import { quickCase } from '@molstar/test-harness/quick';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Camera } from '../../canvas3d/camera.js';
import { createFrame } from '../renderable.js';
import { Renderer } from '../renderer.js';
import { Scene } from '../scene.js';
import { createContext } from '../webgl/context.js';
import { createGl } from './gl.shim.js';
import { createMeshObject, createSpheresObject } from './frame-cases.js';

/** A renderer on the shim and a scene of four meshes and four sphere sets (two programs, alternating) */
function createFrameCase() {
  const ctx = createContext(createGl(64, 64, { preserveDrawingBuffer: true }));
  const camera = new Camera({ position: Vec3.create(0, 0, 50) });
  camera.viewport = { x: 0, y: 0, width: 64, height: 64 };
  camera.update();
  const renderer = Renderer.create(ctx);
  renderer.setViewport(0, 0, 64, 64);
  const scene = Scene.create(ctx);
  for (let i = 0; i < 4; ++i) {
    scene.add(createMeshObject(30, i * 2));
    scene.add(createSpheresObject(8, i * 2));
  }
  scene.commit();
  ctx.resources.finalizePrograms(undefined, true);
  const frame = createFrame();
  renderer.update(camera, scene, frame);
  const opaque = () => {
    ctx.beginRenderPass('drawing-buffer');
    renderer.renderOpaque(scene.primitives, camera);
    ctx.endRenderPass();
  };
  opaque();
  return { ctx, opaque };
}

it('opaque pass at rest', ({ task }) => {
  const { ctx, opaque } = createFrameCase();
  quickCase(task, 'frame: opaque pass at rest (webgl)', () => opaque());
  ctx.destroy();
});
