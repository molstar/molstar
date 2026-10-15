/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import {
  Canvas3DParams,
  type Canvas3DProps,
  DefaultCanvas3DParams,
  type PartialCanvas3DProps,
} from '@molstar/graphics/canvas3d/canvas3d';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Representation } from '@molstar/graphics/repr/representation';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { Color } from '@molstar/core/util/color';
import { now } from '@molstar/core/util/now';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { WebGPUUnsupportedError } from '@molstar/graphics/gl/webgpu/errors';
import { type BackendInstance } from './backend.js';
import { DefaultTestSize, type SceneContext, type VisualTestCase } from './test-case.js';

/** Deterministic baseline; tests opt into effects explicitly */
export const BaseCanvasProps: PartialCanvas3DProps = {
  camera: {
    mode: 'perspective',
    fov: 45,
    helper: { axes: { name: 'off', params: {} } },
    stereo: { name: 'off', params: {} },
    manualReset: true,
  },
  cameraFog: { name: 'off', params: {} },
  transparentBackground: false,
  multiSample: { mode: 'off' },
  postprocessing: {
    enabled: true,
    occlusion: { name: 'off', params: {} },
    shadow: { name: 'off', params: {} },
    outline: { name: 'off', params: {} },
    dof: { name: 'off', params: {} },
    antialiasing: { name: 'off', params: {} },
    sharpening: { name: 'off', params: {} },
    background: { variant: { name: 'off', params: {} } },
    bloom: { name: 'off', params: {} },
  },
  // targetFps 0 makes the tracing quality adaptation independent of frame timing
  illumination: { enabled: false, targetFps: 0 },
  hiZ: { enabled: false },
  renderer: {
    backgroundColor: Color(0xffffff),
    enableAnimation: false,
  },
};

export function getTestCanvasProps(test: VisualTestCase): Canvas3DProps {
  const base = PD.merge(Canvas3DParams, DefaultCanvas3DParams, BaseCanvasProps);
  const merged = PD.merge(Canvas3DParams, base, test.canvasProps ?? {});
  // fills params of mapped values that were switched without giving all params
  return PD.normalizeParams(Canvas3DParams, merged, 'all');
}

export interface RenderOutput {
  image: ImageData;
  errors: string[];
  warnings: string[];
  timeMs: number;
}

function captureConsole() {
  const errors: string[] = [];
  const warnings: string[] = [];
  const { error, warn, log } = console;
  console.error = (...args: any[]) => {
    errors.push(args.map(String).join(' '));
    error.apply(console, args);
  };
  console.warn = (...args: any[]) => {
    warnings.push(args.map(String).join(' '));
    warn.apply(console, args);
  };
  // in debug mode, `checkError` consumes GL errors and only logs them
  console.log = (...args: any[]) => {
    const message = args.map(String).join(' ');
    if (message.startsWith('WebGL error:')) errors.push(message);
    log.apply(console, args);
  };
  return {
    errors,
    warnings,
    restore: () => {
      console.error = error;
      console.warn = warn;
      console.log = log;
    },
  };
}

/** `Canvas3D.boundingSphereVisible` is only recalculated on camera reset, which `manualReset` disables */
function getSceneBoundingSphere(renderObjects: readonly GraphicsRenderObject[]) {
  const helper = new BoundaryHelper('98');
  const spheres = renderObjects
    .filter((ro) => ro.state.visible && ro.values.drawCount.ref.value > 0)
    .map((ro) => ro.values.boundingSphere.ref.value)
    .filter((s) => s.radius > 0);
  for (const s of spheres) helper.includeSphere(s);
  helper.finishedIncludeStep();
  for (const s of spheres) helper.radiusSphere(s);
  return helper.getSphere();
}

function fitCamera(camera: Camera, renderObjects: readonly GraphicsRenderObject[], test: VisualTestCase) {
  const { center, radius } = getSceneBoundingSphere(renderObjects);
  const r = Math.max(radius, 0.01);
  const view = test.view ?? {};

  const direction = Vec3.normalize(Vec3(), Vec3.create(...(view.direction ?? [0, 0, -1])));
  const distance = camera.getTargetDistance(r / (view.zoom ?? 1));
  const position = Vec3.scaleAndAdd(Vec3(), center, direction, -distance);

  camera.setState(
    {
      target: Vec3.clone(center),
      position,
      up: Vec3.create(...(view.up ?? [0, 1, 0])),
      radius: r,
      radiusMax: r,
      ...test.camera,
    },
    0,
  );
}

/** Scene context of non-WebGL backends: cases that build GPU data with the WebGL context are not supported */
function getSceneContext(instance: BackendInstance): SceneContext {
  const { context, canvas3d } = instance;
  if (instance.kind === 'webgl') return { graphics: context.graphics, webgl: context.webgl, canvas3d };
  return {
    graphics: context.graphics,
    get webgl(): never {
      throw new WebGPUUnsupportedError('cases using the WebGL context');
    },
    canvas3d,
  };
}

export async function renderTest(instance: BackendInstance, test: VisualTestCase): Promise<RenderOutput> {
  const { canvas, context, canvas3d, imagePass } = instance;
  const { width, height } = test.size ?? DefaultTestSize;
  const start = performance.now();

  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
    canvas.style.width = `${width}px`;
    canvas.style.height = `${height}px`;
    canvas3d.handleResize();
  }

  context.setProps({ transparency: test.contextProps?.transparency ?? 'blended' });

  const p = getTestCanvasProps(test);
  canvas3d.setProps({
    camera: p.camera,
    cameraFog: p.cameraFog,
    cameraClipping: p.cameraClipping,
    transparentBackground: p.transparentBackground,
    dpoitIterations: p.dpoitIterations,
    multiSample: p.multiSample,
    postprocessing: p.postprocessing,
    marking: p.marking,
    illumination: p.illumination,
    hiZ: p.hiZ,
    renderer: p.renderer,
  });

  instance.takeErrors();
  const log = captureConsole();
  let reprs: Representation.Any[] = [];
  try {
    reprs = await test.scene(getSceneContext(instance));
    for (const r of reprs) canvas3d.add(r);
    canvas3d.commit(true);
    fitCamera(canvas3d.camera, canvas3d.getRenderObjects(), test);

    for (let i = 0; i < (test.frames ?? 0); ++i) {
      canvas3d.requestDraw();
      canvas3d.tick(now(), { isSynchronous: true });
      await context.graphics.waitForGpuCommandsComplete();
    }

    imagePass.setProps({
      transparentBackground: p.transparentBackground,
      dpoitIterations: p.dpoitIterations,
      multiSample: p.multiSample,
      postprocessing: p.postprocessing,
      marking: p.marking,
      illumination: p.illumination,
      renderer: p.renderer,
      cameraHelper: { axes: { name: 'off', params: {} } },
    });
    await imagePass.updateBackground();

    const image = await imagePass.getImageData(SyncRuntimeContext, width, height);
    await context.graphics.waitForGpuCommandsComplete();

    return {
      image,
      errors: [...instance.takeErrors(), ...log.errors],
      warnings: log.warnings,
      timeMs: performance.now() - start,
    };
  } finally {
    log.restore();
    for (const r of reprs) {
      canvas3d.remove(r);
      r.destroy();
    }
    canvas3d.commit(true);
  }
}
