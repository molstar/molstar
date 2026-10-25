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
import { isWebGLContext } from '@molstar/graphics/gl/webgl/context';

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
  /** Loaded from the reference profile's own run, not rendered (see `RunOptions.sharedReferencesUrl`) */
  shared?: boolean;
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
  if (isWebGLContext(context.graphics)) return { graphics: context.graphics, webgl: context.graphics, canvas3d };
  return {
    graphics: context.graphics,
    get webgl(): never {
      throw new WebGPUUnsupportedError('cases using the WebGL context');
    },
    canvas3d,
  };
}

/** The image of `canvas` (its drawing buffer, or the current texture of a WebGPU canvas) */
function readCanvas(canvas: HTMLCanvasElement): ImageData {
  const copy = document.createElement('canvas');
  copy.width = canvas.width;
  copy.height = canvas.height;
  const ctx = copy.getContext('2d')!;
  ctx.drawImage(canvas, 0, 0);
  return ctx.getImageData(0, 0, copy.width, copy.height);
}

/** the error of a release whose context was not restored (`releaseAndRestore`) */
const NotRestored = 'not restored';

/**
 * Releases the context of `instance` and restores it (`Canvas3DContext.release`, `restore`): the release is reported
 * as a loss and not restored on its own, the restore gives a new live context that Canvas3D rebuilds the scene on
 */
async function releaseAndRestore(instance: BackendInstance): Promise<string[]> {
  const { context, canvas3d } = instance;
  const released = context.graphics;
  const objectCount = canvas3d.getRenderObjects().length;
  const lost = new Promise<void>((resolve) => {
    const sub = context.contextLost!.subscribe(() => {
      sub.unsubscribe();
      resolve();
    });
  });
  if (!context.release?.()) return ['release refused'];
  await lost;
  const errors: string[] = [];
  if (!context.isReleased) errors.push('not released once the loss was reported');
  await new Promise((resolve) => setTimeout(resolve, 100));
  if (context.graphics !== released || !released.isContextLost) errors.push('restored without restore');
  // a busy GPU process may hand out the new device late: a released context stays released, restored once more
  if (!(await context.restore!()) && !(await context.restore!())) return [...errors, NotRestored];
  if (context.graphics === released || context.graphics.isContextLost) errors.push('no new live context');
  // the representations update on the new context asynchronously, then are added again
  const deadline = performance.now() + 10_000;
  while (canvas3d.getRenderObjects().length < objectCount && performance.now() < deadline) {
    await new Promise((resolve) => setTimeout(resolve, 10));
  }
  canvas3d.commit(true);
  if (canvas3d.getRenderObjects().length < objectCount) errors.push('scene not rebuilt on the restored context');
  await context.graphics.resources.pending.whenReady();
  return errors;
}

/** How a case renders besides its own props */
export interface RenderOptions {
  /**
   * The renderer's draw lists on (`RendererProps.drawLists`), the case rendered three times and the last captured: the
   * first may skip draws while pipelines compile, the second records, the third replays
   */
  readonly drawLists?: boolean;
}

export async function renderTest(
  instance: BackendInstance,
  test: VisualTestCase,
  options: RenderOptions = {},
): Promise<RenderOutput> {
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
  p.renderer = { ...p.renderer, drawLists: !!options.drawLists };
  const repeats = options.drawLists ? 2 : 0;
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
    viewport: p.viewport,
  });

  instance.takeErrors();
  const log = captureConsole();
  let reprs: Representation.Any[] = [];
  try {
    const sceneContext = getSceneContext(instance);
    reprs = await test.scene(sceneContext);
    for (const r of reprs) canvas3d.add(r);
    canvas3d.commit(true);
    fitCamera(canvas3d.camera, canvas3d.getRenderObjects(), test);
    let interactErrors: string[] = [];
    if (test.interact) {
      // a frame updates the camera to the fitted view, as picking and projecting expect
      await context.graphics.resources.pending.whenReady();
      canvas3d.requestDraw();
      canvas3d.tick(now(), { isSynchronous: true });
      interactErrors = await test.interact(sceneContext, reprs);
      // marks (`canvas3d.mark`) are applied with the next frame
      canvas3d.requestDraw();
      canvas3d.tick(now(), { isSynchronous: true });
    }

    for (let i = 0; i < (test.frames ?? 0); ++i) {
      canvas3d.requestDraw();
      canvas3d.tick(now(), { isSynchronous: true });
      await context.graphics.waitForGpuCommandsComplete();
    }

    const releaseErrors = test.releaseContext ? await releaseAndRestore(instance) : [];
    // capturing would read from the released context: its error instead of the read's
    if (releaseErrors.includes(NotRestored)) throw new Error(releaseErrors.join('; '));

    if (test.capture === 'canvas') {
      for (let i = 0; i <= repeats; ++i) {
        await context.graphics.resources.pending.whenReady();
        canvas3d.requestDraw();
        canvas3d.tick(now(), { isSynchronous: true });
        if (i < repeats) await context.graphics.waitForGpuCommandsComplete();
      }
      // in the task of the frame: a WebGPU canvas keeps its image only until it is presented
      const image = readCanvas(canvas);
      await context.graphics.waitForGpuCommandsComplete();
      return {
        image,
        errors: [...interactErrors, ...releaseErrors, ...instance.takeErrors(), ...log.errors],
        warnings: log.warnings,
        timeMs: performance.now() - start,
      };
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
    // no `updateBackground`: the render waits for the background images itself (pending work)

    for (let i = 0; i < repeats; ++i) {
      await imagePass.getImageData(SyncRuntimeContext, width, height);
      await context.graphics.resources.pending.whenReady();
    }
    const { drawLists } = context.graphics.stats;
    const replayed = drawLists.replayed;
    const image = await imagePass.getImageData(SyncRuntimeContext, width, height);
    await context.graphics.waitForGpuCommandsComplete();
    if (options.drawLists && drawLists.replayed === replayed)
      log.warnings.push('draw lists: the captured frame replayed nothing');

    return {
      image,
      errors: [...interactErrors, ...releaseErrors, ...instance.takeErrors(), ...log.errors],
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
