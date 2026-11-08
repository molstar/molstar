/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StereoCameraParams } from '@molstar/graphics/canvas3d/camera/stereo';
import { Canvas3D, Canvas3DContext, Canvas3DParams, type Canvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { PostprocessingParams } from '@molstar/graphics/canvas3d/passes/postprocessing';
import { IlluminationParams } from '@molstar/graphics/canvas3d/passes/illumination';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { TimerResult } from '@molstar/graphics/gl/timer';
import type { ResourceByteCounts } from '@molstar/graphics/gl/stats';
import { Interval } from '@molstar/core/data/int/interval';
import { Color } from '@molstar/core/util/color';
import { getDistinctBaseColors } from '@molstar/core/util/color/distinct';
import { AutoGpuCullingMinInstances } from '@molstar/graphics/geo/geometry/base';
import { StructureElement } from '@molstar/model/model/structure';
import { Theme } from '@molstar/graphics/theme/theme';
import { ParticleSizeThemeProvider } from '@molstar/graphics/theme/size/particle-size';
import { getThemeRegistries } from '../visual/scenes/structure.js';
import type { Representation, RepresentationProvider } from '@molstar/graphics/repr/representation';
import { BallAndStickRepresentationProvider } from '@molstar/graphics/repr/structure/representation/ball-and-stick';
import { LineRepresentationProvider } from '@molstar/graphics/repr/structure/representation/line';
import { CartoonRepresentationProvider } from '@molstar/graphics/repr/structure/representation/cartoon';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { MolecularSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/molecular-surface';
import { SpacefillRepresentationProvider } from '@molstar/graphics/repr/structure/representation/spacefill';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import type { ValueCell } from '@molstar/core/util/value-cell';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { isWebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { type GraphicsRenderObject, isMergedRenderObject } from '@molstar/graphics/gl/render-object';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { now } from '@molstar/core/util/now';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { OutlineParams } from '@molstar/graphics/canvas3d/passes/outline';
import { SsaoParams } from '@molstar/graphics/canvas3d/passes/ssao';
import { HiZParams } from '@molstar/graphics/canvas3d/passes/hi-z';
import { ClipGhostParams } from '@molstar/graphics/canvas3d/passes/clip-ghost';
import { Structure, Unit } from '@molstar/model/model/structure';
import { CIF } from '@molstar/io/reader/cif';
import { MmcifFormat, trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { getMatrices, operatorGroupsProvider } from '@molstar/model/formats/structure/property/assembly';
import { mergeUnitsWithSameOperator } from '@molstar/model/model/structure/structure/util/unit-merging';
import { Task } from '@molstar/core/task';
import { createParticleListFromMmcifAssembly, getAssemblyIdsFromMmcif } from '@molstar/model/formats/particles/mmcif';
import { Particle, type ParticleList } from '@molstar/model/model/particles/particle-list';
import { EveryLoci } from '@molstar/model/model/loci';
import { OrderedSet } from '@molstar/core/data/int';
import { ParticleTargetRepresentationProvider } from '@molstar/graphics/repr/particles/representation/target/representation';
import { ParticleTargetRepresentationParams } from '@molstar/graphics/repr/particles/representation/target/common';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { setQuadIndexing } from '@molstar/graphics/gl/quad-indices';
import { meshBox, shapeRepr } from '../visual/scenes/primitives.js';
import { structureRepr } from '../visual/scenes/structure.js';
import { type BenchmarkCase, BenchmarkEnv, collectedHeap, type PathRun, waitForGpu } from './framework.js';
import { chainGrid, exampleStructure, globuleStructure } from './structures.js';

/** Drawing buffer of the render cases, at pixel scale 1 */
export const Size = { width: 1280, height: 720 };
/** Camera rotation per frame of rotating cases, about the camera's up axis */
const RotationStep = (2 * Math.PI) / 180;

/** Frames per run, each forced (all passes drawn), issued back to back */
const Frames = 30;

export type RenderInput = {
  readonly scene: (graphics: GraphicsContext) => Promise<Representation.Any[]>;
  readonly props?: Partial<Canvas3DProps>;
  /**
   * multi-sampling `on` (every sample a frame, the samples views of the same culled frame), off by default: progressive
   * multi-sampling would make frames differ
   */
  readonly multiSample?: boolean;
  /** the context's transparency mode (`Canvas3DContext.Props`), the default restored with the scene's disposal */
  readonly transparency?: Canvas3DContext.Props['transparency'];
  /** rotates the camera by `RotationStep` per frame, as an interactive rotation (culling and globals change) */
  readonly rotate?: boolean;
  /** the rotation per frame in degrees with `rotate` (default `RotationStep`'s 2) */
  readonly rotateDegrees?: number;
  /** moves the camera towards the target by this fraction of its distance per frame, from the `zoom` of the case */
  readonly dolly?: number;
  /**
   * yields to the event loop after each frame (a task, as an app's frames do): readbacks complete only then (Hi-Z's;
   * WebGL's fences signal between tasks), the frames in one task would never see them; counted in `wait`
   */
  readonly yieldFrames?: boolean;
  /** drawing buffer, `Size` by default (larger for full-screen passes, so that their GPU times stand out of the noise) */
  readonly size?: { readonly width: number; readonly height: number };
  /** moves the camera towards the target, to this fraction of its distance, after the reset */
  readonly zoom?: number;
  /**
   * the camera's radius (what near and far fit) as this fraction of its radius after the reset: below 1 the near plane
   * cuts the scene's front (solid-interior caps)
   */
  readonly cameraRadius?: number;
  /** false: the scene's quads (spheres, points) drawn unindexed (`setQuadIndexing`), for comparisons */
  readonly quadIndexing?: boolean;
  /** everything highlighted: the marking passes drawn per frame besides the scene's */
  readonly highlight?: boolean;
  /** a hover sweeping across the objects: per frame another representation highlighted whole, the last one not */
  readonly sweep?: boolean;
  /**
   * frames not forced (`requestDraw`): drawn by what changed before them, a marking change alone redrawn over the last
   * image as a hover is (`MarkingPass.redraw`)
   */
  readonly unforced?: boolean;
  /**
   * a pause before each frame, ms, left out of its time: the GPU's clocks fall as between a user's hovers (a redraw after
   * one runs slower than back to back)
   */
  readonly idleMs?: number;
  /** the time the scene took to set up (`SceneSetup`) among the notes */
  readonly setupNotes?: boolean;
  /** the GPU bytes of the scene (`getByteCounts` after its setup) among the notes */
  readonly byteNotes?: boolean;
  /**
   * a change before each frame (a hover, a theme): its time a stage of its own (`update`; with draws it forces, e.g. a
   * commit), its uploads the frame's
   */
  readonly update?: (canvas3d: Canvas3D, reprs: readonly Representation.Any[], frame: number) => void | Promise<void>;
  /** measurements for render bundles (WebGPU), see `RenderBundleCases` */
  readonly bundles?: {
    /** the opaque pass recorded into a render bundle and replayed, every frame (the recording cost) */
    readonly record?: boolean;
    /** the visible set (the objects in the frustum) per frame: its CPU time, and how often it changed */
    readonly visibility?: boolean;
    /** no object culling (the frustum test passes every object) */
    readonly noCulling?: boolean;
  };
};

/** Render bundle recording of the opaque pass (`RenderInput.bundles.record`): recorded and failed recordings */
const bundleStats = { recorded: 0, failed: 0 };

/** `Renderer.create` with the opaque pass recorded into a render bundle and replayed, where the context records */
function withBundleRecording<T>(graphics: GraphicsContext, create: () => T): T {
  if (!isWebGPUContext(graphics)) return create();
  const original = Renderer.create;
  Renderer.create = (...args: Parameters<typeof Renderer.create>) => {
    const renderer = original(...args);
    const renderOpaque = renderer.renderOpaque;
    renderer.renderOpaque = (group, camera) => {
      const bundle = graphics.recordRenderBundle(() => renderOpaque(group, camera));
      if (bundle) {
        bundleStats.recorded += 1;
        graphics.executeRenderBundles([bundle]);
      } else {
        bundleStats.failed += 1;
        renderOpaque(group, camera);
      }
    };
    return renderer;
  };
  try {
    return create();
  } finally {
    Renderer.create = original;
  }
}

/** The objects in the camera's frustum, as a key (one character per render object) */
function getVisibleSet(canvas3d: Canvas3D, frustum: Frustum3D) {
  Frustum3D.fromProjectionMatrix(frustum, canvas3d.camera.projectionView);
  let key = '';
  for (const o of canvas3d.getRenderObjects()) {
    key += o.state.visible && Frustum3D.intersectsSphere3D(frustum, o.values.boundingSphere.ref.value) ? '1' : '0';
  }
  return key;
}

/** Drawing buffer of cases timing full-screen passes */
export const LargeSize = { width: 2560, height: 1440 };

/**
 * The time to set up a case's scene: creating its representations, adding and committing them, then the first frame
 * drawn (pipelines ready, the GPU done; draw lists recorded) in the case's first run
 */
type SceneSetup = {
  readonly reprsMs: number;
  readonly commitMs: number;
  firstFrameMs?: number;
  /** the GPU bytes after the first frame */
  bytes?: ResourceByteCounts;
  /** the CPU bytes of the render objects' values after the first frame */
  cpuBytes?: CpuBytes;
  /** the JS heap before the scene was created (collected, `collectedHeap`) */
  readonly heapBefore?: number;
  /** the JS heap the scene retains after its first frame: its representations, render objects and their resources */
  heapBytes?: number;
  /** the notes of `getRasterSplitNotes` */
  rasterSplit?: string[];
};

type ByteKinds = { geometry: number; instances: number; layers: number };
/** of the render objects, merged objects' members included; `merged` the merged objects' own copies */
type CpuBytes = { objects: ByteKinds; merged: ByteKinds };

const InstanceKeys = /^(aTransform|aInstance|aSegment|instanceGrid)/;
/** the items' data of spheres, points, cylinders, lines and text: geometry, though data arrays (`t…`) */
const ItemGeometryKeys = /^(tPositionGroup|tCylinder|tLine|tText)$/;

/**
 * The CPU bytes of render objects' values, each `ArrayBuffer` counted once: theme layers (`t*`, data texture images),
 * instances (transforms, instance and segment data), the rest geometry
 */
function getCpuBytes(reprs: readonly Representation.Any[]): CpuBytes {
  const objects: ByteKinds = { geometry: 0, instances: 0, layers: 0 };
  const merged: ByteKinds = { geometry: 0, instances: 0, layers: 0 };
  const seen = new Set<ArrayBufferLike>();
  const add = (o: GraphicsRenderObject, into: ByteKinds) => {
    const values = o.values as { [k: string]: { ref?: { value: unknown } } | undefined };
    for (const k of Object.keys(values)) {
      const v = values[k]?.ref?.value as { array?: unknown } | undefined;
      const array = ArrayBuffer.isView(v) ? v : ArrayBuffer.isView(v?.array) ? v.array : undefined;
      if (!array || seen.has(array.buffer)) continue;
      seen.add(array.buffer);
      const kind =
        k.startsWith('t') && !ItemGeometryKeys.test(k) ? 'layers' : InstanceKeys.test(k) ? 'instances' : 'geometry';
      into[kind] += array.buffer.byteLength;
    }
  };
  for (const r of reprs) {
    for (const o of r.renderObjects) {
      if (isMergedRenderObject(o)) for (const m of o.members) add(m, objects);
    }
  }
  for (const r of reprs) {
    for (const o of r.renderObjects) add(o, isMergedRenderObject(o) ? merged : objects);
  }
  return { objects, merged };
}

/** The scene of the case rendered last, kept across its runs */
let current: { id: string; canvas3d: Canvas3D; reprs: Representation.Any[]; setup: SceneSetup } | undefined;

export async function getScene(id: string, graphics: GraphicsContext, input: RenderInput) {
  if (current?.id === id) return current.canvas3d;
  disposeScene();

  const context = BenchmarkEnv.context;
  if (!context?.canvas || context.graphics !== graphics) throw new Error('render cases need the Canvas3D context');
  const canvas = context.canvas;
  const { width, height } = input.size ?? Size;
  canvas.width = width;
  canvas.height = height;
  canvas.style.width = `${width}px`;
  canvas.style.height = `${height}px`;

  if (input.transparency) context.setProps({ transparency: input.transparency });
  // progressive multi-sampling would make frames differ, picking warm-up would run in between
  const create = () =>
    Canvas3D.create(
      context,
      {
        ...input.props,
        ...(BenchmarkEnv.drawLists !== undefined
          ? { renderer: { ...(input.props?.renderer ?? DrawLists.renderer!), drawLists: BenchmarkEnv.drawLists } }
          : {}),
        multiSample: { ...PD.getDefaultValues(Canvas3DParams).multiSample, mode: input.multiSample ? 'on' : 'off' },
      },
      { warmUp: false },
    );
  const canvas3d = input.bundles?.record ? withBundleRecording(graphics, create) : create();
  canvas3d.handleResize();
  const heapBefore = collectedHeap();
  const t0 = performance.now();
  const reprs = await input.scene(graphics);
  const t1 = performance.now();
  // render items decide at creation whether they draw their quads indexed
  setQuadIndexing(input.quadIndexing ?? true);
  try {
    for (const r of reprs) canvas3d.add(r);
    canvas3d.commit(true);
  } finally {
    setQuadIndexing(true);
  }
  const t2 = performance.now();
  if (input.highlight) canvas3d.mark({ loci: { kind: 'every-loci' } }, MarkerAction.Highlight);
  canvas3d.requestCameraReset({ durationMs: 0 });
  canvas3d.tick(now(), { isSynchronous: true });
  const setup: SceneSetup = { reprsMs: t1 - t0, commitMs: t2 - t1, heapBefore };
  if (input.zoom) {
    const { position, target } = canvas3d.camera.getSnapshot();
    const offset = Vec3.scale(Vec3(), Vec3.sub(Vec3(), position, target), input.zoom);
    canvas3d.camera.setState({ position: Vec3.add(Vec3(), target, offset) }, 0);
    canvas3d.tick(now(), { isSynchronous: true });
  }
  if (input.cameraRadius) {
    canvas3d.camera.setState({ radius: canvas3d.camera.getSnapshot().radius * input.cameraRadius }, 0);
    canvas3d.tick(now(), { isSynchronous: true });
  }
  current = { id, canvas3d, reprs, setup };
  return canvas3d;
}

export function disposeScene() {
  if (!current) return;
  current.canvas3d.dispose();
  BenchmarkEnv.context?.setProps({ transparency: Canvas3DContext.DefaultProps.transparency });
  for (const r of current.reprs) r.destroy();
  current = undefined;
}

function drawCalls(graphics: GraphicsContext) {
  const { calls } = graphics.stats;
  return calls.drawInstanced + calls.drawInstancedBase + calls.multiDrawInstancedBase;
}

function resetDrawCalls(graphics: GraphicsContext) {
  const { calls } = graphics.stats;
  calls.drawInstanced = calls.drawInstancedBase = calls.multiDrawInstancedBase = calls.counts = 0;
}

/**
 * `Frames` forced frames issued back to back, then waited for: per frame, the time to issue it (`submit`, the CPU side)
 * and the rest of the time until the GPU completed all (`wait`); their sum is the frame time at full throughput (CPU and
 * GPU overlap, the wait for the last frame is spread over all)
 */
export async function renderFrames(id: string, graphics: GraphicsContext, input: RenderInput): Promise<PathRun> {
  const canvas3d = await getScene(id, graphics, input);
  if (input.update && current) {
    // an update first, outside the timing: the programs it needs (a color type, a marker type) compiled, as after any
    // first change of that kind
    await input.update(canvas3d, current.reprs, 0);
    // its programs are created by the first frames drawing it and compile in the background, draws skipped meanwhile
    // (`whenReady` waits only for what is compiling already), frames held while what they draw with is prepared: frames
    // until a few in a row are drawn without skipping (the first after the update may not draw what it changed yet)
    canvas3d.notifyDidDraw = true;
    for (let k = 0, clean = 0; k < 50 && clean < 3; ++k) {
      const skipped = graphics.stats.skippedDraws;
      let drawn = false;
      const sub = canvas3d.didDraw.subscribe(() => (drawn = true));
      canvas3d.requestDraw();
      canvas3d.tick(now(), { isSynchronous: true });
      sub.unsubscribe();
      await graphics.resources.pending.whenReady();
      clean = drawn && graphics.stats.skippedDraws === skipped ? clean + 1 : 0;
      if (clean === 0) await new Promise((resolve) => setTimeout(resolve, 10));
    }
  }
  const first = performance.now();
  // the first frame that draws everything: programs (WebGL) and pipelines (WebGPU) compiling in the background skip
  // draws, and a frame links what it needs first (`whenReady` waits only for what is compiling already)
  for (let k = 0; k < 200; ++k) {
    await graphics.resources.pending.whenReady();
    const skipped = graphics.stats.skippedDraws;
    canvas3d.requestDraw();
    canvas3d.tick(now(), { isSynchronous: true });
    await waitForGpu(graphics);
    if (graphics.stats.skippedDraws === skipped) break;
    await new Promise((resolve) => setTimeout(resolve, 5));
  }
  if (current && current.setup.firstFrameMs === undefined) {
    current.setup.firstFrameMs = performance.now() - first;
    current.setup.bytes = { ...graphics.resources.getByteCounts() };
    current.setup.cpuBytes = getCpuBytes(current.reprs);
    const heap = current.setup.heapBefore === undefined ? undefined : collectedHeap();
    if (heap !== undefined) current.setup.heapBytes = heap - current.setup.heapBefore!;
  }

  let submit = 0;
  let updateMs = 0;
  let calls = 0;
  let drawn = 0;
  canvas3d.notifyDidDraw = true;
  const sub = canvas3d.didDraw.subscribe(() => drawn++);
  drawn = 0;
  bundleStats.recorded = bundleStats.failed = 0;
  graphics.stats.drawLists.replayed = graphics.stats.drawLists.recorded = 0;
  graphics.stats.drawLists.unplanned = graphics.stats.drawLists.flushed = 0;
  graphics.stats.drawLists.chunksRecorded = graphics.stats.drawLists.chunksDirect = 0;
  const frustum = Frustum3D();
  let visibilityMs = 0;
  let visibleKey: string | undefined;
  let visibleChanges = 0;
  const intersects = Frustum3D.intersectsSphere3D;
  if (input.bundles?.noCulling) Frustum3D.intersectsSphere3D = () => true;
  const skippedBefore = graphics.stats.skippedDraws;
  let idle = 0;
  const start = performance.now();
  try {
    for (let i = 0; i < Frames; ++i) {
      if (input.idleMs) {
        // after the last frame's GPU work, so that the pause is idle; only the pause left out (the wait is the
        // frame's)
        await waitForGpu(graphics);
        const ti = performance.now();
        await new Promise((resolve) => setTimeout(resolve, input.idleMs));
        idle += performance.now() - ti;
      }
      if (input.rotate) rotateCamera(canvas3d, input.rotateDegrees);
      if (input.dolly) dollyCamera(canvas3d, input.dolly, i);
      if (input.sweep && current) sweepHighlight(canvas3d, current.reprs, i);
      if (input.update && current) {
        const tu = performance.now();
        const drawnBefore = drawn;
        await input.update(canvas3d, current.reprs, i);
        // draws the update forced (a commit) are its own, not the frame's
        drawn = drawnBefore;
        updateMs += performance.now() - tu;
      }
      resetDrawCalls(graphics);
      const t0 = performance.now();
      if (!input.unforced) canvas3d.requestDraw();
      canvas3d.tick(now(), { isSynchronous: true });
      submit += performance.now() - t0;
      calls = drawCalls(graphics);
      if (input.yieldFrames) await nextTask();
      if (input.bundles?.visibility) {
        const t1 = performance.now();
        const key = getVisibleSet(canvas3d, frustum);
        visibilityMs += performance.now() - t1;
        if (visibleKey !== undefined && key !== visibleKey) visibleChanges += 1;
        visibleKey = key;
      }
    }
  } finally {
    Frustum3D.intersectsSphere3D = intersects;
  }
  await waitForGpu(graphics);
  const total = performance.now() - start - idle;
  sub.unsubscribe();
  const notes: string[] = [];
  if (input.setupNotes && current) {
    const { reprsMs, commitMs, firstFrameMs } = current.setup;
    notes.push(
      `setup: representations ${reprsMs.toFixed(0)} ms, add and commit ${commitMs.toFixed(0)} ms, first frame ${firstFrameMs?.toFixed(0)} ms`,
    );
  }
  if ((input.byteNotes || BenchmarkEnv.byteNotes) && current?.setup.bytes) {
    const b = current.setup.bytes;
    const mb = (n: number) => (n / 1048576).toFixed(1);
    notes.push(
      `GPU bytes: attributes ${mb(b.attribute)} MB, elements ${mb(b.elements)} MB, textures ${mb(b.texture)} MB, uniforms ${mb(b.uniform)} MB, storage ${mb(b.storageBuffer ?? 0)} MB`,
    );
    const c = current.setup.cpuBytes;
    if (c) {
      const kinds = (k: ByteKinds) =>
        `geometry ${mb(k.geometry)}, instances ${mb(k.instances)}, layers ${mb(k.layers)}`;
      notes.push(`CPU bytes: objects ${kinds(c.objects)} MB; merged copies ${kinds(c.merged)} MB`);
    }
    const heap = current.setup.heapBytes;
    if (heap !== undefined) notes.push(`JS heap retained by the scene: ${mb(heap)} MB`);
  }
  // draws of objects whose pipelines were not ready (WebGPU): frames that did less than the case measures
  const skipped = graphics.stats.skippedDraws - skippedBefore;
  if (skipped > 0) notes.push(`skipped ${skipped} draws while pipelines compiled`);
  const { drawLists } = graphics.stats;
  if (drawLists.replayed || drawLists.recorded) {
    notes.push(
      `draw lists: replayed ${drawLists.replayed} (unplanned ${drawLists.unplanned}, flushed ${drawLists.flushed}), recorded ${drawLists.recorded} (${drawLists.chunksRecorded} chunks), ${drawLists.chunksDirect} chunks drawn directly`,
    );
  }
  // the replay walks spheres in JS: once per scene, and not without a raster (every sphere of the scene to walk)
  if (BenchmarkEnv.rasterSplitNotes && current && input.props?.renderer?.sphereRaster !== 'off') {
    current.setup.rasterSplit ??= getRasterSplitNotes(canvas3d, current.reprs, input);
    notes.push(...current.setup.rasterSplit);
  }
  if (input.bundles?.record) notes.push(`opaque pass recorded ${bundleStats.recorded}, failed ${bundleStats.failed}`);
  if (input.bundles?.visibility) {
    const visible = visibleKey?.split('').filter((v) => v === '1').length ?? 0;
    notes.push(
      `visible set: ${(visibilityMs / Frames).toFixed(3)} ms per frame, changed in ${visibleChanges} of ${Frames - 1} frames, ${visible} of ${visibleKey?.length ?? 0} objects visible at the end`,
    );
  }
  if (drawn !== Frames) throw new Error(`${drawn} of ${Frames} frames drawn`);
  const gpuPasses = BenchmarkEnv.gpuTiming ? await getGpuPasses(graphics) : undefined;
  if (BenchmarkEnv.images && BenchmarkEnv.context?.canvas) {
    // a frame of its own, read in the task that drew it (before the canvas is presented and cleared); rotating cases
    // turned on as in their frames (what lags behind the view shows)
    if (input.rotate) rotateCamera(canvas3d, input.rotateDegrees);
    if (input.dolly) dollyCamera(canvas3d, input.dolly, Frames);
    canvas3d.requestDraw();
    canvas3d.tick(now(), { isSynchronous: true });
    BenchmarkEnv.images[id] = BenchmarkEnv.context.canvas.toDataURL('image/png');
  }
  return {
    stages: [
      ...(input.update ? [{ label: 'update', ms: updateMs / Frames }] : []),
      { label: 'submit', ms: submit / Frames },
      { label: 'wait', ms: (total - submit - updateMs) / Frames },
    ],
    triangles: 0,
    drawCalls: calls,
    gpuPasses,
    notes: notes.length ? notes : undefined,
  };
}

/**
 * The GPU times of the frames' pass groups, averaged per frame by path of labels (e.g. `Canvas3D.render > DrawPass.render
 * > pass opaque`); waits for the timestamp queries of the last frames
 */
async function getGpuPasses(graphics: GraphicsContext) {
  const sums: { [path: string]: number } = {};
  let frames = 0;
  const add = (results: TimerResult[], prefix: string) => {
    for (const r of results) {
      const path = prefix ? `${prefix} > ${r.label}` : r.label;
      if (!isNaN(r.gpuElapsed)) sums[path] = (sums[path] ?? 0) + r.gpuElapsed / 1e6;
      add(r.children, path);
    }
  };
  for (let i = 0; i < 50; ++i) {
    const results = graphics.timer.resolve();
    frames += results.filter((r) => r.label === 'Canvas3D.render').length;
    add(results, '');
    if (frames >= Frames) break;
    await new Promise((resolve) => setTimeout(resolve, 10));
  }
  const avgs: { [path: string]: number } = {};
  for (const k of Object.keys(sums)) avgs[k] = sums[k] / Math.max(frames, 1);
  return avgs;
}

/**
 * The sphere raster's split per instance (`rasterized` in `cell-cull.wgsl.ts`) replayed on the CPU for the merged sphere
 * objects of the scene (`--raster-split-notes`): of their instances in view, those the raster draws and those sent to
 * the impostors, through the near plane or too large at their nearest point; of the latter's spheres in view, those
 * through the near plane, too large, and those the raster could draw one by one. Every sphere at its member's largest
 * radius (`padding`, as the split); no occlusion, no levels (the instances near the camera are at the finest)
 */
function getRasterSplitNotes(canvas3d: Canvas3D, reprs: readonly Representation.Any[], input: RenderInput): string[] {
  const { camera } = canvas3d;
  const { width, height } = input.size ?? Size;
  const P = camera.projection;
  const scale = Math.max(P[0] * width, P[5] * height) / 2;
  const limit = input.props?.renderer?.sphereRaster === 'off' ? 0 : (input.props?.renderer?.sphereRasterRadius ?? 24);
  const near = camera.near;
  const view = camera.view;
  const frustum = Frustum3D.fromProjectionMatrix(Frustum3D(), camera.projectionView);
  // the view-space depth (the distance from the camera plane) of the point at `o` in `a`
  const depth = (a: ArrayLike<number>, o: number, t: Float32Array, to: number) => {
    const x = a[o];
    const y = a[o + 1];
    const z = a[o + 2];
    const wx = t[to] * x + t[to + 4] * y + t[to + 8] * z + t[to + 12];
    const wy = t[to + 1] * x + t[to + 5] * y + t[to + 9] * z + t[to + 13];
    const wz = t[to + 2] * x + t[to + 6] * y + t[to + 10] * z + t[to + 14];
    s.center[0] = wx;
    s.center[1] = wy;
    s.center[2] = wz;
    return -(view[2] * wx + view[6] * wy + view[10] * wz + view[14]);
  };
  const s = Sphere3D();
  const rasterizes = (d: number, r: number) =>
    d - r > near && r * scale <= limit * Math.sqrt(Math.max(d * d - r * r, 0));
  const c = {
    raster: 0,
    rasterSpheres: 0,
    near: 0,
    large: 0,
    spheres: 0,
    spheresNear: 0,
    spheresLarge: 0,
    split: 0,
    area: 0,
  };
  for (const repr of reprs) {
    for (const o of repr.renderObjects) {
      if (!isMergedRenderObject(o) || o.type !== 'spheres') continue;
      for (const m of o.members) {
        const v = m.values as unknown as {
          invariantBoundingSphere: ValueCell<Sphere3D>;
          padding: ValueCell<number>;
          aTransform: ValueCell<Float32Array>;
          instanceCount: ValueCell<number>;
          uVertexCount: ValueCell<number>;
          tPositionGroup: ValueCell<{ array: Float32Array }>;
        };
        const inv = v.invariantBoundingSphere.ref.value;
        const pad = v.padding.ref.value;
        const t = v.aTransform.ref.value;
        const count = v.uVertexCount.ref.value;
        const positions = v.tPositionGroup.ref.value.array;
        for (let i = 0, il = v.instanceCount.ref.value; i < il; ++i) {
          const d = depth(inv.center, 0, t, i * 16) - inv.radius;
          s.radius = inv.radius;
          if (!Frustum3D.intersectsSphere3D(frustum, s)) continue;
          // the instance's nearest point, its largest sphere there (as the split)
          if (rasterizes(d, pad)) {
            c.raster += 1;
            c.rasterSpheres += count;
            continue;
          }
          if (d - pad <= near) c.near += 1;
          else c.large += 1;
          s.radius = pad;
          for (let k = 0; k < count; ++k) {
            const ds = depth(positions, k * 4, t, i * 16);
            if (!Frustum3D.intersectsSphere3D(frustum, s)) continue;
            c.spheres += 1;
            if (ds - pad <= near * 1.001) c.spheresNear += 1;
            else if (!rasterizes(ds, pad)) c.spheresLarge += 1;
            else c.split += 1;
            // the screen area of its disc (the tangents' radius), at most the viewport's
            if (ds - pad > near) {
              const px = (pad * scale) / Math.sqrt(Math.max(ds * ds - pad * pad, 1e-6));
              c.area += Math.min(Math.PI * px * px, width * height);
            }
          }
        }
      }
    }
  }
  const k = (n: number) => (n >= 1e4 ? `${(n / 1e3).toFixed(0)}k` : `${n}`);
  return [
    `raster split (CPU, ${limit} px, no occlusion): instances in view to the raster ${k(c.raster)} (${k(c.rasterSpheres)} spheres), to the impostors ${k(c.near + c.large)} (through near ${k(c.near)}, too large ${k(c.large)}); their spheres in view ${k(c.spheres)}: through near ${k(c.spheresNear)}, too large ${k(c.spheresLarge)}, each rasterizable ${k(c.split)}; their discs ${(c.area / (width * height)).toFixed(1)}× the viewport`,
  ];
}

/** the highlight of frame `i` of a sweep (`RenderInput.sweep`): representations far apart in turn */
function sweepHighlight(canvas3d: Canvas3D, reprs: readonly Representation.Any[], i: number) {
  const at = (k: number) => reprs[(k * 7919) % reprs.length];
  if (i > 0) canvas3d.mark({ loci: { kind: 'every-loci' }, repr: at(i - 1) }, MarkerAction.RemoveHighlight);
  canvas3d.mark({ loci: { kind: 'every-loci' }, repr: at(i) }, MarkerAction.Highlight);
}

/** resolves in a task of its own (a message, not clamped as nested timeouts are) */
function nextTask() {
  return new Promise<void>((resolve) => {
    const channel = new MessageChannel();
    channel.port1.onmessage = () => {
      channel.port1.close();
      resolve();
    };
    channel.port2.postMessage(null);
  });
}

/** the camera's offset from the target when a dolly started, per canvas: frame `i` of a run is at the same place in each */
const dollyStart = new WeakMap<Canvas3D, Vec3>();
function dollyCamera(canvas3d: Canvas3D, step: number, i: number) {
  const { position, target } = canvas3d.camera.getSnapshot();
  let start = dollyStart.get(canvas3d);
  if (!start) dollyStart.set(canvas3d, (start = Vec3.sub(Vec3(), position, target)));
  const offset = Vec3.scale(Vec3(), start, (1 - step) ** i);
  canvas3d.camera.setState({ position: Vec3.add(Vec3(), target, offset) }, 0);
}

export function rotateCamera(canvas3d: Canvas3D, degrees?: number) {
  const { position, target, up } = canvas3d.camera.getSnapshot();
  const rotation = Mat4.fromRotation(Mat4(), degrees === undefined ? RotationStep : (degrees * Math.PI) / 180, up);
  const offset = Vec3.transformMat4(Vec3(), Vec3.sub(Vec3(), position, target), rotation);
  canvas3d.camera.setState({ position: Vec3.add(Vec3(), target, offset) }, 0);
}

/** `count` boxes, each its own representation (a draw call with its own object uniforms) */
async function boxes(count: number) {
  const n = Math.ceil(Math.sqrt(count));
  const reprs: Representation.Any[] = [];
  for (let i = 0; i < count; ++i) {
    const t = Mat4.fromTranslation(Mat4(), Vec3.create((i % n) * 2, Math.floor(i / n) * 2, 0));
    reprs.push(await shapeRepr(`box-${i}`, meshBox(1), { transforms: [t] }));
  }
  return reprs;
}

export function renderCase(
  id: string,
  description: string,
  setup: () => Promise<RenderInput>,
): BenchmarkCase<RenderInput> {
  const path = (graphics: GraphicsContext, input: RenderInput) => renderFrames(id, graphics, input);
  return { id, description, setup, paths: { webgl: path, webgpu: path } };
}

let structure: Promise<Structure> | undefined;
export const getStructure = () => (structure ??= exampleStructure('7qpd.fw2.cif'));
/**
 * 60 copies of 7qpd on a ring of radius 150 Å (\`Structure.instances\`): the copies of a unit share its model and
 * elements, so each unit is drawn instanced, as the units of an assembly
 */
async function assembly60() {
  const s = await getStructure();
  const center = s.boundary.sphere.center;
  const transforms: Mat4[] = [];
  for (let i = 0; i < 60; ++i) {
    const t = Mat4.fromTranslation(Mat4(), Vec3.negate(Vec3(), center));
    const out = Mat4.fromTranslation(Mat4(), Vec3.create(150, 0, 0));
    const turn = Mat4.fromRotation(Mat4(), (i / 60) * 2 * Math.PI, Vec3.create(0, 1, 0));
    transforms.push(Mat4.mul(Mat4(), turn, Mat4.mul(Mat4(), out, t)));
  }
  return Structure.instances(s, transforms);
}
let assembly: Promise<Structure> | undefined;
const getAssembly = () => (assembly ??= assembly60());
let rootAssembly: Promise<Structure> | undefined;
/** `getAssembly` as a root structure (as the plugin's assemblies): themes reading the root's units see the copies */
const getRootAssembly = () => (rootAssembly ??= getAssembly().then((a) => Structure.create([...a.units])));

let globule: Promise<Structure> | undefined;
/** 300k atoms: fragment and vertex work */
export const getGlobule = () => (globule ??= globuleStructure(300_000));

export const RenderCases: BenchmarkCase<any>[] = [
  renderCase(
    'render/boxes-1000',
    '1000 boxes, each its own representation: draw calls and per-object uniforms',
    async () => ({
      scene: () => boxes(1000),
    }),
  ),
  renderCase('render/boxes-100', '100 boxes, each its own representation', async () => ({
    scene: () => boxes(100),
  })),
  renderCase(
    'render/units-1000-spacefill',
    '1000 chains as spacefill: more than 100 units, merged into one render item (renderMerged auto)',
    async () => {
      const s = await chainGrid(1000);
      return { scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, {})]) };
    },
  ),
  renderCase(
    'render/units-1000-spacefill-unmerged',
    '1000 chains as spacefill, a render item per unit (renderMerged off), sharing the material',
    async () => {
      const s = await chainGrid(1000);
      return {
        scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, { renderMerged: false })]),
      };
    },
  ),
  renderCase(
    'render/viewer-7qpd-rotating',
    'PDB entry 7qpd as cartoon, ball-and-stick, molecular and Gaussian surface (a Viewer-like scene), rotating',
    async () => {
      const s = await getStructure();
      return {
        scene: (g) =>
          Promise.all([
            structureRepr(g, CartoonRepresentationProvider, s, {}),
            structureRepr(g, BallAndStickRepresentationProvider, s, {}),
            structureRepr(g, MolecularSurfaceRepresentationProvider, s, { alpha: 0.5 } as any),
            structureRepr(g, GaussianSurfaceRepresentationProvider, s, {}),
          ]),
        rotate: true,
      };
    },
  ),
  renderCase('render/assembly-60-cartoon-rotating', '60 instanced copies of 7qpd as cartoon, rotating', async () => {
    const s = await getAssembly();
    return { scene: (g) => Promise.all([structureRepr(g, CartoonRepresentationProvider, s, {})]), rotate: true };
  }),
  renderCase(
    'render/assembly-60-spacefill-rotating',
    '60 instanced copies of 7qpd as spacefill, rotating',
    async () => {
      const s = await getAssembly();
      return { scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, {})]), rotate: true };
    },
  ),
  renderCase(
    'render/assembly-60-spacefill-ssao',
    '60 instanced copies of 7qpd as spacefill with SSAO, 2560x1440',
    async () => {
      const s = await getAssembly();
      return {
        size: LargeSize,
        scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, {})]),
        props: {
          postprocessing: {
            ...PD.getDefaultValues(PostprocessingParams),
            occlusion: { name: 'on', params: PD.getDefaultValues(SsaoParams) },
          },
        },
      };
    },
  ),
  renderCase(
    'render/assembly-60-spacefill-ssao-multiscale',
    '60 instanced copies of 7qpd as spacefill with multi-scale SSAO (four levels), 2560x1440',
    async () => {
      const s = await getAssembly();
      const ssao = PD.getDefaultValues(SsaoParams);
      return {
        size: LargeSize,
        scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, {})]),
        props: {
          postprocessing: {
            ...PD.getDefaultValues(PostprocessingParams),
            occlusion: {
              name: 'on',
              params: { ...ssao, multiScale: { name: 'on', params: SsaoParams.multiScale.map('on').defaultValue } },
            },
          },
        },
      };
    },
  ),
  renderCase('render/7qpd-cartoon', 'PDB entry 7qpd as cartoon (a mesh per chain)', async () => {
    const s = await getStructure();
    return { scene: (g) => Promise.all([structureRepr(g, CartoonRepresentationProvider, s, {})]) };
  }),
  renderCase(
    'render/7qpd-illumination-rotating',
    'PDB entry 7qpd as cartoon with illumination, rotating: the first tracing iteration per frame (fixed quality)',
    async () => {
      const s = await getStructure();
      return {
        scene: (g) => Promise.all([structureRepr(g, CartoonRepresentationProvider, s, {})]),
        props: {
          illumination: { ...PD.getDefaultValues(IlluminationParams), enabled: true, rendersPerFrame: [4, 4] },
        },
        rotate: true,
      };
    },
  ),
  renderCase('render/7qpd-ball-and-stick', 'PDB entry 7qpd as ball-and-stick (impostors)', async () => {
    const s = await getStructure();
    return { scene: (g) => Promise.all([structureRepr(g, BallAndStickRepresentationProvider, s, {})]) };
  }),
  renderCase('render/7qpd-spacefill', 'PDB entry 7qpd as spacefill (sphere impostors)', async () => {
    const s = await getStructure();
    return { scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, {})]) };
  }),
  renderCase('render/7qpd-gaussian-surface', 'PDB entry 7qpd as Gaussian surface (GPU texture meshes)', async () => {
    const s = await getStructure();
    return { scene: (g) => Promise.all([structureRepr(g, GaussianSurfaceRepresentationProvider, s, {})]) };
  }),
  renderCase('render/globule-300k-spacefill', '300k atoms as spacefill (sphere impostors)', async () => {
    const s = await getGlobule();
    return { scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, {})]) };
  }),
  renderCase(
    'render/globule-300k-spacefill-large',
    '300k atoms as spacefill (sphere impostors), 2560x1440: fragment-bound',
    async () => {
      const s = await getGlobule();
      return { size: LargeSize, scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, {})]) };
    },
  ),
  renderCase('render/globule-300k-spacefill-mesh', '300k atoms as spacefill meshes (detail 1)', async () => {
    const s = await getGlobule();
    return {
      scene: (g) =>
        Promise.all([
          structureRepr(g, SpacefillRepresentationProvider, s, { tryUseImpostor: false, quality: 'custom', detail: 1 }),
        ]),
    };
  }),
  renderCase('render/7qpd-cartoon-postprocessing', 'PDB entry 7qpd as cartoon with SSAO and outlines', async () => {
    const s = await getStructure();
    return {
      scene: (g) => Promise.all([structureRepr(g, CartoonRepresentationProvider, s, {})]),
      props: {
        postprocessing: {
          ...PD.getDefaultValues(PostprocessingParams),
          occlusion: { name: 'on', params: PD.getDefaultValues(SsaoParams) },
          outline: { name: 'on', params: PD.getDefaultValues(OutlineParams) },
        },
      },
    };
  }),
];

/**
 * Measurements for render bundles (WebGPU), not in the default runs (`--filter 'render-bundles/*'`): the recording cost
 * of the opaque pass against the direct cases (`render/boxes-1000`, `render/units-1000-spacefill-unmerged`), how often
 * the visible set changes while rotating or zoomed in, and the GPU time of drawing what culling skips (`--gpu-timing`)
 */
/** The renderer replaying the opaque pass from draw lists (`RendererProps.drawLists`) */
const DrawLists: Partial<Canvas3DProps> = {
  renderer: { ...PD.getDefaultValues(Canvas3DParams).renderer, drawLists: true },
};

export const RenderBundleCases: BenchmarkCase<any>[] = [
  renderCase(
    'render-bundles/units-1356-gaussian-highlighted-rotating',
    '1356 units as Gaussian surfaces, highlighted (passes per frame: opaque, marking depth and mask), direct, rotating',
    async () => {
      const s = await chainGrid(1356);
      return {
        scene: (g) =>
          Promise.all([structureRepr(g, GaussianSurfaceRepresentationProvider, s, { renderMerged: false })]),
        highlight: true,
        rotate: true,
      };
    },
  ),
  renderCase(
    'render-bundles/units-1356-gaussian-highlighted-drawlists-rotating',
    '1356 units as Gaussian surfaces, highlighted (passes per frame: opaque, marking depth and mask), with draw lists, rotating',
    async () => {
      const s = await chainGrid(1356);
      return {
        scene: (g) =>
          Promise.all([structureRepr(g, GaussianSurfaceRepresentationProvider, s, { renderMerged: false })]),
        highlight: true,
        props: DrawLists,
        rotate: true,
      };
    },
  ),
  renderCase(
    'render-bundles/units-1356-gaussian',
    '1356 units as Gaussian surfaces (a render item each, as PDB entry 3j3q), direct',
    async () => {
      const s = await chainGrid(1356);
      return {
        scene: (g) =>
          Promise.all([structureRepr(g, GaussianSurfaceRepresentationProvider, s, { renderMerged: false })]),
      };
    },
  ),
  renderCase(
    'render-bundles/units-1356-gaussian-drawlists',
    '1356 units as Gaussian surfaces (a render item each, as PDB entry 3j3q), with draw lists',
    async () => {
      const s = await chainGrid(1356);
      return {
        scene: (g) =>
          Promise.all([structureRepr(g, GaussianSurfaceRepresentationProvider, s, { renderMerged: false })]),
        props: DrawLists,
      };
    },
  ),
  renderCase(
    'render-bundles/units-1356-gaussian-rotating',
    '1356 units as Gaussian surfaces (a render item each, as PDB entry 3j3q), direct, rotating',
    async () => {
      const s = await chainGrid(1356);
      return {
        scene: (g) =>
          Promise.all([structureRepr(g, GaussianSurfaceRepresentationProvider, s, { renderMerged: false })]),
        rotate: true,
      };
    },
  ),
  renderCase(
    'render-bundles/units-1356-gaussian-drawlists-rotating',
    '1356 units as Gaussian surfaces (a render item each, as PDB entry 3j3q), with draw lists, rotating',
    async () => {
      const s = await chainGrid(1356);
      return {
        scene: (g) =>
          Promise.all([structureRepr(g, GaussianSurfaceRepresentationProvider, s, { renderMerged: false })]),
        props: DrawLists,
        rotate: true,
      };
    },
  ),
  renderCase('render-bundles/boxes-1000-drawlists', '`render/boxes-1000` with draw lists', async () => ({
    scene: () => boxes(1000),
    props: DrawLists,
  })),
  renderCase(
    'render-bundles/boxes-1000-drawlists-rotating',
    '`render/boxes-1000` with draw lists, rotating',
    async () => ({
      scene: () => boxes(1000),
      rotate: true,
      props: DrawLists,
    }),
  ),
  renderCase(
    'render-bundles/boxes-1000-zoom0.3-drawlists-rotating',
    '`render/boxes-1000` zoomed to 0.3 of the distance with draw lists, rotating (the visible set changes)',
    async () => ({ scene: () => boxes(1000), zoom: 0.3, rotate: true, props: DrawLists }),
  ),
  renderCase(
    'render-bundles/units-1000-unmerged-drawlists',
    '`render/units-1000-spacefill-unmerged` with draw lists',
    async () => {
      const s = await chainGrid(1000);
      return {
        scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, { renderMerged: false })]),
        props: DrawLists,
      };
    },
  ),
  renderCase(
    'render-bundles/boxes-1000-record',
    '`render/boxes-1000`, the opaque pass recorded every frame',
    async () => ({
      scene: () => boxes(1000),
      bundles: { record: true },
    }),
  ),
  renderCase(
    'render-bundles/units-1000-unmerged-record',
    '`render/units-1000-spacefill-unmerged`, the opaque pass recorded every frame',
    async () => {
      const s = await chainGrid(1000);
      return {
        scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, { renderMerged: false })]),
        bundles: { record: true },
      };
    },
  ),
  renderCase(
    'render-bundles/boxes-1000-rotating',
    '`render/boxes-1000` rotating: the visible set per frame',
    async () => ({
      scene: () => boxes(1000),
      rotate: true,
      bundles: { visibility: true },
    }),
  ),
  ...[0.3, 0.1].flatMap((zoom) => [
    renderCase(
      `render-bundles/boxes-1000-zoom${zoom}-rotating`,
      `\`render/boxes-1000\` zoomed to ${zoom} of the distance, rotating: the visible set per frame`,
      async () => ({ scene: () => boxes(1000), zoom, rotate: true, bundles: { visibility: true } }),
    ),
    renderCase(
      `render-bundles/boxes-1000-zoom${zoom}-unculled`,
      `\`render/boxes-1000\` zoomed to ${zoom} of the distance, without object culling`,
      async () => ({ scene: () => boxes(1000), zoom, bundles: { noCulling: true, visibility: true } }),
    ),
    renderCase(
      `render-bundles/boxes-1000-zoom${zoom}-culled`,
      `\`render/boxes-1000\` zoomed to ${zoom} of the distance, with object culling`,
      async () => ({ scene: () => boxes(1000), zoom, bundles: { visibility: true } }),
    ),
  ]),
];

/**
 * Many small objects (`draw-table/*`, not in the default runs): 10k boxes, each its own representation, still and
 * rotating, direct and with draw lists, one highlighted (the marking passes); the per-object CPU work of the
 * renderer at the scale of user-generated scenes, and what setting them up costs
 */
export const DrawTableCases: BenchmarkCase<any>[] = (() => {
  const count = 10_000;
  const scene = () => boxes(count);
  const variants: { suffix: string; label: string; input: Partial<RenderInput> }[] = [
    { suffix: '', label: 'direct', input: {} },
    { suffix: '-rotating', label: 'direct, rotating', input: { rotate: true } },
    { suffix: '-drawlists', label: 'with draw lists', input: { props: DrawLists } },
    { suffix: '-drawlists-rotating', label: 'with draw lists, rotating', input: { props: DrawLists, rotate: true } },
    {
      suffix: '-highlighted-rotating',
      label: 'all highlighted, direct, rotating',
      input: { highlight: true, rotate: true },
    },
    {
      suffix: '-highlighted-drawlists-rotating',
      label: 'all highlighted, with draw lists, rotating',
      input: { highlight: true, props: DrawLists, rotate: true },
    },
    { suffix: '-sweep', label: 'a hover sweeping across the boxes, direct', input: { sweep: true } },
    {
      suffix: '-drawlists-sweep',
      label: 'a hover sweeping across the boxes, with draw lists',
      input: { sweep: true, props: DrawLists },
    },
  ];
  return variants.map(({ suffix, label, input }) =>
    renderCase(`draw-table/boxes-10k${suffix}`, `${count} boxes, each its own representation, ${label}`, async () => ({
      scene,
      setupNotes: true,
      ...input,
    })),
  );
})();

/** A unit of `structure` (by index, wrapping) as a loci of all its elements */
function unitLoci(structure: Structure, index: number) {
  const unit = structure.units[index % structure.units.length];
  return StructureElement.Loci(structure, [{ unit, indices: Interval.ofBounds(0, unit.elements.length) }]);
}

const MergingKinds: { name: string; provider: RepresentationProvider<Structure, any> }[] = [
  { name: 'spacefill', provider: SpacefillRepresentationProvider },
  { name: 'cartoon', provider: CartoonRepresentationProvider },
];
const MergingCounts = [10, 30, 100, 300, 1000];
const MergingModes = [
  { suffix: '', label: 'still' },
  { suffix: '-hover', label: 'a unit highlighted per frame (the one before not)' },
  { suffix: '-recolor', label: "a uniform color per frame, through the representation's update" },
  { suffix: '-rotating', label: 'the camera rotating (culled per frame)' },
] as const;

/**
 * Merging (`merged-render-objects.md`, on demand): `count` chains of 10 atoms as the units of one representation,
 * merged into one render object or not (`renderMerged` on or off): frames still, with a hover per frame, or with a color
 * theme change per frame (`update`); with the setup times and the GPU bytes
 */
export const MergingCases: BenchmarkCase<RenderInput>[] = MergingCounts.flatMap((count) =>
  MergingKinds.flatMap((kind) =>
    [true, false].flatMap((merged) =>
      MergingModes.map((mode) =>
        renderCase(
          `merging/units-${count}-${kind.name}-${merged ? 'merged' : 'unmerged'}${mode.suffix}`,
          `${count} chains as ${kind.name}, ${merged ? 'merged' : 'a render object per unit'}: ${mode.label}`,
          async () => {
            const s = await chainGrid(count);
            const props = { renderMerged: merged };
            let update: RenderInput['update'];
            if (mode.suffix === '-hover') {
              update = (canvas3d, _reprs, i) => {
                if (i > 0) canvas3d.mark({ loci: unitLoci(s, (i - 1) * 7919) }, MarkerAction.RemoveHighlight);
                canvas3d.mark({ loci: unitLoci(s, i * 7919) }, MarkerAction.Highlight);
              };
            } else if (mode.suffix === '-recolor') {
              update = async (canvas3d, reprs, i) => {
                const repr = reprs[0];
                const { provider } = kind;
                repr.setTheme(
                  Theme.create(
                    getThemeRegistries(),
                    { structure: s },
                    {
                      colorTheme: { name: 'uniform', params: { value: Color(i % 2 ? 0xdd4422 : 0x2244dd) } },
                      sizeTheme: {
                        name: provider.defaultSizeTheme.name,
                        params: provider.defaultSizeTheme.props ?? {},
                      },
                    },
                  ),
                );
                // committed by the frame's tick (synchronous) and drawn once there: a commit here draws or not by what
                // earlier work left (a forced draw after objects were added), the frames then measure one draw or two
                await repr.createOrUpdate({}, s).run();
              };
            }
            return {
              scene: (g) => Promise.all([structureRepr(g, kind.provider, s, props)]),
              setupNotes: true,
              byteNotes: true,
              update,
              rotate: mode.suffix === '-rotating',
            };
          },
        ),
      ),
    ),
  ),
);

/**
 * Large meshes merged or not: the 60 copies of 7qpd as a root structure (420 units of about 1.6 k atoms) as Gaussian
 * surfaces, many vertices per unit (on WebGL2 without base instances, merged objects read their instance data per
 * vertex from a texture, `dInstanceData`)
 */
export const MergingMeshCases: BenchmarkCase<RenderInput>[] = [true, false].flatMap((merged) =>
  [undefined, { width: 320, height: 180 }, LargeSize].map((size) =>
    renderCase(
      `merging/assembly-gaussian-${merged ? 'merged' : 'unmerged'}${size ? `-${size.width}` : ''}`,
      `60 copies of 7qpd as a root structure (420 units) as Gaussian surfaces, ${merged ? 'merged' : 'a render object per unit'}`,
      async () => {
        const s = await getRootAssembly();
        return {
          size,
          scene: (g) =>
            Promise.all([structureRepr(g, GaussianSurfaceRepresentationProvider, s, { renderMerged: merged })]),
        };
      },
    ),
  ),
);

/** The folder of the large models some cases load (from the repository root; not checked in, `.gitignore`d) */
export const BenchmarkDataDir = '.benchmark';

/** A file of `BenchmarkDataDir`, from the benchmark's server */
export async function fetchBenchmarkFile(file: string) {
  const response = await fetch(`/${BenchmarkDataDir}/${file}`);
  if (!response.ok) throw new Error(`${file} not found in ${BenchmarkDataDir}/ (${response.status})`);
  return new Uint8Array(await response.arrayBuffer());
}

const benchmarkEntries = new Map<string, Promise<Structure>>();
/** A binary CIF file of `BenchmarkDataDir`, its first model */
export const getBenchmarkEntry = (file: string) => {
  let entry = benchmarkEntries.get(file);
  if (!entry) {
    entry = (async () => {
      const parsed = await CIF.parseBinary(await fetchBenchmarkFile(file)).run();
      if (parsed.isError) throw new Error(parsed.message);
      const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0]).run();
      return Structure.ofModel(trajectory.representative);
    })();
    benchmarkEntries.set(file, entry);
  }
  return entry;
};

/**
 * Unique meshes merged or not, no two units sharing a symmetry group (unmerged, a render object per unit): 4V6X, a
 * ribosome (about 89 chains, `4v6x.bcif`), and 3j3q, the HIV-1 capsid (1356 chains, `3j3q.bcif`), as Gaussian surfaces;
 * CPU meshes on every backend (`tryUseGpu: false`: WebGPU builds compute surfaces, texture meshes, never merged); the
 * files go in `.benchmark/` (from models.rcsb.org)
 */
export const MergingUniqueMeshCases: BenchmarkCase<RenderInput>[] = [
  ['4v6x', 'a ribosome, about 89 chains'],
  ['3j3q', 'the HIV-1 capsid, 1356 chains'],
].flatMap(([id, about]) =>
  [true, false].flatMap((merged) =>
    [undefined, { width: 320, height: 180 }, LargeSize].map((size) =>
      renderCase(
        `merging/${id}-gaussian-${merged ? 'merged' : 'unmerged'}${size ? `-${size.width}` : ''}`,
        `${id} (${about}, unique geometry each) as Gaussian surfaces, ${merged ? 'merged' : 'a render object per unit'}`,
        async () => {
          const s = await getBenchmarkEntry(`${id}.bcif`);
          return {
            size,
            setupNotes: true,
            byteNotes: true,
            scene: (g) =>
              Promise.all([
                structureRepr(g, GaussianSurfaceRepresentationProvider, s, { renderMerged: merged, tryUseGpu: false }),
              ]),
          };
        },
      ),
    ),
  ),
);

/**
 * The memory of large scenes, with `--byte-notes` and `--heap-snapshot` (the files in `.benchmark/`): 3j3q, the
 * HIV-1 capsid, as Gaussian surfaces with the defaults (texture meshes, the small units merged CPU meshes) and as
 * cartoon at medium quality
 */
export const MergingMemoryCases: BenchmarkCase<RenderInput>[] = [
  renderCase('merging/3j3q-gaussian-defaults', '3j3q as Gaussian surfaces, the defaults', async () => {
    const s = await getBenchmarkEntry('3j3q.bcif');
    return {
      setupNotes: true,
      byteNotes: true,
      scene: (g) => Promise.all([structureRepr(g, GaussianSurfaceRepresentationProvider, s, {})]),
    };
  }),
  renderCase('merging/3j3q-cartoon-medium', '3j3q as cartoon, medium quality', async () => {
    const s = await getBenchmarkEntry('3j3q.bcif');
    return {
      setupNotes: true,
      byteNotes: true,
      scene: (g) => Promise.all([structureRepr(g, CartoonRepresentationProvider, s, { quality: 'medium' })]),
    };
  }),
];

/**
 * 4V6X's Gaussian surfaces (CPU meshes) with the near plane at half the scene's radius (`cameraRadius`): the
 * solid-interior caps of the near cut (on by default for Gaussian surfaces), merged or a render object per unit
 */
export const MergingNearCapCases: BenchmarkCase<RenderInput>[] = [true, false].map((merged) =>
  renderCase(
    `merging/4v6x-gaussian-${merged ? 'merged' : 'unmerged'}-near`,
    `4v6x as Gaussian surfaces, ${merged ? 'merged' : 'a render object per unit'}, the near plane at half the scene's radius (caps)`,
    async () => {
      const s = await getBenchmarkEntry('4v6x.bcif');
      return {
        cameraRadius: 0.5,
        scene: (g) =>
          Promise.all([
            structureRepr(g, GaussianSurfaceRepresentationProvider, s, { renderMerged: merged, tryUseGpu: false }),
          ]),
      };
    },
  ),
);

/**
 * 4V6X as cartoon and ball-and-stick, merged or a render object per unit: the setup notes (add and commit, the first
 * frame with its programs compiled); run alone with the browser's shader caches off (`--browser-arg
 * --disable-gpu-program-cache --browser-arg --disable-gpu-shader-disk-cache`) for the compile cost
 */
const SetupProviders: [string, RepresentationProvider<any, any, any>][] = [
  ['cartoon', CartoonRepresentationProvider],
  ['ball-and-stick', BallAndStickRepresentationProvider],
];
export const MergingSetupCases: BenchmarkCase<RenderInput>[] = SetupProviders.flatMap(([name, provider]) =>
  [true, false].map((merged) =>
    renderCase(
      `merging/4v6x-${name}-${merged ? 'merged' : 'unmerged'}`,
      `4v6x as ${name}, ${merged ? 'merged' : 'a render object per unit'}: the setup and first frame`,
      async () => {
        const s = await getBenchmarkEntry('4v6x.bcif');
        return {
          setupNotes: true,
          scene: (g) => Promise.all([structureRepr(g, provider, s, { renderMerged: merged })]),
        };
      },
    ),
  ),
);

/**
 * Cylinders: the bonds of 4V6X's ball-and-stick (`4v6x.bcif` in `.benchmark/`), merged or not, with their CPU and
 * GPU bytes: the baseline of cylinders as per-item data
 */
export const CylinderCases: BenchmarkCase<RenderInput>[] = [true, false].map((merged) =>
  renderCase(
    `cylinders/4v6x-bonds-${merged ? 'merged' : 'unmerged'}`,
    `4v6x ball-and-stick bonds only (cylinder impostors), ${merged ? 'merged' : 'a render object per unit'}`,
    async () => {
      const s = await getBenchmarkEntry('4v6x.bcif');
      return {
        setupNotes: true,
        byteNotes: true,
        scene: (g) =>
          Promise.all([
            structureRepr(g, BallAndStickRepresentationProvider, s, {
              visuals: ['intra-bond', 'inter-bond'],
              renderMerged: merged,
            }),
          ]),
      };
    },
  ),
);

/** Lines: 4V6X as the line representation (bonds as lines), merged or not, with their CPU and GPU bytes */
export const LineCases: BenchmarkCase<RenderInput>[] = [true, false].map((merged) =>
  renderCase(
    `lines/4v6x-${merged ? 'merged' : 'unmerged'}`,
    `4v6x as lines (bonds), ${merged ? 'merged' : 'a render object per unit'}`,
    async () => {
      const s = await getBenchmarkEntry('4v6x.bcif');
      return {
        setupNotes: true,
        byteNotes: true,
        scene: (g) => Promise.all([structureRepr(g, LineRepresentationProvider, s, { renderMerged: merged })]),
      };
    },
  ),
);

/** `units` chains of 10 atoms (`chainGrid`), `copies` instances of each on a grid of copies: a root structure */
async function instancedChainGrid(units: number, copies: number) {
  const s = await chainGrid(units);
  const extent = s.boundary.sphere.radius * 2.5;
  const n = Math.ceil(Math.cbrt(copies));
  const transforms: Mat4[] = [];
  for (let i = 0; i < copies; ++i) {
    const offset = Vec3.create((i % n) * extent, (Math.floor(i / n) % n) * extent, Math.floor(i / (n * n)) * extent);
    transforms.push(Mat4.fromTranslation(Mat4(), offset));
  }
  return Structure.create([...Structure.instances(s, transforms).units]);
}

/**
 * Merging with many instances (`merged-render-objects.md` step 6, on demand): 100 units of `copies` instances each as
 * spacefill, the camera rotating: the CPU culls every unit's instances per frame (instance grids); 300 and 500 copies
 * between, where culling on the GPU starts to pay (`AutoGpuCullingMinInstances`)
 */
export const MergingInstanceCases: BenchmarkCase<RenderInput>[] = [100, 300, 500, 1000].flatMap((copies) =>
  (['merged', 'merged-gpu', 'unmerged'] as const).flatMap((kind) =>
    [true, false].map((rotate) =>
      renderCase(
        `merging/instances-${copies}-spacefill-${kind}${rotate ? '-rotating' : ''}`,
        `100 chains of ${copies} copies each as spacefill, ${kind === 'unmerged' ? 'a render object per unit' : kind === 'merged' ? 'merged' : 'merged, culled on the GPU'}${rotate ? ', rotating' : ''}`,
        async () => {
          const s = await instancedChainGrid(100, copies);
          const props = { renderMerged: kind !== 'unmerged', gpuCulling: kind === 'merged-gpu' };
          return {
            scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, props)]),
            setupNotes: true,
            byteNotes: true,
            rotate,
          };
        },
      ),
    ),
  ),
);

/**
 * As `MergingInstanceCases` as cartoon, merged and culled on the CPU or the GPU: the instance threshold of
 * `gpuCulling: 'auto'` (`AutoGpuCullingMinInstances`) without the sphere raster's (spacefill takes merged objects with
 * many spheres to the GPU anyway, `AutoGpuCullingMinRasterSpheres`)
 */
export const MergingInstanceCartoonCases: BenchmarkCase<RenderInput>[] = [100, 300, 500, 1000].flatMap((copies) =>
  (['merged', 'merged-gpu'] as const).flatMap((kind) =>
    [true, false].map((rotate) =>
      renderCase(
        `merging/instances-${copies}-cartoon-${kind}${rotate ? '-rotating' : ''}`,
        `100 chains of ${copies} copies each as cartoon, merged${kind === 'merged-gpu' ? ', culled on the GPU' : ''}${rotate ? ', rotating' : ''}`,
        async () => {
          const s = await instancedChainGrid(100, copies);
          const props = { renderMerged: true, gpuCulling: kind === 'merged-gpu' };
          return {
            scene: (g) => Promise.all([structureRepr(g, CartoonRepresentationProvider, s, props)]),
            rotate,
          };
        },
      ),
    ),
  ),
);

/**
 * The merged instances (100 chains of 1000 copies, rotating) with what reads their segment's sphere: tumble (the
 * animation's noise, by the instance's bounding sphere) and instance clipping (an inverted sphere around everything:
 * the test per instance, cutting nothing); merged, and culled on the GPU; on demand
 */
export const MergingSegmentDataCases: BenchmarkCase<RenderInput>[] = (['merged', 'merged-gpu'] as const).flatMap(
  (kind) =>
    (['tumble', 'clip'] as const).map((use) =>
      renderCase(
        `merging/instances-1000-spacefill-${kind}-${use}-rotating`,
        `100 chains of 1000 copies each as spacefill, ${kind === 'merged' ? 'merged' : 'merged, culled on the GPU'}, ${use === 'tumble' ? 'tumbling' : 'clipped per instance (cutting nothing)'}, rotating`,
        async () => {
          const s = await instancedChainGrid(100, 1000);
          const props = {
            renderMerged: true,
            gpuCulling: kind === 'merged-gpu',
            ...(use === 'tumble'
              ? {
                  // every field: a partial group replaces the defaults
                  animation: {
                    wiggleMode: 'position' as const,
                    wiggleSpeed: 7,
                    wiggleAmplitude: 0,
                    wiggleFrequency: 0.2,
                    tumbleSpeed: 1,
                    tumbleAmplitude: 10,
                    tumbleFrequency: 0.2,
                  },
                }
              : {
                  clip: {
                    variant: 'instance' as const,
                    objects: [
                      {
                        type: 'sphere' as const,
                        invert: true,
                        position: Vec3(),
                        rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
                        scale: Vec3.create(1e5, 1e5, 1e5),
                        transform: Mat4.identity(),
                      },
                    ],
                  },
                }),
          };
          return {
            scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, props)]),
            setupNotes: true,
            rotate: true,
          };
        },
      ),
    ),
);

/**
 * One unit of `copies` instances as spacefill, rotating, culled on the CPU or the GPU: the fixed cost of a GPU-culled
 * object against the CPU's per instance (`gpuCulling: 'auto'`, `AutoGpuCullingMinInstances`)
 */
export const MergingSingleInstanceCases: BenchmarkCase<RenderInput>[] = [100, 1000, 10000].flatMap((copies) =>
  ([false, true] as const).map((gpu) =>
    renderCase(
      `merging/instances-1x${copies}-spacefill-${gpu ? 'gpu' : 'cpu'}-rotating`,
      `one chain of ${copies} copies as spacefill, culled on the ${gpu ? 'GPU' : 'CPU'}, rotating`,
      async () => {
        const s = await instancedChainGrid(1, copies);
        return {
          scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, { gpuCulling: gpu })]),
          rotate: true,
        };
      },
    ),
  ),
);

let postsynModels: Promise<Structure[]> | undefined;
/**
 * `postsynstubby_2-no_bonds.bcif` from `.benchmark/` as the mesoscale explorer builds it (`StructureFromPetworld`): a
 * structure per model, the model's units merged into one, instanced by the operators of the model's group in the first
 * assembly (the explorer splits a model's membrane entities into a structure of their own, left out here)
 */
const getPostsynModels = () =>
  (postsynModels ??= (async () => {
    const parsed = await CIF.parseBinary(await fetchBenchmarkFile('postsynstubby_2-no_bonds.bcif')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0]).run();
    const first = trajectory.representative;
    if (!MmcifFormat.is(first.sourceData)) throw new Error('not mmCIF');
    const { db } = first.sourceData.data;
    const gen = db.pdbx_struct_assembly_gen;
    const assemblyId = db.pdbx_struct_assembly.id.value(0);
    const generators = [];
    for (let i = 0; i < gen._rowCount; ++i) {
      if (gen.assembly_id.value(i) !== assemblyId) continue;
      generators.push({ assemblyId, expression: gen.oper_expression.value(i), asymIds: gen.asym_id_list.value(i) });
    }
    const groups = operatorGroupsProvider(generators, getMatrices(db.pdbx_struct_oper_list))();
    const structures: Structure[] = [];
    for (let i = 0; i < trajectory.frameCount; ++i) {
      const model = await Task.resolveInContext(trajectory.getFrameAtIndex(i));
      const unit = mergeUnitsWithSameOperator(Structure.ofModel(model).units, 0);
      const builder = Structure.Builder();
      for (const oper of groups[i].operators) {
        builder.addUnit(
          unit.kind,
          unit.model,
          oper,
          unit.elements,
          unit.traits | Unit.Trait.FastBoundary,
          unit.invariantId,
        );
      }
      structures.push(builder.getStructure());
    }
    return structures;
  })());

/**
 * PetWorld as the mesoscale explorer draws it: a spacefill representation per model (`getPostsynModels`, an object
 * each, merged and culled on the GPU), the explorer's `quality` parameters (levels by size on screen, instance clipping,
 * no light), a uniform color each, Hi-Z on; whole and zoomed to 0.3, still and rotating. `auto`: as on a device with
 * `atomicMinMax64` (`gpuCulling: 'auto'` marks objects of many spheres for the raster, the sphere raster `auto`: the
 * `sphereRasterObjects` largest culled on the GPU and drawn by it); `auto-objectsN`: `auto` with N of them (`all`: 256);
 * `off`: the same objects, the raster off (those for the raster culled on the CPU); `noraster`: as on a device without
 * (the raster off, culled on the GPU what `auto` culls there: 50k instances and more); on demand
 */
export const MesoscaleSplitCases: BenchmarkCase<RenderInput>[] = [1, 0.3].flatMap((zoom) =>
  [false, true].flatMap((rotate) =>
    (['noraster', 'off', 'auto', 'auto-objects0', 'auto-objects4', 'auto-objects16', 'auto-objectsall'] as const).map(
      (variant) =>
        renderCase(
          `mesoscale-split/postsyn-${zoom === 1 ? 'whole' : `zoom${zoom}`}${rotate ? '-rotating' : ''}-${variant}`,
          `postsynstubby_2 as the mesoscale explorer draws it (a representation per model), ${zoom === 1 ? 'whole' : `zoomed to ${zoom} of the distance`}${rotate ? ', rotating' : ''}; ${variant === 'noraster' ? 'as without atomicMinMax64' : `sphere raster ${variant}`}`,
          async () => {
            const structures = await getPostsynModels();
            const lod = Spheres.LodPixelPresets.quality;
            const props = {
              ...SpacefillRepresentationProvider.defaultValues,
              ignoreHydrogens: true,
              instanceGranularity: true,
              ignoreLight: true,
              lodLevels: lod.lodLevels,
              lodPixelSize: lod.lodPixelSize,
              quality: 'lowest',
              clip: { variant: 'instance', objects: [] },
              clipPrimitive: true,
              approximate: false,
              alphaThickness: 12,
            };
            return {
              rotate,
              zoom: zoom === 1 ? undefined : zoom,
              yieldFrames: true,
              setupNotes: true,
              scene: (g) => {
                const colors = getDistinctBaseColors(structures.length, 0);
                return Promise.all(
                  structures.map((s, i) => {
                    // without the raster `auto` culls on the GPU from many instances only (a unit an instance here)
                    const gpuCulling = variant === 'noraster' ? s.units.length >= AutoGpuCullingMinInstances : 'auto';
                    return structureRepr(g, SpacefillRepresentationProvider, s, { ...props, gpuCulling } as never, {
                      name: 'uniform',
                      params: { value: colors[i] },
                    });
                  }),
                );
              },
              props: {
                hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true },
                renderer: {
                  ...PD.getDefaultValues(Canvas3DParams).renderer,
                  sphereRaster: variant.startsWith('auto') ? 'auto' : 'off',
                  ...(variant.startsWith('auto-objects') && {
                    sphereRasterObjects: variant === 'auto-objectsall' ? 256 : Number(variant.slice(12)),
                  }),
                },
              },
            };
          },
        ),
    ),
  ),
);

let postsyn: Promise<ParticleList> | undefined;
/**
 * `postsynstubby_2-no_bonds.bcif` from `.benchmark/` (not checked in: a PetWorld model, 3.68M atoms, 517k
 * operators) as the particle list the viewer loads (`mmcif_particles`): the molecule types the targets, instanced at
 * their particles
 */
export const getPostsyn = () =>
  (postsyn ??= (async () => {
    const parsed = await CIF.parseBinary(await fetchBenchmarkFile('postsynstubby_2-no_bonds.bcif')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const cif = parsed.result;
    const assemblyId = getAssemblyIdsFromMmcif(cif)[0];
    return Task.create('particles', (ctx) =>
      createParticleListFromMmcifAssembly(ctx, cif, { assemblyId, variant: 'petworld' }),
    ).run();
  })());

/** The particle targets as the viewer shows them, with the mesoscale preset's sphere settings */
type LodLevels = (typeof Spheres.LodLevelsPresets)['quality'];

export async function postsynRepr(
  graphics: GraphicsContext,
  list: ParticleList,
  gpuCulling: boolean,
  cellSize?: number,
  lod?: keyof typeof Spheres.LodLevelsPresets | LodLevels,
  structure?: { readonly ignoreLight?: boolean; readonly approximate?: boolean },
  lodPixelSize?: number,
) {
  const registries = getThemeRegistries();
  const provider = ParticleTargetRepresentationProvider;
  const repr = provider.factory({ graphics, ...registries }, provider.getParams);
  repr.setTheme(
    Theme.create(
      registries,
      { particles: list },
      { colorTheme: { name: 'particle-entity', params: {} }, sizeTheme: { name: 'uniform', params: {} } },
    ),
  );
  const defaults = PD.getDefaultValues(ParticleTargetRepresentationParams);
  const props = {
    ...defaults,
    gpuCulling,
    ...(cellSize ? { cellSize } : {}),
    // the distances of the `quality` preset unless asked otherwise (the particles' default is by size on screen), as
    // the cases measured them first
    lodLevels: typeof lod === 'string' ? Spheres.LodLevelsPresets[lod] : (lod ?? Spheres.LodLevelsPresets.quality),
    lodPixelSize: lodPixelSize ?? 0,
    structure: { ...defaults.structure, approximate: false, alphaThickness: 12, clipPrimitive: true, ...structure },
  };
  await repr.createOrUpdate(props, list).run();
  return repr as unknown as Representation.Any;
}

/**
 * The PetWorld particles with the particle target's defaults as the viewer adds them (sized by `particle-size`: the
 * targets' elements by their physical radii), culled on the GPU
 */
async function postsynViewerRepr(graphics: GraphicsContext, list: ParticleList, props: { [k: string]: unknown } = {}) {
  const registries = getThemeRegistries();
  if (!registries.sizeThemeRegistry.has(ParticleSizeThemeProvider))
    registries.sizeThemeRegistry.add(ParticleSizeThemeProvider);
  const provider = ParticleTargetRepresentationProvider;
  const repr = provider.factory({ graphics, ...registries }, provider.getParams);
  repr.setTheme(
    Theme.create(
      registries,
      { particles: list },
      {
        colorTheme: { name: provider.defaultColorTheme.name, params: {} },
        sizeTheme: { name: ParticleSizeThemeProvider.name, params: { scale: 1, defaultSize: 0 } },
      },
    ),
  );
  await repr
    .createOrUpdate({ ...PD.getDefaultValues(ParticleTargetRepresentationParams), gpuCulling: true, ...props }, list)
    .run();
  return repr as unknown as Representation.Any;
}

/**
 * GPU culling on a large PetWorld model (on demand, needs `postsynstubby_2-no_bonds.bcif` in `.benchmark/`):
 * its particle targets as spacefill (merged by `auto`), culled on the CPU or the GPU, still or rotating
 */
/**
 * What the one-phase occlusion test of GPU culling misses while the camera moves (`two-phase-culling.md`): PetWorld
 * whole and zoomed, rotating 2, 10, 30 and 90 (a jump) degrees a frame, with the GPU's occlusion and without, the
 * impostors and with the sphere raster on (`-raster`, `sphereRaster` 'on' at 6 px: its instances those of the same
 * culling, `sphere-compute-raster.md` step 5); with `--images` the last frame drawn after one more step, the images of
 * `-occlusion` and `-none` differ where it popped
 */
export const OcclusionLagCases: BenchmarkCase<RenderInput>[] = [2, 10, 30, 90].flatMap((degrees) =>
  [undefined, 0.3].flatMap((zoom) =>
    [false, true].flatMap((raster) =>
      (['occlusion', 'one-phase', 'none'] as const).map((occlusion) =>
        renderCase(
          `occlusion-lag/postsyn${zoom ? '-zoom' : ''}-${degrees}deg${raster ? '-raster' : ''}-${occlusion}`,
          `postsynstubby_2 culled on the GPU${zoom ? `, zoomed to ${zoom} of the distance` : ''}, rotating ${degrees}° a frame, occlusion on the GPU: ${occlusion}${raster ? ', the sphere raster on' : ''}`,
          async () => {
            const list = await getPostsyn();
            return {
              scene: (g) => Promise.all([postsynRepr(g, list, true)]),
              zoom,
              rotate: true,
              rotateDegrees: degrees,
              props: {
                hiZ: {
                  ...PD.getDefaultValues(HiZParams),
                  gpu: occlusion !== 'none',
                  twoPhase: occlusion === 'occlusion',
                },
                ...(raster && {
                  renderer: {
                    ...PD.getDefaultValues(Canvas3DParams).renderer,
                    sphereRaster: 'on',
                    sphereRasterRadius: 6,
                  },
                }),
              },
            };
          },
        ),
      ),
    ),
  ),
);

/**
 * What the one-phase occlusion test misses while zooming in: the models as the viewer draws them, the camera moving
 * towards the target by 1, 3 and 8 % of its distance a frame, from the whole view and from 0.3 of the distance, with the
 * GPU's occlusion in two phases, in one and without; with `--images` the frame after the last, the images of
 * `-one-phase` and `-none` differ where it lost something, those of `-occlusion` and `-none` not
 */
export const OcclusionDollyCases: BenchmarkCase<RenderInput>[] = (['postsyn', 'mg'] as const).flatMap((model) =>
  [1, 3, 8].flatMap((percent) =>
    [undefined, 0.3].flatMap((zoom) =>
      (['occlusion', 'one-phase', 'none'] as const).map((occlusion) =>
        renderCase(
          `occlusion-dolly/${model}${zoom ? '-zoom' : ''}-${percent}pct-${occlusion}`,
          `${model} as the viewer draws it${zoom ? `, from ${zoom} of the distance` : ''}, moving in by ${percent} % a frame; occlusion on the GPU: ${occlusion}`,
          async () => {
            const list = await (model === 'mg' ? getMg() : getPostsyn());
            return {
              scene: (g) => Promise.all([postsynViewerRepr(g, list)]),
              zoom,
              dolly: percent / 100,
              yieldFrames: true,
              props: {
                hiZ: {
                  ...PD.getDefaultValues(HiZParams),
                  gpu: occlusion !== 'none',
                  twoPhase: occlusion === 'occlusion',
                },
              },
            };
          },
        ),
      ),
    ),
  ),
);

/**
 * The GPU's occlusion in the views the frame's pyramid is not built for or not alone: stereo's eyes, still and moving,
 * and multi-sampling's samples moving; with `--images`, as `OcclusionDollyCases`, `-occlusion` and `-none` differ where it
 * lost something
 */
export const OcclusionViewCases: BenchmarkCase<RenderInput>[] = (['stereo', 'multi-sample'] as const).flatMap((view) =>
  (view === 'stereo' ? [0, 0.03] : [0.03]).flatMap((dolly) =>
    (['occlusion', 'none'] as const).map((occlusion) =>
      renderCase(
        `occlusion-view/postsyn-zoom-${view}-${dolly ? `${dolly * 100}pct` : 'still'}-${occlusion}`,
        `postsyn as the viewer draws it, from 0.3 of the distance, ${view}, ${dolly ? `moving in by ${dolly * 100} % a frame` : 'still'}; occlusion on the GPU: ${occlusion}`,
        async () => {
          const list = await getPostsyn();
          const { camera } = PD.getDefaultValues(Canvas3DParams);
          return {
            scene: (g) => Promise.all([postsynViewerRepr(g, list)]),
            zoom: 0.3,
            dolly: dolly || undefined,
            yieldFrames: true,
            multiSample: view === 'multi-sample',
            props: {
              hiZ: { ...PD.getDefaultValues(HiZParams), gpu: occlusion === 'occlusion' },
              ...(view === 'stereo' && {
                camera: { ...camera, stereo: { name: 'on' as const, params: PD.getDefaultValues(StereoCameraParams) } },
              }),
            },
          };
        },
      ),
    ),
  ),
);

export const MergingPostsynCases: BenchmarkCase<RenderInput>[] = (
  ['cpu', 'cpu-hiz', 'gpu', 'gpu-hiz', 'gpu-hiz-cell50', 'gpu-hiz-cell20'] as const
).flatMap((culling) =>
  [false, true].flatMap((rotate) =>
    [undefined, 0.3].map((zoom) =>
      renderCase(
        `merging/postsyn-${culling}${zoom ? '-zoom' : ''}${rotate ? '-rotating' : ''}`,
        `postsynstubby_2 particles as spacefill (mesoscale settings), culled on the ${culling.toUpperCase()}${zoom ? `, zoomed to ${zoom} of the distance` : ''}${rotate ? ', rotating' : ''}`,
        async () => {
          const list = await getPostsyn();
          return {
            scene: (g) =>
              Promise.all([
                postsynRepr(
                  g,
                  list,
                  culling.startsWith('gpu'),
                  culling.endsWith('cell50') ? 50 : culling.endsWith('cell20') ? 20 : undefined,
                ),
              ]),
            setupNotes: true,
            byteNotes: true,
            rotate,
            zoom,
            // occlusion culling of the CPU (Hi-Z, read back with a lag of frames)
            props: culling.includes('hiz') ? { hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true } } : undefined,
          };
        },
      ),
    ),
  ),
);

/**
 * The zoomed postsyn view at a quarter and four times the pixels, with other sphere LOD presets and cheaper impostors:
 * what of the opaque pass scales with the pixels (fragments), what with the spheres drawn (vertices, instances) and
 * what the impostors' shading and projection cost
 */
export const MergingPostsynSizeCases: BenchmarkCase<RenderInput>[] = (
  [
    { label: 'small', note: 'at 640×360', size: { width: 640, height: 360 } },
    { label: 'large', note: 'at 2560×1440', size: LargeSize },
    { label: 'lod-performance', note: 'LOD preset performance', lod: 'performance' },
    { label: 'lod-ultra', note: 'LOD preset ultra', lod: 'ultra' },
    { label: 'ignore-light', note: 'without lighting', structure: { ignoreLight: true } },
    { label: 'approximate', note: 'approximate impostors', structure: { approximate: true } },
  ] as const
).map(
  ({
    label,
    note,
    ...v
  }: {
    label: string;
    note: string;
    size?: { width: number; height: number };
    lod?: keyof typeof Spheres.LodLevelsPresets;
    structure?: { ignoreLight?: boolean; approximate?: boolean };
  }) =>
    renderCase(
      `merging/postsyn-gpu-hiz-zoom-${label}`,
      `postsynstubby_2 particles as spacefill (mesoscale settings), culled on the GPU with Hi-Z, zoomed to 0.3 of the distance, ${note}`,
      async () => {
        const list = await getPostsyn();
        return {
          size: v.size,
          scene: (g) => Promise.all([postsynRepr(g, list, true, undefined, v.lod, v.structure)]),
          setupNotes: true,
          byteNotes: true,
          zoom: 0.3,
          props: { hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true } },
        };
      },
    ),
);

/**
 * The `quality` preset's levels with distances from their size on screen (screen-space LOD, step 1): level `j` where a
 * pixel covers from `2 r s_(j-1)^(1 / b_(j-1)) / D` to `2 r s_j^(1 / b_j) / D` Å (its spheres at least `D` pixels
 * across, those of the finer level not), as distances for a perspective camera of field of view `fov` and `height`
 * pixels; `r` 1.7 Å
 */
function pixelLodLevels(
  pixelSize: number,
  fov: number,
  height: number,
  levels: LodLevels = Spheres.LodLevelsPresets.quality,
): LodLevels {
  const k = (2 * Math.tan(fov / 2)) / height;
  let minDistance = levels[0].minDistance;
  return levels.map((l, j) => {
    const scale = Math.pow(l.stride, 1 / l.scaleBias);
    const maxDistance = j === levels.length - 1 ? l.maxDistance : (2 * 1.7 * scale) / pixelSize / k;
    const level = { ...l, minDistance, maxDistance };
    minDistance = maxDistance;
    return level;
  });
}

/**
 * The postsyn views (whole and zoomed, at 1280x720 and 2560x1440) with levels chosen by their size on screen at a
 * smallest sphere of `D` pixels, against the `quality` preset (`-quality`): times and, with `--images`, the pictures
 */
export const MergingPostsynPixelCases: BenchmarkCase<RenderInput>[] = ([undefined, 1.5, 2, 3, 4] as const).flatMap(
  (pixelSize) =>
    [undefined, 0.3].flatMap((zoom) =>
      [Size, LargeSize].map((size) =>
        renderCase(
          `merging/postsyn-lod-${pixelSize ? `px${pixelSize}` : 'quality'}${zoom ? '-zoom' : ''}${size === LargeSize ? '-large' : ''}`,
          `postsynstubby_2 particles as spacefill (mesoscale settings), culled on the GPU with Hi-Z, ${zoom ? `zoomed to ${zoom} of the distance, ` : ''}at ${size.width}×${size.height}, ${pixelSize ? `levels for spheres of at least ${pixelSize} px` : 'LOD preset quality'}`,
          async () => {
            const list = await getPostsyn();
            const fov = (PD.getDefaultValues(Canvas3DParams).camera.fov * Math.PI) / 180;
            const lod = pixelSize ? pixelLodLevels(pixelSize, fov, size.height) : undefined;
            return {
              size,
              scene: (g) => Promise.all([postsynRepr(g, list, true, undefined, lod)]),
              zoom,
              props: { hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true } },
            };
          },
        ),
      ),
    ),
);

/** Coarser chains of 4 levels (strides, scale biases) for the pixel rule: what reaches the whole view */
const CoarseChains: { [name: string]: [stride: number, scaleBias: number][] } = {
  'x10-b2.3': [
    [1, 1],
    [10, 3],
    [100, 2.7],
    [1000, 2.3],
  ],
  'x10-b2': [
    [1, 1],
    [10, 2.5],
    [100, 2.2],
    [1000, 2],
  ],
  'x20-b2': [
    [1, 1],
    [20, 2.5],
    [400, 2.2],
    [8000, 2],
  ],
};

/** The postsyn views at 1280x720 with coarser chains at `D` 3 (`merging/postsyn-lod-px3[-zoom]` the `quality` chain) */
export const MergingPostsynChainCases: BenchmarkCase<RenderInput>[] = Object.entries(CoarseChains).flatMap(
  ([name, chain]) =>
    [undefined, 0.3].map((zoom) =>
      renderCase(
        `merging/postsyn-lod-chain-${name}${zoom ? '-zoom' : ''}`,
        `postsynstubby_2 particles as spacefill (mesoscale settings), culled on the GPU with Hi-Z, ${zoom ? `zoomed to ${zoom} of the distance, ` : ''}levels of strides ${chain.map((c) => c[0]).join(', ')} for spheres of at least 3 px`,
        async () => {
          const list = await getPostsyn();
          const fov = (PD.getDefaultValues(Canvas3DParams).camera.fov * Math.PI) / 180;
          const base = Spheres.LodLevelsPresets.quality.map((l, j) => ({
            ...l,
            stride: chain[j][0],
            scaleBias: chain[j][1],
          }));
          const lod = pixelLodLevels(3, fov, Size.height, base);
          return {
            scene: (g) => Promise.all([postsynRepr(g, list, true, undefined, lod)]),
            zoom,
            props: { hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true } },
          };
        },
      ),
    ),
);

/**
 * Levels by size on screen in the renderer (`lodPixelSize` 3) with the first coarser chain, whole and zoomed, at 720p and
 * 1440p: as `merging/postsyn-lod-chain-x10-b2.3[-zoom]` (distances computed in the benchmark) at 720p
 */
type PixelModeCase = {
  size: typeof Size;
  ortho: boolean;
  pixelSize: number;
  /** Hi-Z (the CPU's test) off: the GPU's own, as by default */
  noHiZ?: boolean;
  /** no occlusion test at all (`hiZ.gpu` off as well) */
  noOcclusion?: boolean;
  /** culled on the CPU (`gpuCulling` off; no occlusion test with Hi-Z off) */
  cpuCulling?: boolean;
  rotate?: boolean;
};

export const MergingPostsynPixelModeCases: BenchmarkCase<RenderInput>[] = [undefined, 0.3].flatMap((zoom) =>
  (
    [
      { size: Size, ortho: false, pixelSize: 3 },
      { size: LargeSize, ortho: false, pixelSize: 3 },
      { size: Size, ortho: true, pixelSize: 3 },
      ...[1.5, 2, 4, 6].map((pixelSize) => ({ size: Size, ortho: false, pixelSize })),
      // the GPU's occlusion test alone (as by default), none: what it gains
      { size: Size, ortho: false, pixelSize: 3, noHiZ: true },
      { size: Size, ortho: false, pixelSize: 3, noHiZ: true, noOcclusion: true },
      // rotating, with it and without: what it misses while the view moves (the previous frame's depth); the CPU's
      // test, whose readback lags a few frames
      { size: Size, ortho: false, pixelSize: 3, rotate: true },
      { size: Size, ortho: false, pixelSize: 3, noHiZ: true, rotate: true },
      { size: Size, ortho: false, pixelSize: 3, noHiZ: true, noOcclusion: true, rotate: true },
      // culled on the CPU: the GPU's culling pass against the CPU's (on a GPU-bound device)
      { size: Size, ortho: false, pixelSize: 3, noHiZ: true, cpuCulling: true },
      { size: Size, ortho: false, pixelSize: 3, noHiZ: true, cpuCulling: true, rotate: true },
    ] as PixelModeCase[]
  ).map(({ size, ortho, pixelSize, noHiZ, noOcclusion, cpuCulling, rotate }) =>
    renderCase(
      `merging/postsyn-lod-pixel${pixelSize}${zoom ? '-zoom' : ''}${size === LargeSize ? '-large' : ''}${ortho ? '-ortho' : ''}${cpuCulling ? '-cpu' : noOcclusion ? '-noocclusion' : noHiZ ? '-nohiz' : ''}${rotate ? '-rotating' : ''}`,
      `postsynstubby_2 particles as spacefill (mesoscale settings), culled on the ${cpuCulling ? 'CPU' : 'GPU'}${cpuCulling || noOcclusion ? '' : noHiZ ? ' with its occlusion test' : ' with Hi-Z'}, ${zoom ? `zoomed to ${zoom} of the distance, ` : ''}at ${size.width}×${size.height}${ortho ? ', orthographic' : ''}${rotate ? ', rotating' : ''}, levels of strides 1, 10, 100, 1000 by size on screen (${pixelSize} px)`,
      async () => {
        const list = await getPostsyn();
        const chain = CoarseChains['x10-b2.3'];
        const lod = Spheres.LodLevelsPresets.quality.map((l, j) => ({
          ...l,
          stride: chain[j][0],
          scaleBias: chain[j][1],
        }));
        const camera = PD.getDefaultValues(Canvas3DParams).camera;
        return {
          size,
          scene: (g) => Promise.all([postsynRepr(g, list, !cpuCulling, undefined, lod, undefined, pixelSize)]),
          setupNotes: true,
          byteNotes: true,
          zoom,
          rotate,
          // the CPU's test reads the levels back
          yieldFrames: !noHiZ,
          props: {
            hiZ: { ...PD.getDefaultValues(HiZParams), enabled: !noHiZ, gpu: !noOcclusion },
            ...(ortho ? { camera: { ...camera, mode: 'orthographic' as const } } : {}),
          },
        };
      },
    ),
  ),
);

/**
 * The sphere raster (`sphere-compute-raster.md`) through the renderer (`sphereRaster`) on the PetWorld particles as
 * `merging/postsyn-lod-pixel1` draws them (levels by size on screen at 1 px: more spheres, and larger, than at the
 * default 3), whole and zoomed, at 1280×720 and 2560×1440, and multi-sampled (`-ms`: four samples a frame, views of
 * a frame culled already, as temporal multi-sampling draws them at rest): off (the impostors), on at
 * `sphereRasterRadius` 4, 6 and 8 CSS pixels (larger spheres as impostors), `only` (every sphere by the raster, no
 * impostor); needs `atomic-vec2u-min-max` (Chrome 156: `--browser-arg --enable-webgpu-developer-features`)
 */
export const SphereRasterRendererCases: BenchmarkCase<RenderInput>[] = [undefined, 0.3].flatMap((zoom) =>
  [Size, LargeSize, 'multi-sample' as const, 'viewer' as const].flatMap((variant) =>
    (
      [
        { label: 'off', sphereRaster: 'off', radius: 6 },
        ...[4, 6, 8].map((radius) => ({ label: `on${radius}`, sphereRaster: 'on' as const, radius })),
        { label: 'auto', sphereRaster: 'auto', radius: 6 },
        { label: 'only', sphereRaster: 'only', radius: 6 },
      ] as const
    ).map(({ label, sphereRaster, radius }) => {
      const size = variant === 'multi-sample' || variant === 'viewer' ? Size : variant;
      const name =
        variant === LargeSize ? '-large' : variant === 'multi-sample' ? '-ms' : variant === 'viewer' ? '-viewer' : '';
      return renderCase(
        `sphere-raster-renderer/postsyn${zoom ? '-zoom' : ''}${name}-${label}`,
        `postsynstubby_2 particles as spacefill (mesoscale settings), culled on the GPU with Hi-Z, ${zoom ? `zoomed to ${zoom} of the distance, ` : ''}at ${size.width}×${size.height}${variant === 'multi-sample' ? ', multi-sampled (4 samples a frame)' : ''}, levels by size on screen (1 px); sphere raster ${sphereRaster}${sphereRaster === 'on' ? ` up to ${radius} px` : ''}`,
        async () => {
          const list = await getPostsyn();
          const chain = CoarseChains['x10-b2.3'];
          const lod = Spheres.LodLevelsPresets.quality.map((l, j) => ({
            ...l,
            stride: chain[j][0],
            scaleBias: chain[j][1],
          }));
          // the particles' defaults as the viewer draws them (`viewer`): their levels by size on screen at 3 px
          const viewer = variant === 'viewer';
          return {
            size,
            scene: (g) =>
              Promise.all([
                viewer ? postsynViewerRepr(g, list) : postsynRepr(g, list, true, undefined, lod, undefined, 1),
              ]),
            setupNotes: true,
            zoom,
            yieldFrames: true,
            multiSample: variant === 'multi-sample',
            props: {
              hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true },
              renderer: {
                ...PD.getDefaultValues(Canvas3DParams).renderer,
                sphereRaster,
                sphereRasterRadius: radius,
              },
            },
          };
        },
      );
    }),
  ),
);

/**
 * What still goes to the impostors as spheres grow past the raster's limit (`sphere-compute-raster.md` option A, tiles
 * for larger spheres): PetWorld as the viewer draws it, whole and zoomed to 0.3, 0.1 and 0.03 of the distance; `off`,
 * `auto` with the limit at 6, 24, 64 and 256 px (larger spheres to the raster, the large ones by its workgroups),
 * `only` (every sphere); on demand
 */
export const SphereRasterZoomCases: BenchmarkCase<RenderInput>[] = [1, 0.3, 0.1, 0.03].flatMap((zoom) =>
  (
    [
      { label: 'off', sphereRaster: 'off', radius: 6 },
      ...[6, 24, 64, 256].map((radius) => ({ label: `auto${radius}`, sphereRaster: 'auto' as const, radius })),
      { label: 'only', sphereRaster: 'only', radius: 6 },
    ] as const
  ).map(({ label, sphereRaster, radius }) =>
    renderCase(
      `sphere-raster-zoom/postsyn-viewer-${zoom === 1 ? 'whole' : `zoom${zoom}`}-${label}`,
      `postsynstubby_2 as the viewer draws it, ${zoom === 1 ? 'whole' : `zoomed to ${zoom} of the distance`}; sphere raster ${sphereRaster}${sphereRaster === 'auto' ? ` up to ${radius} px` : ''}`,
      async () => {
        const list = await getPostsyn();
        return {
          scene: (g) => Promise.all([postsynViewerRepr(g, list)]),
          zoom: zoom === 1 ? undefined : zoom,
          yieldFrames: true,
          props: {
            hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true },
            renderer: { ...PD.getDefaultValues(Canvas3DParams).renderer, sphereRaster, sphereRasterRadius: radius },
          },
        };
      },
    ),
  ),
);

type MarkingVariant = {
  readonly marking: 'none' | 'particle' | 'hover' | 'all';
  /** forced frames by default; `redraw`: the camera still, the marking redrawn over the last image; `rotating`: the
   * camera rotating, the scene and marking drawn per frame */
  readonly mode?: 'redraw' | 'rotating';
  /** a pause of 300 ms before each frame (`RenderInput.idleMs`): the GPU's clocks fall as between hovers */
  readonly idle?: boolean;
  /** dimming on; an unmarked transparent box (alpha 0.5) in the model's middle; multi-sampling */
  readonly extra?: 'dim' | 'transparent' | 'multi-sample';
};

const MarkingVariants: readonly MarkingVariant[] = [
  { marking: 'none' },
  { marking: 'particle' },
  { marking: 'hover' },
  { marking: 'all' },
  { marking: 'hover', mode: 'redraw' },
  { marking: 'hover', mode: 'redraw', idle: true },
  { marking: 'none', mode: 'rotating' },
  { marking: 'hover', mode: 'rotating' },
];

const MarkingExtraVariants: readonly MarkingVariant[] = (['dim', 'transparent', 'multi-sample'] as const).flatMap(
  (extra) => [
    { marking: 'hover', mode: 'redraw', extra },
    { marking: 'none', mode: 'rotating', extra },
    { marking: 'hover', mode: 'rotating', extra },
  ],
);

/**
 * Marking on postsynstubby_2 as the viewer draws it, whole and zoomed to 0.3: nothing marked, one particle highlighted,
 * a hover (per frame another particle highlighted, the one before not) and everything highlighted, in forced frames;
 * the hover as a still camera draws it (`-redraw`: the marking alone, redrawn over the last image; `-redraw-idle`: 300
 * ms apart, at the GPU's idle clocks) and as a moving one does (`-rotating`, against nothing marked rotating); zoomed,
 * those also with dimming, an unmarked transparent box and multi-sampling; on demand
 */
export const MarkingPostsynCases: BenchmarkCase<RenderInput>[] = [undefined, 0.3].flatMap((zoom) =>
  [...MarkingVariants, ...(zoom ? MarkingExtraVariants : [])].map(({ marking, mode, idle, extra }) =>
    renderCase(
      `marking/postsyn-${zoom ? 'zoom' : 'whole'}-${marking}${mode ? `-${mode}` : ''}${idle ? '-idle' : ''}${extra ? `-${extra}` : ''}`,
      `postsynstubby_2 as the viewer draws it${zoom ? `, zoomed to ${zoom} of the distance` : ''}; marked: ${marking}${mode === 'redraw' ? ', redrawn over the last image' : mode === 'rotating' ? ', rotating' : ''}${idle ? ', 300 ms apart' : ''}${extra ? `, ${extra}` : ''}`,
      async () => {
        const list = await getPostsyn();
        const particle = (i: number) => Particle.Loci(list, OrderedSet.ofSingleton((i * 7919) % list.count));
        const { center, radius } = Particle.getBoundary(list).sphere;
        const box = () =>
          shapeRepr('transparent-box', meshBox(radius * 0.6), {
            transforms: [Mat4.fromTranslation(Mat4(), center)],
            props: { alpha: 0.5 },
          });
        const { renderer } = PD.getDefaultValues(Canvas3DParams);
        return {
          scene: async (g) => [await postsynViewerRepr(g, list), ...(extra === 'transparent' ? [await box()] : [])],
          zoom,
          yieldFrames: true,
          unforced: mode === 'redraw',
          idleMs: idle ? 300 : undefined,
          rotate: mode === 'rotating',
          multiSample: extra === 'multi-sample',
          highlight: marking === 'all',
          update:
            marking === 'particle'
              ? (canvas3d, _, i) => {
                  if (i === 0) canvas3d.mark({ loci: particle(0) }, MarkerAction.Highlight);
                }
              : marking === 'hover'
                ? (canvas3d, _, i) => {
                    canvas3d.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
                    canvas3d.mark({ loci: particle(i + 1) }, MarkerAction.Highlight);
                  }
                : undefined,
          props: {
            hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true },
            ...(extra === 'dim' && { renderer: { ...renderer, dimStrength: 0.5 } }),
          },
        };
      },
    ),
  ),
);

/**
 * Marking on 60 instanced copies of 7qpd as spacefill, a render object per unit (not merged: instanced, each culled by
 * its instance grid): a hover over the chains of the copies (per frame another copy's chain highlighted, the one before
 * not), redrawn over the last image and rotating
 */
export const MarkingAssemblyCases: BenchmarkCase<RenderInput>[] = (['redraw', 'rotating'] as const).map((mode) =>
  renderCase(
    `marking/assembly-60-hover-${mode}`,
    `60 instanced copies of 7qpd as spacefill, unmerged; a chain of a copy hovered, ${mode === 'redraw' ? 'redrawn over the last image' : 'rotating'}`,
    async () => {
      const s = await getAssembly();
      const unit = (i: number) => {
        const u = s.units[(i * 7) % s.units.length];
        return StructureElement.Loci(s, [
          { unit: u, indices: OrderedSet.ofBounds(0, u.elements.length) as StructureElement.UnitIndex[] & never },
        ]);
      };
      return {
        scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, { renderMerged: false })]),
        yieldFrames: true,
        unforced: mode === 'redraw',
        rotate: mode === 'rotating',
        update: (canvas3d, _, i) => {
          canvas3d.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
          canvas3d.mark({ loci: unit(i + 1) }, MarkerAction.Highlight);
        },
      };
    },
  ),
);

let mg: Promise<ParticleList> | undefined;
/** `MG_1184_cellpack_atom_instances.bcif` in `.benchmark/` as the particle list the viewer loads */
const getMg = () =>
  (mg ??= (async () => {
    const parsed = await CIF.parseBinary(await fetchBenchmarkFile('MG_1184_cellpack_atom_instances.bcif')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const cif = parsed.result;
    const assemblyId = getAssemblyIdsFromMmcif(cif)[0];
    return Task.create('particles', (ctx) => createParticleListFromMmcifAssembly(ctx, cif, { assemblyId })).run();
  })());

/**
 * The impostors against the sphere raster on a cellPACK model unlike PetWorld (848 units, 2.1 M instances of about 76
 * spheres, a scene a fifth of the size), as the viewer draws it: whole and zoomed, 1280×720; `off` (culled on the GPU),
 * `auto`, `only` and `cpu` (the impostors culled on the CPU); whole and zoomed to 0.3 at 2560×1440 and rotating; on demand
 */
export const SphereRasterMgCases: BenchmarkCase<RenderInput>[] = (
  [
    ...[1, 0.3, 0.1, 0.03].flatMap((zoom) =>
      (['off', 'auto', 'only', 'cpu'] as const).map((variant) => ({ zoom, variant, view: '' })),
    ),
    ...[1, 0.3].flatMap((zoom) =>
      (['off', 'auto'] as const).flatMap((variant) => [
        { zoom, variant, view: 'large' },
        { zoom, variant, view: 'rotate' },
      ]),
    ),
  ] as { zoom: number; variant: 'off' | 'auto' | 'only' | 'cpu'; view: string }[]
).map(({ zoom, variant, view }) =>
  renderCase(
    `sphere-raster-mg/${zoom === 1 ? 'whole' : `zoom${zoom}`}${view ? `-${view}` : ''}-${variant}`,
    `MG_1184 as the viewer draws it, ${zoom === 1 ? 'whole' : `zoomed to ${zoom} of the distance`}${view === 'large' ? ', 2560×1440' : view === 'rotate' ? ', rotating' : ''}; ${variant === 'cpu' ? 'culled on the CPU, the impostors' : `sphere raster ${variant}`}`,
    async () => {
      const list = await getMg();
      return {
        size: view === 'large' ? LargeSize : Size,
        rotate: view === 'rotate',
        setupNotes: true,
        byteNotes: true,
        zoom: zoom === 1 ? undefined : zoom,
        yieldFrames: true,
        scene: (g) => Promise.all([postsynViewerRepr(g, list, variant === 'cpu' ? { gpuCulling: false } : {})]),
        props: {
          hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true },
          renderer: {
            ...PD.getDefaultValues(Canvas3DParams).renderer,
            sphereRaster: variant === 'cpu' ? 'off' : variant,
          },
        },
      };
    },
  ),
);

/**
 * Illumination on the mesoscale models as the viewer draws them (PetWorld, MG_1184; Hi-Z and the sphere raster at their
 * defaults, the ray distance of the mesoscale explorer): rotating, the first tracing iteration per frame, its input
 * (tracing and the back depth) drawn each frame; `still`, the camera still and the frames not forced, the iterations
 * going on (no last one) at the quality adapted up to the props' maximum, the input drawn once; whole and zoomed to 0.3,
 * 1280×720 and 2560×1440; the thickness by the back depth (`auto`, the default) or fixed (no back depth) at 4 Å (the
 * default, `-fixed`), 16, 64 and 256 Å (`-fixed16` …); on demand (`--gpu-timing`: the trace and its inputs apart)
 */
export const IlluminationMesoscaleCases: BenchmarkCase<RenderInput>[] = (
  [
    ['postsyn', () => getPostsyn()],
    ['mg', () => getMg()],
  ] as const
).flatMap(([name, list]) =>
  [1, 0.3].flatMap((zoom) =>
    [false, true].flatMap((large) =>
      [false, true].flatMap((still) =>
        ([undefined, 4, 16, 64, 256] as const).map((thickness) =>
          renderCase(
            `illumination-mesoscale/${name}-${zoom === 1 ? 'whole' : `zoom${zoom}`}${large ? '-large' : ''}${still ? '-still' : ''}${thickness === undefined ? '' : thickness === 4 ? '-fixed' : `-fixed${thickness}`}`,
            `${name === 'postsyn' ? 'postsynstubby_2' : 'MG_1184'} as the viewer draws it with illumination, ${zoom === 1 ? 'whole' : `zoomed to ${zoom} of the distance`}, ${still ? 'still' : 'rotating'}${large ? ', 2560×1440' : ''}, thickness ${thickness === undefined ? 'auto' : `fixed at ${thickness} Å`}`,
            async () => {
              const l = await list();
              const illumination = {
                ...PD.getDefaultValues(IlluminationParams),
                enabled: true,
                rayDistance: 1024,
                ...(thickness === undefined ? {} : { thicknessMode: 'fixed' as const, thickness }),
              };
              return {
                size: large ? LargeSize : Size,
                rotate: !still,
                unforced: still,
                zoom: zoom === 1 ? undefined : zoom,
                yieldFrames: true,
                scene: (g) => Promise.all([postsynViewerRepr(g, l)]),
                props: {
                  hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true },
                  illumination: still ? { ...illumination, maxIterations: 16 } : illumination,
                },
              };
            },
          ),
        ),
      ),
    ),
  ),
);

/**
 * The sphere raster's prop on scenes other than PetWorld (`sphere-compute-raster.md` step 6, what `auto` decides from):
 * spacefill with the default props, Hi-Z on, of a crambin (327 atoms), a globule of 6k and one of 300k atoms (an
 * unmerged object each: the raster's split per sphere), 60 copies of 7qpd as an object per unit (instanced, unmerged)
 * and merged and culled on the GPU (the split per instance); whole and zoomed, at 1280×720 and 2560×1440; `off`,
 * `auto`, `on` and `only` at the default limit (`sphereRasterRadius`)
 */
export const SphereRasterSceneCases: BenchmarkCase<RenderInput>[] = (
  [
    ['crambin', () => exampleStructure('1crn.cif'), {}],
    ['globule-6k', () => globuleStructure(6000), {}],
    ['globule-300k', getGlobule, {}],
    ['assembly-60', getRootAssembly, {}],
    ['assembly-60-merged', getRootAssembly, { renderMerged: true, gpuCulling: true }],
  ] as const
).flatMap(([scene, structure, props]) =>
  [undefined, 0.3].flatMap((zoom) =>
    [Size, LargeSize].flatMap((size) =>
      (['off', 'auto', 'on', 'only'] as const).map((sphereRaster) =>
        renderCase(
          `sphere-raster-scene/${scene}${zoom ? '-zoom' : ''}${size === LargeSize ? '-large' : ''}-${sphereRaster}`,
          `${scene} as spacefill${zoom ? `, zoomed to ${zoom} of the distance` : ''}, at ${size.width}×${size.height}; sphere raster ${sphereRaster}`,
          async () => {
            const s = await structure();
            return {
              size,
              zoom,
              yieldFrames: true,
              scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, props)]),
              props: {
                hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true },
                renderer: { ...PD.getDefaultValues(Canvas3DParams).renderer, sphereRaster },
              },
            };
          },
        ),
      ),
    ),
  ),
);

/**
 * The sphere raster's prop where PetWorld is too large (an integrated GPU: `--browser-arg --force_low_power_gpu`):
 * 3j3q, the HIV-1 capsid (`3j3q.bcif` in `.benchmark/`, about 2.4 M atoms, 1356 chains), as spacefill merged and
 * culled on the GPU (what `auto` draws), Hi-Z on, the default limit; whole and zoomed at 1280×720, whole at 2560×1440,
 * and zoomed cut through its center by a plane by pixel; `off`, `auto` and `only`
 */
export const SphereRasterIGpuCases: BenchmarkCase<RenderInput>[] = (
  [
    { label: 'whole', size: Size },
    { label: 'whole-large', size: LargeSize },
    { label: 'zoom', size: Size, zoom: 0.3 },
    { label: 'zoom-plane', size: Size, zoom: 0.3, plane: true },
  ] as const
).flatMap(({ label, size, zoom, plane }: { label: string; size: typeof Size; zoom?: number; plane?: boolean }) =>
  (['off', 'auto', 'only'] as const).map((sphereRaster) =>
    renderCase(
      `sphere-raster-igpu/3j3q-${label}-${sphereRaster}`,
      `3j3q as spacefill, merged and culled on the GPU${zoom ? `, zoomed to ${zoom} of the distance` : ''}${plane ? ', cut through its center by a plane' : ''}, at ${size.width}×${size.height}; sphere raster ${sphereRaster}`,
      async () => {
        const s = await getBenchmarkEntry('3j3q.bcif');
        const clip = plane
          ? {
              variant: 'pixel' as const,
              objects: [
                {
                  type: 'plane' as const,
                  invert: false,
                  position: Vec3.copy(Vec3(), s.boundary.sphere.center),
                  rotation: { axis: Vec3.create(1, 0, 0), angle: 30 },
                  scale: Vec3.create(1, 1, 1),
                  transform: Mat4.identity(),
                },
              ],
            }
          : undefined;
        return {
          size,
          zoom,
          yieldFrames: true,
          scene: (g) =>
            Promise.all([
              structureRepr(g, SpacefillRepresentationProvider, s, {
                renderMerged: true,
                gpuCulling: true,
                ...(clip && { clip }),
              }),
            ]),
          props: {
            hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true },
            renderer: { ...PD.getDefaultValues(Canvas3DParams).renderer, sphereRaster },
          },
        };
      },
    ),
  ),
);

/**
 * A recolor per frame of the PetWorld particles through the representation's update (a theme change): a uniform color
 * per molecule type, or the entity colors by name and by function in turn (a color per particle)
 */
export const MergingPostsynUpdateCases: BenchmarkCase<RenderInput>[] = (['uniform', 'entity'] as const).map((theme) =>
  renderCase(
    `merging/postsyn-lod-pixel3-recolor${theme === 'entity' ? '-entity' : ''}`,
    `postsynstubby_2 particles as spacefill (mesoscale settings, levels by size on screen), recolored every frame (${theme === 'entity' ? 'entity colors by name and by function' : 'a uniform color'})`,
    async () => {
      const list = await getPostsyn();
      const chain = CoarseChains['x10-b2.3'];
      const lod = Spheres.LodLevelsPresets.quality.map((l, j) => ({
        ...l,
        stride: chain[j][0],
        scaleBias: chain[j][1],
      }));
      const registries = getThemeRegistries();
      return {
        scene: (g) => Promise.all([postsynRepr(g, list, true, undefined, lod, undefined, 3)]),
        update: async (_canvas3d, reprs, i) => {
          const repr = reprs[0];
          repr.setTheme(
            Theme.create(
              registries,
              { particles: list },
              {
                colorTheme:
                  theme === 'entity'
                    ? { name: 'particle-entity', params: { by: i % 2 ? 'function' : 'name' } }
                    : { name: 'uniform', params: { value: Color(i % 2 ? 0xdd4422 : 0x2244dd) } },
                sizeTheme: { name: 'uniform', params: {} },
              },
            ),
          );
          await repr.createOrUpdate({}, list).run();
        },
      };
    },
  ),
);

/**
 * Quads (sphere impostors) drawn unindexed (`setQuadIndexing(false)`): against the indexed draws of the same cases
 * (`render/7qpd-spacefill`, `render/globule-300k-spacefill`, `merging/postsyn-gpu-hiz-zoom`, `merging/postsyn-lod-pixel3`)
 */

export const UnindexedQuadCases: BenchmarkCase<RenderInput>[] = [
  renderCase('render/7qpd-spacefill-unindexed', 'PDB entry 7qpd as spacefill, quads unindexed', async () => {
    const s = await getStructure();

    return {
      quadIndexing: false,

      scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, {})]),
    };
  }),

  renderCase('render/globule-300k-spacefill-unindexed', '300k atoms as spacefill, quads unindexed', async () => {
    const s = await getGlobule();

    return {
      quadIndexing: false,

      scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, {})]),
    };
  }),

  ...[undefined, 0.3].map((zoom) =>
    renderCase(
      `merging/postsyn-lod-pixel3${zoom ? '-zoom' : ''}-unindexed`,

      `postsynstubby_2 particles as spacefill, GPU culled with Hi-Z, levels by size on screen (3 px), ${zoom ? 'zoomed, ' : ''}quads unindexed`,

      async () => {
        const list = await getPostsyn();

        const lod = Spheres.LodPixelPresets.quality.lodLevels;

        return {
          quadIndexing: false,

          scene: (g) => Promise.all([postsynRepr(g, list, true, undefined, lod, undefined, 3)]),

          zoom,

          props: { hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true } },
        };
      },
    ),
  ),
];

/** An `element-index` color theme (per group and instance) of a two-color palette, alternating per frame */
function elementIndexColor(i: number) {
  const colors = i % 2 ? [Color(0xdd4422), Color(0x2244dd)] : [Color(0x22dd44), Color(0xdddd22)];
  return {
    name: 'element-index',
    params: { palette: { name: 'colors', params: { list: { kind: 'interpolate', colors } } } },
  };
}

/**
 * Merging with instances (`merged-render-objects.md`, on demand): the 60 copies of 7qpd (its chains the units, 60
 * instances each) colored by `element-index`, per group and instance: the color layers have groups times instances
 * texels, copied into the merged object. Frames still or with the palette changed per frame
 */
export const MergingAssemblyCases: BenchmarkCase<RenderInput>[] = MergingKinds.flatMap((kind) =>
  [true, false].flatMap((merged) =>
    ['', '-recolor'].map((suffix) =>
      renderCase(
        `merging/assembly-60-${kind.name}-${merged ? 'merged' : 'unmerged'}${suffix}`,
        `60 copies of 7qpd as ${kind.name} colored by element index, ${merged ? 'merged' : 'a render object per unit'}${suffix ? ': a palette change per frame' : ''}`,
        async () => {
          const s = await getRootAssembly();
          const update: RenderInput['update'] = suffix
            ? async (canvas3d, reprs, i) => {
                const { provider } = kind;
                const sizeTheme = {
                  name: provider.defaultSizeTheme.name,
                  params: provider.defaultSizeTheme.props ?? {},
                };
                reprs[0].setTheme(
                  Theme.create(getThemeRegistries(), { structure: s }, { colorTheme: elementIndexColor(i), sizeTheme }),
                );
                // committed and drawn by the frame's tick, as `MergingCases`
                await reprs[0].createOrUpdate({}, s).run();
              }
            : undefined;
          return {
            scene: (g) =>
              Promise.all([structureRepr(g, kind.provider, s, { renderMerged: merged }, elementIndexColor(1))]),
            setupNotes: true,
            byteNotes: true,
            update,
          };
        },
      ),
    ),
  ),
);

/**
 * Clip objects that cut nothing (inverted spheres around everything) on the globule at 2560x1440 (fragment-bound): the
 * cost of the per-pixel clip tests, their loops bounded by a uniform (`uClipObjectCount`), on demand; with the sphere
 * raster (`-on`, `-only`: its tests per covered pixel of every sphere, the impostors' per fragment)
 */
export const ClipCases: BenchmarkCase<RenderInput>[] = [1, 3].flatMap((count) =>
  (['off', 'on', 'only'] as const).map((sphereRaster) =>
    renderCase(
      `clip/globule-300k-spacefill-large-${count}${sphereRaster === 'off' ? '' : `-${sphereRaster}`}`,
      `300k atoms as spacefill, 2560x1440, ${count} clip object${count > 1 ? 's' : ''} cutting nothing (pixel)${sphereRaster === 'off' ? '' : `, sphere raster ${sphereRaster}`}`,
      async () => {
        const s = await getGlobule();
        const object = {
          type: 'sphere' as const,
          invert: true,
          position: Vec3(),
          rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
          scale: Vec3.create(1e5, 1e5, 1e5),
          transform: Mat4.identity(),
        };
        const clip = { variant: 'pixel' as const, objects: Array.from({ length: count }, () => ({ ...object })) };
        return {
          size: LargeSize,
          scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, { clip })]),
          props: { renderer: { ...PD.getDefaultValues(Canvas3DParams).renderer, sphereRaster } },
        };
      },
    ),
  ),
);

/**
 * PetWorld as the viewer draws it, zoomed, cut through its center by pixel (the caps of the cut spheres' solid
 * interiors) by a plane, an ellipsoid (the inside removed) or an infinite cone: the impostors (`off`, culled on the GPU;
 * `cpu`, on the CPU; `cpu-instance`, clipped by instance, culled on the CPU) and the raster (`auto`; its spheres
 * classified against the clip objects by distance bounds, per pixel only where they reach a surface); `-ghost`: with
 * the ghost of what the clip object removes (`clipGhost`), on demand
 */
const ClipPostsynObjects = {
  plane: { type: 'plane' as const, scale: Vec3.create(1, 1, 1) },
  ellipsoid: { type: 'sphere' as const, scale: Vec3.create(600, 300, 400) },
  cone: { type: 'infiniteCone' as const, scale: Vec3.create(1, 1, 1) },
};
export const ClipPostsynCases: BenchmarkCase<RenderInput>[] = (['plane', 'ellipsoid', 'cone'] as const).flatMap(
  (shape) =>
    (['off', 'auto', 'cpu', 'cpu-instance', 'off-ghost', 'auto-ghost', 'cpu-ghost'] as const).map((name) => {
      const ghost = name.endsWith('-ghost');
      const variant = ghost ? (name.slice(0, -'-ghost'.length) as 'off' | 'auto' | 'cpu') : name;
      return renderCase(
        `clip/postsyn-zoom-${shape}-${name}`,
        `postsynstubby_2 as the viewer draws it, zoomed to 0.3 of the distance, cut by ${shape === 'ellipsoid' ? 'an' : 'a'} ${shape} ${variant === 'cpu-instance' ? 'by instance' : 'by pixel'}; ${variant === 'cpu' || variant === 'cpu-instance' ? 'culled on the CPU, the impostors' : `sphere raster ${variant}`}${ghost ? ', with the ghost of what it removes' : ''}`,
        async () => {
          const list = await getPostsyn();
          const { type, scale } = ClipPostsynObjects[shape];
          // the particles' center (their bounds'): the model is not about the origin
          const min = Vec3.create(Infinity, Infinity, Infinity);
          const max = Vec3.create(-Infinity, -Infinity, -Infinity);
          for (let i = 0; i < list.count; ++i) {
            const p = Vec3.fromArray(Vec3(), list.coordinates, i * 3);
            Vec3.min(min, min, p);
            Vec3.max(max, max, p);
          }
          const center = Vec3.scale(Vec3(), Vec3.add(Vec3(), min, max), 0.5);
          const clip = {
            variant: variant === 'cpu-instance' ? ('instance' as const) : ('pixel' as const),
            objects: [
              {
                type,
                invert: false,
                position: center,
                rotation: { axis: Vec3.create(1, 0, 0), angle: 30 },
                scale,
                transform: Mat4.identity(),
              },
            ],
          };
          return {
            scene: (g) =>
              Promise.all([
                postsynViewerRepr(g, list, {
                  clip,
                  ...((variant === 'cpu' || variant === 'cpu-instance') && { gpuCulling: false }),
                }),
              ]),
            zoom: 0.3,
            yieldFrames: true,
            props: {
              hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true },
              renderer: {
                ...PD.getDefaultValues(Canvas3DParams).renderer,
                sphereRaster: variant === 'auto' ? 'auto' : 'off',
              },
              ...(ghost && { clipGhost: { name: 'on' as const, params: PD.getDefaultValues(ClipGhostParams) } }),
            },
          };
        },
      );
    }),
);

/**
 * 60 copies of 7qpd as cartoon and as Gaussian surfaces (meshes), zoomed to 0.5 and cut through the center by a plane
 * by pixel: instanced, culled on the CPU, and merged, culled on the GPU; the culling of what the plane clips whole on
 * geometry other than spheres (`clip-culling.md`), on demand
 */
export const ClipMeshCases: BenchmarkCase<RenderInput>[] = (['cartoon', 'gaussian'] as const).flatMap((repr) =>
  ([false, true] as const).map((merged) =>
    renderCase(
      `clip/assembly-60-${repr}-plane${merged ? '-merged' : ''}`,
      `60 copies of 7qpd as ${repr === 'cartoon' ? 'cartoon' : 'Gaussian surfaces'}, zoomed to 0.5, cut through the center by a plane by pixel; ${merged ? 'merged, culled on the GPU' : 'instanced, culled on the CPU'}`,
      async () => {
        const s = await getRootAssembly();
        const clip = {
          variant: 'pixel' as const,
          objects: [
            {
              type: 'plane' as const,
              invert: false,
              position: Vec3.copy(Vec3(), s.boundary.sphere.center),
              rotation: { axis: Vec3.create(1, 0, 0), angle: 30 },
              scale: Vec3.create(1, 1, 1),
              transform: Mat4.identity(),
            },
          ],
        };
        const props = { clip, ...(merged && { renderMerged: true, gpuCulling: true }) };
        return {
          zoom: 0.5,
          yieldFrames: true,
          scene: (g) =>
            Promise.all([
              repr === 'cartoon'
                ? structureRepr(g, CartoonRepresentationProvider, s, props)
                : structureRepr(g, GaussianSurfaceRepresentationProvider, s, props),
            ]),
          props: { hiZ: { ...PD.getDefaultValues(HiZParams), enabled: true } },
        };
      },
    ),
  ),
);

/** Cel shading on (`uCelShaded`, a runtime branch of the lighting), fragment-bound and small, on demand */
export const CelCases: BenchmarkCase<RenderInput>[] = [
  renderCase('cel/globule-300k-spacefill-large', '300k atoms as spacefill, 2560x1440, cel shaded', async () => {
    const s = await getGlobule();
    return {
      size: LargeSize,
      scene: (g) => Promise.all([structureRepr(g, SpacefillRepresentationProvider, s, { celShaded: true })]),
    };
  }),
  renderCase('cel/7qpd-cartoon', 'PDB entry 7qpd as cartoon, cel shaded', async () => {
    const s = await getStructure();
    return { scene: (g) => Promise.all([structureRepr(g, CartoonRepresentationProvider, s, { celShaded: true })]) };
  }),
];
