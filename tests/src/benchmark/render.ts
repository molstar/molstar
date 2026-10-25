/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Canvas3D, Canvas3DContext, Canvas3DParams, type Canvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { PostprocessingParams } from '@molstar/graphics/canvas3d/passes/postprocessing';
import { IlluminationParams } from '@molstar/graphics/canvas3d/passes/illumination';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { TimerResult } from '@molstar/graphics/gl/timer';
import type { ResourceByteCounts } from '@molstar/graphics/gl/stats';
import { Interval } from '@molstar/core/data/int/interval';
import { Color } from '@molstar/core/util/color';
import { StructureElement } from '@molstar/model/model/structure';
import { Theme } from '@molstar/graphics/theme/theme';
import { getThemeRegistries } from '../visual/scenes/structure.js';
import type { Representation, RepresentationProvider } from '@molstar/graphics/repr/representation';
import { BallAndStickRepresentationProvider } from '@molstar/graphics/repr/structure/representation/ball-and-stick';
import { CartoonRepresentationProvider } from '@molstar/graphics/repr/structure/representation/cartoon';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { MolecularSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/molecular-surface';
import { SpacefillRepresentationProvider } from '@molstar/graphics/repr/structure/representation/spacefill';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { isWebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { type GraphicsRenderObject, isMergedRenderObject } from '@molstar/graphics/gl/render-object';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { now } from '@molstar/core/util/now';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { OutlineParams } from '@molstar/graphics/canvas3d/passes/outline';
import { SsaoParams } from '@molstar/graphics/canvas3d/passes/ssao';
import { HiZParams } from '@molstar/graphics/canvas3d/passes/hi-z';
import { Structure } from '@molstar/model/model/structure';
import { CIF } from '@molstar/io/reader/cif';
import { Task } from '@molstar/core/task';
import { createParticleListFromMmcifAssembly, getAssemblyIdsFromMmcif } from '@molstar/model/formats/particles/mmcif';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';
import { ParticleTargetRepresentationProvider } from '@molstar/graphics/repr/particles/representation/target/representation';
import { ParticleTargetRepresentationParams } from '@molstar/graphics/repr/particles/representation/target/common';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { setQuadIndexing } from '@molstar/graphics/gl/quad-indices';
import { meshBox, shapeRepr } from '../visual/scenes/primitives.js';
import { atomStructure, structureRepr } from '../visual/scenes/structure.js';
import { type BenchmarkCase, BenchmarkEnv, type PathRun, waitForGpu } from './framework.js';
import { exampleStructure, globuleStructure } from './structures.js';

/** Drawing buffer of the render cases, at pixel scale 1 */
const Size = { width: 1280, height: 720 };
/** Camera rotation per frame of rotating cases, about the camera's up axis */
const RotationStep = (2 * Math.PI) / 180;

/** Frames per run, each forced (all passes drawn), issued back to back */
const Frames = 30;

export type RenderInput = {
  readonly scene: (graphics: GraphicsContext) => Promise<Representation.Any[]>;
  readonly props?: Partial<Canvas3DProps>;
  /** the context's transparency mode (`Canvas3DContext.Props`), the default restored with the scene's disposal */
  readonly transparency?: Canvas3DContext.Props['transparency'];
  /** rotates the camera by `RotationStep` per frame, as an interactive rotation (culling and globals change) */
  readonly rotate?: boolean;
  /**
   * yields to the event loop after each frame (a task, as an app's frames do): readbacks complete only then (Hi-Z's;
   * WebGL's fences signal between tasks), the frames in one task would never see them; counted in `wait`
   */
  readonly yieldFrames?: boolean;
  /** drawing buffer, `Size` by default (larger for full-screen passes, so that their GPU times stand out of the noise) */
  readonly size?: { readonly width: number; readonly height: number };
  /** moves the camera towards the target, to this fraction of its distance, after the reset */
  readonly zoom?: number;
  /** false: the scene's quads (spheres, points) drawn unindexed (`setQuadIndexing`), for comparisons */
  readonly quadIndexing?: boolean;
  /** everything highlighted: the marking passes drawn per frame besides the scene's */
  readonly highlight?: boolean;
  /** a hover sweeping across the objects: per frame another representation highlighted whole, the last one not */
  readonly sweep?: boolean;
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
};

type ByteKinds = { geometry: number; instances: number; layers: number };
/** of the render objects, merged objects' members included; `merged` the merged objects' own copies */
type CpuBytes = { objects: ByteKinds; merged: ByteKinds };

const InstanceKeys = /^(aTransform|aInstance|aSegment|instanceGrid)/;

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
      const kind = k.startsWith('t') ? 'layers' : InstanceKeys.test(k) ? 'instances' : 'geometry';
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

async function getScene(id: string, graphics: GraphicsContext, input: RenderInput) {
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
        multiSample: { ...PD.getDefaultValues(Canvas3DParams).multiSample, mode: 'off' },
      },
      { warmUp: false },
    );
  const canvas3d = input.bundles?.record ? withBundleRecording(graphics, create) : create();
  canvas3d.handleResize();
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
  const setup: SceneSetup = { reprsMs: t1 - t0, commitMs: t2 - t1 };
  if (input.zoom) {
    const { position, target } = canvas3d.camera.getSnapshot();
    const offset = Vec3.scale(Vec3(), Vec3.sub(Vec3(), position, target), input.zoom);
    canvas3d.camera.setState({ position: Vec3.add(Vec3(), target, offset) }, 0);
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
async function renderFrames(id: string, graphics: GraphicsContext, input: RenderInput): Promise<PathRun> {
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
  // pipelines compiling in the background (WebGPU) would skip draws
  await graphics.resources.pending.whenReady();
  canvas3d.requestDraw();
  canvas3d.tick(now(), { isSynchronous: true });
  await waitForGpu(graphics);
  if (current && current.setup.firstFrameMs === undefined) {
    current.setup.firstFrameMs = performance.now() - first;
    current.setup.bytes = { ...graphics.resources.getByteCounts() };
    current.setup.cpuBytes = getCpuBytes(current.reprs);
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
  const start = performance.now();
  try {
    for (let i = 0; i < Frames; ++i) {
      if (input.rotate) rotateCamera(canvas3d);
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
      canvas3d.requestDraw();
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
  const total = performance.now() - start;
  sub.unsubscribe();
  const notes: string[] = [];
  if (input.setupNotes && current) {
    const { reprsMs, commitMs, firstFrameMs } = current.setup;
    notes.push(
      `setup: representations ${reprsMs.toFixed(0)} ms, add and commit ${commitMs.toFixed(0)} ms, first frame ${firstFrameMs?.toFixed(0)} ms`,
    );
  }
  if (input.byteNotes && current?.setup.bytes) {
    const b = current.setup.bytes;
    const mb = (n: number) => (n / 1048576).toFixed(1);
    notes.push(
      `GPU bytes: attributes ${mb(b.attribute)} MB, elements ${mb(b.elements)} MB, textures ${mb(b.texture)} MB, uniforms ${mb(b.uniform)} MB`,
    );
    const c = current.setup.cpuBytes;
    if (c) {
      const kinds = (k: ByteKinds) =>
        `geometry ${mb(k.geometry)}, instances ${mb(k.instances)}, layers ${mb(k.layers)}`;
      notes.push(`CPU bytes: objects ${kinds(c.objects)} MB; merged copies ${kinds(c.merged)} MB`);
    }
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
    if (input.rotate) rotateCamera(canvas3d);
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

function rotateCamera(canvas3d: Canvas3D) {
  const { position, target, up } = canvas3d.camera.getSnapshot();
  const rotation = Mat4.fromRotation(Mat4(), RotationStep, up);
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
/** `count` chains of 10 atoms on a grid: a unit each, the visuals of one representation share its material */
function chainGrid(count: number) {
  const atoms = 10;
  const n = Math.ceil(Math.sqrt(count));
  const x = new Float32Array(count * atoms);
  const y = new Float32Array(count * atoms);
  const z = new Float32Array(count * atoms);
  const chainOf = new Int32Array(count * atoms);
  for (let c = 0; c < count; ++c) {
    for (let a = 0; a < atoms; ++a) {
      const i = c * atoms + a;
      x[i] = (c % n) * 8 + a * 0.5;
      y[i] = Math.floor(c / n) * 8;
      z[i] = (a % 2) * 1.5;
      chainOf[i] = c;
    }
  }
  const chains = Array.from({ length: count }, (_, c) => `C${c}`);
  return atomStructure('chain-grid', x, y, z, chains, chainOf);
}

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
const getGlobule = () => (globule ??= globuleStructure(300_000));

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
 * spacefill, the camera rotating: the CPU culls every unit's instances per frame (instance grids)
 */
export const MergingInstanceCases: BenchmarkCase<RenderInput>[] = [100, 1000].flatMap((copies) =>
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
            rotate,
          };
        },
      ),
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

let postsyn: Promise<ParticleList> | undefined;
/**
 * `postsynstubby_2-no_bonds.bcif` from the repository root (not checked in: a PetWorld model, 3.68M atoms, 517k
 * operators) as the particle list the viewer loads (`mmcif_particles`): the molecule types the targets, instanced at
 * their particles
 */
const getPostsyn = () =>
  (postsyn ??= (async () => {
    const response = await fetch('/postsynstubby_2-no_bonds.bcif');
    if (!response.ok) {
      throw new Error(`postsynstubby_2-no_bonds.bcif not found in the repository root (${response.status})`);
    }
    const parsed = await CIF.parseBinary(new Uint8Array(await response.arrayBuffer())).run();
    if (parsed.isError) throw new Error(parsed.message);
    const cif = parsed.result;
    const assemblyId = getAssemblyIdsFromMmcif(cif)[0];
    return Task.create('particles', (ctx) =>
      createParticleListFromMmcifAssembly(ctx, cif, { assemblyId, variant: 'petworld' }),
    ).run();
  })());

/** The particle targets as the viewer shows them, with the mesoscale preset's sphere settings */
type LodLevels = (typeof Spheres.LodLevelsPresets)['quality'];

async function postsynRepr(
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
 * GPU culling on a large PetWorld model (on demand, needs `postsynstubby_2-no_bonds.bcif` in the repository root):
 * its particle targets as spacefill (merged by `auto`), culled on the CPU or the GPU, still or rotating
 */
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
 * cost of the per-pixel clip tests, their loops bounded by a uniform (`uClipObjectCount`), on demand
 */
export const ClipCases: BenchmarkCase<RenderInput>[] = [1, 3].map((count) =>
  renderCase(
    `clip/globule-300k-spacefill-large-${count}`,
    `300k atoms as spacefill, 2560x1440, ${count} clip object${count > 1 ? 's' : ''} cutting nothing (pixel)`,
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
      };
    },
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
