/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Canvas3D, Canvas3DParams, type Canvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { PostprocessingParams } from '@molstar/graphics/canvas3d/passes/postprocessing';
import { IlluminationParams } from '@molstar/graphics/canvas3d/passes/illumination';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { TimerResult } from '@molstar/graphics/gl/timer';
import type { Representation } from '@molstar/graphics/repr/representation';
import { BallAndStickRepresentationProvider } from '@molstar/graphics/repr/structure/representation/ball-and-stick';
import { CartoonRepresentationProvider } from '@molstar/graphics/repr/structure/representation/cartoon';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { MolecularSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/molecular-surface';
import { SpacefillRepresentationProvider } from '@molstar/graphics/repr/structure/representation/spacefill';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { isWebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { now } from '@molstar/core/util/now';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { OutlineParams } from '@molstar/graphics/canvas3d/passes/outline';
import { SsaoParams } from '@molstar/graphics/canvas3d/passes/ssao';
import { Structure } from '@molstar/model/model/structure';
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

type RenderInput = {
  readonly scene: (graphics: GraphicsContext) => Promise<Representation.Any[]>;
  readonly props?: Partial<Canvas3DProps>;
  /** rotates the camera by `RotationStep` per frame, as an interactive rotation (culling and globals change) */
  readonly rotate?: boolean;
  /** drawing buffer, `Size` by default (larger for full-screen passes, so that their GPU times stand out of the noise) */
  readonly size?: { readonly width: number; readonly height: number };
  /** moves the camera towards the target, to this fraction of its distance, after the reset */
  readonly zoom?: number;
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
const LargeSize = { width: 2560, height: 1440 };

/** The scene of the case rendered last, kept across its runs */
let current: { id: string; canvas3d: Canvas3D; reprs: Representation.Any[] } | undefined;

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

  // progressive multi-sampling would make frames differ, picking warm-up would run in between
  const create = () =>
    Canvas3D.create(
      context,
      { ...input.props, multiSample: { ...PD.getDefaultValues(Canvas3DParams).multiSample, mode: 'off' } },
      { warmUp: false },
    );
  const canvas3d = input.bundles?.record ? withBundleRecording(graphics, create) : create();
  canvas3d.handleResize();
  const reprs = await input.scene(graphics);
  for (const r of reprs) canvas3d.add(r);
  canvas3d.commit(true);
  canvas3d.requestCameraReset({ durationMs: 0 });
  canvas3d.tick(now(), { isSynchronous: true });
  if (input.zoom) {
    const { position, target } = canvas3d.camera.getSnapshot();
    const offset = Vec3.scale(Vec3(), Vec3.sub(Vec3(), position, target), input.zoom);
    canvas3d.camera.setState({ position: Vec3.add(Vec3(), target, offset) }, 0);
    canvas3d.tick(now(), { isSynchronous: true });
  }
  current = { id, canvas3d, reprs };
  return canvas3d;
}

export function disposeScene() {
  if (!current) return;
  current.canvas3d.dispose();
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
  // pipelines compiling in the background (WebGPU) would skip draws
  await graphics.resources.pending.whenReady();
  canvas3d.requestDraw();
  canvas3d.tick(now(), { isSynchronous: true });
  await waitForGpu(graphics);

  let submit = 0;
  let calls = 0;
  let drawn = 0;
  canvas3d.notifyDidDraw = true;
  const sub = canvas3d.didDraw.subscribe(() => drawn++);
  drawn = 0;
  bundleStats.recorded = bundleStats.failed = 0;
  const frustum = Frustum3D();
  let visibilityMs = 0;
  let visibleKey: string | undefined;
  let visibleChanges = 0;
  const intersects = Frustum3D.intersectsSphere3D;
  if (input.bundles?.noCulling) Frustum3D.intersectsSphere3D = () => true;
  const start = performance.now();
  try {
    for (let i = 0; i < Frames; ++i) {
      if (input.rotate) rotateCamera(canvas3d);
      resetDrawCalls(graphics);
      const t0 = performance.now();
      canvas3d.requestDraw();
      canvas3d.tick(now(), { isSynchronous: true });
      submit += performance.now() - t0;
      calls = drawCalls(graphics);
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
  if (input.bundles?.record) notes.push(`opaque pass recorded ${bundleStats.recorded}, failed ${bundleStats.failed}`);
  if (input.bundles?.visibility) {
    const visible = visibleKey?.split('').filter((v) => v === '1').length ?? 0;
    notes.push(
      `visible set: ${(visibilityMs / Frames).toFixed(3)} ms per frame, changed in ${visibleChanges} of ${Frames - 1} frames, ${visible} of ${visibleKey?.length ?? 0} objects visible at the end`,
    );
  }
  if (drawn !== Frames) throw new Error(`${drawn} of ${Frames} frames drawn`);
  const gpuPasses = BenchmarkEnv.gpuTiming ? await getGpuPasses(graphics) : undefined;
  return {
    stages: [
      { label: 'submit', ms: submit / Frames },
      { label: 'wait', ms: (total - submit) / Frames },
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

function renderCase(id: string, description: string, setup: () => Promise<RenderInput>): BenchmarkCase<RenderInput> {
  const path = (graphics: GraphicsContext, input: RenderInput) => renderFrames(id, graphics, input);
  return { id, description, setup, paths: { webgl: path, webgpu: path } };
}

let structure: Promise<Structure> | undefined;
const getStructure = () => (structure ??= exampleStructure('7qpd.fw2.cif'));
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
export const RenderBundleCases: BenchmarkCase<any>[] = [
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
