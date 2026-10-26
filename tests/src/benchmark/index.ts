/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { setTimingMode } from '@molstar/core/util/debug';
import { Canvas3DContext } from '@molstar/graphics/canvas3d/canvas3d';
import { createWebGLBackend, WebGL2Backend, WebGLBackend } from '@molstar/graphics/gl/webgl/backend';
import { isWebGLContext } from '@molstar/graphics/gl/webgl/context';
import { createWebGPUBackend } from '@molstar/graphics/gl/webgpu/backend';
import { isWebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { matchesAny } from '../visual/framework/test-case.js';
import { ColorSmoothingCases } from './color-smoothing.js';
import { type BenchmarkCase, BenchmarkEnv, type BenchmarkPath, type PathResult, runPath } from './framework.js';
import { GaussianSurfaceCases, GaussianSurfaceCutoffCases } from './gaussian-surface.js';
import { MarchingCubesCases } from './marching-cubes.js';
import {
  DrawTableCases,
  MergingCases,
  MergingMeshCases,
  MergingUniqueMeshCases,
  MergingAssemblyCases,
  MergingInstanceCases,
  MergingSingleInstanceCases,
  MergingPostsynCases,
  MergingPostsynSizeCases,
  MergingPostsynPixelCases,
  MergingPostsynChainCases,
  MergingPostsynPixelModeCases,
  MergingPostsynUpdateCases,
  UnindexedQuadCases,
  ClipCases,
  CelCases,
  RenderBundleCases,
  RenderCases,
  disposeScene,
} from './render.js';
import { DirectVolumeCases } from './direct-volume.js';
import { DrawOverheadCases } from './draw-overhead.js';
import { PassCases } from './passes.js';
import { RepresentationCases } from './representations.js';

/**
 * Benchmark page: `backend` (`webgl2`, which also runs the CPU paths, `webgl2-no-ubo`, `webgl1` or `webgpu`), `paths` (default: all of the
 * backend), `filter` (case id globs), `runs` and `warmup`, `dataBinding` (WebGPU: the theme layers' binding limit, bytes).
 * Cases run sequentially; results in `window.__benchmark`.
 */

const AllCases: BenchmarkCase<any>[] = [
  ...MarchingCubesCases,
  ...GaussianSurfaceCases,
  ...RepresentationCases,
  ...ColorSmoothingCases,
  ...RenderCases,
];

/** Cases run only when a filter selects them (measurements for a decision, not the default run) */
const OnDemandCases: BenchmarkCase<any>[] = [
  ...GaussianSurfaceCutoffCases,
  ...RenderBundleCases,
  ...PassCases,
  ...DrawTableCases,
  ...DrawOverheadCases,
  ...MergingCases,
  ...MergingMeshCases,
  ...MergingUniqueMeshCases,
  ...MergingAssemblyCases,
  ...MergingInstanceCases,
  ...MergingSingleInstanceCases,
  ...MergingPostsynCases,
  ...MergingPostsynSizeCases,
  ...MergingPostsynPixelCases,
  ...MergingPostsynChainCases,
  ...MergingPostsynPixelModeCases,
  ...MergingPostsynUpdateCases,
  ...UnindexedQuadCases,
  ...ClipCases,
  ...CelCases,
  ...DirectVolumeCases,
];

const BackendPaths: { [backend: string]: BenchmarkPath[] } = {
  webgl2: ['cpu', 'webgl'],
  'webgl2-no-ubo': ['webgl'],
  webgl1: ['webgl'],
  webgpu: ['webgpu'],
};

const Backends = {
  webgl2: WebGL2Backend,
  'webgl2-no-ubo': createWebGLBackend('webgl2', { uniformBlocks: false }),
  webgl1: WebGLBackend,
  webgpu: createWebGPUBackend(),
};

type BenchmarkState = {
  done: boolean;
  info: { [k: string]: string };
  results: PathResult[];
  images?: { [id: string]: string };
  error?: string;
};
declare global {
  interface Window {
    __benchmark: BenchmarkState;
  }
}

const state: BenchmarkState = { done: false, info: {}, results: [] };
window.__benchmark = state;

const log = (line: string) => {
  document.getElementById('log')!.textContent += `${line}\n`;
};

async function createContext(backend: string, dataBinding?: number) {
  const canvas = document.createElement('canvas');
  canvas.width = canvas.height = 64;
  document.body.appendChild(canvas);
  const descriptor =
    backend === 'webgpu' && dataBinding
      ? createWebGPUBackend({ maxDataBindingSize: dataBinding })
      : Backends[backend as keyof typeof Backends];
  const context = await Canvas3DContext.create(
    canvas,
    new DefaultAssetManager(),
    [descriptor],
    {},
    { pixelScale: 1, resolutionMode: 'native' },
  );
  BenchmarkEnv.context = context;
  return context.graphics;
}

async function run() {
  const params = new URL(window.location.href).searchParams;
  const backend = params.get('backend') ?? 'webgl2';
  const allPaths = BackendPaths[backend];
  if (!allPaths) throw new Error(`unknown backend '${backend}'`);
  const paths =
    (params.get('paths')?.split(',') as BenchmarkPath[] | undefined)?.filter((p) => allPaths.includes(p)) ?? allPaths;
  const filter = params.get('filter')?.split(',') ?? [];
  const options = { runs: Number(params.get('runs') ?? 7), warmup: Number(params.get('warmup') ?? 2) };

  if (params.get('gpuTiming')) {
    // GPU times of the passes of render cases (adds timer queries and CPU work to the frames)
    setTimingMode(true);
    BenchmarkEnv.gpuTiming = true;
  }
  if (params.get('images')) state.images = BenchmarkEnv.images = {};
  if (params.get('drawLists')) BenchmarkEnv.drawLists = params.get('drawLists') === '1';
  const graphics = await createContext(backend, Number(params.get('dataBinding') ?? 0));
  if (isWebGLContext(graphics)) {
    const ext = graphics.gl.getExtension('WEBGL_debug_renderer_info');
    state.info.renderer = ext ? graphics.gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : 'unknown';
  } else if (isWebGPUContext(graphics)) {
    const info = (graphics.device as GPUDevice & { adapterInfo?: GPUAdapterInfo }).adapterInfo;
    state.info.renderer = info ? `${info.vendor} ${info.architecture} ${info.description}`.trim() : 'unknown';
  }
  state.info.backend = backend;
  state.info.userAgent = navigator.userAgent;
  log(`${backend}: ${state.info.renderer}`);

  for (const c of filter.length > 0 ? [...AllCases, ...OnDemandCases] : AllCases) {
    if (filter.length > 0 && !matchesAny(c.id, filter)) continue;
    const casePaths = paths.filter((p) => c.paths[p]);
    if (casePaths.length === 0) continue;
    // yield between cases (CPU paths await nothing asynchronous)
    await new Promise((r) => setTimeout(r, 0));
    const input = await c.setup();
    for (const path of casePaths) {
      const result = await runPath(c, path, graphics, input, options);
      if (path === casePaths[casePaths.length - 1]) disposeScene();
      state.results.push(result);
      log(
        `${c.id} ${path}: ${result.status === 'ok' ? `median ${result.median!.toFixed(1)} ms` : `${result.status} ${result.reason}`}`,
      );
    }
  }
}

// after the page has loaded, the runs block the event loop for long stretches
window.addEventListener('load', () => {
  run()
    .catch((e) => {
      state.error = e instanceof Error ? `${e.message}\n${e.stack}` : String(e);
      log(state.error);
    })
    .finally(() => {
      state.done = true;
    });
});
