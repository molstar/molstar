/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getRegularGrid3dDelta, type RegularGrid3d } from '@molstar/core/math/geometry/common';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { RuntimeContext } from '@molstar/core/task';
import { UUID, ValueCell } from '@molstar/core/util';
import { objectForEach } from '@molstar/core/util/object';
import {
  type ComputeContext,
  KernelCode,
  type KernelRequest,
  type KernelSchema,
  StorageBufferSpec,
} from '../compute.js';
import type { GraphicsContext } from '../context.js';
import {
  DefineSpec,
  type RenderableSchema,
  TextureSpec,
  type UnboxedValues,
  UniformSpec,
  type Values,
} from '../renderable/schema.js';
import type { Texture } from '../texture.js';
import { getDispatchSize } from './prefix-sum.js';

/**
 * Whether `computable` (from `createGrid3dComputeRenderable`) runs on the GPU of `graphics`: compute kernels for specs
 * with WGSL code, or fragment shaders (WebGL); without `computable`, whether fragment shaders are available
 */
export function canComputeGrid3dOnGPU(graphics?: GraphicsContext, computable?: Grid3dComputable<any>): boolean {
  if (!graphics) return false;
  return computable ? computable.isSupported(graphics) : !!getFragmentImplementation(graphics);
}

/** The computation of a spec with fragment shaders */
export type Grid3dFragmentCompute<P> = (
  ctx: RuntimeContext,
  graphics: GraphicsContext,
  grid: RegularGrid3d,
  params: P,
) => Promise<Float32Array>;

/** grid3d with fragment shaders, an implementation a backend registers (WebGL, `gl/webgl/compute/grid3d`) */
export interface Grid3dFragmentImplementation {
  /** whether it runs on `graphics` (its backend) */
  isSupported(graphics: GraphicsContext): boolean;
  create<S extends RenderableSchema, P, CS>(
    id: string,
    spec: Grid3DComputeRenderableSpec<S, P, CS>,
  ): Grid3dFragmentCompute<P>;
}

const fragmentImplementations: Grid3dFragmentImplementation[] = [];

/** Called by backends that implement it, when they create a context */
export function registerGrid3dFragmentImplementation(implementation: Grid3dFragmentImplementation) {
  if (!fragmentImplementations.includes(implementation)) fragmentImplementations.push(implementation);
}

function getFragmentImplementation(graphics: GraphicsContext) {
  return fragmentImplementations.find((i) => i.isSupported(graphics));
}

/** Computes the grid values (axis order fast to slow z, y, x) */
export interface Grid3dComputable<P> {
  (ctx: RuntimeContext, graphics: GraphicsContext, grid: RegularGrid3d, params: P): Promise<Float32Array>;
  isSupported(graphics: GraphicsContext): boolean;
}

export interface Grid3DComputeRenderableSpec<S extends RenderableSchema, P, CS> {
  schema: S;
  // indicate which params are loop bounds for WebGL1 compat
  loopBounds?: (keyof S)[];
  utilCode?: string;
  mainCode: string;
  returnCode: string;
  /**
   * The code in WGSL, for compute kernels (WebGPU): uniforms are `uParams.<name>`, textures are read with `textureLoad`,
   * `xyz` is the cell position and `current` the accumulated value (cumulative)
   */
  wgsl?: { utilCode?: string; mainCode: string; returnCode: string };

  values(params: P, grid: RegularGrid3d): UnboxedValues<S>;

  cumulative?: {
    states(params: P): CS[];
    update(params: P, state: CS, values: Values<S>): void;
    // call gl.readPixes every 'yieldPeriod' states to split the computation
    // into multiple parts, if not set, the computation will be synchronous
    yieldPeriod?: number;
  };
}

export function createGrid3dComputeRenderable<S extends RenderableSchema, P, CS>(
  spec: Grid3DComputeRenderableSpec<S, P, CS>,
): Grid3dComputable<P> {
  const id = UUID.create22();

  objectForEach(spec.schema, (u) => {
    if (u.type === 'uniform' && u.kind.endsWith('[]')) throw new Error('array uniforms are not supported');
  });

  const kernel = spec.wgsl ? createGrid3dKernel(id, spec) : undefined;
  const useKernel = (graphics: GraphicsContext) => !!kernel && !!graphics.compute?.isKernelSupported(kernel.code);
  // the fragment computations are created on first use (backends register their implementation with a context)
  const fragmentComputes = new Map<Grid3dFragmentImplementation, Grid3dFragmentCompute<P>>();

  const compute = async (ctx: RuntimeContext, graphics: GraphicsContext, grid: RegularGrid3d, params: P) => {
    if (useKernel(graphics)) {
      return kernel!.run(graphics as GraphicsContext & { compute: ComputeContext }, grid, params);
    }
    const implementation = getFragmentImplementation(graphics);
    if (!implementation) throw new Error('grid3d needs compute kernels (for specs with WGSL code) or fragment shaders');
    let fragmentCompute = fragmentComputes.get(implementation);
    if (!fragmentCompute) {
      fragmentCompute = implementation.create(id, spec);
      fragmentComputes.set(implementation, fragmentCompute);
    }
    return fragmentCompute(ctx, graphics, grid, params);
  };
  return Object.assign(compute, {
    isSupported: (graphics: GraphicsContext) => useKernel(graphics) || !!getFragmentImplementation(graphics),
  });
}

/** The WGSL kernel of a spec: a thread per cell, the value (or the sum, cumulative) in a storage buffer */
function getGrid3dKernelSource(wgsl: NonNullable<Grid3DComputeRenderableSpec<any, any, any>['wgsl']>) {
  return /* wgsl */ `
${wgsl.utilCode ?? ''}

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) workgroupCount: vec3<u32>) {
    let offset = id.x + id.y * workgroupCount.x * WorkgroupSizeX;
    let dim = vec3<u32>(uParams.uDimensions);
    if (offset >= dim.x * dim.y * dim.z) {
        return;
    }
    // axis order fast to slow z, y, x, as on WebGL
    let k = offset % dim.z;
    let kk = offset / dim.z;
    let j = kk % dim.y;
    let i = kk / dim.y;
    let xyz = uParams.uMin + uParams.uDelta * vec3<f32>(f32(i), f32(j), f32(k));
    #if defined(dCumulative)
        let current = outValues[offset];
    #endif

    ${wgsl.mainCode}

    outValues[offset] = ${wgsl.returnCode};
}
`;
}

const KernelBaseSchema = {
  uDimensions: UniformSpec('v3'),
  uMin: UniformSpec('v3'),
  uDelta: UniformSpec('v3'),
  outValues: StorageBufferSpec('f32', 'read-write'),
  dCumulative: DefineSpec('boolean'),
};

/** Per-dimension dispatch limit of the WebGPU core */
const MaxWorkgroups = 65535;

function createGrid3dKernel<S extends RenderableSchema, P, CS>(
  id: string,
  spec: Grid3DComputeRenderableSpec<S, P, CS>,
) {
  const code = KernelCode(id, [64], { wgsl: getGrid3dKernelSource(spec.wgsl!) });

  // image values are uploaded into textures, kernels read textures
  const schema: { [k: string]: KernelSchema[string] } = { ...KernelBaseSchema };
  objectForEach(spec.schema, (s, k) => {
    if (s.type === 'uniform') schema[k as string] = s;
    else if (s.type === 'texture') schema[k as string] = TextureSpec('texture', s.format, s.dataType, s.filter);
    else throw new Error(`grid3d kernels support uniforms and textures, not '${s.type}' ('${String(k)}')`);
  });
  /** the kernel of `run`, see `ComputeContext.prepareKernels` */
  const kernels: readonly KernelRequest[] = [{ code, schema, defines: { dCumulative: !!spec.cumulative } }];

  async function run(graphics: GraphicsContext & { compute: ComputeContext }, grid: RegularGrid3d, params: P) {
    const { compute, resources } = graphics;
    // compiled in the background, not blocking the queue (and frames) at the first dispatch
    await compute.prepareKernels(kernels);
    const [nx, ny, nz] = grid.dimensions;
    const cellCount = nx * ny * nz;
    const output = compute.createBuffer(cellCount * 4, 'grid3d values');

    // the values as the spec has them (images for textures), and the kernel's (textures for images)
    const specValues: { [k: string]: ValueCell<any> } = {};
    objectForEach(spec.values(params, grid), (v, k) => (specValues[k as string] = ValueCell.create(v)));
    const textures: { [k: string]: { texture: Texture; version: number } } = {};
    const kernelValues: { [k: string]: ValueCell<any> } = {
      uDimensions: ValueCell.create(Vec3.clone(grid.dimensions)),
      uMin: ValueCell.create(Vec3.clone(grid.box.min)),
      uDelta: ValueCell.create(getRegularGrid3dDelta(grid)),
      outValues: ValueCell.create(output),
      dCumulative: ValueCell.create(!!spec.cumulative),
    };
    objectForEach(spec.schema, (s, k) => {
      const key = k as string;
      if (s.type !== 'texture') {
        kernelValues[key] = specValues[key];
        return;
      }
      const texture = resources.texture(s.kind as never, s.format, s.dataType, s.filter);
      texture.load(specValues[key].ref.value);
      textures[key] = { texture, version: specValues[key].ref.version };
      kernelValues[key] = ValueCell.create(texture);
    });
    const kernel = compute.createKernel(code, schema, kernelValues as never);
    const { x, y } = getDispatchSize(Math.max(1, Math.ceil(cellCount / 64)), MaxWorkgroups);
    const dispatch = () => {
      compute.beginComputePass({ label: 'grid3d' });
      kernel.dispatch(x, y);
      compute.endComputePass();
    };

    if (spec.cumulative) {
      for (const state of spec.cumulative.states(params)) {
        spec.cumulative.update(params, state, specValues as unknown as Values<S>);
        // images replaced by the state
        for (const key of Object.keys(textures)) {
          const t = textures[key];
          if (specValues[key].ref.version === t.version) continue;
          t.texture.load(specValues[key].ref.value);
          t.version = specValues[key].ref.version;
        }
        dispatch();
      }
    } else {
      dispatch();
    }

    const values = new Float32Array(await output.readAsync(0, cellCount * 4));
    kernel.destroy();
    output.destroy();
    for (const key of Object.keys(textures)) textures[key].texture.destroy();
    return values;
  }

  return { code, kernels, run };
}
