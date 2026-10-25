/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { createWebGPUContext, isWebGPUContext, requestWebGPUDevice } from '@molstar/graphics/gl/webgpu/context';
import { getWgslProgram } from '@molstar/graphics/gl/webgpu/shader/programs';
import type { Representation } from '@molstar/graphics/repr/representation';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { GaussianVolumeRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-volume';
import { IsosurfaceRepresentationProvider } from '@molstar/graphics/repr/volume/isosurface';
import { Volume } from '@molstar/model/model/volume';
import { type VisualTestCase } from '../framework/test-case.js';
import { ObjExporter } from '@molstar/geo-export-extension/obj-exporter';
import type { TextureKind } from '@molstar/graphics/gl/texture';
import { Box3D, Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { CopySchema, CopyShaderCode, QuadValues } from '@molstar/graphics/gl/compute/util';
import { createComputeRenderable } from '@molstar/graphics/gl/renderable';
import type { Values } from '@molstar/graphics/gl/renderable/schema';
import { RenderState } from '@molstar/graphics/gl/render-state';
import { RuntimeContext } from '@molstar/core/task';
import { sphericalCollocation } from '@molstar/alpha-orbitals-extension/collocation';
import type { AlphaOrbital, Basis, CubeGridInfo } from '@molstar/alpha-orbitals-extension/data-model';
import {
  canComputeAlphaOrbitalsOnGPU,
  gpuComputeAlphaOrbitalsDensityGridValues,
  gpuComputeAlphaOrbitalsGridValues,
} from '@molstar/alpha-orbitals-extension/gpu/compute';
import { exampleStructure, helixStructure, structureRepr } from '../scenes/structure.js';
import { blobVolume, segmentedBlobVolume, volumeRepr } from '../scenes/volume.js';
import { SegmentRepresentationProvider } from '@molstar/graphics/repr/volume/segment';
import type { Camera } from '@molstar/graphics/canvas3d/camera';
import type { DrawPass } from '@molstar/graphics/canvas3d/passes/draw';
import { HiZPass } from '@molstar/graphics/canvas3d/passes/hi-z';

/**
 * The context as the backend reports support without strict mode: the strict WebGPU backend of the visual tests
 * reports every program as supported (so that unported ones fail visibly), which hides a GPU path falling back to the
 * CPU because a support check fails
 */
function withoutStrictSupport(graphics: GraphicsContext): GraphicsContext {
  if (!isWebGPUContext(graphics)) return graphics;
  const wrapped = Object.create(graphics) as GraphicsContext;
  Object.defineProperty(wrapped, 'isProgramSupported', {
    value: (shaderCode: { name: string }) => getWgslProgram(shaderCode.name) !== undefined,
  });
  return wrapped;
}

function expectGeometry(name: string, repr: Representation.Any, type: string) {
  const types = repr.renderObjects.map((ro) => ro.type);
  repr.destroy();
  if (types.length === 0) return [`${name}: no render objects`];
  const other = types.filter((t) => t !== type);
  return other.length ? [`${name}: '${other[0]}' render objects, expected '${type}' (the GPU path)`] : [];
}

async function checkGpuPaths(graphics: GraphicsContext) {
  const g = withoutStrictSupport(graphics);
  const structure = await helixStructure(24);
  return [
    ...expectGeometry(
      'volume isosurface',
      await volumeRepr(g, IsosurfaceRepresentationProvider, blobVolume(32), {
        isoValue: Volume.IsoValue.absolute(0.5),
        tryUseGpu: true,
      }),
      'texture-mesh',
    ),
    ...expectGeometry(
      'volume segment',
      await volumeRepr(g, SegmentRepresentationProvider, segmentedBlobVolume(32), { segments: [1], tryUseGpu: true }),
      'texture-mesh',
    ),
    ...expectGeometry(
      'gaussian surface',
      await structureRepr(g, GaussianSurfaceRepresentationProvider, structure, { resolution: 0.5, tryUseGpu: true }),
      'texture-mesh',
    ),
    ...expectGeometry(
      'gaussian volume',
      await structureRepr(g, GaussianVolumeRepresentationProvider, structure, { resolution: 0.5 }),
      'direct-volume',
    ),
  ];
}

/** Vertices per unit (by group count), of indirectly drawn texture meshes the count read back */
function unitVertexCounts(repr: Representation.Any) {
  const counts = new Map<number, number>();
  for (const { values } of repr.renderObjects) {
    const v = values as unknown as {
      drawCount: { ref: { value: number } };
      uGroupCount: { ref: { value: number } };
      meta?: { ref: { value: { readVertexCount?: number } } };
    };
    counts.set(v.uGroupCount.ref.value, v.meta?.ref.value.readVertexCount ?? v.drawCount.ref.value);
  }
  return counts;
}

/**
 * The GPU Gaussian surface of a PDB entry against the CPU one, per unit: the vertex counts within 2% (resolutions or
 * densities that differ between the paths change them more, e.g. a resolution coarsened for a 2d texture limit)
 */
async function checkGaussianSurfaceUnits(graphics: GraphicsContext) {
  const structure = await exampleStructure('7qpd.fw2.cif');
  const props = { resolution: 1, tryUseGpu: true };
  const gpu = await structureRepr(
    withoutStrictSupport(graphics),
    GaussianSurfaceRepresentationProvider,
    structure,
    props,
  );
  await graphics.resources.pending.whenReady();
  const cpu = await structureRepr(graphics, GaussianSurfaceRepresentationProvider, structure, {
    ...props,
    tryUseGpu: false,
  });
  const gpuCounts = unitVertexCounts(gpu);
  const cpuCounts = unitVertexCounts(cpu);
  gpu.destroy();
  cpu.destroy();
  const errors: string[] = [];
  for (const [groups, expected] of cpuCounts) {
    const actual = gpuCounts.get(groups);
    if (actual === undefined) errors.push(`unit of ${groups} atoms: no GPU surface`);
    else if (Math.abs(actual - expected) > 0.02 * expected) {
      errors.push(`unit of ${groups} atoms: ${actual} vertices on the GPU, ${expected} on the CPU`);
    }
  }
  return errors;
}

/** Three atoms with a shell of each angular momentum 0–4, so all spherical harmonics of the GPU code are used */
function testBasis(): Basis {
  const centers: [number, number, number][] = [
    [0, 0, 0],
    [1.2, 0.4, -0.3],
    [-0.5, 1.1, 0.8],
  ];
  return {
    atoms: centers.map((center) => ({
      center,
      shells: [0, 1, 2, 3, 4].map((L) => ({
        angularMomentum: [L],
        coefficients: [[0.6, 0.4]],
        exponents: [1.5, 0.4],
      })),
    })),
  };
}

/** Coefficients of the 25 basis functions per atom, deterministic */
function testOrbital(seed: number, occupancy: number): AlphaOrbital {
  const alpha: number[] = [];
  let s = seed;
  for (let i = 0; i < 3 * 25; ++i) {
    s = (s * 1664525 + 1013904223) >>> 0;
    alpha.push((s / 4294967296) * 2 - 1);
  }
  return { energy: 0, occupancy, alpha };
}

function testGrid(basis: Basis, n: number): CubeGridInfo {
  const box = Box3D.create(Vec3.create(-2, -2, -2), Vec3.create(3, 3, 3));
  const size = Box3D.size(Vec3(), box);
  return {
    params: { basis, cutoffThreshold: 0, sphericalOrder: 'cca-reverse', boxExpand: 0, gridSpacing: [] },
    box,
    size,
    dimensions: Vec3.create(n, n, n),
    npoints: n * n * n,
    delta: Vec3.scale(Vec3(), size, 1 / (n - 1)),
  };
}

function compareValues(name: string, actual: Float32Array, expected: ArrayLike<number>) {
  if (actual.length !== expected.length) return [`${name}: ${actual.length} values, expected ${expected.length}`];
  let max = 0;
  for (let i = 0; i < expected.length; ++i) max = Math.max(max, Math.abs(expected[i]));
  for (let i = 0; i < expected.length; ++i) {
    if (Math.abs(actual[i] - expected[i]) > 1e-4 * max) {
      return [`${name}: value ${i} is ${actual[i]}, expected ${expected[i]} (of at most ${max})`];
    }
  }
  return [];
}

/**
 * Alpha orbitals and their density on the GPU (grid3d: WebGL fragment shader, or a compute kernel) against the CPU
 * collocation, for all angular momenta
 */
async function checkAlphaOrbitals(graphics: GraphicsContext) {
  if (!canComputeAlphaOrbitalsOnGPU(graphics)) return ['alpha orbitals not computed on the GPU'];
  const basis = testBasis();
  const grid = testGrid(basis, 12);
  const orbitals = [testOrbital(1, 2), testOrbital(2, 1)];
  const sync = RuntimeContext.Synchronous;

  const cpu = await Promise.all(orbitals.map((o) => sphericalCollocation(grid, o, sync)));
  const density = new Float32Array(grid.npoints);
  orbitals.forEach((o, k) => {
    for (let i = 0; i < density.length; ++i) density[i] += o.occupancy * cpu[k][i] * cpu[k][i];
  });

  return [
    ...compareValues('orbital', await gpuComputeAlphaOrbitalsGridValues(sync, graphics, grid, orbitals[0]), cpu[0]),
    ...compareValues(
      'density',
      await gpuComputeAlphaOrbitalsDensityGridValues(sync, graphics, grid, orbitals),
      density,
    ),
  ];
}

/**
 * `Texture.loadRows` on every backend: rows of an image changed in place after its upload (a font atlas's new glyphs)
 * written in place, the other rows kept, the image the texture's source after (copied from on a backend switch); an
 * image of another size loaded whole
 */
async function checkLoadRows(graphics: GraphicsContext) {
  // an odd width: rows are not aligned to anything
  const width = 37;
  const height = 24;
  const errors: string[] = [];
  for (const c of [
    { format: 'r', channels: 1 },
    { format: 'rgba', channels: 4 },
  ] as const) {
    const texture = graphics.resources.texture('image-uint8', c.format, 'ubyte', 'nearest');
    if (!texture.loadRows) {
      texture.destroy();
      errors.push(`${c.format}: no loadRows`);
      continue;
    }
    const rowLength = width * c.channels;
    const array = Uint8Array.from({ length: rowLength * height }, (_, i) => (i * 37) % 256);
    const image = { array, width, height };
    texture.load(image);
    // three ranges of rows changed, the last row among them
    const rows = [3, 7, 20, 21, 23, 24];
    for (let k = 0; k < rows.length; k += 2) {
      for (let i = rows[k] * rowLength; i < rows[k + 1] * rowLength; ++i) array[i] = (i * 11 + 101) % 256;
    }
    texture.loadRows(image, rows);
    const compare = async (name: string, expected: { array: Uint8Array; width: number; height: number }) => {
      const read = await graphics.readTextureData(texture);
      if (read.width !== expected.width || read.height !== expected.height) {
        errors.push(`${name}: read ${read.width}×${read.height}`);
        return;
      }
      for (let i = 0; i < expected.width * expected.height; ++i) {
        for (let k = 0; k < c.channels; ++k) {
          const value = expected.array[i * c.channels + k];
          if (read.array[i * 4 + k] !== value) {
            const row = Math.floor(i / expected.width);
            errors.push(`${name}: row ${row} texel ${i} channel ${k} is ${read.array[i * 4 + k]}, not ${value}`);
            return;
          }
        }
      }
    };
    await compare(`${c.format} rows`, image);
    if (texture.getSource?.().data !== image) errors.push(`${c.format}: the image is not the texture's source`);

    // another size: loaded whole
    const other = {
      array: Uint8Array.from({ length: rowLength * 2 * 8 }, (_, i) => (i * 7) % 256),
      width: 2 * width,
      height: 8,
    };
    texture.loadRows(other, [0, 2]);
    await compare(`${c.format} other size`, other);
    texture.destroy();
  }
  return errors;
}

/** Formats `readTextureData` reads on every backend; float textures where they are render targets (WebGL reads them) */
async function checkReadTextureData(graphics: GraphicsContext) {
  const width = 3;
  const height = 2;
  const cases: { kind: TextureKind; format: 'rgba' | 'r'; type: 'ubyte' | 'float'; channels: number }[] = [
    { kind: 'image-uint8', format: 'rgba', type: 'ubyte', channels: 4 },
    { kind: 'image-uint8', format: 'r', type: 'ubyte', channels: 1 },
  ];
  if (graphics.capabilities.floatRenderTarget) {
    cases.push({ kind: 'image-float32', format: 'rgba', type: 'float', channels: 4 });
  }
  const errors: string[] = [];
  for (const c of cases) {
    const n = width * height * c.channels;
    const data =
      c.type === 'float'
        ? Float32Array.from({ length: n }, (_, i) => i * 0.5 - 1)
        : Uint8Array.from({ length: n }, (_, i) => (i * 37) % 256);
    const texture = graphics.resources.texture(c.kind, c.format, c.type, 'nearest');
    texture.load({ array: data, width, height });
    const read = await graphics.readTextureData(texture);
    texture.destroy();
    const name = `${c.format} ${c.type}`;
    if (read.width !== width || read.height !== height) {
      errors.push(`${name}: read ${read.width}×${read.height}`);
      continue;
    }
    // data order, four values per texel (single-channel formats in r)
    for (let i = 0; i < width * height; ++i) {
      for (let k = 0; k < c.channels; ++k) {
        if (read.array[i * 4 + k] !== data[i * c.channels + k]) {
          errors.push(
            `${name}: texel ${i} channel ${k} is ${read.array[i * 4 + k]}, uploaded ${data[i * c.channels + k]}`,
          );
          break;
        }
      }
      if (errors.length) break;
    }
  }
  return errors;
}

/**
 * A viewport partly outside the target is cropped, as in GL and WebGPU core, not squashed into the target: a quad over
 * the bottom-left quarter of the viewport is drawn where the viewport places it
 */
async function checkViewportCropping(graphics: GraphicsContext) {
  const size = 64;
  const white = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
  white.load({ array: new Uint8Array(size * size * 4).fill(255), width: size, height: size });
  const target = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
  target.define(size, size);
  const framebuffer = graphics.resources.framebuffer();
  target.attachFramebuffer(framebuffer, 'color0');
  // a copy with a quad scale of its own (`QuadValues.uQuadScale` is shared by the quads of all passes)
  const values: Values<typeof CopySchema> = {
    ...QuadValues,
    uQuadScale: ValueCell.create(Vec2.create(0.5, 0.5)),
    tColor: ValueCell.create(white),
    uTexSize: ValueCell.create(Vec2.create(size, size)),
    dPackDepth: ValueCell.create(false),
  };
  const copy = createComputeRenderable(
    graphics.createComputeRenderItem('triangles', CopyShaderCode, { ...CopySchema }, values),
    values,
  );

  // pixels (top-left origin, as the viewports) the quad covers when cropped, and ones it covers only when squashed
  const cases = [
    { name: 'below and left of the target', viewport: [-16, 16, 64, 64], drawn: [8, 55], empty: [20, 43] },
    { name: 'above and right of the target', viewport: [32, -32, 64, 64], drawn: [56, 7], empty: [20, 43] },
  ] as const;
  const errors: string[] = [];
  const draw = (viewport: readonly number[]) => {
    graphics.disableScissor();
    graphics.beginRenderPass(framebuffer, { label: 'viewport cropping', clear: { color: [0, 0, 0, 1] } });
    graphics.setRenderState(RenderState.Default);
    graphics.setViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    copy.render();
    graphics.endRenderPass();
  };
  for (const c of cases) {
    // the first draw requests the pipeline, which compiles in the background (draws skipped meanwhile)
    draw(c.viewport);
    await graphics.resources.pending.whenReady();
    draw(c.viewport);
    // rows top to bottom (rendered content is stored flipped on WebGPU, so not `readTextureData`)
    const array = await graphics.readTextureAsync(target, 0, 0, size, size);
    const red = ([x, y]: readonly [number, number]) => array[(y * size + x) * 4];
    if (red(c.drawn) !== 255) errors.push(`${c.name}: pixel ${c.drawn} not drawn`);
    if (red(c.empty) !== 0) errors.push(`${c.name}: pixel ${c.empty} drawn (viewport squashed into the target)`);
  }
  copy.dispose();
  framebuffer.destroy();
  target.destroy();
  white.destroy();
  return errors;
}

/** OBJ triangles of the render objects of `repr` */
async function exportedTriangles(graphics: GraphicsContext, repr: Representation.Any) {
  const exporter = new ObjExporter('test', Box3D.create(Vec3.create(-1, -1, -1), Vec3.create(1, 1, 1)));
  for (const ro of repr.renderObjects) await exporter.add(ro, graphics, RuntimeContext.Synchronous);
  const { obj } = await exporter.getData();
  return obj.split('\n').filter((line) => line.startsWith('f ')).length;
}

/**
 * A GPU isosurface updated to another iso value (into the capacity of the first extraction, drawn indirectly on WebGPU,
 * with a degenerate tail) exported against the CPU one: triangle counts within 2%
 */
async function checkExportGpuSurface(graphics: GraphicsContext) {
  const volume = blobVolume(32);
  const iso = (v: number) => Volume.IsoValue.absolute(v);
  const gpu = await volumeRepr(withoutStrictSupport(graphics), IsosurfaceRepresentationProvider, volume, {
    isoValue: iso(0.5),
    tryUseGpu: true,
  });
  await gpu.createOrUpdate({ isoValue: iso(0.6) }, volume).run();
  const cpu = await volumeRepr(graphics, IsosurfaceRepresentationProvider, volume, {
    isoValue: iso(0.6),
    tryUseGpu: false,
  });
  const types = gpu.renderObjects.map((ro) => ro.type);
  const [gpuTriangles, cpuTriangles] = [await exportedTriangles(graphics, gpu), await exportedTriangles(graphics, cpu)];
  gpu.destroy();
  cpu.destroy();
  if (!types.includes('texture-mesh')) return [`the GPU surface is '${types[0]}', not a texture mesh`];
  if (gpuTriangles === 0) return ['no triangles exported from the GPU surface'];
  if (Math.abs(gpuTriangles - cpuTriangles) > 0.02 * cpuTriangles) {
    return [`${gpuTriangles} triangles exported from the GPU surface, ${cpuTriangles} from the CPU one`];
  }
  return [];
}

/**
 * A GPU isosurface colored by the volume value at its vertices (vertex granularity): the vertices are read back for the
 * position iterator, so the vertex colors vary (without positions, every vertex gets the default color)
 */
async function checkVertexColorsGpuSurface(graphics: GraphicsContext) {
  const repr = await volumeRepr(
    withoutStrictSupport(graphics),
    IsosurfaceRepresentationProvider,
    blobVolume(32),
    { isoValue: Volume.IsoValue.absolute(0.5), tryUseGpu: true },
    { name: 'volume-value', params: { domain: { name: 'custom', params: [0.45, 0.55] } } },
  );
  const ro = repr.renderObjects[0];
  const errors: string[] = [];
  if (ro?.type !== 'texture-mesh') errors.push(`the GPU surface is '${ro?.type}', not a texture mesh`);
  else {
    const values = ro.values as unknown as {
      tColor: { ref: { value: { array: Uint8Array; count: number } } };
      uVertexCount: { ref: { value: number } };
    };
    // rgba per vertex (the data array may be longer)
    const stride = 4;
    const { array } = values.tColor.ref.value;
    const colors = new Set<number>();
    for (let i = 0; i < values.uVertexCount.ref.value * stride; i += stride) {
      colors.add((array[i] << 16) | (array[i + 1] << 8) | array[i + 2]);
    }
    if (colors.size < 2) errors.push('all vertices have the same color (no vertex positions)');
  }
  repr.destroy();
  return errors;
}

/**
 * The Hi-Z pass over a depth target with a near bottom-left quadrant (depth 0.8), far elsewhere (depth 0; depth is
 * reversed), seen by an orthographic camera: a sphere behind the quadrant is occluded, ones above or right of it are not
 */
async function checkHiZOcclusion(graphics: GraphicsContext) {
  if (!graphics.createTextureReader().asynchronous) return []; // Hi-Z needs asynchronous reads (WebGL2, WebGPU)
  const errors: string[] = [];
  const { width, height } = graphics.getDrawingBufferSize();

  const depth = graphics.resources.texture('image-float32', 'r', 'float', 'nearest');
  depth.define(width, height);
  const framebuffer = graphics.resources.framebuffer();
  depth.attachFramebuffer(framebuffer, 'color0');
  graphics.beginRenderPass(framebuffer, { label: 'hi-z depth' });
  graphics.setViewport(0, 0, width, height);
  // clears cover the scissor, which earlier work may have left set (e.g. the GPU Gaussian density)
  graphics.disableScissor();
  graphics.clear({ color: [0, 0, 0, 0] });
  // the bottom-left quadrant (rectangles have a top-left origin)
  graphics.setScissor(0, height / 2, width / 2, height / 2);
  graphics.clear({ color: [0.8, 0.8, 0.8, 0.8] });
  graphics.disableScissor();

  const hiZ = new HiZPass(graphics, { depthTextureOpaque: depth } as DrawPass, undefined, {
    enabled: true,
    minLevel: 1,
    maxFrameLag: 30,
  });
  hiZ.setViewport(0, 0, width, height);
  // view space z = 10 * depth - 11: -3 in the quadrant, -11 elsewhere
  const near = 1;
  const far = 11;
  const camera = {
    near,
    far,
    view: Mat4.identity(),
    projection: Mat4.ortho(Mat4(), -100, 100, 100, -100, near, far),
  } as unknown as Camera;
  const behind = Sphere3D.create(Vec3.create(-50, -50, -6), 1);
  const above = Sphere3D.create(Vec3.create(-50, 50, -6), 1);
  const right = Sphere3D.create(Vec3.create(50, -50, -6), 1);
  // renders once no work is pending (the pass skips frames with pipelines compiling in the background, e.g. of other
  // cases) until a read resolved
  for (let i = 0; i < 100 && !hiZ.isOccluded(behind); ++i) {
    await graphics.resources.pending.whenReady();
    hiZ.render(camera);
    await new Promise((resolve) => setTimeout(resolve, 20));
    hiZ.tick();
  }
  if (!hiZ.isOccluded(behind)) errors.push('a sphere behind the near quadrant is not occluded');
  if (hiZ.isOccluded(above)) errors.push('a sphere above the near quadrant is occluded');
  if (hiZ.isOccluded(right)) errors.push('a sphere right of the near quadrant is occluded');

  hiZ.dispose();
  framebuffer.destroy();
  depth.destroy();
  return errors;
}

/**
 * A WebGPU context on a canvas without size (as before the page is laid out): frames go to a placeholder until the canvas
 * has a size (Firefox rejects configuring a zero-sized canvas), then to the canvas. Not covered: resizing a configured
 * canvas to zero, which Firefox answers by reconfiguring it itself (a validation error the context cannot prevent);
 * Canvas3D keeps the drawing buffer at least 1×1
 */
async function checkZeroSizeCanvas(graphics: GraphicsContext) {
  if (!isWebGPUContext(graphics)) return [];
  const errors: string[] = [];
  const canvas = document.createElement('canvas');
  canvas.width = 0;
  canvas.height = 0;
  // Chrome accepts configuring a zero-sized canvas, Firefox does not: record the sizes `configure` sees
  const configure = GPUCanvasContext.prototype.configure;
  const configured: string[] = [];
  GPUCanvasContext.prototype.configure = function (this: GPUCanvasContext, config: GPUCanvasConfiguration) {
    if (this.canvas === canvas) configured.push(`${canvas.width}×${canvas.height}`);
    return configure.call(this, config);
  };
  let ctx: ReturnType<typeof createWebGPUContext>;
  try {
    ctx = createWebGPUContext(graphics.device, { width: 0, height: 0, canvas });
  } finally {
    GPUCanvasContext.prototype.configure = configure;
  }
  const clear = (color: [number, number, number, number]) => {
    ctx.beginRenderPass('drawing-buffer');
    ctx.clear({ color });
    ctx.flush();
  };
  const formats = new Set<GPUTextureFormat>();
  GPUCanvasContext.prototype.configure = function (this: GPUCanvasContext, config: GPUCanvasConfiguration) {
    if (this.canvas === canvas) {
      configured.push(`${canvas.width}×${canvas.height}`);
      formats.add(config.format);
    }
    return configure.call(this, config);
  };
  try {
    clear([0, 1, 0, 1]);
    canvas.width = 16;
    canvas.height = 16;
    clear([1, 0, 0, 1]);
  } finally {
    GPUCanvasContext.prototype.configure = configure;
  }
  if (configured.some((size) => size.startsWith('0×') || size.endsWith('×0'))) {
    errors.push(`configured at the sizes ${configured.join(', ')}`);
  }
  // the browser's preferred format (no copy per frame), read back as rgba
  const preferred = navigator.gpu.getPreferredCanvasFormat();
  if ([...formats].some((f) => f !== preferred))
    errors.push(`configured as ${[...formats].join(', ')}, not ${preferred}`);
  // read in the task that rendered the frame
  const image = await ctx.readDrawingBuffer();
  if (image.width !== 16 || image.height !== 16) errors.push(`drawing buffer of ${image.width}×${image.height}`);
  const [r, g, b, a] = image.data;
  if (r !== 255 || g !== 0 || b !== 0 || a !== 255) errors.push(`cleared to ${[r, g, b, a].join(',')}, not red`);
  errors.push(...ctx.takeErrors());
  ctx.destroy();
  return errors;
}

/**
 * Devices are requested only from adapters with the required features (`AutoWebGPUBackend` requires
 * `float32-blendable`, which e.g. Firefox lacks), else none, so that a list of backends falls back to WebGL
 */
async function checkRequiredFeatures(graphics: GraphicsContext) {
  if (!isWebGPUContext(graphics)) return [];
  const errors: string[] = [];
  const missing = await requestWebGPUDevice('default', navigator.gpu, 'core', 'none', [
    'not-a-feature' as GPUFeatureName,
  ]);
  if (missing) {
    errors.push('a device for an adapter without a required feature');
    missing.destroy();
  }
  const adapter = await navigator.gpu.requestAdapter();
  if (adapter?.features.has('float32-blendable')) {
    const device = await requestWebGPUDevice('default', navigator.gpu, 'core', 'none', ['float32-blendable']);
    if (!device) errors.push('no device although the adapter has float32-blendable');
    else if (!device.features.has('float32-blendable')) errors.push('the required feature was not requested');
    device?.destroy();
  }
  return errors;
}

export const GpuPathCases: VisualTestCase[] = [
  {
    id: 'gpu-paths/representations',
    description:
      'Representations with GPU paths take them (texture meshes, the GPU density) with support as reported without strict mode.',
    requires: ['floatTextureLinear', 'floatRenderTarget', 'multipleRenderTargets', 'halfFloatRenderTarget'],
    scene: async () => [],
    check: checkGpuPaths,
  },
  {
    id: 'gpu-paths/gaussian-surface-units',
    description: 'The GPU Gaussian surface of PDB entry 7qpd per unit against the CPU one (vertex counts within 2%).',
    // without flat varyings the CPU surface splits triangles of mixed groups (`Mesh.uniformTriangleGroup`)
    requires: ['gaussianSurfaceGpu', 'halfFloatRenderTarget', 'flatInterpolation'],
    scene: async () => [],
    check: checkGaussianSurfaceUnits,
  },
  {
    id: 'gpu-paths/alpha-orbitals',
    description: 'Alpha orbitals and density with grid3d on the GPU against the CPU collocation (angular momenta 0–4).',
    scene: async () => [],
    check: checkAlphaOrbitals,
  },
  {
    id: 'gpu-paths/viewport-cropping',
    description: 'A viewport partly outside the target is cropped, not squashed into the target.',
    scene: async () => [],
    check: checkViewportCropping,
  },
  {
    id: 'gpu-paths/read-texture-data',
    description: 'Texture readback in data order on every backend: rgba and r bytes, rgba floats where renderable.',
    scene: async () => [],
    check: checkReadTextureData,
  },
  {
    id: 'gpu-paths/load-rows',
    description:
      'Rows of an image changed after its upload written in place (a font atlas), the rest kept; another size whole.',
    scene: async () => [],
    check: checkLoadRows,
  },
  {
    id: 'gpu-paths/export-gpu-surface',
    description:
      'A GPU isosurface, updated into its capacity, exported (OBJ) against the CPU one: triangles within 2%.',
    requires: ['floatRenderTarget', 'multipleRenderTargets'],
    scene: async () => [],
    check: checkExportGpuSurface,
  },
  {
    id: 'gpu-paths/vertex-colors-gpu-surface',
    description: 'A GPU isosurface colored by the volume value at its vertices (vertex positions read back).',
    requires: ['floatRenderTarget', 'multipleRenderTargets'],
    scene: async () => [],
    check: checkVertexColorsGpuSurface,
  },
  {
    id: 'gpu-paths/hi-z-occlusion',
    description:
      'Hi-Z levels of a depth target read back asynchronously: a sphere behind its near bottom-left quadrant is occluded, ones above or right of it are not.',
    requires: ['floatRenderTarget'],
    scene: async () => [],
    check: checkHiZOcclusion,
  },
  {
    id: 'gpu-paths/zero-size-canvas',
    description: 'A WebGPU context created on a canvas without size renders once the canvas is sized (WebGPU only).',
    scene: async () => [],
    check: checkZeroSizeCanvas,
  },
  {
    id: 'gpu-paths/required-features',
    description: 'WebGPU devices only from adapters with the required features, else none (WebGPU only).',
    scene: async () => [],
    check: checkRequiredFeatures,
  },
];
