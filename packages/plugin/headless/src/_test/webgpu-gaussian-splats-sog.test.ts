/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { gaussianSplatsFromSog, getSogShLayout } from '@molstar/model/formats/shape/gaussian-splat/sog';
import { GaussianSplats } from '@molstar/graphics/geo/geometry/gaussian-splats/gaussian-splats';
import { int32ToFloatBits } from '@molstar/graphics/geo/geometry/center-group';
import { buildGaussianSplats, createSplatShape } from '@molstar/graphics/formats/shape/gaussian-splat';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { ShapeRepresentation } from '@molstar/graphics/repr/shape/representation';
import type { PostprocessingProps } from '@molstar/graphics/canvas3d/passes/postprocessing';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { createSplatSorter, SplatSortKernels } from '@molstar/graphics/gl/compute/splat-sort';
import { WebGPU } from './webgpu-module.js';

type SogFile = Parameters<typeof gaussianSplatsFromSog>[0];
type SogImage = SogFile['sh0'];

const Size = 96;

const NoPostprocessing: Partial<PostprocessingProps> = {
  occlusion: { name: 'off', params: {} },
  outline: { name: 'off', params: {} },
  antialiasing: { name: 'off', params: {} },
  shadow: { name: 'off', params: {} },
};

/** the codebook of the values used (log scales, base color, view-dependent color): a byte each */
const Codebook = [-1.5, -1, -0.5, 0, 0.25, 0.5, 0.75, 1, 1.5];
const code = (v: number) => Codebook.indexOf(v);

interface TestSplat {
  /** the center in the log domain (`sign(x) ln(|x| + 1)`), at most 2 an axis */
  logPosition: number[];
  logScale: number[];
  rotation: Quat;
  opacity: number;
  dc: number[];
}

/** overlapping splats of three colors, turned, of several opacities, a palette label each */
const Splats: TestSplat[] = [
  {
    logPosition: [0, 0, 0],
    logScale: [1, 0.5, 0],
    rotation: Quat.create(0.2, -0.5, 0.1, 0.8),
    opacity: 0.9,
    dc: [1, -1, -1],
  },
  { logPosition: [1, 0.5, 0.5], logScale: [0.5, 0.5, 0.5], rotation: Quat.identity(), opacity: 0.6, dc: [-1, 1, -0.5] },
  {
    logPosition: [-1, -0.5, -1],
    logScale: [1, 0, 0.5],
    rotation: Quat.create(-0.6, 0.3, 0.2, 0.7),
    opacity: 0.75,
    dc: [-0.5, 0, 1.5],
  },
];
/** two palette entries of degree 1 (3 coefficients, rgb codebook values each), the splats' labels */
const Palette = [
  [
    [0.5, -0.5, 0],
    [0, 0.25, -0.5],
    [-0.5, 0.75, 0.25],
  ],
  [
    [-0.5, 0, 0.5],
    [0.25, -0.5, 0],
    [0, 0.5, -1],
  ],
];
const Labels = [1, 0, 1];
const MeansRange = 2;

function image(width: number, height: number, fill: (i: number, out: Uint8Array) => void): SogImage {
  const data = new Uint8Array(width * height * 4);
  const px = new Uint8Array(4);
  for (let i = 0; i < width * height; ++i) {
    px.fill(0);
    fill(i, px);
    data.set(px, i * 4);
  }
  return { width, height, data };
}

/** version 3's exponential map of `q` (`w, x, y, z`, w made positive), as splat-transform writes it */
function expMapBytes(q: Quat) {
  let [x, y, z, w] = Quat.normalize(Quat(), q);
  if (w < 0) [x, y, z, w] = [-x, -y, -z, -w];
  const sv = Math.sqrt(x * x + y * y + z * z);
  const k = sv > 1e-12 ? Math.atan2(sv, w) / sv / (Math.PI / 2) : 1 / (Math.PI / 2);
  return [x, y, z].map((v) => Math.min(255, Math.max(0, Math.round(v * k * 127 + 128))));
}

/** a version 3 SOG of `Splats` in images of 2x2 (one pixel unused) with the palette */
function sogFile(): SogFile {
  const count = Splats.length;
  const q16 = (i: number, c: number) =>
    Math.round(((Splats[i].logPosition[c] + MeansRange) / (2 * MeansRange)) * 65535);
  const per = (fill: (s: TestSplat, i: number, out: Uint8Array) => void) =>
    image(2, 2, (i, out) => {
      if (i < count) fill(Splats[i], i, out);
    });
  return {
    meta: {
      version: 3,
      count,
      means: { mins: [-MeansRange, -MeansRange, -MeansRange], maxs: [MeansRange, MeansRange, MeansRange], files: [] },
      scales: { codebook: Codebook, files: [] },
      quats: { files: [] },
      sh0: { codebook: Codebook, files: [] },
      shN: { codebook: Codebook, files: [], count: Palette.length, bands: 1 },
    },
    version: 3,
    count,
    meansLow: per((_, i, out) => out.set([0, 1, 2].map((c) => q16(i, c) & 0xff))),
    meansHigh: per((_, i, out) => out.set([0, 1, 2].map((c) => q16(i, c) >> 8))),
    scales: per((s, _, out) => out.set(s.logScale.map(code))),
    quats: per((s, _, out) => out.set(expMapBytes(s.rotation))),
    sh0: per((s, _, out) => {
      out.set(s.dc.map(code));
      out[3] = Math.round(Math.sqrt(s.opacity) * 204);
    }),
    shN: {
      // 64 entries a row, 3 coefficients (pixels) each
      centroids: image(64 * 3, 1, (p, out) => {
        const entry = Palette[Math.floor(p / 3)];
        if (entry) out.set(entry[p % 3].map(code));
      }),
      labels: per((_, i, out) => out.set([Labels[i] & 0xff, Labels[i] >> 8, 0])),
    },
  };
}

/** `file` read in place on `graphics` (the `sog` encoding), its images uploaded as they are */
function sogSplats(graphics: GraphicsContext, file: SogFile, shDegree: number) {
  const upload = (img: SogImage) => {
    const texture = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    texture.load({ array: img.data, width: img.width, height: img.height });
    return texture;
  };
  // the centers as the CPU decodes them (`gaussianSplatsFromSog`)
  const data = gaussianSplatsFromSog(file, { shDegree: 0 });
  const positionGroup = new Float32Array(file.count * 4);
  const s = { x: 0, y: 0, z: 0 } as Parameters<typeof data.getSplat>[1];
  for (let i = 0; i < file.count; ++i) {
    data.getSplat(i, s);
    positionGroup.set([s.x, s.y, s.z, int32ToFloatBits(i)], i * 4);
  }
  const codebook = new Float32Array(256 * 4);
  Codebook.forEach((v, i) => codebook.set([v, v, v, 0], i * 4));
  const { centroids, labels } = file.shN!;
  const layout = getSogShLayout(file.meta, centroids.width, centroids.height);
  return GaussianSplats.createSog({
    graphics,
    count: file.count,
    positionGroup,
    // the means' ranges out of the log domain
    centersBounds: Sphere3D.create(Vec3(), (Math.exp(MeansRange) - 1) * Math.sqrt(3)),
    textures: {
      meansLow: upload(file.meansLow),
      meansHigh: upload(file.meansHigh),
      scales: upload(file.scales),
      quats: upload(file.quats),
      sh0: upload(file.sh0),
      shLabels: upload(labels),
      shCentroids: upload(centroids),
    },
    meansMin: Vec3.create(-MeansRange, -MeansRange, -MeansRange),
    meansMax: Vec3.create(MeansRange, MeansRange, MeansRange),
    codebook,
    version: 3,
    sh: { degree: Math.min(shDegree, layout.bands), ...layout },
    maxScale: Math.exp(Math.max(...Codebook)),
  });
}

describe.skipIf(!WebGPU.module)('SOG read in place on WebGPU', () => {
  /** the splats of `make` (with `colors` by group: the packed path's base colors, else the shader's) seen from `eye` */
  async function render(make: (graphics: GraphicsContext) => GaussianSplats, eye: Vec3, colors?: Uint8Array) {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const canvas = plugin.canvas3d!;
      canvas.setProps({ renderer: { backgroundColor: Color(0x000000) }, cameraFog: { name: 'off', params: {} } });
      const shape = createSplatShape('splats', {}, undefined, {
        splats: make(canvas.graphics),
        colors,
        upright: false,
      });
      const repr = ShapeRepresentation(() => shape, GaussianSplats.Utils, { graphics: canvas.graphics });
      await repr.createOrUpdate(PD.getDefaultValues(GaussianSplats.Params), {}).run();
      canvas.add(repr);
      canvas.commit(true);
      canvas.camera.setState({ position: eye, target: Vec3(), up: Vec3.unitY }, 0);
      const { pending } = canvas.graphics.resources;
      let image = await plugin.getImageRaw({ width: Size, height: Size }, NoPostprocessing);
      for (let i = 0; i < 10; i++) {
        await pending.whenReady();
        image = await plugin.getImageRaw({ width: Size, height: Size }, NoPostprocessing);
      }
      expect(canvas.graphics.takeErrors?.() ?? []).toEqual([]);
      return image.data;
    } finally {
      plugin.dispose();
    }
  }

  it('draws as the CPU decode in the packed layout: positions, scales, rotations, colors, palette', async () => {
    const file = sogFile();
    // a runtime that never yields
    const runtime = { shouldUpdate: false, update: async () => {} } as never;
    const built = await buildGaussianSplats(runtime, gaussianSplatsFromSog(file, { shDegree: 1 }));
    for (const eye of [Vec3.create(0, 0, 12), Vec3.create(9, 4, -6)]) {
      const packed = await render(() => built.splats, eye, built.colors);
      const sog = await render((graphics) => sogSplats(graphics, file, 1), eye);
      let max = 0;
      let lit = 0;
      for (let i = 0; i < packed.length; i += 4) {
        for (let c = 0; c < 3; ++c) max = Math.max(max, Math.abs(packed[i + c] - sog[i + c]));
        if (packed[i] + packed[i + 1] + packed[i + 2] > 30) lit++;
      }
      // the splats in view; the packed layout quantizes once more (rotations, scales, colors as bytes)
      expect(lit).toBeGreaterThan(500);
      expect(max).toBeLessThanOrEqual(8);
    }
  }, 120_000);

  // the encoding's programs (its defines, the direct colors) prepared from a list recorded where it was drawn
  it('a precompile list of SOG read in place prepares what its frames draw it with, before it exists', async () => {
    const create = async () => {
      const plugin = await HeadlessPluginContext.create(
        { webgpu: WebGPU.module! },
        DefaultPluginSpec(),
        { width: 32, height: 32 },
        { backend: 'webgpu' },
      );
      await plugin.init();
      return plugin;
    };
    const add = async (plugin: HeadlessPluginContext) => {
      const canvas = plugin.canvas3d!;
      const shape = createSplatShape('splats', {}, undefined, {
        splats: sogSplats(canvas.graphics, sogFile(), 1),
        upright: false,
      });
      const repr = ShapeRepresentation(() => shape, GaussianSplats.Utils, { graphics: canvas.graphics });
      await repr.createOrUpdate(PD.getDefaultValues(GaussianSplats.Params), {}).run();
      canvas.add(repr);
      canvas.commit(true);
    };
    const a = await create();
    const b = await create();
    try {
      a.canvas3d!.resetPrecompileList();
      await add(a);
      const recorded = JSON.parse(JSON.stringify(a.canvas3d!.getPrecompileList()));
      const splats = recorded.objects.find((o: { type: string }) => o.type === 'gaussian-splats');
      expect(splats?.defines.dSplatEncoding).toBe('sog');

      const canvas = b.canvas3d!;
      const { pending } = canvas.graphics.resources;
      await canvas.prepare({ interaction: true });
      await pending.whenReady();
      await canvas.precompile(recorded);
      // added: their frames and the first interaction compile nothing more
      const before = pending.trackedCount;
      await add(b);
      await canvas.prepare({ interaction: true });
      await pending.whenReady();
      expect(pending.trackedSince(before)).toEqual([]);
      expect(canvas.graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      a.dispose();
      b.dispose();
    }
  }, 60_000);

  // the GPU sorts SOG by the centers it decodes itself (the CPU's read later): the same as the CPU's
  it("decodes the centers on the GPU as the CPU does, for the GPU's sort", async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const compute = graphics.compute!;
      await compute.prepareKernels(SplatSortKernels);
      const splats = sogSplats(graphics, sogFile(), 0);
      const sorter = createSplatSorter(compute);
      const positions = compute.createBuffer(splats.splatCount * 16, 'centers');
      const { sog } = splats;
      sorter.decodeSogCenters(
        sog.meansLow.ref.value,
        sog.meansHigh.ref.value,
        sog.meansMin.ref.value,
        sog.meansMax.ref.value,
        splats.splatCount,
        positions,
      );
      const gpu = new Float32Array(await positions.readAsync());
      const cpu = splats.positionGroup.ref.value.array;
      for (let i = 0; i < splats.splatCount; ++i) {
        for (let c = 0; c < 3; ++c) expect(gpu[i * 4 + c]).toBeCloseTo(cpu[i * 4 + c], 4);
      }
      sorter.destroy();
      positions.destroy();
      GaussianSplats.dispose(splats);
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
