/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { writeParticlePack } from '@molstar/model/formats/particles/pack';
import type { ExternalModules } from '@molstar/plugin-headless/screenshot';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { ParticlePackProvider } from '@molstar/plugin/state/formats/particles/pack';
import { MmcifParticlesProvider } from '@molstar/plugin/state/formats/particles/mmcif-assembly';
import type { HeadlessPluginContext as Plugin } from '@molstar/plugin-headless/context';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { Asset } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Particle } from '@molstar/model/model/particles/particle-list';
import { RelionStarParticlesProvider } from '@molstar/plugin/state/formats/particles/star';
import { ParticleListWithTargets } from '@molstar/plugin/state/transforms/particles/ops';
import { WebGPU } from './webgpu-module.js';

const readExample = (name: string) => readFileSync(new URL(`../../../../../data/examples/${name}`, import.meta.url));

const require = createRequire(import.meta.url);
let glModule: typeof import('gl') | undefined;
try {
  glModule = require('gl') as typeof import('gl');
  if (!glModule(1, 1)) glModule = undefined;
} catch {
  glModule = undefined;
}

/** a target of three oxygens about the origin, a particle at the identity */
const Cif = `data_site_colors
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 O O HOH A 1 1 0 0 0 1
ATOM 2 O O HOH A 1 1 0.4 0 0 1
ATOM 3 O O HOH A 1 1 -0.4 0 0 1

loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 (1) A

loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 'identity operation' 1 0 0 0 0 1 0 0 0 0 1 0
`;

/** The assembly's list as read by the plugin (`MmcifParticlesProvider`), written as a pack: its sites quantized */
async function createPack(plugin: Plugin) {
  const cif = await plugin.builders.data.rawData({ data: Cif, label: 'cif' });
  const { list } = await MmcifParticlesProvider.parse(plugin, cif, { assemblyId: '1', variant: 'standard' });
  return writeParticlePack(list.cell!.obj!.data);
}

const Size = 64;

/**
 * The center pixel of the target drawn green (its particles' color) from each of `distances` Å (a field of view of
 * 45°): an oxygen's radius on screen about 20 px from 6 Å, 3 px from 40 Å (the second color from 4 px, whole from 12 px)
 */
async function renderCenter(
  modules: ExternalModules,
  backend: 'webgl' | 'webgpu',
  distances: number[],
  raster?: 'off' | 'only',
) {
  const spec = DefaultPluginSpec();
  const plugin = await HeadlessPluginContext.create(
    modules,
    { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
    { width: Size, height: Size },
    { backend },
  );
  try {
    await plugin.init();
    if (raster === 'only' && !plugin.canvas3d!.graphics.capabilities.atomicMinMax64) return undefined;
    // the second color is the sphere impostors' (a sphere mesh has none): none without frag depth (headless-gl)
    if (!plugin.canvas3d!.graphics.capabilities.fragDepth) return undefined;
    plugin.canvas3d!.setProps({
      renderer: raster ? { sphereRaster: raster, sphereRasterRadius: 32 } : {},
      camera: { mode: 'perspective' } as never,
      cameraFog: { name: 'off', params: {} },
      multiSample: { mode: 'off' } as never,
      postprocessing: {
        occlusion: { name: 'off', params: {} },
        outline: { name: 'off', params: {} },
      },
    });
    const data = await plugin.builders.data.rawData({ data: await createPack(plugin), label: 'pack' });
    const { list } = await ParticlePackProvider.parse(plugin, data);
    await plugin
      .build()
      .to(list)
      .apply(ParticlesRepresentation3D, {
        type: { name: 'target', params: {} },
        colorTheme: { name: 'uniform', params: { value: 0x00ff00 } },
      } as never)
      .commit();
    const o = ((Size / 2) * Size + Size / 2) * 4;
    return await renderAt(plugin, distances, Vec3(), (pixels) => [pixels[o], pixels[o + 1], pixels[o + 2]]);
  } finally {
    plugin.dispose();
  }
}

/**
 * The pixels redder than green of a structure target (1crn at the particles of a RELION list) drawn green from
 * each of `distances` Å of the first particle (where the structure's centroid is placed; a field of view of 45°): its atoms by
 * element near the camera (oxygens red), by the particles' color far from it
 */
async function renderStructureTarget(distances: number[], raster: 'off' | 'only') {
  const spec = DefaultPluginSpec();
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
    { width: Size, height: Size },
    { backend: 'webgpu' },
  );
  try {
    await plugin.init();
    if (raster === 'only' && !plugin.canvas3d!.graphics.capabilities.atomicMinMax64) return undefined;
    plugin.canvas3d!.setProps({
      renderer: { sphereRaster: raster, sphereRasterRadius: 32 },
      camera: { mode: 'perspective' } as never,
      cameraFog: { name: 'off', params: {} },
      multiSample: { mode: 'off' } as never,
      postprocessing: { occlusion: { name: 'off', params: {} }, outline: { name: 'off', params: {} } },
    });
    const star = await plugin.builders.data.rawData({
      data: readExample('particles/relion-subtomo.star').toString(),
      label: 'star',
    });
    const { list } = await RelionStarParticlesProvider.parse(plugin, star);
    const particles = list.cell!.obj!.data;
    const targetId = particles.targetInfo.keys().next().value!;
    const file = Asset.File(new File([readExample('1crn.cif')], '1crn.cif'));
    await plugin
      .build()
      .to(list)
      .apply(ParticleListWithTargets, {
        ...PD.getDefaultValues(ParticleListWithTargets.definition.params!(undefined, plugin)),
        targets: [{ targetId, source: { name: 'file', params: { file, format: 'auto' } } }],
      } as never)
      .apply(ParticlesRepresentation3D, {
        type: { name: 'target', params: {} },
        colorTheme: { name: 'uniform', params: { value: 0x00ff00 } },
      } as never)
      .commit();
    const center = Particle.position(Vec3(), Particle.Location(particles, 0));
    return await renderAt(plugin, distances, center, (pixels) => {
      let red = 0;
      for (let i = 0; i < pixels.length; i += 4) if (pixels[i] > pixels[i + 1] + 40) red += 1;
      return red;
    });
  } finally {
    plugin.dispose();
  }
}

/**
 * `read` of the images looking at `target` from each of `distances` Å along z: pipelines and kernels compile in the
 * background (frames drawn before they are ready may draw nothing); the camera set after the first frame (which
 * focuses the scene)
 */
async function renderAt<T>(plugin: Plugin, distances: number[], target: Vec3, read: (pixels: Uint8ClampedArray) => T) {
  const image = () => plugin.getImageRaw({ width: Size, height: Size });
  await image();
  await plugin.canvas3d!.graphics.resources.pending.whenReady();
  const results: T[] = [];
  for (const distance of distances) {
    const position = Vec3.add(Vec3(), target, Vec3.create(0, 0, distance));
    plugin.canvas3d!.camera.setState({ target, position, up: Vec3.create(0, 1, 0), fov: Math.PI / 4 }, 0);
    for (let i = 0; i < 2; ++i) {
      await image();
      await plugin.canvas3d!.graphics.resources.pending.whenReady();
    }
    results.push(read((await image()).data));
  }
  expect(plugin.canvas3d!.graphics.takeErrors?.() ?? []).toEqual([]);
  return results;
}

/** near: the oxygens' red; far: the particles' green */
function expectNearAndFar(near: number[], far: number[]) {
  expect(near[0]).toBeGreaterThan(near[1] + 40);
  expect(far[1]).toBeGreaterThan(far[0] + 40);
}

// the second color of a target's sites: by element near the camera, by the particle's theme
// far from it
describe.skipIf(!WebGPU.module)('particle site colors on WebGPU', () => {
  for (const raster of ['off', 'only'] as const) {
    // the sphere raster's case medium (a second each)
    const lane = raster === 'only' ? { tags: ['medium'] } : {};
    it(`blend in by size on screen (sphereRaster ${raster})`, lane, async () => {
      const modules = { webgpu: WebGPU.module! };
      const colors = await renderCenter(modules, 'webgpu', [6, 40], raster);
      if (!colors) return;
      expectNearAndFar(colors[0], colors[1]);
    });
  }
});

describe.skipIf(!WebGPU.module)('structure target colors on WebGPU', () => {
  it.each(['off', 'only'] as const)(
    'blend in by size on screen, a texel an element (sphereRaster %s)',
    { tags: ['medium'] },
    async (raster) => {
      const reds = await renderStructureTarget([25, 400], raster);
      if (reds === undefined) return;
      expect(reds[0]).toBeGreaterThan(20);
      expect(reds[1]).toBe(0);
    },
  );
});

describe('particle site colors on WebGL', () => {
  (glModule ? it : it.skip)('blend in by size on screen, a site each', { tags: ['medium'] }, async (ctx) => {
    const modules = { gl: glModule! };
    const colors = await renderCenter(modules, 'webgl', [6, 40]);
    // headless-gl has no EXT_frag_depth: the spheres are meshes, which have no second color
    if (!colors) return ctx.skip();
    expectNearAndFar(colors[0], colors[1]);
  });
});
