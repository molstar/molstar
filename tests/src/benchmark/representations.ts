/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ValueCell } from '@molstar/core/util';
import type { TextureMesh } from '@molstar/graphics/geo/geometry/texture-mesh/texture-mesh';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Representation, RepresentationProvider } from '@molstar/graphics/repr/representation';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { IsosurfaceRepresentationProvider } from '@molstar/graphics/repr/volume/isosurface';
import { Theme } from '@molstar/graphics/theme/theme';
import type { Structure } from '@molstar/model/model/structure';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { Volume } from '@molstar/model/model/volume';
import { blobVolume } from '../visual/scenes/volume.js';
import { type BenchmarkCase, type PathRun, StageTimer, withGpuPasses } from './framework.js';
import { exampleStructure, getThemeRegistries, globuleStructure } from './structures.js';

type ReprInput = {
  provider: RepresentationProvider<any, any, any>;
  data: Structure | Volume;
  props: { [k: string]: unknown };
  sizeTheme: string;
};

/**
 * Triangles: indices of meshes, vertices of texture meshes, three each; of indirectly drawn ones the count read back
 * (after pending work settled), not the capacity drawn
 */
function drawnTriangles(repr: Representation.Any) {
  let triangles = 0;
  for (const { values } of repr.renderObjects) {
    const meta = (values as { meta?: ValueCell<TextureMesh['meta']> }).meta?.ref.value;
    const indirect = 'drawIndirect' in values && !!values.drawIndirect.ref.value;
    const read: number | undefined = indirect ? (meta?.readVertexCount as number | undefined) : undefined;
    const count: number = read ?? (values.drawCount.ref.value as number);
    triangles += count / 3;
  }
  return triangles;
}

/**
 * Creates the representation from scratch, until its geometry is created and the GPU work completed, then until
 * pending work settled (the vertex count of indirectly drawn surfaces read back, and extracted again on overflow)
 */
async function createRepresentation(graphics: GraphicsContext, input: ReprInput, tryUseGpu: boolean): Promise<PathRun> {
  const { provider, props, sizeTheme } = input;
  // a volume caches its GPU texture: a fresh volume (sharing the grid) per run, so each run uploads as a new one would
  const data =
    'grid' in input.data ? { ...input.data, _propertyData: {}, customProperties: new CustomProperties() } : input.data;
  const registries = getThemeRegistries();
  const timer = new StageTimer(graphics);
  const repr: Representation.Any = provider.factory({ graphics, ...registries }, provider.getParams);
  const themeData = 'grid' in data ? { volume: data } : { structure: data };
  repr.setTheme(
    Theme.create(registries, themeData, { colorTheme: { name: 'uniform' }, sizeTheme: { name: sizeTheme } }),
  );
  const params = PD.getDefaultValues(provider.getParams(registries as never, data as never));
  const { gpuPasses } = await withGpuPasses(graphics, async () => {
    await repr.createOrUpdate({ ...params, ...props, tryUseGpu }, data).run();
    await timer.stage('representation');
    await graphics.resources.pending.whenReady();
    await timer.stage('settled');
  });
  const triangles = drawnTriangles(repr);
  repr.destroy();
  if ('grid' in data) data.customProperties.dispose();
  return { stages: timer.stages, triangles, gpuPasses };
}

type UpdateInput = ReprInput & {
  /** two values of a prop that recreates the geometry, alternated between runs */
  update: { name: string; values: [unknown, unknown] };
};

/** The representation updated by runs of a path, created on its first run (not timed) */
const updatedReprs = new WeakMap<UpdateInput, Map<string, { repr: Representation.Any; next: number }>>();

/**
 * Updates the representation (a prop that recreates the geometry, alternating between two values) until its geometry is
 * updated and the GPU work completed, then until pending work settled
 */
async function updateRepresentation(
  graphics: GraphicsContext,
  input: UpdateInput,
  path: string,
  tryUseGpu: boolean,
): Promise<PathRun> {
  const { provider, data, props, sizeTheme, update } = input;
  const registries = getThemeRegistries();
  let perPath = updatedReprs.get(input);
  if (!perPath) {
    perPath = new Map();
    updatedReprs.set(input, perPath);
  }
  let state = perPath.get(path);
  if (!state) {
    const repr: Representation.Any = provider.factory({ graphics, ...registries }, provider.getParams);
    const themeData = 'grid' in data ? { volume: data } : { structure: data };
    repr.setTheme(
      Theme.create(registries, themeData, { colorTheme: { name: 'uniform' }, sizeTheme: { name: sizeTheme } }),
    );
    const params = PD.getDefaultValues(provider.getParams(registries as never, data as never));
    const first = { ...params, ...props, [update.name]: update.values[0], tryUseGpu };
    await repr.createOrUpdate(first, data).run();
    await graphics.resources.pending.whenReady();
    state = { repr, next: 1 };
    perPath.set(path, state);
  }
  const timer = new StageTimer(graphics);
  const { gpuPasses } = await withGpuPasses(graphics, async () => {
    await state.repr.createOrUpdate({ [update.name]: update.values[state.next] }, data).run();
    state.next = 1 - state.next;
    await timer.stage('update');
    await graphics.resources.pending.whenReady();
    await timer.stage('settled');
  });
  return { stages: timer.stages, triangles: drawnTriangles(state.repr), gpuPasses };
}

function updateCase(
  id: string,
  description: string,
  setup: () => Promise<UpdateInput> | UpdateInput,
): BenchmarkCase<UpdateInput> {
  return {
    id: `mc/repr/${id}`,
    description,
    setup,
    paths: {
      cpu: (graphics, input) => updateRepresentation(graphics, input, 'cpu', false),
      webgl: (graphics, input) => updateRepresentation(graphics, input, 'webgl', true),
      webgpu: (graphics, input) => updateRepresentation(graphics, input, 'webgpu', true),
    },
  };
}

function reprCase(
  id: string,
  description: string,
  setup: () => Promise<ReprInput> | ReprInput,
): BenchmarkCase<ReprInput> {
  return {
    id: `mc/repr/${id}`,
    description,
    setup,
    paths: {
      cpu: (graphics, input) => createRepresentation(graphics, input, false),
      webgl: (graphics, input) => createRepresentation(graphics, input, true),
      webgpu: (graphics, input) => createRepresentation(graphics, input, true),
    },
  };
}

export const RepresentationCases: BenchmarkCase<any>[] = [
  reprCase('isosurface-blobs-128', 'Volume isosurface representation of three blobs on 128³ cells, iso 0.5.', () => ({
    provider: IsosurfaceRepresentationProvider,
    data: blobVolume(128),
    props: { isoValue: Volume.IsoValue.absolute(0.5) },
    sizeTheme: 'uniform',
  })),
  reprCase(
    'gaussian-surface-globule-10k',
    'Gaussian surface representation of a ball of 10000 atoms at a resolution of 1 Å.',
    async () => ({
      provider: GaussianSurfaceRepresentationProvider,
      data: await globuleStructure(10000),
      props: { resolution: 1 },
      sizeTheme: 'physical',
    }),
  ),
  reprCase(
    'gaussian-surface-7qpd',
    'Gaussian surface representation of PDB entry 7qpd (about 11.5 k atoms, a surface per chain) at 1 Å.',
    async () => ({
      provider: GaussianSurfaceRepresentationProvider,
      data: await exampleStructure('7qpd.fw2.cif'),
      props: { resolution: 1 },
      sizeTheme: 'physical',
    }),
  ),
  updateCase(
    'isosurface-blobs-128-update',
    'Volume isosurface representation of three blobs on 128³ cells, the iso value changed between 0.5 and 0.55.',
    () => ({
      provider: IsosurfaceRepresentationProvider,
      data: blobVolume(128),
      props: {},
      sizeTheme: 'uniform',
      update: { name: 'isoValue', values: [Volume.IsoValue.absolute(0.5), Volume.IsoValue.absolute(0.55)] },
    }),
  ),
  updateCase(
    'gaussian-surface-globule-10k-update',
    'Gaussian surface representation of a ball of 10000 atoms at 1 Å, the smoothness changed between 1.5 and 1.6.',
    async () => ({
      provider: GaussianSurfaceRepresentationProvider,
      data: await globuleStructure(10000),
      props: { resolution: 1 },
      sizeTheme: 'physical',
      update: { name: 'smoothness', values: [1.5, 1.6] },
    }),
  ),
  updateCase(
    'gaussian-surface-7qpd-update',
    'Gaussian surface representation of PDB entry 7qpd at 1 Å, the smoothness changed between 1.5 and 1.6.',
    async () => ({
      provider: GaussianSurfaceRepresentationProvider,
      data: await exampleStructure('7qpd.fw2.cif'),
      props: { resolution: 1 },
      sizeTheme: 'physical',
      update: { name: 'smoothness', values: [1.5, 1.6] },
    }),
  ),
];
