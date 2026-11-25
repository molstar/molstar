/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Quat, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { Color } from '@molstar/core/util/color';
import type { Clip } from '@molstar/core/util/clip';
import { Clipping } from '@molstar/graphics/theme/clipping';
import { isExpanded } from '@molstar/graphics/gl/instance-expansion';
import { OrderedSet } from '@molstar/core/data/int';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import {
  ParticleTargetRepresentationParams,
  ParticleTargetRepresentationProvider,
} from '@molstar/graphics/repr/particles/representation/target/representation';
import type { Representation } from '@molstar/graphics/repr/representation';
import { Theme } from '@molstar/graphics/theme/theme';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { Particle, type ParticleList, type ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import { DefaultTestSize, type SceneContext, type VisualTestCase } from '../framework/test-case.js';
import { getThemeRegistries } from '../scenes/structure.js';
import { Impostor } from './geometry.js';

// A list with templates: drawn by its copies' poses on WebGL2 (the reference: the vertex stage composes a template's
// particles with each copy's pose, `instance-expansion.ts`), flattened on WebGL1, by a page store on WebGPU (culled and
// composed on the GPU); the same image

/** three sites a line along `axis` */
const line = (axis: number): ParticleTarget => ({
  kind: 'sites',
  sites: {
    count: 3,
    positions: Float32Array.from({ length: 9 }, (_, k) => (k % 3 === axis ? (Math.floor(k / 3) - 1) * 3 : 0)),
    radii: new Float32Array(3).fill(1.5),
  },
});

function createList(fields: Partial<ParticleList> & Pick<ParticleList, 'count' | 'coordinates' | 'targets'>) {
  return {
    targetInfo: new Map(),
    getParticleLabel: (i: number) => `particle ${i}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    ...fields,
  } as ParticleList;
}

/** 16 particles on a 4² grid of 12 Å, turned, two targets (lines along x and y), an entity a target */
function createTemplate() {
  const n = 16;
  const coordinates = new Float32Array(n * 3);
  const rotations = new Float32Array(n * 4);
  const q = Quat();
  for (let i = 0; i < n; ++i) {
    coordinates.set([(i % 4) * 12, Math.floor(i / 4) * 12, 0], i * 3);
    Quat.setAxisAngle(q, Vec3.create(0, 0, 1), i * 0.5);
    rotations.set([q[0], q[1], q[2], q[3]], i * 4);
  }
  const targets = Int32Array.from({ length: n }, (_, i) => (i < 8 ? 0 : 1));
  return createList({
    label: 'template',
    count: n,
    coordinates,
    rotations,
    targets,
    targetInfo: new Map([
      [0, { name: 'x' }],
      [1, { name: 'y' }],
    ]),
    targetMapping: new Map([
      [0, line(0)],
      [1, line(1)],
    ]),
    entities: Int32Array.from(targets),
    entityInfo: new Map([
      [0, { name: 'A' }],
      [1, { name: 'B' }],
    ]),
  });
}

/**
 * three copies of the template about its center (the second turned 90° about z, the third tilted out of the plane),
 * and a leaf of its own
 */
function createNested() {
  const template = createTemplate();
  const turn = Quat.setAxisAngle(Quat(), Vec3.create(0, 0, 1), Math.PI / 2);
  const tilt = Quat.setAxisAngle(Quat(), Vec3.create(1, 1, 0), 0.6);
  return createList({
    label: 'nested',
    count: 4,
    coordinates: Float32Array.of(0, 0, 0, 60, 0, 0, 0, 60, 0, 60, 60, 0),
    rotations: Float32Array.of(0, 0, 0, 1, ...turn, ...tilt, 0, 0, 0, 1),
    targets: Int32Array.of(0, 0, 0, 1),
    targetInfo: new Map([
      [0, { name: 'copy' }],
      [1, { name: 'leaf' }],
    ]),
    targetMapping: new Map<number, ParticleTarget>([
      [0, { kind: 'particles', particles: template, center: Vec3.create(18, 18, 0) }],
      [1, line(0)],
    ]),
    entities: Int32Array.of(0, 0, 0, 0),
    entityInfo: new Map([[0, { name: 'T' }]]),
  });
}

/** a site of radius 3 */
const dot: ParticleTarget = {
  kind: 'sites',
  sites: { count: 1, positions: new Float32Array(3), radii: Float32Array.of(3) },
};

/**
 * 8³ particles 8 Å apart in a cube, two targets (a half each, an entity a target): a template of a coarse level (every
 * eighth of a target, `getTemplateLevels`)
 */
function createCube() {
  const n = 512;
  const coordinates = new Float32Array(n * 3);
  for (let i = 0; i < n; ++i)
    coordinates.set([(i % 8) * 8, (Math.floor(i / 8) % 8) * 8, Math.floor(i / 64) * 8], i * 3);
  const targets = Int32Array.from({ length: n }, (_, i) => (i % 2 ? 1 : 0));
  return createList({
    label: 'cube',
    count: n,
    coordinates,
    targets,
    targetInfo: new Map([
      [0, { name: 'a' }],
      [1, { name: 'b' }],
    ]),
    targetMapping: new Map([
      [0, dot],
      [1, dot],
    ]),
    entities: Int32Array.from(targets),
    entityInfo: new Map([
      [0, { name: 'A' }],
      [1, { name: 'B' }],
    ]),
  });
}

/** copies of the cube receding from the camera (`FarView`): the near ones whole, the far ones at its coarse level */
function createFar() {
  const cube = createCube();
  const depths = [0, -250, -500, -900, -1400, -2000];
  const coordinates = new Float32Array(depths.length * 3);
  depths.forEach((z, i) => coordinates.set([(i % 2 ? 1 : -1) * (60 + i * 50), (i - 2.5) * 40, z], i * 3));
  const turn = Quat.setAxisAngle(Quat(), Vec3.create(1, 1, 0), 0.5);
  const rotations = new Float32Array(depths.length * 4);
  for (let i = 0; i < depths.length; ++i) rotations.set(i % 2 ? turn : [0, 0, 0, 1], i * 4);
  return createList({
    label: 'far',
    count: depths.length,
    coordinates,
    rotations,
    targets: new Int32Array(depths.length),
    targetInfo: new Map([[0, { name: 'copy' }]]),
    targetMapping: new Map<number, ParticleTarget>([
      [0, { kind: 'particles', particles: cube, center: Vec3.create(28, 28, 28) }],
    ]),
    entities: new Int32Array(depths.length),
    entityInfo: new Map([[0, { name: 'T' }]]),
  });
}

/** the wiggle by position the layers' cases animate (at time 0: a fixed displacement), its amplitude the layer's */
const LayerAnimation = { wiggleMode: 'position', wiggleSpeed: 1, wiggleAmplitude: 0, wiggleFrequency: 0.3 } as const;

async function nestedRepr(graphics: GraphicsContext, list: ParticleList, impostor: boolean, clip?: Clip.Props) {
  const registries = getThemeRegistries();
  const provider = ParticleTargetRepresentationProvider;
  const repr = provider.factory({ graphics, ...registries }, provider.getParams);
  repr.setTheme(
    Theme.create(
      registries,
      // by its leaves by level (no copy expanded), as the plugin themes it
      { particles: ParticleNesting.getThemeList(list) },
      { colorTheme: { name: 'particle-entity', params: {} }, sizeTheme: { name: 'uniform', params: {} } },
    ),
  );
  const defaults = PD.getDefaultValues(ParticleTargetRepresentationParams);
  await repr
    .createOrUpdate(
      {
        ...defaults,
        structure: { ...defaults.structure, tryUseImpostor: impostor },
        ...(clip ? { clip, animation: LayerAnimation } : {}),
      },
      list,
    )
    .run();
  return repr as unknown as Representation.Any;
}

/** the view given, not fitted: the bounds differ by the form (an object's sphere a template's leaf target's) */
const NestedView = {
  target: Vec3.create(30, 30, 0),
  position: Vec3.create(30, 30, 160),
  up: Vec3.create(0, 1, 0),
  radius: 70,
  radiusMax: 70,
};

/** down the receding copies of `createFar` */
const FarView = {
  target: Vec3.create(0, 0, -1000),
  position: Vec3.create(0, 0, 400),
  up: Vec3.create(0, 1, 0),
  radius: 1400,
  radiusMax: 1400,
};

function nestedCase(
  name: string,
  impostor: boolean,
  description: string,
  list = createNested,
  camera = NestedView,
  /** the top particles picked at least, how far a pick from its instance (a line of three sites 3 Å apart, of 1.5) */
  picks = { tops: 4, reach: 6 },
  /** whether far copies draw the template's coarse level (checked where drawn by copies) */
  coarse = false,
): VisualTestCase {
  return {
    id: `particles/nested/${name}`,
    description,
    requires: impostor ? Impostor : undefined,
    camera,
    scene: async ({ graphics }) => {
      const repr = await nestedRepr(graphics, list(), impostor);
      // the reference's form: a template's leaf targets drawn by their copies' poses, not flattened
      if (graphics.capabilities.instanceExpansion && !repr.renderObjects.some((ro) => isExpanded(ro.values))) {
        throw new Error('no object drawn by its copies: the list was flattened');
      }
      return [repr];
    },
    interact: async (ctx, reprs) => {
      const errors = await pickInstances(ctx, reprs, picks.tops, picks.reach);
      // the view's last cull: a coarse level's object saw copies, the template whole's too
      if (coarse && ctx.graphics.capabilities.instanceExpansion) {
        const seen = new Set<number>();
        for (const ro of reprs[0].renderObjects) {
          if (!isExpanded(ro.values)) continue;
          const e = ro.values.instanceExpansion.ref.value;
          if (e.seen.count > 0) seen.add(e.level);
        }
        if (!seen.has(0) || seen.size < 2)
          errors.push(`copies seen at levels ${[...seen].join(', ')}, not 0 and coarse`);
      }
      return errors;
    },
  };
}

/**
 * Picked on a grid over the canvas: each pick one instance by its path (a copy's template particle, or the leaf), its
 * position the flattened list's instance's within `reach` of the picked point; at least `tops` top particles picked
 */
async function pickInstances(
  ctx: SceneContext,
  reprs: Representation.Any[],
  tops: number,
  reach: number,
): Promise<string[]> {
  // WebGPU's: a page store's picks resolved on the GPU once its kernel compiled (frames after the first), tested where
  // they wait for it (`webgpu-particle-nested.test.ts`)
  if (ctx.graphics.compute) return [];
  const errors: string[] = [];
  const copies = new Set<number>();
  let picks = 0;
  const p = Vec3();
  const { width, height } = DefaultTestSize;
  for (let y = 2; y < height; y += 4) {
    for (let x = 2; x < width; x += 4) {
      const picked = await ctx.canvas3d.asyncIdentify(Vec2.create(x, y))?.wait();
      if (!picked) continue;
      const loci = reprs[0].getLoci(picked.id);
      if (!Particle.isLoci(loci)) {
        if (errors.length < 3) errors.push(`pick at ${x}, ${y}: no particle loci (${loci.kind})`);
        continue;
      }
      const instances: [number, number][] = [];
      ParticleNesting.forEachInstance(loci, (top, local) => instances.push([top, local]));
      if (instances.length !== 1) {
        if (errors.length < 3) errors.push(`pick at ${x}, ${y}: ${instances.length} instances`);
        continue;
      }
      picks++;
      const [top, local] = instances[0];
      copies.add(top);
      const list = loci.particles;
      const flat = ParticleNesting.flatten(list);
      Particle.position(p, Particle.Location(flat, ParticleNesting.getExpandedIndex(list, top, local)));
      if (Vec3.distance(p, picked.position) > reach) {
        const d = Vec3.distance(p, picked.position).toFixed(1);
        if (errors.length < 3) errors.push(`pick at ${x}, ${y}: instance ${top}/${local} ${d} Å away`);
      }
    }
  }
  if (picks === 0) errors.push('nothing picked');
  else if (copies.size < tops) errors.push(`picked top particles ${[...copies].join(', ')}, not ${tops}`);
  return errors;
}

/**
 * Marked by path: one instance of the first copy selected, the third copy highlighted whole, a template particle
 * selected in every copy, the leaf highlighted; the others not
 */
function markedCase(name: string, impostor: boolean, description: string): VisualTestCase {
  return {
    id: `particles/nested/marked-${name}`,
    description,
    requires: impostor ? Impostor : undefined,
    camera: NestedView,
    scene: async ({ graphics }) => [await nestedRepr(graphics, createNested(), impostor)],
    interact: async (ctx, reprs) => {
      const repr = reprs[0];
      const list = (repr.getAllLoci()[0] as Particle.Loci).particles;
      const template = ParticleNesting.getTemplate(list, 0)!.particles;
      const errors: string[] = [];
      const marks: [Particle.Loci, MarkerAction][] = [
        [ParticleNesting.ofInstance(list, 0, 3), MarkerAction.Select],
        [ParticleNesting.ofCopy(list, 2), MarkerAction.Highlight],
        [ParticleNesting.ofTemplateParticles(list, template, OrderedSet.ofSingleton(9)), MarkerAction.Select],
        [ParticleNesting.ofCopy(list, 3), MarkerAction.Highlight],
      ];
      // through the canvas (its scene's marker state updated)
      for (const [loci, action] of marks) ctx.canvas3d.mark({ loci, repr }, action);
      return errors;
    },
  };
}

/** a plane through `position` facing `axis`, by instance: what lies on its one side culled whole */
function planeClip(position: Vec3, axis: Vec3, angle: number): Clip.Props {
  return {
    variant: 'instance',
    objects: [
      {
        type: 'plane',
        invert: false,
        position,
        rotation: { axis, angle },
        scale: Vec3.create(1, 1, 1),
        transform: Mat4.identity(),
      },
    ],
  };
}

/** the layers of `marks` (an overpaint, a transparency, a clipping exempting from the clip object) */
function getLayers(layers: {
  red?: Particle.Loci[];
  green?: Particle.Loci[];
  transparent?: Particle.Loci[];
  unclipped?: Particle.Loci[];
  wiggled?: Particle.Loci[];
}): Partial<Representation.State> {
  const kind = 'particle-loci' as const;
  return {
    overpaint: {
      kind,
      layers: [
        ...(layers.red ?? []).map((loci) => ({ loci, color: Color(0xee2211), clear: false })),
        ...(layers.green ?? []).map((loci) => ({ loci, color: Color(0x22cc33), clear: false })),
      ],
    },
    transparency: { kind, layers: (layers.transparent ?? []).map((loci) => ({ loci, value: 0.7 })) },
    wiggle: { kind, layers: (layers.wiggled ?? []).map((loci) => ({ loci, value: 1 })) },
    // exempt from the clip object (the first: its flag the objects it is exempt from)
    clipping: { kind, layers: (layers.unclipped ?? []).map((loci) => ({ loci, groups: Clipping.Groups.Flag.One })) },
  };
}

/**
 * A flat list's layers by particle: two of the first row overpainted, the third row transparent, the half past the
 * plane culled but two of them
 */
const flatLayersCase: VisualTestCase = {
  id: 'particles/layers',
  description:
    'A flat list of 16 particles, its layers by particle: two of its first row overpainted red, its third row ' +
    'transparent, the half past a plane clipped but two exempt (a clipping layer), two wiggled (at time 0).',
  requires: Impostor,
  canvasProps: { renderer: { enableAnimation: true } },
  camera: {
    target: Vec3.create(18, 18, 0),
    position: Vec3.create(18, 18, 110),
    up: Vec3.create(0, 1, 0),
    radius: 40,
    radiusMax: 40,
  },
  scene: async ({ graphics }) => {
    const list = createTemplate();
    const repr = await nestedRepr(graphics, list, true, planeClip(Vec3.create(20, 0, 0), Vec3.create(0, 0, 1), 90));
    const of = (...i: number[]) => Particle.Loci(list, OrderedSet.ofSortedArray(Int32Array.from(i)));
    repr.setState(
      getLayers({ red: [of(2, 3)], transparent: [of(8, 9, 10, 11)], unclipped: [of(4, 12)], wiggled: [of(14, 15)] }),
    );
    return [repr];
  },
};

/**
 * A list with templates' layers by path: a copy whole and a template particle in every copy overpainted, the leaf
 * transparent, the copies past a plane clipped but two instances of one
 */
function layersCase(name: string, impostor: boolean, description: string): VisualTestCase {
  return {
    id: `particles/nested/layers-${name}`,
    description,
    requires: impostor ? Impostor : undefined,
    canvasProps: { renderer: { enableAnimation: true } },
    camera: NestedView,
    scene: async ({ graphics }) => {
      const list = createNested();
      const clip = planeClip(Vec3.create(0, 0, 0), Vec3.create(0, 0, 1), 90);
      const repr = await nestedRepr(graphics, list, impostor, clip);
      const template = ParticleNesting.getTemplate(list, 0)!.particles;
      repr.setState(
        getLayers({
          red: [ParticleNesting.ofCopy(list, 1)],
          green: [ParticleNesting.ofTemplateParticles(list, template, OrderedSet.ofSingleton(6))],
          transparent: [Particle.Loci(list, OrderedSet.ofSingleton(3))],
          unclipped: [ParticleNesting.ofInstance(list, 2, 0), ParticleNesting.ofInstance(list, 2, 12)],
          wiggled: [ParticleNesting.ofCopy(list, 0)],
        }),
      );
      return [repr];
    },
  };
}

export const ParticleCases: VisualTestCase[] = [
  nestedCase(
    'spheres',
    true,
    'A list with templates (three copies, two turned, and a leaf), its sites as sphere impostors: by its copies on ' +
      'WebGL2, flattened on WebGL1, a page store on WebGPU.',
  ),
  nestedCase(
    'mesh',
    false,
    'A list with templates (three copies, two turned, and a leaf), its sites as meshes: by its copies on WebGL2, ' +
      'flattened on WebGL1 and WebGPU.',
  ),
  nestedCase(
    'far',
    true,
    'Copies of a cube of 512 particles receding from the camera: the near ones whole, the far ones at its coarse ' +
      'level (every eighth particle, drawn larger), on WebGL2 by its copies, on WebGPU by a page store.',
    createFar,
    FarView,
    // the far copies' sample particles drawn larger (radius 6 at the coarse level), the farthest a pixel or two
    { tops: 4, reach: 9 },
    true,
  ),
  markedCase(
    'spheres',
    true,
    'Marked by path (an instance of one copy, a copy whole, a template particle in every copy, the leaf): only those ' +
      'drawn marked; on WebGL2 by its copies (their markers apart), on WebGPU a page store.',
  ),
  markedCase(
    'mesh',
    false,
    'Marked by path (an instance of one copy, a copy whole, a template particle in every copy, the leaf), its sites ' +
      'as meshes: on WebGL2 by its copies (their markers apart), flattened on WebGPU.',
  ),
  flatLayersCase,
  layersCase(
    'spheres',
    true,
    'Layers by path (a copy overpainted whole, a template particle in every copy, the leaf transparent, two ' +
      'instances of a clipped copy exempt): on WebGL2 by its copies, on WebGPU a page store.',
  ),
  layersCase(
    'mesh',
    false,
    'Layers by path (a copy overpainted whole, a template particle in every copy, the leaf transparent, two ' +
      'instances of a clipped copy exempt), its sites as meshes: on WebGL2 by its copies, flattened on WebGPU.',
  ),
];
