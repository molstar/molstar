/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Quat, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { OrderedSet } from '@molstar/core/data/int';
import { now } from '@molstar/core/util/now';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { Particle, type ParticleList, type ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ParticleTargetRepresentationParams } from '@molstar/graphics/repr/particles/representation/target/common';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import type { MergedGraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { RendererParams } from '@molstar/graphics/gl/renderer';
import { WebGPU } from './webgpu-module.js';

// A list with templates drawn by a page store's object (`paged.ts`, `addNested`): a template's particles stored once,
// each copy's instances culled and their poses composed on the GPU (`paged-cull.wgsl.ts`); its frames its flattened
// list's (`ParticleNesting.flatten`) drawn flat.

// Medium: scenarios of many frames each, a page store's culling settling
// @module-tag medium

const Size = 96;

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

/**
 * 16 particles on a 4² grid of 12 Å, turned, two targets (lines along x and y), an entity a target: none overlapping
 * (spheres of one depth overlapping drawn in another order elsewhere)
 */
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
 * three copies of the template (the second turned 90° about z), about its center, and a leaf of target 1; all moved
 * by `offset` along x
 */
function createNested(offset = 0) {
  const template = createTemplate();
  const turn = Quat.setAxisAngle(Quat(), Vec3.create(0, 0, 1), Math.PI / 2);
  return createList({
    label: 'nested',
    count: 4,
    coordinates: Float32Array.of(offset, 0, 0, offset + 60, 0, 0, offset, 60, 0, offset + 60, 60, 0),
    rotations: Float32Array.of(0, 0, 0, 1, ...turn, 0, 0, 0, 1, 0, 0, 0, 1),
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

const ListOf = PluginStateTransform.BuiltIn({
  name: 'test-particle-nested-list-of',
  from: SO.Root,
  to: SO.Particle.List,
  params: { list: PD.Value<ParticleList>(undefined as any, { isHidden: true }) },
})({
  apply: ({ params }) => new SO.Particle.List(params.list, { label: 'nested' }),
});

describe.skipIf(!WebGPU.module)('a list with templates as a page store on WebGPU', () => {
  it("draws its flattened list's frames: a template's particles once, its copies' culled and composed on the GPU", async () => {
    const spec = DefaultPluginSpec();
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    // Node has no `ImageData` (`readDrawingBuffer` returns one); the setup restores stubbed globals after the file
    if (typeof ImageData === 'undefined') {
      vi.stubGlobal(
        'ImageData',
        class {
          constructor(
            readonly data: Uint8ClampedArray,
            readonly width: number,
            readonly height: number,
          ) {}
        },
      );
    }
    try {
      await plugin.init();
      const canvas = plugin.canvas3d!;
      canvas.setProps({
        camera: { mode: 'orthographic' } as never,
        cameraFog: { name: 'off', params: {} },
        multiSample: { mode: 'off' } as never,
        postprocessing: { occlusion: { name: 'off', params: {} }, outline: { name: 'off', params: {} } },
      });
      const graphics = canvas.graphics as WebGPUContext;
      let time = now();
      /**
       * the canvas' frame at `focus` once its work settled (the culling's kernels compiled: frames before draw none;
       * pages written as the cull wants them): three frames alike, nothing started meanwhile; its near plane `cut` in
       * front of the focus where given
       */
      const frame = async (focus: Vec3, radius: number, cut?: number) => {
        const { pending } = graphics.resources;
        let image: number[] = [];
        let alike = 0;
        for (let i = 0; i < 12 && alike < 2; ++i) {
          const tracked = pending.trackedCount;
          canvas.commit(true);
          const state = canvas.camera.getFocus(focus, radius);
          canvas.camera.setState(cut === undefined ? state : { ...state, radius: cut }, 0);
          canvas.requestDraw();
          canvas.tick((time += 100) as now.Timestamp, { isSynchronous: true });
          await pending.whenReady();
          const next = Array.from((await graphics.readDrawingBuffer()).data);
          const same = pending.trackedCount === tracked && next.every((v, k) => v === image[k]);
          alike = same ? alike + 1 : 0;
          image = next;
        }
        return image;
      };
      /** `list` drawn as targets, its structure targets as impostors (the store's) or meshes (`impostors` false) */
      const draw = async (list: ParticleList, impostors = true) => {
        const data = await plugin.build().toRoot().apply(ListOf, { list }).commit();
        const structure = { ...PD.getDefaultValues(ParticleTargetRepresentationParams).structure };
        structure.tryUseImpostor = impostors;
        const repr = await plugin
          .build()
          .to(data.ref)
          .apply(ParticlesRepresentation3D, {
            type: { name: 'target', params: { lodPixelSize: 0, structure } },
            colorTheme: { name: 'particle-entity', params: {} },
          } as never)
          .commit();
        return { data, repr };
      };
      /** pixels differing (a composed pose's quaternion snorm16: an edge pixel off at most) */
      const differing = (a: number[], b: number[]) => {
        let n = 0;
        for (let i = 0; i < a.length; i += 4) {
          if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 6) n += 1;
        }
        return n;
      };
      const whole: [Vec3, number] = [Vec3.create(30, 30, 0), 60];
      // zoomed on the second copy: the others out of view whole, culled by their groups
      const zoomed: [Vec3, number] = [Vec3.create(60, 0, 0), 12];

      // the reference: the flattened list of another (the same) nested list, drawn flat
      const flat = await draw(ParticleNesting.flatten(createNested()));
      const flatWhole = await frame(...whole);
      const flatZoomed = await frame(...zoomed);
      expect(new Set(flatWhole).size).toBeGreaterThan(3);
      expect(differing(flatZoomed, flatWhole)).toBeGreaterThan(50);
      await plugin.build().delete(flat.data.ref).commit();

      // drawn as it is: never flattened (its copies expanded on the GPU alone), themed by its leaves by level
      const flatten = vi.spyOn(ParticleNesting, 'flatten');
      const nested = createNested();
      const { data, repr } = await draw(nested);
      const r = repr.data!.repr;
      expect(r.renderObjects).toHaveLength(1);
      const object = r.renderObjects[0] as MergedGraphicsRenderObject;
      expect(object.merged.source?.paged).toBeDefined();
      const nestedWhole = await frame(...whole);
      expect(differing(nestedWhole, flatWhole)).toBeLessThanOrEqual(8);
      expect(differing(await frame(...zoomed), flatZoomed)).toBeLessThanOrEqual(8);
      expect(flatten.mock.calls.some(([list]) => list === nested)).toBe(false);

      // picks through the canvas: a copy's instance by its path (its copy the pick's, gathered in its draw), the
      // flattened list's instance of that path the nearest the hit; a leaf by its index
      const reference = ParticleNesting.flatten(createNested());
      const picks = async (focus: Vec3, radius: number, step: number) => {
        await frame(focus, radius);
        const tops: number[] = [];
        for (let y = 1; y < Size; y += step) {
          for (let x = 1; x < Size; x += step) {
            let hit = canvas.identify(Vec2.create(x, y));
            for (let i = 0; i < 50 && hit === 'pending'; ++i) {
              await new Promise((resolve) => setTimeout(resolve, 10));
              hit = canvas.identify(Vec2.create(x, y));
            }
            if (!hit || hit === 'pending' || hit.id.objectId !== object.id) continue;
            const loci = r.getLoci(hit.id) as Particle.Loci;
            expect(Particle.isLoci(loci)).toBe(true);
            expect(loci.particles).toBe(nested);
            expect(Particle.lociSize(loci)).toBe(1);
            const top = loci.partial ? loci.partial.tops[0] : OrderedSet.start(loci.indices);
            const local = loci.partial ? OrderedSet.start(loci.partial.locals[0]) : 0;
            // the group the geometry's (a site of the line of three), the copy apart (the copies the top particles 0-2)
            expect(hit.id.groupId).toBeLessThan(3);
            expect(hit.id.copyId).toBe(loci.partial ? top : undefined);
            tops.push(top);
            // the instance nearest the hit (its line of sites 3 Å each way, 1.5 Å spheres; the instances 12 Å apart,
            // the copies 60: the hit's position by the pick's depth, coarser than the frame's)
            const index = ParticleNesting.getExpandedIndex(nested, top, local);
            const at = Vec3.fromArray(Vec3(), reference.coordinates as Float32Array, index * 3);
            let nearest = Infinity;
            for (let i = 0; i < reference.count; ++i) {
              const p = Vec3.fromArray(Vec3(), reference.coordinates as Float32Array, i * 3);
              nearest = Math.min(nearest, Vec3.distance(p, hit.position));
            }
            expect(Vec3.distance(at, hit.position)).toBeLessThan(nearest + 1);
            expect(Vec3.distance(at, hit.position)).toBeLessThan(8);
          }
        }
        return tops;
      };
      // every copy's instances picked by their copy
      expect([...new Set(await picks(...whole, 2))].sort()).toEqual([0, 1, 2]);
      // the leaf near (a line two pixels across seen whole: smaller than the pick's texels, at a quarter of the size)
      expect(new Set(await picks(Vec3.create(60, 60, 0), 8, 4))).toEqual(new Set([3]));
      await frame(...whole);

      // its loci by path: all of it its top particles whole; its leaf marked by its store index
      const [all] = r.getAllLoci() as Particle.Loci[];
      expect(all.particles).toBe(nested);
      expect(Particle.lociSize(all)).toBe(3 * 16 + 1);
      const markers = object.values as unknown as { markerAverage: { ref: { value: number } } };
      expect(r.mark(ParticleNesting.ofCopy(nested, 3), MarkerAction.Highlight)).toBe(true);
      expect(markers.markerAverage.ref.value).toBeGreaterThan(0);
      expect(r.mark(ParticleNesting.ofCopy(nested, 3), MarkerAction.RemoveHighlight)).toBe(true);
      expect(markers.markerAverage.ref.value).toBe(0);

      // a copy's instances marked by their copy (a template particle drawn once a copy from one store index): one
      // instance by its path its copy's alone, the same particle of another copy elsewhere, a copy whole its own
      const unmarked = await frame(...whole);
      /** the pixels of `image` other than unmarked */
      const changedOf = (image: number[]) => {
        const changed = new Set<number>();
        for (let i = 0; i < image.length; i += 4) {
          const d =
            Math.abs(image[i] - unmarked[i]) +
            Math.abs(image[i + 1] - unmarked[i + 1]) +
            Math.abs(image[i + 2] - unmarked[i + 2]);
          if (d > 6) changed.add(i / 4);
        }
        return changed;
      };
      /** the pixels `loci` highlighted changes, cleared after */
      const highlighted = async (loci: Particle.Loci) => {
        // through the canvas (its marks update the objects they changed)
        canvas.mark({ loci, repr: r }, MarkerAction.Highlight);
        const image = await frame(...whole);
        expect(markers.markerAverage.ref.value).toBeGreaterThan(0);
        canvas.mark({ loci, repr: r }, MarkerAction.RemoveHighlight);
        await frame(...whole);
        expect(markers.markerAverage.ref.value).toBe(0);
        return changedOf(image);
      };
      const one = await highlighted(ParticleNesting.ofInstance(nested, 1, 5));
      const other = await highlighted(ParticleNesting.ofInstance(nested, 0, 5));
      const copy = await highlighted(ParticleNesting.ofCopy(nested, 1));
      expect(one.size).toBeGreaterThan(0);
      expect(other.size).toBeGreaterThan(0);
      expect([...one].some((i) => other.has(i))).toBe(false);
      // the copy's sixteen instances, its one among them
      expect(copy.size).toBeGreaterThan(one.size * 8);
      expect([...one].every((i) => copy.has(i))).toBe(true);
      expect([...other].some((i) => copy.has(i))).toBe(false);
      // the same particle of every copy (an entity's: `ofTemplateParticles`) the template's, in each copy; one instance
      // of it taken away its copy's alone (the template's marks made each copy's first)
      const template = (nested.targetMapping!.get(0) as { particles: ParticleList }).particles;
      const every = ParticleNesting.ofTemplateParticles(nested, template, OrderedSet.ofSingleton(5));
      canvas.mark({ loci: every, repr: r }, MarkerAction.Highlight);
      const everyCopy = changedOf(await frame(...whole));
      expect([...one, ...other].every((i) => everyCopy.has(i))).toBe(true);
      // no copy's block (a copy's first word its flag where it has one) nor marker of its own: the template's alone
      const words = object.merged.source!.paged.copyMarkers!.words;
      expect(words.every((w) => w === 0)).toBe(true);
      canvas.mark({ loci: ParticleNesting.ofInstance(nested, 1, 5), repr: r }, MarkerAction.RemoveHighlight);
      const rest = changedOf(await frame(...whole));
      expect([...one].some((i) => rest.has(i))).toBe(false);
      expect([...other].every((i) => rest.has(i))).toBe(true);
      canvas.mark({ loci: every, repr: r }, MarkerAction.RemoveHighlight);
      await frame(...whole);
      expect(markers.markerAverage.ref.value).toBe(0);
      // cleared: as before
      expect(differing(await frame(...whole), unmarked)).toBe(0);
      await plugin.build().delete(data.ref).commit();

      // its copies far from the template's own frame (where its pages' boxes are), in perspective with the sphere
      // raster taking what it can: drawn as its flattened list is (a member's sphere where its copies are, not its
      // template's pages': a level's impostors not predicted where the raster would take that sphere's all)
      canvas.setProps({
        camera: { mode: 'perspective' } as never,
        renderer: { sphereRaster: 'on', sphereRasterRadius: 64 } as never,
      });
      // its near plane through the copies' spheres (the raster draws none crossing it)
      const far: [Vec3, number, number] = [Vec3.create(1030, 30, 0), 60, 0.5];
      const flatFar = await draw(ParticleNesting.flatten(createNested(1000)));
      const flatFarImage = await frame(...far);
      expect(new Set(flatFarImage).size).toBeGreaterThan(3);
      await plugin.build().delete(flatFar.data.ref).commit();
      const nestedFarList = createNested(1000);
      const nestedFar = await draw(nestedFarList);
      expect(differing(await frame(...far), flatFarImage)).toBeLessThanOrEqual(8);
      // a copy marked whole, its spheres the raster's (the resolve reads its copy's markers too: tinted in the material's
      // color, no marking pass drawing its edges by impostors)
      canvas.setProps({ marking: { enabled: false } as never });
      const farRepr = nestedFar.repr.data!.repr;
      const farWhole: [Vec3, number] = [Vec3.create(1030, 30, 0), 60];
      const farUnmarked = await frame(...farWhole);
      canvas.mark({ loci: ParticleNesting.ofCopy(nestedFarList, 0), repr: farRepr }, MarkerAction.Highlight);
      expect(differing(await frame(...farWhole), farUnmarked)).toBeGreaterThan(20);
      canvas.setProps({ marking: { enabled: true } as never });
      await plugin.build().delete(nestedFar.data.ref).commit();
      canvas.setProps({
        camera: { mode: 'orthographic' } as never,
        renderer: {
          sphereRaster: 'auto',
          sphereRasterRadius: PD.getDefaultValues(RendererParams).sphereRasterRadius,
        } as never,
      });

      // its targets as meshes (not the store's): flattened, themed by its leaves by level all the same
      flatten.mockClear();
      const meshes = createNested();
      const flatMeshes = await draw(ParticleNesting.flatten(createNested()), false);
      const flatMeshImage = await frame(...whole);
      await plugin.build().delete(flatMeshes.data.ref).commit();
      const mesh = await draw(meshes, false);
      expect(
        mesh.repr.data!.repr.renderObjects.some((o) => (o as MergedGraphicsRenderObject).merged?.source?.paged),
      ).toBe(false);
      expect(flatten.mock.calls.some(([list]) => list === meshes)).toBe(true);
      expect(differing(await frame(...whole), flatMeshImage)).toBeLessThanOrEqual(8);

      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("replays a hover's marking redraws from its draw list: the store's draws many, its markers flushed", async () => {
    const spec = DefaultPluginSpec();
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const canvas = plugin.canvas3d!;
      canvas.setProps({
        camera: { mode: 'orthographic' } as never,
        multiSample: { mode: 'off' } as never,
        renderer: { ...canvas.props.renderer, drawLists: true },
      });
      const graphics = canvas.graphics as WebGPUContext;
      // a template of 100 targets, 10 particles each (a site each target, of its own radius) placed twice: the store's
      // draws a target's (their geometry its own), more than a chunk of few rows draws directly
      const targetCount = 100;
      const n = targetCount * 10;
      const coordinates = new Float32Array(n * 3);
      for (let i = 0; i < n; ++i)
        coordinates.set([(i % 10) * 4, (Math.floor(i / 10) % 10) * 4, Math.floor(i / 100) * 4], i * 3);
      const targetMapping = new Map<number, ParticleTarget>();
      const targetInfo = new Map<number, { name: string }>();
      for (let t = 0; t < targetCount; ++t) {
        targetMapping.set(t, {
          kind: 'sites',
          sites: { count: 1, positions: new Float32Array(3), radii: Float32Array.of(1 + t / targetCount) },
        });
        targetInfo.set(t, { name: `site ${t}` });
      }
      const template = createList({
        label: 'template',
        count: n,
        coordinates,
        targets: Int32Array.from({ length: n }, (_, i) => Math.floor(i / 10)),
        targetInfo,
        targetMapping,
        entities: new Int32Array(n),
        entityInfo: new Map([[0, { name: 'A' }]]),
      });
      const copies = 2;
      const nested = createList({
        label: 'nested',
        count: copies,
        coordinates: Float32Array.of(0, 0, 0, 60, 0, 0),
        rotations: Float32Array.of(0, 0, 0, 1, 0, 0, 0, 1),
        targets: new Int32Array(copies),
        targetInfo: new Map([[0, { name: 'copy' }]]),
        targetMapping: new Map<number, ParticleTarget>([
          [0, { kind: 'particles', particles: template, center: Vec3.create(18, 18, 18) }],
        ]),
        entities: new Int32Array(copies),
        entityInfo: new Map([[0, { name: 'T' }]]),
      });
      const data = await plugin.build().toRoot().apply(ListOf, { list: nested }).commit();
      const repr = await plugin
        .build()
        .to(data.ref)
        .apply(ParticlesRepresentation3D, {
          type: { name: 'target', params: { lodPixelSize: 0 } },
          colorTheme: { name: 'particle-entity', params: {} },
        } as never)
        .commit();
      const r = repr.data!.repr;
      let time = now();
      for (let i = 0; i < 12; ++i) {
        canvas.commit(true);
        canvas.camera.setState(canvas.camera.getFocus(Vec3.create(30, 0, 0), 60), 0);
        canvas.requestDraw();
        canvas.tick((time += 100) as now.Timestamp, { isSynchronous: true });
        await graphics.resources.pending.whenReady();
      }
      const { stats } = graphics;
      /** another copy highlighted, the one before not: a marking redraw (no frame requested), its draws counted */
      const hover = async (copy: number) => {
        canvas.mark(
          { loci: ParticleNesting.ofCopy(nested, (copy + copies - 1) % copies), repr: r },
          MarkerAction.RemoveHighlight,
        );
        canvas.mark({ loci: ParticleNesting.ofCopy(nested, copy % copies), repr: r }, MarkerAction.Highlight);
        const calls = stats.calls.drawInstanced + stats.calls.drawInstancedBase;
        canvas.tick((time += 100) as now.Timestamp, { isSynchronous: true });
        await graphics.resources.pending.whenReady();
        return stats.calls.drawInstanced + stats.calls.drawInstancedBase - calls;
      };
      for (let i = 0; i < 4; ++i) await hover(i);
      const before = { ...stats.drawLists };
      const draws: number[] = [];
      for (let i = 4; i < 12; ++i) draws.push(await hover(i));
      const after = stats.drawLists;
      expect(Math.min(...draws)).toBeGreaterThan(64);
      // the store's chunk changing every frame, one row of many draws: walked, its markers flushed, replayed, not
      // drawn directly (its draws' encoding a frame the driver's: more than a flush and a replay)
      expect(after.chunksDirect - before.chunksDirect).toBe(0);
      expect(after.flushed - before.flushed).toBe(8);
      expect(after.recorded - before.recorded).toBe(0);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it("keeps in a view's still frames what its second phase found, the sphere raster's regions too (a zoom stopped)", async () => {
    const spec = DefaultPluginSpec();
    const size = 256;
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
      // large enough for the pyramid to hide copies behind copies
      { width: size, height: size },
      { backend: 'webgpu' },
    );
    if (typeof ImageData === 'undefined') {
      vi.stubGlobal(
        'ImageData',
        class {
          constructor(
            readonly data: Uint8ClampedArray,
            readonly width: number,
            readonly height: number,
          ) {}
        },
      );
    }
    try {
      await plugin.init();
      const canvas = plugin.canvas3d!;
      canvas.setProps({
        camera: { mode: 'perspective' } as never,
        cameraFog: { name: 'off', params: {} },
        multiSample: { mode: 'off' } as never,
        postprocessing: { occlusion: { name: 'off', params: {} }, outline: { name: 'off', params: {} } },
        hiZ: { ...canvas.props.hiZ, enabled: false, gpu: true, twoPhase: true },
      });
      canvas.handleResize();
      const graphics = canvas.graphics as WebGPUContext;
      // a cube of 10^3 sites 3 apart placed by 4x4x4 copies 40 apart: copies behind copies
      const m = 10;
      const n = m * m * m;
      const coordinates = new Float32Array(n * 3);
      for (let i = 0; i < n; ++i) {
        coordinates.set([(i % m) * 3, (Math.floor(i / m) % m) * 3, Math.floor(i / (m * m)) * 3], i * 3);
      }
      const site: ParticleTarget = {
        kind: 'sites',
        sites: { count: 1, positions: new Float32Array(3), radii: Float32Array.of(1.6) },
      };
      const template = createList({
        label: 'template',
        count: n,
        coordinates,
        targets: new Int32Array(n),
        targetMapping: new Map([[0, site]]),
        entities: Int32Array.from({ length: n }, (_, i) => i % 7),
        entityInfo: new Map(Array.from({ length: 7 }, (_, e) => [e, { name: `E${e}` }])),
      });
      const g = 4;
      const copies = g * g * g;
      const nested = createList({
        label: 'nested',
        count: copies,
        coordinates: Float32Array.from({ length: copies * 3 }, (_, k) => {
          const c = Math.floor(k / 3);
          return [c % g, Math.floor(c / g) % g, Math.floor(c / (g * g))][k % 3] * 40;
        }),
        rotations: Float32Array.from({ length: copies * 4 }, (_, k) => (k % 4 === 3 ? 1 : 0)),
        targets: new Int32Array(copies),
        targetMapping: new Map<number, ParticleTarget>([
          [0, { kind: 'particles', particles: template, center: Vec3.create(13.5, 13.5, 13.5) }],
        ]),
        entities: new Int32Array(copies),
        entityInfo: new Map([[0, { name: 'T' }]]),
      });
      const data = await plugin.build().toRoot().apply(ListOf, { list: nested }).commit();
      await plugin
        .build()
        .to(data.ref)
        .apply(ParticlesRepresentation3D, {
          type: { name: 'target', params: { lodPixelSize: 0 } },
          colorTheme: { name: 'particle-entity', params: {} },
        } as never)
        .commit();
      let time = now();
      const tick = async () => {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 100) as now.Timestamp, { isSynchronous: true });
        await graphics.resources.pending.whenReady();
        return Array.from((await graphics.readDrawingBuffer()).data);
      };
      const target = Vec3.create(73.5, 73.5, 73.5);
      const look = (d: number) =>
        canvas.camera.setState(
          {
            target,
            position: Vec3.add(Vec3(), target, Vec3.create(d * 0.3, d * 0.2, d)),
            up: Vec3.create(0, 1, 0),
            radius: 150,
          },
          0,
        );
      const occlusion = (gpu: boolean) => canvas.setProps({ hiZ: { ...canvas.props.hiZ, gpu } });
      const differing = (a: number[], b: number[]) => {
        let k = 0;
        for (let i = 0; i < a.length; i += 4) {
          if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 6) k += 1;
        }
        return k;
      };
      look(600);
      for (let i = 0; i < 12; ++i) await tick();
      // zoomed in and out a step at a time: a step's first frame culled anew (its second phase catching what the last
      // pyramid hid), the next ones the same view (not culled again), as drawn without occlusion culling
      const missing: string[] = [];
      for (const d of [300, 160, 90, 160, 300]) {
        look(d);
        const first = await tick();
        const still = await tick();
        occlusion(false);
        await tick();
        const reference = await tick();
        occlusion(true);
        await tick();
        await tick();
        const a = differing(first, reference);
        const b = differing(still, reference);
        if (a > 8 || b > 8) missing.push(`${d}: first ${a}, still ${b}`);
      }
      expect(missing).toEqual([]);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);
});
