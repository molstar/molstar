/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Mat4, Quat, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { ClipGhostParams } from '@molstar/graphics/canvas3d/passes/clip-ghost';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { createFrame } from '@molstar/graphics/gl/renderable';
import { Renderer, type RendererProps } from '@molstar/graphics/gl/renderer';
import { Scene } from '@molstar/graphics/gl/scene';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { now } from '@molstar/core/util/now';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomProperties } from '@molstar/model/model/custom-property';
import {
  Particle,
  type ParticleList,
  type ParticleTarget,
  ParticleValues,
} from '@molstar/model/model/particles/particle-list';
import { ParticleSlots } from '@molstar/model/model/particles/particle-slots';
import { ValueCell } from '@molstar/core/util/value-cell';

import { OrderedSet } from '@molstar/core/data/int';
import type { MergedGraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { ParticleTileSetFromList, ParticleTilesView } from '@molstar/plugin/state/transforms/particles/tiles';
import { WebGPU } from './webgpu-module.js';

// A slot view of particle tiles as a page store's object (`paged.ts`): its instances culled
// and composed on the GPU (`merged-paged-culler.ts`), drawn by pose region; its frames the flat view's of the same
// tiles, as tiles are let go and added again; a pick its particle.

const Size = 64;
const Count = 64;
const OcclusionSize = 256;

/** 64 particles on a 4³ lattice over 0-90 Å, rotated, two targets of three sites a line (along x, along y) */
function createList(): ParticleList {
  const coordinates = new Float32Array(Count * 3);
  const rotations = new Float32Array(Count * 4);
  const targets = new Int32Array(Count);
  const q = Quat();
  for (let i = 0; i < Count; ++i) {
    const x = i % 4;
    coordinates.set([x * 30, (Math.floor(i / 4) % 4) * 30, Math.floor(i / 16) * 30], i * 3);
    targets[i] = Math.floor(i / 16) < 2 ? 1 : 2;
    Quat.setAxisAngle(q, Vec3.create(0, 0, 1), i * 0.4);
    rotations.set([q[0], q[1], q[2], q[3]], i * 4);
  }
  const line = (axis: number): ParticleTarget => ({
    kind: 'sites',
    sites: {
      count: 3,
      positions: Float32Array.from({ length: 9 }, (_, k) => (k % 3 === axis ? (Math.floor(k / 3) - 1) * 6 : 0)),
      radii: new Float32Array(3).fill(3),
    },
  });
  return {
    label: 'lattice',
    count: Count,
    targets,
    targetInfo: new Map([
      [1, { name: 'one' }],
      [2, { name: 'two' }],
    ]),
    targetMapping: new Map([
      [1, line(0)],
      [2, line(1)],
    ]),
    coordinates,
    rotations,
    getParticleLabel: (i) => `particle ${i}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

/**
 * 16 large spheres (14 Å) at z 90 over a 4² grid of 30 Å, before 48 small ones (2 Å) at z 0, 30 and 60: seen along z
 * the large hide the small, along x the small are in sight
 */
function createOccluderList(): ParticleList {
  const coordinates = new Float32Array(Count * 3);
  const targets = new Int32Array(Count);
  for (let i = 0; i < Count; ++i) {
    const layer = Math.floor(i / 16);
    coordinates.set([(i % 4) * 30, (Math.floor(i / 4) % 4) * 30, layer === 3 ? 90 : layer * 30], i * 3);
    targets[i] = layer === 3 ? 2 : 1;
  }
  const sphere = (r: number): ParticleTarget => ({
    kind: 'sites',
    sites: { count: 1, positions: new Float32Array(3), radii: Float32Array.of(r) },
  });
  return {
    ...createList(),
    targets,
    coordinates,
    rotations: undefined,
    targetMapping: new Map([
      [1, sphere(2)],
      [2, sphere(14)],
    ]),
  };
}

const ListOf = PluginStateTransform.BuiltIn({
  name: 'test-particle-paged-list-of',
  from: SO.Root,
  to: SO.Particle.List,
  params: { list: PD.Value<ParticleList>(undefined as any, { isHidden: true }) },
})({
  apply: ({ params }) => new SO.Particle.List(params.list, { label: 'lattice' }),
});

describe.skipIf(!WebGPU.module)('a slot view of particle tiles as a page store on WebGPU', () => {
  it("draws the flat view's frames, tiles let go and added again; picks a particle by its store index", async () => {
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
      const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();

      let time = now();
      /** the canvas' frame once its work settled (the culling's kernels compiled: frames before draw none) */
      const frame = async () => {
        for (let i = 0; i < 12; ++i) {
          canvas.commit(true);
          canvas.camera.setState(canvas.camera.getFocus(Vec3.create(45, 45, 45), 70), 0);
          canvas.requestDraw();
          canvas.tick((time += 100) as now.Timestamp, { isSynchronous: true });
          await graphics.resources.pending.whenReady();
        }
        return Array.from((await graphics.readDrawingBuffer()).data);
      };
      const view = async (ids: string, slots: boolean) => {
        const v = await plugin
          .build()
          .to(tiles.ref)
          .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids } }, frame: 0, slots })
          .commit();
        const repr = await plugin
          .build()
          .to(v.ref)
          .apply(ParticlesRepresentation3D, {
            type: { name: 'target', params: { lodPixelSize: 0 } },
            colorTheme: { name: 'particle-index', params: {} },
          } as never)
          .commit();
        return { view: v, repr };
      };
      const flatFrame = async (ids: string) => {
        const { view: v } = await view(ids, false);
        const image = await frame();
        await plugin.build().delete(v.ref).commit();
        return image;
      };
      /** pixels differing (none on Dawn): a composed pose's quaternion snorm16, a pixel off at an edge on others */
      const differing = (a: number[], b: number[]) => {
        let n = 0;
        for (let i = 0; i < a.length; i += 4) {
          if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 6) n += 1;
        }
        return n;
      };

      const allFlat = await flatFrame('0-7');
      const someFlat = await flatFrame('0, 2, 5');
      expect(new Set(allFlat).size).toBeGreaterThan(3);
      expect(differing(someFlat, allFlat)).toBeGreaterThan(20);

      const { view: v, repr } = await view('0-7', true);
      const r = repr.data!.repr;
      expect(r.renderObjects).toHaveLength(1);
      const object = r.renderObjects[0] as MergedGraphicsRenderObject;
      expect(object.merged.source?.paged).toBeDefined();
      const allPaged = await frame();
      expect(differing(allPaged, allFlat)).toBeLessThanOrEqual(8);

      const setIds = async (ids: string) => {
        await plugin
          .build()
          .to(v.ref)
          .update(ParticleTilesView, (old) => ({ ...old, tiles: { name: 'ids', params: { ids } } }))
          .commit();
        return frame();
      };
      // tiles let go: their pages free; the object kept
      expect(differing(await setIds('0, 2, 5'), someFlat)).toBeLessThanOrEqual(8);
      expect(r.renderObjects[0]).toBe(object);
      // added again
      expect(differing(await setIds('0-7'), allFlat)).toBeLessThanOrEqual(8);
      expect(r.renderObjects[0]).toBe(object);

      // a pick: an instance's store index its particle's slot
      const slots = v.cell!.obj!.data;
      const loci = r.getLoci({ objectId: object.id, instanceId: 0, groupId: 0 }) as Particle.Loci;
      expect(Particle.isLoci(loci)).toBe(true);
      expect(loci.particles).toBe(slots);
      const slot = OrderedSet.getAt(loci.indices, 0);
      expect(ParticleSlots.findPage(slots, slot)).toBeGreaterThanOrEqual(0);
      // a store index past the pages' particles none
      expect(Particle.isLoci(r.getLoci({ objectId: object.id, instanceId: 1 << 30, groupId: 0 }))).toBe(false);

      // a pick through the canvas: an instance drawn, its particle
      let picked: Particle.Loci | undefined;
      for (let y = 2; y < Size && !picked; y += 4) {
        for (let x = 2; x < Size && !picked; x += 4) {
          let hit = canvas.identify(Vec2.create(x, y));
          for (let i = 0; i < 50 && hit === 'pending'; ++i) {
            await new Promise((r) => setTimeout(r, 10));
            hit = canvas.identify(Vec2.create(x, y));
          }
          if (!hit || hit === 'pending' || hit.id.objectId !== object.id) continue;
          const l = r.getLoci(hit.id);
          if (Particle.isLoci(l)) picked = l;
        }
      }
      expect(picked).toBeDefined();
      expect(ParticleSlots.findPage(slots, OrderedSet.getAt(picked!.indices, 0))).toBeGreaterThanOrEqual(0);

      // marked by store index: the object's markers (a highlight not drawn in a headless frame, the tile store's alike)
      const markers = object.values as unknown as { markerAverage: { ref: { value: number } } };
      expect(r.mark(picked!, MarkerAction.Highlight)).toBe(true);
      expect(markers.markerAverage.ref.value).toBeGreaterThan(0);
      // its tile let go and added again (other pages, or the same reused): its markers clear
      const tileOf = (slot: number) =>
        tiles.data!.set.tiles.values().find((t) => slot >= t.base && slot < t.base + t.count)!.id;
      const tile = tileOf(OrderedSet.getAt(picked!.indices, 0));
      const others = [0, 1, 2, 3, 4, 5, 6, 7].filter((id) => id !== tile).join(',');
      await setIds(others);
      await setIds('0-7');
      expect(markers.markerAverage.ref.value).toBe(0);
      expect(differing(await frame(), allFlat)).toBeLessThanOrEqual(8);

      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("culls by occlusion in two phases as the flat view's culler: the frames of a jump to a view along the lattice's rows the same", async () => {
    const spec = DefaultPluginSpec();
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
      // large enough for the pyramid to resolve a small sphere behind a large one
      { width: OcclusionSize, height: OcclusionSize },
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
        // the GPU's pyramid, two phases (the second the frame's depth so far)
        hiZ: { ...canvas.props.hiZ, enabled: false, gpu: true, twoPhase: true },
      });
      canvas.handleResize();
      const graphics = canvas.graphics as WebGPUContext;
      const list = await plugin.build().toRoot().apply(ListOf, { list: createOccluderList() }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();
      let time = now();
      const tick = async () => {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 100) as now.Timestamp, { isSynchronous: true });
        await graphics.resources.pending.whenReady();
        return Array.from((await graphics.readDrawingBuffer()).data);
      };
      /** a view along z settled, then a jump to one along x: its first frame (the second phase's), and settled */
      const jump = async (slots: boolean) => {
        const v = await plugin
          .build()
          .to(tiles.ref)
          .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids: '0-7' } }, frame: 0, slots })
          .commit();
        await plugin
          .build()
          .to(v.ref)
          .apply(ParticlesRepresentation3D, {
            type: { name: 'target', params: { lodPixelSize: 0 } },
            colorTheme: { name: 'particle-index', params: {} },
          } as never)
          .commit();
        // along z: each column's back particles hidden by its front one (the pyramid's)
        const look = (position: Vec3) =>
          canvas.camera.setState(
            { target: Vec3.create(45, 45, 45), position, up: Vec3.create(0, 1, 0), radius: 80 },
            0,
          );
        look(Vec3.create(45, 45, 400));
        for (let i = 0; i < 12; ++i) await tick();
        // along x: those in sight, culled by the first phase (the last pyramid's), found by the second
        look(Vec3.create(400, 45, 45));
        const first = await tick();
        for (let i = 0; i < 6; ++i) await tick();
        const settled = await tick();
        await plugin.build().delete(v.ref).commit();
        return { first, settled };
      };
      const differing = (a: number[], b: number[]) => {
        let n = 0;
        for (let i = 0; i < a.length; i += 4) {
          if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 6) n += 1;
        }
        return n;
      };
      const flat = await jump(false);
      const paged = await jump(true);
      expect(new Set(flat.settled).size).toBeGreaterThan(3);
      expect(differing(paged.first, flat.first)).toBeLessThanOrEqual(8);
      expect(differing(paged.settled, flat.settled)).toBeLessThanOrEqual(8);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("clips as the flat view's culler, drawing the ghost of what the clip objects remove: by pixel and by instance", async () => {
    const spec = DefaultPluginSpec();
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
      { width: Size, height: Size },
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
        camera: { mode: 'orthographic' } as never,
        cameraFog: { name: 'off', params: {} },
        multiSample: { mode: 'off' } as never,
        postprocessing: { occlusion: { name: 'off', params: {} }, outline: { name: 'off', params: {} } },
      });
      const graphics = canvas.graphics as WebGPUContext;
      const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();
      let time = now();
      const frame = async () => {
        for (let i = 0; i < 12; ++i) {
          canvas.commit(true);
          canvas.camera.setState(canvas.camera.getFocus(Vec3.create(45, 45, 45), 70), 0);
          canvas.requestDraw();
          canvas.tick((time += 100) as now.Timestamp, { isSynchronous: true });
          await graphics.resources.pending.whenReady();
        }
        return Array.from((await graphics.readDrawingBuffer()).data);
      };
      /** a plane through the lattice's middle, tilted: a particle cut by it, others whole on either side */
      const clip = (variant: 'pixel' | 'instance') => ({
        variant,
        objects: [
          {
            type: 'plane' as const,
            invert: false,
            position: Vec3.create(44, 45, 45),
            rotation: { axis: Vec3.create(0, 1, 0.3), angle: 70 },
            scale: Vec3.create(1, 1, 1),
            transform: Mat4.identity(),
          },
        ],
      });
      const image = async (slots: boolean, variant: 'pixel' | 'instance') => {
        const v = await plugin
          .build()
          .to(tiles.ref)
          .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids: '0-7' } }, frame: 0, slots })
          .commit();
        await plugin
          .build()
          .to(v.ref)
          .apply(ParticlesRepresentation3D, {
            type: { name: 'target', params: { lodPixelSize: 0, clip: clip(variant) } },
            colorTheme: { name: 'particle-index', params: {} },
          } as never)
          .commit();
        const ghosted = await frame();
        canvas.setProps({ clipGhost: { name: 'off', params: {} } });
        const plain = await frame();
        canvas.setProps({ clipGhost: { name: 'on', params: { ...PD.getDefaultValues(ClipGhostParams), alpha: 0.5 } } });
        await plugin.build().delete(v.ref).commit();
        return { ghosted, plain };
      };
      const differing = (a: number[], b: number[]) => {
        let n = 0;
        for (let i = 0; i < a.length; i += 4) {
          if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 6) n += 1;
        }
        return n;
      };
      canvas.setProps({ clipGhost: { name: 'on', params: { ...PD.getDefaultValues(ClipGhostParams), alpha: 0.5 } } });
      const whole = await (async () => {
        const v = await plugin
          .build()
          .to(tiles.ref)
          .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids: '0-7' } }, frame: 0, slots: false })
          .commit();
        await plugin
          .build()
          .to(v.ref)
          .apply(ParticlesRepresentation3D, {
            type: { name: 'target', params: { lodPixelSize: 0 } },
            colorTheme: { name: 'particle-index', params: {} },
          } as never)
          .commit();
        const f = await frame();
        await plugin.build().delete(v.ref).commit();
        return f;
      })();
      for (const variant of ['pixel', 'instance'] as const) {
        const flat = await image(false, variant);
        const paged = await image(true, variant);
        // something clipped, and its ghost drawn
        expect(differing(flat.plain, whole)).toBeGreaterThan(20);
        expect(differing(flat.ghosted, flat.plain)).toBeGreaterThan(20);
        expect(differing(paged.plain, flat.plain)).toBeLessThanOrEqual(8);
        expect(differing(paged.ghosted, flat.ghosted)).toBeLessThanOrEqual(8);
      }
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("splits for the sphere raster as the flat view's culler: its regions' instances drawn by the raster, the same image", async () => {
    const spec = DefaultPluginSpec();
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics as WebGPUContext;
      if (!graphics.capabilities.atomicMinMax64) return;
      const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();
      /** a view's objects: the flat view's (an object a target) or the page store's */
      const objectOf = async (slots: boolean) => {
        const v = await plugin
          .build()
          .to(tiles.ref)
          .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids: '0-7' } }, frame: 0, slots })
          .commit();
        const repr = await plugin
          .build()
          .to(v.ref)
          .apply(ParticlesRepresentation3D, {
            type: { name: 'target', params: { lodPixelSize: 0 } },
            colorTheme: { name: 'particle-index', params: {} },
          } as never)
          .commit();
        return repr.data!.repr.renderObjects as MergedGraphicsRenderObject[];
      };
      const flatObjects = await objectOf(false);
      const pagedObjects = await objectOf(true);
      expect(pagedObjects).toHaveLength(1);
      expect(pagedObjects[0].merged.source?.paged).toBeDefined();
      expect(flatObjects.some((o) => o.merged?.source?.paged)).toBe(false);

      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, Size, Size);
      camera.setState({ ...camera.getFocus(Vec3.create(45, 45, 45), 70), radiusMax: 200 }, 0);
      camera.update();
      const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
      /** the objects' frame by a renderer of `props` in a scene of their own; their raster's inputs and resolves */
      const render = async (objects: MergedGraphicsRenderObject[], props: Partial<RendererProps>) => {
        const scene = Scene.create(graphics);
        for (const o of objects) scene.add(o);
        scene.commit();
        type Input = ReturnType<NonNullable<(typeof scene.renderables)[0]['getRasterInput']>>;
        const inputs = new Map<number, Input>();
        const restores: (() => void)[] = [];
        let resolved = 0;
        for (const [k, renderable] of scene.renderables.entries()) {
          const getRasterInput = renderable.getRasterInput;
          if (getRasterInput) {
            renderable.getRasterInput = (...args) => {
              const input = getRasterInput(...args) ?? inputs.get(k);
              inputs.set(k, input);
              return input;
            };
          }
          const item = renderable.renderItem;
          const renderFullscreen = item.renderFullscreen!;
          item.renderFullscreen = (...args) => {
            const drawn = renderFullscreen(...args);
            if (drawn) resolved += 1;
            return drawn;
          };
          restores.push(() => {
            renderable.getRasterInput = getRasterInput;
            item.renderFullscreen = renderFullscreen;
          });
        }
        const renderer = Renderer.create(graphics, props);
        await graphics.resources.preparePipelines(
          renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
        );
        const frame = async () => {
          resolved = 0;
          renderer.setViewport(0, 0, Size, Size);
          renderer.update(camera, scene, createFrame());
          renderer.prepareCull(scene);
          renderer.prepareSphereRaster(scene, camera);
          graphics.beginRenderPass(target);
          renderer.clear(true);
          renderer.renderOpaque(scene, camera);
          graphics.flush();
          return Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
        };
        let pixels: number[] = [];
        for (let i = 0; i < 30; ++i) {
          pixels = await frame();
          await graphics.resources.pending.whenReady();
          if (props.sphereRaster === 'off' ? i > 3 : resolved > 0 && i > 3) break;
        }
        pixels = await frame();
        // the raster's regions' instances, by its records
        let rastered = 0;
        for (const input of inputs.values()) {
          if (!input) continue;
          const records = new Uint32Array(await input.records.readAsync());
          for (let r = input.recordBase; r < 2 * input.recordBase; ++r) rastered += records[r * input.recordWords + 1];
        }
        for (const restore of restores) restore();
        renderer.dispose();
        scene.clear();
        return { pixels, resolved, rastered };
      };
      const differing = (a: number[], b: number[]) => {
        let n = 0;
        for (let i = 0; i < a.length; i += 4) {
          if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 24) n += 1;
        }
        return n;
      };
      // every sphere small enough for the raster
      const on: Partial<RendererProps> = { sphereRaster: 'on', sphereRasterRadius: 64 };
      const flat = await render(flatObjects, on);
      const paged = await render(pagedObjects, on);
      // a later scene's culler of the same store too (its first sync every page: the changes taken by the first)
      const impostors = await render(pagedObjects, { sphereRaster: 'off' });
      expect(new Set(impostors.pixels).size).toBeGreaterThan(4);
      expect(flat.resolved).toBeGreaterThan(0);
      expect(paged.resolved).toBeGreaterThan(0);
      // every particle in view the raster's (the flat view's objects split by the raster per instance)
      expect(paged.rastered).toBe(Count);
      expect(differing(paged.pixels, flat.pixels)).toBeLessThanOrEqual(8);
      // the raster's image the impostors', but where its arithmetic and the fragment stage's decide an edge differently
      expect(differing(paged.pixels, impostors.pixels)).toBeLessThanOrEqual(Size * Size * 0.02);
      target.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it('culls on the GPU for views of their own: a frame and camera of its own, a second camera of a frame', async () => {
    const spec = DefaultPluginSpec();
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics as WebGPUContext;
      const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();
      const v = await plugin
        .build()
        .to(tiles.ref)
        .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids: '0-7' } }, frame: 0, slots: true })
        .commit();
      const repr = await plugin
        .build()
        .to(v.ref)
        .apply(ParticlesRepresentation3D, {
          type: { name: 'target', params: { lodPixelSize: 0 } },
          colorTheme: { name: 'particle-index', params: {} },
        } as never)
        .commit();
      const object = repr.data!.repr.renderObjects[0] as MergedGraphicsRenderObject;
      expect(object.merged.source?.paged).toBeDefined();
      const scene = Scene.create(graphics);
      scene.add(object);
      scene.commit();

      const cameraAt = (target: Vec3, radius: number, position: Vec3) => {
        const camera = new Camera();
        Viewport.set(camera.viewport, 0, 0, Size, Size);
        camera.setState({ ...camera.getFocus(target, radius), position, up: Vec3.create(0, 1, 0), radiusMax: 200 }, 0);
        camera.update();
        return camera;
      };
      // the main camera at a corner of the lattice (most particles out of its view), the other at all of it
      const main = cameraAt(Vec3.create(0, 0, 0), 12, Vec3.create(0, 0, 40));
      const other = cameraAt(Vec3.create(45, 45, 45), 70, Vec3.create(240, 80, 120));
      const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
      const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
      await graphics.resources.preparePipelines(
        renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
      );
      /** a frame's view of `camera`: culled by `prepare` (the frame's, a view of its own's, none) */
      const draw = async (
        camera: Camera,
        frame: ReturnType<typeof createFrame>,
        prepare: 'frame' | 'view' | 'none',
      ) => {
        renderer.setViewport(0, 0, Size, Size);
        renderer.update(camera, scene, frame);
        if (prepare === 'frame') renderer.prepareCull(scene);
        else if (prepare === 'view') renderer.prepareCullView(scene);
        graphics.beginRenderPass(target);
        renderer.clear(true);
        renderer.renderOpaque(scene, camera);
        graphics.flush();
        await graphics.resources.pending.whenReady();
        return Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
      };
      // the culling's kernels compile in the background
      let mainImage: number[] = [];
      for (let i = 0; i < 20; ++i) {
        mainImage = await draw(main, createFrame(), 'frame');
        if (new Set(mainImage).size > 4 && i > 3) break;
      }
      expect(new Set(mainImage).size).toBeGreaterThan(4);
      const otherImage = await draw(other, createFrame(), 'frame');
      expect(new Set(otherImage).size).toBeGreaterThan(4);
      const background = await (async () => {
        const s = Scene.create(graphics);
        renderer.update(main, s, createFrame());
        graphics.beginRenderPass(target);
        renderer.clear(true);
        graphics.flush();
        return Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
      })();
      const differing = (a: number[], b: number[]) => {
        let n = 0;
        for (let i = 0; i < a.length; i += 4) {
          if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 6) n += 1;
        }
        return n;
      };
      expect(differing(mainImage, otherImage)).toBeGreaterThan(20);

      await draw(main, createFrame(), 'frame');
      // a frame of its own and another camera: nothing without its culling (no CPU culling), its image with
      expect(await draw(other, createFrame(), 'none')).toEqual(background);
      expect(differing(await draw(other, createFrame(), 'view'), otherImage)).toBe(0);
      // the main view's next frame culled again
      expect(differing(await draw(main, createFrame(), 'frame'), mainImage)).toBe(0);
      // the same camera's frame of its own (the thickness pass) after the frame's: drawn, the frame's passes after too
      const frame = createFrame();
      await draw(main, frame, 'frame');
      expect(differing(await draw(main, createFrame(), 'view'), mainImage)).toBe(0);
      expect(differing(await draw(main, frame, 'none'), mainImage)).toBe(0);
      // a second camera of a frame (an XR eye): culled for it, the first's again by the next frame
      const eyes = createFrame();
      expect(differing(await draw(main, eyes, 'frame'), mainImage)).toBe(0);
      expect(differing(await draw(other, eyes, 'frame'), otherImage)).toBe(0);
      expect(differing(await draw(main, createFrame(), 'frame'), mainImage)).toBe(0);

      // moved by its transient transform (a drag's preview) to where only the moved is in view: culled as drawn, the
      // image of the lattice at rest seen from as far
      const offset = Vec3.create(600, 0, 0);
      const moved = cameraAt(Vec3.create(645, 45, 45), 70, Vec3.create(840, 80, 120));
      const transform = (object.values as { uObjectTransform: ValueCell<Mat4> }).uObjectTransform;
      ValueCell.update(transform, Mat4.fromTranslation(Mat4(), offset));
      let movedImage = await draw(moved, createFrame(), 'frame');
      movedImage = await draw(moved, createFrame(), 'frame');
      expect(differing(movedImage, otherImage)).toBeLessThan(8);
      ValueCell.update(transform, Mat4.identity());

      renderer.dispose();
      scene.clear();
      target.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("colors by member where the theme gives a target's color (by its entity): no color layer, the flat view's frames", async () => {
    const spec = DefaultPluginSpec();
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
      { width: Size, height: Size },
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
        camera: { mode: 'orthographic' } as never,
        cameraFog: { name: 'off', params: {} },
        multiSample: { mode: 'off' } as never,
        postprocessing: { occlusion: { name: 'off', params: {} }, outline: { name: 'off', params: {} } },
      });
      const graphics = canvas.graphics as WebGPUContext;
      // the lattice's targets of an entity each (every particle of a target its entity: `ParticleTargetInfo.entity`)
      const base = createList();
      const entities = Int32Array.from({ length: Count }, (_, k) => ParticleValues.get(base.targets, k) - 1);
      const withEntities: ParticleList = {
        ...base,
        entities,
        entityInfo: new Map([
          [0, { name: 'zero' }],
          [1, { name: 'one' }],
        ]),
        targetInfo: new Map([
          [1, { name: 'one', entity: 0 }],
          [2, { name: 'two', entity: 1 }],
        ]),
      };
      const list = await plugin.build().toRoot().apply(ListOf, { list: withEntities }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();
      let time = now();
      const frame = async () => {
        for (let i = 0; i < 12; ++i) {
          canvas.commit(true);
          canvas.camera.setState(canvas.camera.getFocus(Vec3.create(45, 45, 45), 70), 0);
          canvas.requestDraw();
          canvas.tick((time += 100) as now.Timestamp, { isSynchronous: true });
          await graphics.resources.pending.whenReady();
        }
        return Array.from((await graphics.readDrawingBuffer()).data);
      };
      const view = async (slots: boolean, colorTheme: string) => {
        const v = await plugin
          .build()
          .to(tiles.ref)
          .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids: '0-7' } }, frame: 0, slots })
          .commit();
        const repr = await plugin
          .build()
          .to(v.ref)
          .apply(ParticlesRepresentation3D, {
            type: { name: 'target', params: { lodPixelSize: 0 } },
            colorTheme: { name: colorTheme, params: {} },
          } as never)
          .commit();
        return { view: v, repr };
      };
      const differing = (a: number[], b: number[]) => {
        let n = 0;
        for (let i = 0; i < a.length; i += 4) {
          if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 6) n += 1;
        }
        return n;
      };
      const flat = await view(false, 'particle-entity');
      const flatImage = await frame();
      await plugin.build().delete(flat.view.ref).commit();

      const paged = await view(true, 'particle-entity');
      const object = paged.repr.data!.repr.renderObjects[0] as MergedGraphicsRenderObject;
      expect(object.merged.source?.paged).toBeDefined();
      const colorType = (object.values as unknown as { colorType: { ref: { value: string } } }).colorType.ref.value;
      // a member's color, not a particle's: no layer filled
      expect(colorType).toBe('object');
      const pagedImage = await frame();
      expect(new Set(flatImage).size).toBeGreaterThan(3);
      expect(differing(pagedImage, flatImage)).toBeLessThanOrEqual(8);

      // another theme a particle's color (by index): the object made again with a layer, the flat view's frames too
      await plugin.build().delete(paged.view.ref).commit();
      const flatIndex = await view(false, 'particle-index');
      const flatIndexImage = await frame();
      await plugin.build().delete(flatIndex.view.ref).commit();
      const pagedIndex = await view(true, 'particle-index');
      const indexObject = pagedIndex.repr.data!.repr.renderObjects[0] as MergedGraphicsRenderObject;
      const indexType = (indexObject.values as unknown as { colorType: { ref: { value: string } } }).colorType.ref
        .value;
      expect(indexType).toBe('instance');
      expect(differing(await frame(), flatIndexImage)).toBeLessThanOrEqual(8);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);
});
