/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { now } from '@molstar/core/util/now';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { OrderedSet } from '@molstar/core/data/int';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { Particle, type ParticleList, type ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { MarkerAction } from '@molstar/core/util/marker-action';
import type { MergedGraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import {
  ParticleTileSetFromList,
  ParticleTileSetFromPack,
  ParticleTilesView,
} from '@molstar/plugin/state/transforms/particles/tiles';
import { ParticleTileView } from '@molstar/model/model/particles/particle-tile-view';
import { writeParticlePack } from '@molstar/model/formats/particles/pack';
import { Asset } from '@molstar/core/util/assets';

import { WebGPU } from './webgpu-module.js';

// A slot view of a tile set (`ParticleTilesView`'s `slots`) drawn by the target representation's store-backed object
// (`paged.ts`): the image of the flat view of the same tiles, as tiles are let go and added again; a pick its
// particle, through a canvas multi-sampled over time.

const Size = 64;
const Count = 64;

/**
 * 64 particles on a 4³ lattice over 0-90 Å, targets 1 and 2 alternating along x (sites each, 6 Å); `shaped`: three
 * sites of 3 Å each, target 1's along x, target 2's along y (an instance drawn as the other target's shows), target 1
 * below z 45 (the tiles 0-3), target 2 above
 */
function createList(shaped = false): ParticleList {
  const coordinates = new Float32Array(Count * 3);
  const targets = new Int32Array(Count);
  for (let i = 0; i < Count; ++i) {
    const x = i % 4;
    coordinates.set([x * 30, (Math.floor(i / 4) % 4) * 30, Math.floor(i / 16) * 30], i * 3);
    targets[i] = shaped ? (Math.floor(i / 16) < 2 ? 1 : 2) : x % 2 ? 2 : 1;
  }
  const line = (axis: number): ParticleTarget => ({
    kind: 'sites',
    sites: {
      count: 3,
      positions: Float32Array.from({ length: 9 }, (_, k) => (k % 3 === axis ? (Math.floor(k / 3) - 1) * 6 : 0)),
      radii: new Float32Array(3).fill(3),
    },
  });
  const sites = (n: number): ParticleTarget =>
    shaped
      ? line(n === 2 ? 0 : 1)
      : {
          kind: 'sites',
          sites: {
            count: n,
            positions: Float32Array.from({ length: n * 3 }, (_, k) => (k % 3) * 2),
            radii: new Float32Array(n).fill(6),
          },
        };
  return {
    label: 'lattice',
    count: Count,
    targets,
    targetInfo: new Map([
      [1, { name: 'one' }],
      [2, { name: 'two' }],
    ]),
    targetMapping: new Map([
      [1, sites(2)],
      [2, sites(3)],
    ]),
    coordinates,
    getParticleLabel: (i) => `particle ${i}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

const ListOf = PluginStateTransform.BuiltIn({
  name: 'test-particle-slots-list-of',
  from: SO.Root,
  to: SO.Particle.List,
  params: { list: PD.Value<ParticleList>(undefined as any, { isHidden: true }) },
})({
  apply: ({ params }) => new SO.Particle.List(params.list, { label: 'lattice' }),
});

describe.skipIf(!WebGPU.module)('a slot view of particle tiles on WebGPU', () => {
  it('draws as the flat view of the same tiles, tiles let go and added again; picks a particle by its slot', async () => {
    const spec = DefaultPluginSpec();
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      plugin.canvas3d!.setProps({
        camera: { mode: 'orthographic' } as never,
        cameraFog: { name: 'off', params: {} },
        multiSample: { mode: 'off' } as never,
        postprocessing: { occlusion: { name: 'off', params: {} }, outline: { name: 'off', params: {} } },
        // impostors alone: the sphere raster (`auto`) takes the slot store, culled on the GPU, not the flat view's
        // objects, and decides a sphere's edge pixels unlike the impostors
        renderer: { ...plugin.canvas3d!.props.renderer, sphereRaster: 'off' },
      });
      const list = await plugin.build().toRoot().apply(ListOf, { list: createList() }).commit();
      const tiles = await plugin.build().to(list.ref).apply(ParticleTileSetFromList, { cells: 2 }).commit();

      const image = async () => {
        for (let i = 0; i < 4; ++i) {
          await plugin.getImageRaw({ width: Size, height: Size });
          // the cull's read-backs landed (a room overflowed asks for another cull, tracked once it is read: on a slow
          // device, a software Vulkan's, after the frame), then the work they asked for
          await plugin.canvas3d!.graphics.waitForGpuCommandsComplete();
          await new Promise((r) => setTimeout(r, 0));
          await plugin.canvas3d!.graphics.resources.pending.whenReady();
          const camera = plugin.canvas3d!.camera;
          camera.setState(camera.getFocus(Vec3.create(45, 45, 45), 70), 0);
        }
        return Array.from((await plugin.getImageRaw({ width: Size, height: Size })).data);
      };

      /** the image of a view of `ids`, its representation colored by particle index (by slot: the same colors) */
      const render = async (ids: string, slots: boolean) => {
        const view = await plugin
          .build()
          .to(tiles.ref)
          .apply(ParticleTilesView, { tiles: { name: 'ids', params: { ids } }, frame: 0, slots })
          .commit();
        const repr = await plugin
          .build()
          .to(view.ref)
          .apply(ParticlesRepresentation3D, {
            // no levels: an instance's level is its cell's, and the cells differ (the store's one grid, the flat view's a grid a member)
            type: { name: 'target', params: { lodPixelSize: 0 } },
            colorTheme: { name: 'particle-index', params: {} },
          } as never)
          .commit();
        return { view, repr, image: await image() };
      };
      const flatImage = async (ids: string) => {
        const { view, image } = await render(ids, false);
        await plugin.build().delete(view.ref).commit();
        return image;
      };

      const allFlat = await flatImage('0-7');
      const someFlat = await flatImage('0, 2, 5');
      expect(new Set(allFlat).size).toBeGreaterThan(3);
      expect(someFlat).not.toEqual(allFlat);

      const { view, repr, image: allSlots } = await render('0-7', true);
      const r = repr.data!.repr;
      // one object of the store's instances
      expect(r.renderObjects).toHaveLength(1);
      const object = r.renderObjects[0] as MergedGraphicsRenderObject;
      expect(object.merged.source).toBeDefined();
      expect(allSlots).toEqual(allFlat);

      const setIds = async (ids: string) => {
        await plugin
          .build()
          .to(view.ref)
          .update(ParticleTilesView, (old) => ({ ...old, tiles: { name: 'ids', params: { ids } } }))
          .commit();
        return image();
      };
      // tiles let go: their instances gaps, drawn by none; the object kept
      expect(await setIds('0, 2, 5')).toEqual(someFlat);
      expect(r.renderObjects[0]).toBe(object);
      // added again: in the gaps' place or after them
      expect(await setIds('0-7')).toEqual(allFlat);
      expect(r.renderObjects[0]).toBe(object);

      // a pick: an instance's particle its slot
      const slots = view.cell!.obj!.data;
      const loci = r.getLoci({ objectId: object.id, instanceId: 0, groupId: 0 }) as Particle.Loci;
      expect(Particle.isLoci(loci)).toBe(true);
      expect(loci.particles).toBe(slots);
      const slot = OrderedSet.getAt(loci.indices, 0);
      expect(slot).toBeGreaterThanOrEqual(0);
      expect(slot).toBeLessThan(slots.count);
      // marked by slot: the store's layer, the object's markers
      expect(r.mark(loci, MarkerAction.Highlight)).toBe(true);
      expect((object.values as any).markerAverage.ref.value).toBeGreaterThan(0);
      expect(r.mark(loci, MarkerAction.RemoveHighlight)).toBe(true);
      expect((object.values as any).markerAverage.ref.value).toBe(0);

      // the same slot list colored otherwise: what is marked kept, the object too (not laid out again)
      const recolor = (colorTheme: object) =>
        plugin
          .build()
          .to(repr.ref)
          .update(ParticlesRepresentation3D, (old) => ({ ...old, colorTheme }) as never)
          .commit();
      expect(r.mark(loci, MarkerAction.Highlight)).toBe(true);
      await recolor({ name: 'uniform', params: { value: 0xff0000 } });
      expect(r.renderObjects[0]).toBe(object);
      expect((object.values as any).markerAverage.ref.value).toBeGreaterThan(0);
      await recolor({ name: 'particle-index', params: {} });
      expect(r.mark(loci, MarkerAction.RemoveHighlight)).toBe(true);
      expect(await image()).toEqual(allFlat);

      // a layer by slot: the store's texels (its define on), the object kept; none again, as before
      // as the plugin applies a state: the canvas told
      const overpaint = (layers: { loci: Particle.Loci; color: number; clear: boolean }[]) => {
        r.setState({ overpaint: { kind: 'particle-loci', layers } as never });
        plugin.canvas3d!.update(r);
        plugin.canvas3d!.commit(true);
      };
      // every slot by its index (not every loci: the store indices of the slots)
      const all = Particle.Loci(slots, OrderedSet.ofBounds(0, slots.count));
      overpaint([{ loci: all, color: 0xff0000, clear: false }]);
      expect((object.values as any).dOverpaint.ref.value).toBe(true);
      const painted = await image();
      expect(painted).not.toEqual(allFlat);
      expect(r.renderObjects[0]).toBe(object);
      // the view's tiles let go and added again: the pages written their texels by the layers again, the others kept
      expect(await setIds('0, 2, 5')).not.toEqual(someFlat);
      expect(await setIds('0-7')).toEqual(painted);
      expect(r.renderObjects[0]).toBe(object);
      overpaint([]);
      expect((object.values as any).dOverpaint.ref.value).toBe(false);
      expect(await image()).toEqual(allFlat);

      // picks through the canvas, its frames multi-sampled over time: the canvas' frame the culled one. The first an
      // async one drawn incomplete (a draw skipped, as while a program compiles), then again as work settles: of the
      // canvas' frame then, not one of its own
      const canvas = plugin.canvas3d!;
      canvas.setProps({ multiSample: { mode: 'temporal' } as never });
      canvas.requestDraw();
      for (let i = 0; i < 3; ++i) canvas.tick(now(), { isSynchronous: true });
      // the pick's render reads the skipped draws before and after: one more after
      const stats = canvas.graphics.stats as { skippedDraws: number };
      let skipped = stats.skippedDraws;
      let reads = 0;
      Object.defineProperty(stats, 'skippedDraws', {
        configurable: true,
        get: () => (++reads === 2 ? ++skipped : skipped),
        set: (v: number) => {
          skipped = v;
        },
      });
      const pick = canvas.asyncIdentify(Vec2.create(Size / 2, Size / 2))!;
      Object.defineProperty(stats, 'skippedDraws', { configurable: true, writable: true, value: skipped });
      expect(pick.tryGet()).toBe('pending');
      // work settled since: drawn again
      await canvas.graphics.resources.pending.track(Promise.resolve(), 'test settle');
      let picked = pick.tryGet();
      for (let i = 0; i < 200 && picked === 'pending'; ++i) {
        canvas.tick(now(), { isSynchronous: true });
        await canvas.graphics.resources.pending.whenReady();
        await new Promise((r) => setTimeout(r, 5));
        picked = pick.tryGet();
      }
      expect(picked).not.toBe('pending');
      // the first pixel of a coarse grid that hits (the lattice's particles apart)
      let objectId: number | undefined;
      for (let y = 4; y < Size && objectId === undefined; y += 8) {
        for (let x = 4; x < Size && objectId === undefined; x += 8) {
          let hit = canvas.identify(Vec2.create(x, y));
          for (let i = 0; i < 50 && hit === 'pending'; ++i) {
            await new Promise((r) => setTimeout(r, 10));
            hit = canvas.identify(Vec2.create(x, y));
          }
          if (hit && hit !== 'pending') objectId = hit.id.objectId;
        }
      }
      expect(objectId).toBe(object.id);

      expect(plugin.canvas3d!.graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("draws a coarse node's spheres larger by its pose's scale, as sprites as rasterized", async () => {
    // the lattice a tile a particle (4³ tiles), its levels: the root a sample of it (one particle) drawn 4 times larger
    const bytes = writeParticlePack(createList(), { tiles: 4, levels: true });
    const file = Asset.File(new File([bytes], 'lattice.mspk'));
    /**
     * the pixels drawn (`sphereRaster` off: sprites; only: the raster) of the root node alone, or of the lattice's
     * tiles (`base`: 16 particles seen along z, the others behind them)
     */
    const drawn = async (raster: 'off' | 'only', base: boolean) => {
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
          renderer: { sphereRaster: raster, sphereRasterRadius: 64, backgroundColor: 0x000000 as never },
          camera: { mode: 'orthographic' } as never,
          cameraFog: { name: 'off', params: {} },
          multiSample: { mode: 'off' } as never,
          postprocessing: { occlusion: { name: 'off', params: {} }, outline: { name: 'off', params: {} } },
        });
        const tiles = await plugin
          .build()
          .toRoot()
          .apply(ParticleTileSetFromPack, { source: { name: 'file', params: file }, levels: true } as never)
          .commit();
        const levels = tiles.data!.set.levels!;
        // a camera given that sees the set: the tiles themselves at the default detail, the root alone where any
        // detail is too fine
        const camera = {
          snapshot: {
            mode: 'perspective' as const,
            fov: Math.PI / 4,
            position: Vec3.create(45, 45, 400),
            target: Vec3.create(45, 45, 0),
            up: Vec3.create(0, 1, 0),
            minNear: 1,
          },
          width: 512,
          height: 512,
        };
        const detailPixels = base ? 4 : 1e9;
        const view = await plugin
          .build()
          .to(tiles.ref)
          .apply(ParticleTilesView, {
            tiles: { name: 'camera', params: { minPixels: 0, occlusion: false, follow: false, detailPixels, camera } },
            frame: 0,
            slots: true,
          } as never)
          .commit();
        const shown = Array.from(ParticleTileView.of(view.data!)!.tiles);
        if (base) expect(shown.length).toBe(64);
        else expect(shown).toEqual([levels.getNode(levels.count, 0)]);
        await plugin
          .build()
          .to(view.ref)
          .apply(ParticlesRepresentation3D, {
            type: { name: 'target', params: { lodPixelSize: 0 } },
            colorTheme: { name: 'uniform', params: { value: 0x00ff00 } },
          } as never)
          .commit();
        for (let i = 0; i < 4; ++i) {
          await plugin.getImageRaw({ width: Size, height: Size });
          // the cull's read-backs landed (a room overflowed asks for another cull, tracked once it is read: on a slow
          // device, a software Vulkan's, after the frame), then the work they asked for
          await plugin.canvas3d!.graphics.waitForGpuCommandsComplete();
          await new Promise((r) => setTimeout(r, 0));
          await plugin.canvas3d!.graphics.resources.pending.whenReady();
          const camera = plugin.canvas3d!.camera;
          camera.setState(camera.getFocus(Vec3.create(45, 45, 45), 90), 0);
        }
        const { data } = await plugin.getImageRaw({ width: Size, height: Size });
        expect(plugin.canvas3d!.graphics.takeErrors?.() ?? []).toEqual([]);
        let n = 0;
        for (let i = 0; i < data.length; i += 4) if (data[i + 1] > 40) n += 1;
        return n;
      } finally {
        plugin.dispose();
      }
    };
    for (const raster of ['off', 'only'] as const) {
      const lattice = await drawn(raster, true);
      if (lattice === undefined) continue;
      const root = (await drawn(raster, false))!;
      // a particle of the lattice drawn as it is, the root's 4 times larger: about 16 times its area (its sites
      // overlap more, the discs of whole pixels)
      const one = lattice / 16;
      expect(root, raster).toBeGreaterThan(8 * one);
      expect(root, raster).toBeLessThan(24 * one);
    }
  }, 120_000);
});
