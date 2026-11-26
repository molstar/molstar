/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { OrderedSet } from '@molstar/core/data/int';
import { now } from '@molstar/core/util/now';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { CustomProperties } from '@molstar/model/model/custom-property';
import {
  Particle,
  type ParticleList,
  type ParticleSites,
  type ParticleTarget,
} from '@molstar/model/model/particles/particle-list';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import { ParticleTargetRepresentationParams } from '@molstar/graphics/repr/particles/representation/target/common';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { createWebGPUPlugin } from './webgpu-plugin.js';
// contexts of their own (`own`): a page store drawn before on a shared context changes their frames
import { WebGPU } from './webgpu-module.js';

// A list with templates' far copies (`ParticleNesting.getTemplateLevels`): a copy whose template's particles are a few
// pixels apart drawn at a coarse level (a sample, larger), chosen on the GPU by its distance; near, its template whole.

const Size = 96;

/** one site, overlapping its neighbors (a dense packing: a cover, not a lattice of dots) */
const site: ParticleTarget = {
  kind: 'sites',
  sites: { count: 1, positions: Float32Array.of(0, 0, 0), radii: Float32Array.of(3) },
};

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

/** an 8³ grid of 4 Å, two targets (alternating slabs along z), an entity a target: coarse levels of its own */
function createTemplate(of = site) {
  const n = 512;
  const coordinates = new Float32Array(n * 3);
  for (let i = 0; i < n; ++i)
    coordinates.set([(i % 8) * 4, (Math.floor(i / 8) % 8) * 4, Math.floor(i / 64) * 4], i * 3);
  const targets = Int32Array.from({ length: n }, (_, i) => Math.floor(i / 64) % 2);
  return createList({
    label: 'template',
    count: n,
    coordinates,
    targets,
    targetInfo: new Map([
      [0, { name: 'a' }],
      [1, { name: 'b' }],
    ]),
    targetMapping: new Map([
      [0, of],
      [1, of],
    ]),
    entities: Int32Array.from(targets),
    entityInfo: new Map([
      [0, { name: 'A' }],
      [1, { name: 'B' }],
    ]),
  });
}

/** four copies of the template on a line 100 Å apart */
function createNested(of = site) {
  return createList({
    label: 'nested',
    count: 4,
    coordinates: Float32Array.of(0, 0, 0, 100, 0, 0, 200, 0, 0, 300, 0, 0),
    targets: new Int32Array(4),
    targetInfo: new Map([[0, { name: 'copy' }]]),
    targetMapping: new Map<number, ParticleTarget>([
      [0, { kind: 'particles', particles: createTemplate(of), center: Vec3.create(14, 14, 14) }],
    ]),
  });
}

const ListOf = PluginStateTransform.BuiltIn({
  name: 'test-particle-nested-far-list-of',
  from: SO.Root,
  to: SO.Particle.List,
  params: { list: PD.Value<ParticleList>(undefined as any, { isHidden: true }) },
})({
  apply: ({ params }) => new SO.Particle.List(params.list, { label: 'nested' }),
});

// Medium: views of twelve frames each (the cull's read-backs landed, then the work they ask for), a list each
describe.skipIf(!WebGPU.module)("a list with templates' far copies on WebGPU", () => {
  it('draws a far copy at its coarse level (a sample, larger), a near one whole: four copies of 512 particles', async () => {
    const spec = DefaultPluginSpec();
    const plugin = await createWebGPUPlugin(
      { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
      { width: Size, height: Size },
      { own: true },
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
        marking: { enabled: false } as never,
        // fewer pipelines to compile: no sphere raster or pyramid
        renderer: { ...canvas.props.renderer, sphereRaster: 'off' },
        hiZ: { ...canvas.props.hiZ, enabled: false, gpu: false },
      });
      const graphics = canvas.graphics as WebGPUContext;
      let time = now();
      /** the canvas' frame at `focus` once its work settled (the cull's read-backs landed): two alike */
      const frame = async (focus: Vec3, radius: number) => {
        const { pending } = graphics.resources;
        let image: number[] = [];
        for (let i = 0; i < 12; ++i) {
          const tracked = pending.trackedCount;
          canvas.commit(true);
          canvas.camera.setState(canvas.camera.getFocus(focus, radius), 0);
          canvas.requestDraw();
          canvas.tick((time += 100) as now.Timestamp, { isSynchronous: true });
          await graphics.waitForGpuCommandsComplete();
          await new Promise((r) => setTimeout(r, 0));
          await pending.whenReady();
          const next = Array.from((await graphics.readDrawingBuffer()).data);
          const same = pending.trackedCount === tracked && next.every((v, k) => v === image[k]);
          image = next;
          if (same) break;
        }
        return image;
      };
      const draw = async (list: ParticleList) => {
        const data = await plugin.build().toRoot().apply(ListOf, { list }).commit();
        const structure = { ...PD.getDefaultValues(ParticleTargetRepresentationParams).structure };
        await plugin
          .build()
          .to(data.ref)
          .apply(ParticlesRepresentation3D, {
            type: { name: 'target', params: { lodPixelSize: 0, structure } },
            colorTheme: { name: 'particle-entity', params: {} },
          } as never)
          .commit();
        return data;
      };
      const differing = (a: number[], b: number[]) => {
        let n = 0;
        for (let i = 0; i < a.length; i += 4) {
          if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 6) n += 1;
        }
        return n;
      };
      /** the pixels of the background's color (a corner's) */
      const background = (image: number[]) => {
        let n = 0;
        for (let i = 0; i < image.length; i += 4) {
          if (image[i] === image[0] && image[i + 1] === image[1] && image[i + 2] === image[2]) n += 1;
        }
        return n;
      };
      const near: [Vec3, number] = [Vec3.create(0, 0, 0), 30];
      const far: [Vec3, number] = [Vec3.create(150, 0, 0), 180];

      const flat = await draw(ParticleNesting.flatten(createNested()));
      const flatNear = await frame(...near);
      const flatFar = await frame(...far);
      await plugin.build().delete(flat.ref).commit();

      const nested = await draw(createNested());
      // near: its template whole, as the flattened list; far: a sample drawn larger, about its cover
      expect(differing(await frame(...near), flatNear)).toBeLessThanOrEqual(8);
      const nestedFar = await frame(...far);
      expect(differing(nestedFar, flatFar)).toBeGreaterThan(8);
      const cover = Size * Size - background(flatFar);
      expect(Math.abs(Size * Size - background(nestedFar) - cover)).toBeLessThan(cover * 0.5);
      await plugin.build().delete(nested.ref).commit();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it(
    'draws a far copy at its coarse level (a sample, larger), a near one whole; picks and marks by path',
    { tags: ['medium'] },
    async () => {
      const spec = DefaultPluginSpec();
      const plugin = await createWebGPUPlugin(
        { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
        { width: Size, height: Size },
        { own: true },
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
          marking: { enabled: false } as never,
        });
        const graphics = canvas.graphics as WebGPUContext;
        let time = now();
        const frame = async (focus: Vec3, radius: number) => {
          for (let i = 0; i < 12; ++i) {
            canvas.commit(true);
            canvas.camera.setState(canvas.camera.getFocus(focus, radius), 0);
            canvas.requestDraw();
            canvas.tick((time += 100) as now.Timestamp, { isSynchronous: true });
            // the cull's read-backs landed (the leaf pages it wants, a room overflowed: tracked once read, on a slow
            // device, a software Vulkan's, after the frame), then the work they asked for
            await graphics.waitForGpuCommandsComplete();
            await new Promise((r) => setTimeout(r, 0));
            await graphics.resources.pending.whenReady();
          }
          return Array.from((await graphics.readDrawingBuffer()).data);
        };
        const draw = async (list: ParticleList) => {
          const data = await plugin.build().toRoot().apply(ListOf, { list }).commit();
          const structure = { ...PD.getDefaultValues(ParticleTargetRepresentationParams).structure };
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
        const differing = (a: number[], b: number[]) => {
          let n = 0;
          for (let i = 0; i < a.length; i += 4) {
            if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 6) n += 1;
          }
          return n;
        };
        /** the pixels of the background's color (a corner's) */
        const background = (image: number[]) => {
          let n = 0;
          for (let i = 0; i < image.length; i += 4) {
            if (image[i] === image[0] && image[i + 1] === image[1] && image[i + 2] === image[2]) n += 1;
          }
          return n;
        };
        /** the mean color of the pixels not the background's */
        const meanColor = (image: number[]) => {
          const sum = [0, 0, 0];
          let n = 0;
          for (let i = 0; i < image.length; i += 4) {
            if (image[i] === image[0] && image[i + 1] === image[1] && image[i + 2] === image[2]) continue;
            for (let k = 0; k < 3; ++k) sum[k] += image[i + k];
            n += 1;
          }
          return sum.map((v) => v / Math.max(1, n));
        };
        // near: a copy's particles 4 Å apart about 6 px; far: all four, about a pixel apart
        const near: [Vec3, number] = [Vec3.create(0, 0, 0), 30];
        const far: [Vec3, number] = [Vec3.create(150, 0, 0), 180];

        const flat = await draw(ParticleNesting.flatten(createNested()));
        const flatNear = await frame(...near);
        const flatFar = await frame(...far);
        await plugin.build().delete(flat.data.ref).commit();

        const nested = createNested();
        const template = (nested.targetMapping!.get(0) as { particles: ParticleList }).particles;
        const levels = ParticleNesting.getTemplateLevels(nested, template);
        expect(levels.length).toBeGreaterThan(0);
        const { data, repr } = await draw(nested);
        const r = repr.data!.repr;
        // near: its template whole, as the flattened list
        expect(differing(await frame(...near), flatNear)).toBeLessThanOrEqual(8);
        // far: a sample drawn larger, not the flattened list's image, but about its color and its cover
        const nestedFar = await frame(...far);
        expect(differing(nestedFar, flatFar)).toBeGreaterThan(8);
        const flatMean = meanColor(flatFar);
        const nestedMean = meanColor(nestedFar);
        for (let k = 0; k < 3; ++k) expect(Math.abs(nestedMean[k] - flatMean[k])).toBeLessThan(24);
        const cover = Size * Size - background(flatFar);
        expect(Math.abs(Size * Size - background(nestedFar) - cover)).toBeLessThan(cover * 0.5);

        // a far copy's pick: its path's particle one of a level's sample (its template particle's index), near the hit
        const sampled = new Set<number>();
        for (const l of levels) for (const i of l.indices) sampled.add(i);
        const reference = ParticleNesting.flatten(createNested());
        let picked: { top: number; local: number } | undefined;
        for (let y = 2; y < Size && !picked; y += 3) {
          for (let x = 2; x < Size && !picked; x += 3) {
            let hit = canvas.identify(Vec2.create(x, y));
            for (let i = 0; i < 50 && hit === 'pending'; ++i) {
              await new Promise((resolve) => setTimeout(resolve, 10));
              hit = canvas.identify(Vec2.create(x, y));
            }
            if (!hit || hit === 'pending') continue;
            const loci = r.getLoci(hit.id) as Particle.Loci;
            if (!Particle.isLoci(loci) || !loci.partial) continue;
            const top = loci.partial.tops[0];
            const local = OrderedSet.start(loci.partial.locals[0]);
            expect(sampled.has(local)).toBe(true);
            const index = ParticleNesting.getExpandedIndex(nested, top, local);
            const at = Vec3.fromArray(Vec3(), reference.coordinates as Float32Array, index * 3);
            // the hit on its sphere, drawn larger (by the level's scale, a few times the 1.5 Å)
            expect(Vec3.distance(at, hit.position)).toBeLessThan(20);
            picked = { top, local };
          }
        }
        expect(picked).toBeDefined();

        // marked by its path: its copy's sampled particle tinted (the copy's markers at its level)
        const unmarked = await frame(...far);
        canvas.mark(
          { loci: ParticleNesting.ofInstance(nested, picked!.top, picked!.local), repr: r },
          MarkerAction.Select,
        );
        expect(differing(await frame(...far), unmarked)).toBeGreaterThan(0);
        await plugin.build().delete(data.ref).commit();

        expect(graphics.takeErrors?.() ?? []).toEqual([]);
      } finally {
        plugin.dispose();
      }
    },
    120_000,
  );
  it(
    "blends a far copy's sites' element color by their own size on screen, not its coarse level's scale",
    { tags: ['medium'] },
    async () => {
      const spec = DefaultPluginSpec();
      const plugin = await createWebGPUPlugin(
        { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
        { width: Size, height: Size },
        { own: true },
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
          marking: { enabled: false } as never,
        });
        const graphics = canvas.graphics as WebGPUContext;
        let time = now();
        // the sites oxygens: their second color red
        const oxygen: ParticleTarget = {
          kind: 'sites',
          sites: { ...(site as { sites: ParticleSites }).sites, elements: Uint8Array.of(8) },
        };
        /** the red pixels (an oxygen's color blended in) of the four copies far, the second color from `radius` px */
        const red = async (radius: [number, number]) => {
          const data = await plugin
            .build()
            .toRoot()
            .apply(ListOf, { list: createNested(oxygen) })
            .commit();
          const structure = {
            ...PD.getDefaultValues(ParticleTargetRepresentationParams).structure,
            siteColor: { name: 'element', params: { carbonColor: 'particle', radius } },
          };
          await plugin
            .build()
            .to(data.ref)
            .apply(ParticlesRepresentation3D, {
              type: {
                name: 'target',
                params: { lodPixelSize: 0, structure },
              },
              colorTheme: { name: 'particle-entity', params: {} },
            } as never)
            .commit();
          for (let i = 0; i < 12; ++i) {
            canvas.commit(true);
            canvas.camera.setState(canvas.camera.getFocus(Vec3.create(150, 0, 0), 180), 0);
            canvas.requestDraw();
            canvas.tick((time += 100) as now.Timestamp, { isSynchronous: true });
            // the cull's read-backs landed (the leaf pages it wants, a room overflowed: tracked once read, on a slow
            // device, a software Vulkan's, after the frame), then the work they asked for
            await graphics.waitForGpuCommandsComplete();
            await new Promise((r) => setTimeout(r, 0));
            await graphics.resources.pending.whenReady();
          }
          const image = (await graphics.readDrawingBuffer()).data;
          await plugin.build().delete(data.ref).commit();
          let n = 0;
          for (let i = 0; i < image.length; i += 4)
            if (image[i] > 150 && image[i + 1] < 60 && image[i + 2] < 60) n += 1;
          return n;
        };
        // far, the particles about a pixel apart: a site of 3 Å about 0.75 px in radius, its coarse level's drawn larger
        // (2 or 4 times); blended from 0.2 px red, from 1.2 px not (its own size below it, whatever its level's scale)
        expect(await red([0.2, 0.4])).toBeGreaterThan(100);
        expect(await red([1.2, 1.4])).toBe(0);
        expect(graphics.takeErrors?.() ?? []).toEqual([]);
      } finally {
        plugin.dispose();
      }
    },
    120_000,
  );
});
