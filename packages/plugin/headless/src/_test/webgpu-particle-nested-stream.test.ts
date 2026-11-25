/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { now } from '@molstar/core/util/now';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { Particle, type ParticleList, type ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import { createParticleListFromPack, readParticlePack, writeParticlePack } from '@molstar/model/formats/particles/pack';
import {
  openParticleListFromPackStream,
  openParticlePackStream,
  ParticlePackSource,
} from '@molstar/model/formats/particles/pack-stream';
import { ParticleTargetRepresentationParams } from '@molstar/graphics/repr/particles/representation/target/common';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { WebGPU } from './webgpu-module.js';

// A nested pack read by range (`openParticleListFromPackStream`): its templates' particles read by page as the GPU's
// cull wants them (`ParticlePageStore.requestValues`), drawn as the pack read whole draws them.

const Size = 96;

/** one site, overlapping its neighbors (a dense packing: a cover) */
const site: ParticleTarget = {
  kind: 'sites',
  sites: { count: 1, positions: Float32Array.of(0, 0, 0), radii: Float32Array.of(2.5) },
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

/**
 * a template of 80 x 80 x 32 particles 4 Å apart (204,800: tile-ordered, four tiles an axis, 25 store pages, its levels
 * stored), two targets in alternating slabs, an entity a target; two copies of it side by side (`small`: 64 x 64 x 32,
 * 131,072 particles, two tiles an axis, 16 store pages; the copies 300 Å apart)
 */
function createNested(small = false) {
  const [nx, ny, nz] = small ? [64, 64, 32] : [80, 80, 32];
  const n = nx * ny * nz;
  const coordinates = new Float32Array(n * 3);
  for (let i = 0; i < n; ++i) {
    coordinates.set([(i % nx) * 4, (Math.floor(i / nx) % ny) * 4, Math.floor(i / (nx * ny)) * 4], i * 3);
  }
  const targets = Int32Array.from({ length: n }, (_, i) => Math.floor(i / (nx * ny * 5)) % 2);
  const template = createList({
    label: 'template',
    count: n,
    coordinates,
    targets,
    targetInfo: new Map([
      [0, { name: 'a' }],
      [1, { name: 'b' }],
    ]),
    targetMapping: new Map([
      [0, site],
      [1, site],
    ]),
    entities: Int32Array.from(targets),
    entityInfo: new Map([
      [0, { name: 'A' }],
      [1, { name: 'B' }],
    ]),
  });
  return createList({
    label: 'nested',
    count: 2,
    coordinates: Float32Array.of(0, 0, 0, small ? 300 : 400, 0, 0),
    targets: new Int32Array(2),
    targetInfo: new Map([[0, { name: 'copy' }]]),
    targetMapping: new Map<number, ParticleTarget>([
      [
        0,
        {
          kind: 'particles',
          particles: template,
          center: small ? Vec3.create(126, 126, 62) : Vec3.create(158, 158, 62),
        },
      ],
    ]),
  });
}

const ListOf = PluginStateTransform.BuiltIn({
  name: 'test-particle-nested-stream-list-of',
  from: SO.Root,
  to: SO.Particle.List,
  params: { list: PD.Value<ParticleList>(undefined as any, { isHidden: true }) },
})({
  apply: ({ params }) => new SO.Particle.List(params.list, { label: 'nested' }),
});

// Medium: a nested pack of 204,800 particles a template (25 store pages), drawn read whole and by range

describe.skipIf(!WebGPU.module)('a nested pack read by range on WebGPU', () => {
  it("draws as the pack read whole, near: its template's pages there read alone (16 pages, a small canvas)", async () => {
    const size = 32;
    const spec = DefaultPluginSpec();
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
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
        camera: { mode: 'orthographic', manualReset: true } as never,
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
      /** the canvas' frame at `focus` once its work settled (the cull's read-backs landed, the pages read): two alike */
      const frame = async (focus: Vec3, radius: number) => {
        const { pending } = graphics.resources;
        let image: number[] = [];
        for (let i = 0; i < 16; ++i) {
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
      // a corner of the first copy (its particles about its center): the pages there alone
      const near: [Vec3, number] = [Vec3.create(-100, -100, 62), 10];

      const bytes = writeParticlePack(createNested(true));
      const whole = await draw(createParticleListFromPack(readParticlePack(bytes)));
      const wholeNear = await frame(...near);
      expect(new Set(wholeNear).size).toBeGreaterThan(3);
      await plugin.build().delete(whole.ref).commit();

      // by range, the bytes read counted
      const source = ParticlePackSource.fromBytes(bytes);
      let read = 0;
      const counted: ParticlePackSource = {
        byteLength: source.byteLength,
        read: (offset, length, signal) => {
          read += length;
          return source.read(offset, length, signal);
        },
      };
      const list = await openParticleListFromPackStream(await openParticlePackStream(counted));
      const template = ParticleNesting.getTemplate(list, 0)!.particles;
      expect(Particle.getDeferred(template)).toBeDefined();
      const opened = read;
      const streamed = await draw(list);
      // the pages its view wants read (some, not all of the template's particles' bytes), drawn as whole
      expect(differing(await frame(...near), wholeNear)).toBeLessThanOrEqual(4);
      expect(read - opened).toBeGreaterThan(0);
      expect(read - opened).toBeLessThan(template.count * 8);
      await plugin.build().delete(streamed.ref).commit();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it(
    "draws as the pack read whole: its templates' pages read as the cull wants them, not all of them",
    { tags: ['medium'] },
    async () => {
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
          // the camera the test's alone (none fit to a scene added after the last let go: the pages a view wants)
          camera: { mode: 'orthographic', manualReset: true } as never,
          cameraFog: { name: 'off', params: {} },
          multiSample: { mode: 'off' } as never,
          postprocessing: { occlusion: { name: 'off', params: {} }, outline: { name: 'off', params: {} } },
          marking: { enabled: false } as never,
        });
        const graphics = canvas.graphics as WebGPUContext;
        let time = now();
        /** the canvas' frame at `focus` once its work settled (the pages read and written: pending work) */
        const frame = async (focus: Vec3, radius: number) => {
          for (let i = 0; i < 16; ++i) {
            canvas.commit(true);
            canvas.camera.setState(canvas.camera.getFocus(focus, radius), 0);
            canvas.requestDraw();
            canvas.tick((time += 100) as now.Timestamp, { isSynchronous: true });
            await graphics.resources.pending.whenReady();
          }
          return Array.from((await graphics.readDrawingBuffer()).data);
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
        // far: both copies whole; near: a corner of the first copy, the pages there alone
        const far: [Vec3, number] = [Vec3.create(358, 158, 62), 480];
        const near: [Vec3, number] = [Vec3.create(20, 20, 124), 14];

        const bytes = writeParticlePack(createNested());
        const whole = await draw(createParticleListFromPack(readParticlePack(bytes)));
        const wholeFar = await frame(...far);
        const wholeNear = await frame(...near);
        await plugin.build().delete(whole.ref).commit();

        // by range, the bytes read counted
        const source = ParticlePackSource.fromBytes(bytes);
        let read = 0;
        const counted: ParticlePackSource = {
          byteLength: source.byteLength,
          read: (offset, length, signal) => {
            read += length;
            return source.read(offset, length, signal);
          },
        };
        const list = await openParticleListFromPackStream(await openParticlePackStream(counted));
        const template = ParticleNesting.getTemplate(list, 0)!.particles;
        expect(Particle.getDeferred(template)).toBeDefined();
        const opened = read;
        const streamed = await draw(list);
        // near first: the pages its view wants read (some, not all of the template's particles' bytes), drawn as whole
        expect(differing(await frame(...near), wholeNear)).toBeLessThanOrEqual(8);
        const afterNear = read;
        const particleBytes = template.count * 8;
        expect(afterNear - opened).toBeGreaterThan(0);
        expect(afterNear - opened).toBeLessThan(particleBytes);
        // far: the rest read as wanted, drawn as whole
        expect(differing(await frame(...far), wholeFar)).toBeLessThanOrEqual(8);
        await plugin.build().delete(streamed.ref).commit();
        expect(graphics.takeErrors?.() ?? []).toEqual([]);
      } finally {
        plugin.dispose();
      }
    },
    120_000,
  );
});
