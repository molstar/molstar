/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { Clip } from '@molstar/core/util/clip';
import { ValueCell } from '@molstar/core/util/value-cell';
import { ClipStorageVec4s, packStorageClipObjects } from '@molstar/graphics/gl/clip-storage';
import { CustomProperties } from '@molstar/model/model/custom-property';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';
import { createPagedCuller, PagedPoseWords } from '@molstar/graphics/gl/compute/paged-cull';
import type { OcclusionPyramid } from '@molstar/graphics/gl/compute/cell-cull';
import type { StorageBuffer } from '@molstar/graphics/gl/compute';
import {
  PageRecordWords,
  ParticlePageStore,
  StorePageBits,
  StorePageSize,
  type StorePlacement,
} from '@molstar/graphics/repr/particles/representation/target/page-store';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { createDataArray } from '@molstar/graphics/gl/renderable/util';
import { createMergedValues, type MergeableValues } from '@molstar/graphics/gl/renderable/merged';
import { createPagedMergedCuller } from '@molstar/graphics/gl/renderable/merged-paged-culler';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

// A page store's particles culled and their poses composed on the GPU (`paged-cull.ts`): the regions' counts and the
// instances in them (their pose slots, the poses composed as the store's pages are written) as the CPU classifies the
// store's composed poses (`ParticlePageStore.compose`).

/** the regions' pose slots and the poses they read: per region its slots (from its first, its count's) */
async function readRegions(culler: ReturnType<typeof createPagedCuller>, counts: StorageBuffer) {
  const regionCount = culler.regionCount;
  const c = new Uint32Array(await counts.readAsync(0, regionCount * 4));
  const bases = new Uint32Array(await culler.bases.readAsync(0, regionCount * 4));
  const indices = new Uint32Array(await culler.indices.readAsync(0, culler.indexRoom * 4));
  const poses = new Uint32Array(await culler.poses.readAsync(0, culler.poses.byteLength));
  const regions: number[][] = [];
  for (let region = 0; region < regionCount; ++region) {
    const slots: number[] = [];
    for (let k = bases[region]; k < bases[region] + c[region]; ++k) slots.push(indices[k]);
    regions.push(slots);
  }
  return { regions, poses, total: c.reduce((a, b) => a + b, 0) };
}

/** `n` particles in a 100 Å box (a deterministic scatter), targets by `targetOf`, rotated, radii a particle */
function createList(n: number, targetOf: (i: number) => number, grid: boolean): ParticleList {
  const coordinates = new Float32Array(n * 3);
  const rotations = new Float32Array(n * 4);
  const q = Quat();
  for (let i = 0; i < n; ++i) {
    coordinates.set([(i * 37) % 100, (i * 61) % 97, (i * 13) % 89], i * 3);
    Quat.setAxisAngle(q, Vec3.normalize(Vec3(), Vec3.create(1, i % 3, 2)), i * 0.3);
    rotations.set([q[0], q[1], q[2], q[3]], i * 4);
  }
  return {
    count: n,
    targets: Int32Array.from({ length: n }, (_, i) => targetOf(i)),
    targetInfo: new Map(),
    coordinates: grid ? GridPositions.fromFloat(coordinates) : coordinates,
    rotations: grid ? Int16Array.from(rotations, (x) => Math.round(x * 32767)) : rotations,
    radii: Float32Array.from({ length: n }, (_, i) => 1 + (i % 5) * 0.5),
    getParticleLabel: (i) => `${i}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

const centers = new Map<number, Vec3>([
  [1, Vec3.create(1, 0, 0)],
  [2, Vec3.create(0, 2, 0)],
  [3, Vec3.create(0, 0, -1)],
  [4, Vec3.create(0.5, 0.5, 0.5)],
]);
const placement: StorePlacement = {
  centerOf: (t) => centers.get(t)!,
  scaleByRadius: (t) => t === 2,
  accepts: (t) => t !== 4,
};

/** a camera at `eye` looking at `target`, a narrow frustum: its planes, its camera plane, two levels by distance */
function createView(eye: Vec3, target: Vec3) {
  const view = Mat4.lookAt(Mat4(), eye, target, Vec3.create(0, 1, 0));
  const t = 0.08;
  const proj = Mat4.perspective(Mat4(), -t, t, t, -t, 1, 400);
  const frustum = Frustum3D.fromProjectionMatrix(Frustum3D(), Mat4.mul(Mat4(), proj, view));
  const dir = Vec3.normalize(Vec3(), Vec3.sub(Vec3(), target, eye));
  const cameraPlane = Plane3D.create(dir, -Vec3.dot(dir, eye));
  const split = Vec3.distance(eye, target);
  const levels = Mat4.zero();
  levels[0] = 0;
  levels[1] = split;
  levels[4] = split;
  levels[5] = 1e9;
  return { frustum, cameraPlane, levels, levelCount: 2 };
}

describe.skipIf(!WebGPU.module)('a page store culled on the GPU', () => {
  it('counts the visible instances by member and level and writes their poses into their regions', async () => {
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
      const store = new ParticlePageStore({ pageCapacity: 2 });
      // a grid list of two pages (its runs across their end), a float one converted, alternating targets (a member a
      // particle), a list let go
      store.add(
        1,
        createList(70_000, (i) => (i < 30_000 ? 1 : i < 50_000 ? 2 : i < 52_000 ? 4 : 3), true),
        0,
        placement,
      );
      store.add(
        2,
        createList(5_000, (i) => (i < 2_500 ? 3 : 1), false),
        70_000,
        placement,
      );
      store.add(
        3,
        createList(64, (i) => [1, 2, 3][i % 3], true),
        75_000,
        placement,
      );
      store.add(
        4,
        createList(500, () => 1, true),
        76_000,
        placement,
      );
      store.remove(4);
      const reaches = [2, 3, 1.5, 4];
      const culler = createPagedCuller(compute, graphics.resources);
      culler.sync(store, store.takeChanges(), reaches);
      const view = createView(Vec3.create(30, 40, 250), Vec3.create(50, 50, 40));

      /** the store culled with `view`: the GPU's regions as the CPU classifies the composed poses; the slots seen */
      const check = async () => {
        culler.cull(view);
        const { regions, poses, total } = await readRegions(culler, culler.counts);
        expect(total).toBeGreaterThan(100);
        const poseFloats = new Float32Array(poses.buffer);

        // the CPU's: each instance's particle from its composed pose (p = t + s R c), classified as the kernel does;
        // within a hair of a plane or a level's bound undecided
        const certain = new Map<number, number>();
        const undecided = new Set<number>();
        const cpuPoses = new Map<number, Float32Array>();
        const p = Vec3();
        const c = Vec3();
        for (const block of store.blockList) {
          const composed = store.compose(block);
          for (let i = 0; i < composed.slots.length; ++i) {
            const m = composed.members[i];
            const pose = composed.poses.subarray(i * 8, i * 8 + 8);
            cpuPoses.set(composed.slots[i], pose);
            const q = Quat.create(pose[4], pose[5], pose[6], pose[7]);
            Vec3.set(c, store.pool.centers[m * 3], store.pool.centers[m * 3 + 1], store.pool.centers[m * 3 + 2]);
            Vec3.transformQuat(p, c, q);
            Vec3.scaleAndAdd(p, Vec3.create(pose[0], pose[1], pose[2]), p, pose[3]);
            const r = reaches[m] * pose[3];
            let margin = Infinity;
            let inside = true;
            for (let k = 0; k < 6; ++k) {
              const plane = view.frustum[k as 0];
              const dist = Vec3.dot(plane.normal, p) + plane.constant + r;
              margin = Math.min(margin, Math.abs(dist));
              if (dist < 0) inside = false;
            }
            const d = Vec3.dot(view.cameraPlane.normal, p) + view.cameraPlane.constant;
            margin = Math.min(margin, Math.abs(d - view.levels[1]));
            const level = d < view.levels[1] ? 0 : 1;
            const slot = composed.slots[i];
            if (margin < 1e-2) undecided.add(slot);
            else if (inside) certain.set(slot, m * 2 + level);
          }
        }

        // the GPU's: each region's instances (their poses by their pose slots), by store slot
        const seen = new Map<number, number>();
        for (let region = 0; region < regions.length; ++region) {
          for (const k of regions[region]) {
            const o = k * PagedPoseWords;
            // its store index, its slot by its page's record
            const index = poses[o + 6];
            const page = index >> StorePageBits;
            const slot = store.records[page * PageRecordWords + 7] + (index & (StorePageSize - 1));
            expect(seen.has(slot)).toBe(false);
            seen.set(slot, region);
            expect(poses[o + 7]).toBe(Math.floor(region / 2));
            const cpu = cpuPoses.get(slot)!;
            expect(cpu).toBeDefined();
            // t and the scale as f32, the quaternion as snorm16
            for (let a = 0; a < 4; ++a) expect(Math.abs(poseFloats[o + a] - cpu[a])).toBeLessThan(1e-3);
            const q = [poses[o + 4] & 0xffff, poses[o + 4] >>> 16, poses[o + 5] & 0xffff, poses[o + 5] >>> 16].map(
              (v) => Math.max(-1, ((v << 16) >> 16) / 32767),
            );
            for (let a = 0; a < 4; ++a) expect(Math.abs(q[a] - cpu[4 + a])).toBeLessThan(1e-4);
          }
        }
        for (const [slot, region] of certain) expect(seen.get(slot), `slot ${slot}`).toBe(region);
        for (const slot of seen.keys()) expect(certain.has(slot) || undecided.has(slot), `slot ${slot}`).toBe(true);
        return seen;
      };

      let seen = await check();
      // a let go block's and a left out target's none
      expect([...seen.keys()].some((s) => s >= 76_000)).toBe(false);
      expect([...seen.keys()].some((s) => s >= 50_000 && s < 52_000)).toBe(false);
      expect([...seen.keys()].some((s) => s >= 70_000 && s < 75_000)).toBe(true);

      // another sync: a list of five pages (the values' buffers grown, what was there copied before the new pages are
      // written), the float list let go (its pages free, reused)
      const pages = store.pageCapacity;
      store.add(
        5,
        createList(300_000, (i) => (i < 150_000 ? 2 : 3), true),
        100_000,
        placement,
      );
      store.remove(2);
      expect(store.pageCapacity).toBeGreaterThan(pages);
      culler.sync(store, store.takeChanges(), reaches);
      seen = await check();
      expect([...seen.keys()].some((s) => s >= 70_000 && s < 75_000)).toBe(false);
      expect([...seen.keys()].some((s) => s >= 100_000)).toBe(true);
      // the first list's instances as before the growth
      expect([...seen.keys()].some((s) => s < 50_000)).toBe(true);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
      culler.destroy();
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it('culls by a pyramid: the occluded marked, the second phase finding those it shows', async () => {
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
      const store = new ParticlePageStore();
      store.add(
        1,
        createList(20_000, (i) => (i % 2 ? 1 : 3), true),
        0,
        placement,
      );
      const culler = createPagedCuller(compute, graphics.resources);
      culler.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      const eye = Vec3.create(30, 40, 250);
      const target = Vec3.create(50, 50, 40);
      const view = createView(eye, target);
      const total = async (counts: StorageBuffer) => {
        let n = 0;
        for (const c of new Uint32Array(await counts.readAsync(0, culler.regionCount * 4))) n += c;
        return n;
      };
      /** a pyramid of one texel of `depth` (reversed: 1 at the near plane) over a viewport of a pixel */
      const pyramidOf = (depth: number): OcclusionPyramid => {
        const texture = graphics.resources.texture('image-float32', 'r', 'float', 'nearest');
        texture.load({ array: Float32Array.of(depth), width: 1, height: 1 });
        const viewMatrix = Mat4.lookAt(Mat4(), eye, target, Vec3.create(0, 1, 0));
        const t = 0.08;
        return {
          texture,
          levelCount: 1,
          minLevel: 0,
          width: 1,
          height: 1,
          view: viewMatrix,
          projection: Mat4.perspective(Mat4(), -t, t, t, -t, 1, 400),
          near: 1,
          far: 400,
          version: 1,
        } as OcclusionPyramid;
      };
      culler.cull(view);
      const visible = await total(culler.counts);
      expect(visible).toBeGreaterThan(100);
      // all at the near plane: a sphere in the view behind it (one reaching past the view's edge not tested: kept)
      culler.cull(view, pyramidOf(1));
      const unoccluded = await total(culler.counts);
      expect(unoccluded).toBeLessThan(visible - 500);
      // the second phase against a pyramid at the far plane: those the first occluded, all in sight
      expect(culler.cullSecond(pyramidOf(0))).toBe(true);
      expect(await total(culler.secondCounts)).toBe(visible - unoccluded);
      // none again without a cull of a pyramid
      expect(culler.cullSecond(pyramidOf(0))).toBe(false);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
      culler.destroy();
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it('culls what clip objects clip whole and writes the ghost of what they reach into its regions, by pixel and by instance', async () => {
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
      const store = new ParticlePageStore();
      store.add(
        1,
        createList(20_000, (i) => (i % 2 ? 1 : 3), true),
        0,
        placement,
      );
      const reaches = [2, 3, 1.5, 4];
      const culler = createPagedCuller(compute, graphics.resources);
      culler.sync(store, store.takeChanges(), reaches);
      const view = { ...createView(Vec3.create(30, 40, 250), Vec3.create(50, 50, 40)), ghost: true };
      // a tilted plane through the view's target
      const clip = Clip.getClip({
        variant: 'pixel',
        objects: [
          {
            type: 'plane',
            invert: false,
            position: Vec3.create(50, 50, 40),
            rotation: { axis: Vec3.create(0, 1, 0.4), angle: 60 },
            scale: Vec3.create(1, 1, 1),
            transform: Mat4.identity(),
          },
        ],
      });
      const data = new Float32Array(ClipStorageVec4s * 4);
      packStorageClipObjects(
        {
          dClipObjects: ValueCell.create(true),
          dClipVariant: ValueCell.create('pixel'),
          dClipping: ValueCell.create(false),
          uClipObjectCount: ValueCell.create(1),
          uClipObjectFrame: ValueCell.create(clip.objects.frame),
          uClipObjectShape: ValueCell.create(clip.objects.shape),
        },
        data,
      );

      /** the slots in the levels' regions and in the ghost's, from the poses' store indices */
      const cullSlots = async () => {
        culler.cull(view);
        const ghostBase = culler.memberCount * view.levelCount;
        expect(culler.regionCount).toBe(ghostBase + culler.memberCount);
        const { regions, poses } = await readRegions(culler, culler.counts);
        const levels = new Set<number>();
        const ghost = new Set<number>();
        for (let region = 0; region < regions.length; ++region) {
          const into = region < ghostBase ? levels : ghost;
          for (const k of regions[region]) {
            const index = poses[k * PagedPoseWords + 6];
            const page = index >> StorePageBits;
            const slot = store.records[page * PageRecordWords + 7] + (index & (StorePageSize - 1));
            expect(into.has(slot)).toBe(false);
            into.add(slot);
          }
        }
        return { levels, ghost };
      };

      /**
       * the CPU's, as the kernel classifies (`Clip.classifySphere`, a center clipped where it is `Some` at no radius):
       * per slot whether it is drawn at a level and whether the ghost's; within a hair of a plane or the clip plane
       * undecided
       */
      const expected = (instance: boolean, clipped: boolean) => {
        const levels = new Set<number>();
        const ghost = new Set<number>();
        const undecided = new Set<number>();
        const p = Vec3();
        const c = Vec3();
        const sphere = Sphere3D();
        const classOf = (r: number) => {
          Sphere3D.set(sphere, p, r);
          return Clip.classifySphere(clip.objects, sphere);
        };
        for (const block of store.blockList) {
          const composed = store.compose(block);
          for (let i = 0; i < composed.slots.length; ++i) {
            const m = composed.members[i];
            const pose = composed.poses.subarray(i * 8, i * 8 + 8);
            const q = Quat.create(pose[4], pose[5], pose[6], pose[7]);
            Vec3.set(c, store.pool.centers[m * 3], store.pool.centers[m * 3 + 1], store.pool.centers[m * 3 + 2]);
            Vec3.transformQuat(p, c, q);
            Vec3.scaleAndAdd(p, Vec3.create(pose[0], pose[1], pose[2]), p, pose[3]);
            const r = reaches[m] * pose[3];
            const slot = composed.slots[i];
            let inside = true;
            let margin = Infinity;
            for (let k = 0; k < 6; ++k) {
              const plane = view.frustum[k as 0];
              const dist = Vec3.dot(plane.normal, p) + plane.constant + r;
              margin = Math.min(margin, Math.abs(dist));
              if (dist < 0) inside = false;
            }
            const near =
              margin < 1e-2 ||
              (clipped && (classOf(r + 1e-2) !== classOf(Math.max(0, r - 1e-2)) || classOf(1e-2) === 1));
            if (near) {
              undecided.add(slot);
              continue;
            }
            if (!inside) continue;
            const sphereClass = clipped ? classOf(r) : 0;
            const centerClipped = instance && clipped && classOf(0) !== 0;
            if (sphereClass !== 2 && !centerClipped) levels.add(slot);
            if ((!clipped || sphereClass !== 0) && (!instance || centerClipped)) ghost.add(slot);
          }
        }
        return { levels, ghost, undecided };
      };
      const check = async (instance: boolean, clipped: boolean) => {
        const gpu = await cullSlots();
        const cpu = expected(instance, clipped);
        for (const [name, g, c] of [
          ['levels', gpu.levels, cpu.levels],
          ['ghost', gpu.ghost, cpu.ghost],
        ] as const) {
          for (const slot of c) expect(g.has(slot), `${name} ${slot}`).toBe(true);
          for (const slot of g) expect(c.has(slot) || cpu.undecided.has(slot), `${name} ${slot}`).toBe(true);
        }
        return gpu;
      };

      // no clip objects: the ghost's every instance in view
      const none = await check(false, false);
      expect(none.ghost.size).toBe(none.levels.size);
      culler.setClip({ data, count: 1, instance: false });
      const pixel = await check(false, true);
      expect(pixel.levels.size).toBeLessThan(none.levels.size - 100);
      expect(pixel.ghost.size).toBeGreaterThan(100);
      expect(pixel.ghost.size).toBeLessThan(none.ghost.size - 100);
      culler.setClip({ data, count: 1, instance: true });
      const instance = await check(true, true);
      expect(instance.levels.size).toBeLessThan(pixel.levels.size);
      expect(instance.ghost.size).toBeLessThan(pixel.ghost.size);
      expect(instance.ghost.size).toBeGreaterThan(100);
      // kept through the tables written whole (their sections grown)
      store.add(
        2,
        createList(200_000, () => 1, true),
        20_000,
        placement,
      );
      culler.sync(store, store.takeChanges(), reaches);
      await check(true, true);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
      culler.destroy();
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("keeps its poses and indices within the device's buffer limit: the particles past the poses' room not drawn", async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      // a limit the store's values fit (its 3 pages' positions), not its 20k instances' poses (32 B each); not a multiple
      // of 16, as the devices' (2^31 - 4)
      const limit = 3 * 8192 * 8 + 4;
      const compute = Object.create(graphics.compute!, { maxStorageBufferSize: { value: limit } });
      const store = new ParticlePageStore();
      store.add(
        1,
        createList(20_000, (i) => (i % 2 ? 1 : 3), true),
        0,
        placement,
      );
      const culler = createPagedCuller(compute, graphics.resources);
      culler.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      // a view of many of them (the far plane at 400)
      const view = createView(Vec3.create(50, 50, 395), Vec3.create(50, 50, 0));
      culler.cull(view);
      expect(culler.poses.byteLength).toBeLessThanOrEqual(limit);
      expect(culler.indices.byteLength).toBeLessThanOrEqual(limit);
      // the pages with poses drawn (some), those past the poses' room not: fewer than the instances
      const room = Math.floor(limit / (PagedPoseWords * 4));
      const { regions, poses, total } = await readRegions(culler, culler.counts);
      expect(total).toBeGreaterThan(100);
      expect(total).toBeLessThanOrEqual(room);
      // every slot written a pose of the store's: a store index of its pages
      const pages = new Set(store.getBlock(1)!.pages);
      for (const slots of regions) {
        for (const k of slots) {
          expect(k).toBeLessThan(room);
          expect(pages.has(poses[k * PagedPoseWords + 6] >> StorePageBits), `pose ${k}`).toBe(true);
        }
      }
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
      culler.destroy();
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("grows its poses' buffer to the device's limit, not past it (a growth's half again capped, its size rounded)", async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      // not a multiple of 16, as the devices' (2^31 - 4): a growth capped at it, rounded up to 16, past it
      const limit = 4 * 8192 * 8 + 4;
      const compute = Object.create(graphics.compute!, { maxStorageBufferSize: { value: limit } });
      const poseBytes = PagedPoseWords * 4;
      // a pose slot an instance: about 80% of the limit, then 85%
      const first = Math.floor((0.8 * limit) / poseBytes);
      const second = Math.floor((0.85 * limit) / poseBytes) - first;
      const store = new ParticlePageStore();
      store.add(
        1,
        createList(first, (i) => (i % 2 ? 1 : 3), true),
        0,
        placement,
      );
      const culler = createPagedCuller(compute, graphics.resources);
      culler.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      const view = createView(Vec3.create(50, 50, 395), Vec3.create(50, 50, 0));
      culler.cull(view);
      const before = culler.poses.byteLength;
      store.add(
        2,
        createList(second, (i) => (i % 2 ? 1 : 3), true),
        10_000,
        placement,
      );
      culler.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      culler.cull(view);
      // grown, to the limit at most
      expect(culler.poses.byteLength).toBeGreaterThan(before);
      expect(culler.poses.byteLength).toBeLessThanOrEqual(limit);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
      culler.destroy();
    } finally {
      plugin.dispose();
    }
  }, 120_000);
  it('keeps its data version where a sync changed nothing it reads (the markers alone), else advances it', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const store = new ParticlePageStore();
      store.add(
        1,
        createList(100, (i) => (i % 2 ? 1 : 3), true),
        0,
        placement,
      );
      // the store's members (targets 1 and 3) and the gaps' member, as a target representation's pool
      const members: MergeableValues[] = [];
      for (let i = 0; i < 3; ++i) {
        const values = Spheres.Utils.createValuesSimple(
          Spheres.createEmpty(),
          PD.getDefaultValues(Spheres.Params),
          ColorNames.orange,
          1,
        ) as unknown as MergeableValues;
        ValueCell.update(values.colorType, 'instance');
        ValueCell.update(values.markerType, 'instance');
        members.push(values);
      }
      const markers = ValueCell.create(createDataArray(store.pageCapacity * StorePageSize, 1, Uint8Array));
      let pending = store.takeChanges();
      const merged = createMergedValues('spheres', members, true, false, {
        layers: new Map([
          ['tColor', ValueCell.create(createDataArray(store.pageCapacity * StorePageSize, 4, Uint8Array))],
          ['tMarker', markers],
        ]),
        paged: {
          source: store,
          takeChanges: () => {
            const taken = pending;
            pending = { written: [], released: [], runs: [], extras: [], members: new Set(), freed: [], grown: false };
            return taken;
          },
          reaches: [2, 3, 1.5],
        },
      });
      const culler = createPagedMergedCuller(graphics.compute!, graphics.resources, merged, members, false);
      culler.sync();
      const synced = culler.dataVersion;
      // a marker set (a hover): the frame's cull stands, the marking pass drawing its lists
      markers.ref.value.array[3] = 1;
      ValueCell.update(markers, markers.ref.value);
      culler.sync();
      expect(culler.dataVersion).toBe(synced);
      // a tile added: what it culls changed
      store.add(
        2,
        createList(50, (i) => (i % 2 ? 3 : 1), true),
        100,
        placement,
      );
      pending = store.takeChanges();
      culler.sync();
      expect(culler.dataVersion).toBeGreaterThan(synced);
      // a member's geometry replaced alone (another target's, its reach with it): no page written, culled again
      const added = culler.dataVersion;
      pending = { written: [], released: [], runs: [], extras: [], members: new Set([1]), freed: [], grown: false };
      culler.sync();
      expect(culler.dataVersion).toBeGreaterThan(added);
      culler.destroy();
    } finally {
      plugin.dispose();
    }
  }, 120_000);
});
