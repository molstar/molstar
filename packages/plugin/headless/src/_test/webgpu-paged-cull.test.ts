/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
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
import {
  createPagedCuller,
  PagedMaxCoarseBias,
  PagedPoseWords,
  PagedRegionMinRoom,
  PagedRegionRoomShare,
} from '@molstar/graphics/gl/compute/paged-cull';
import type { OcclusionPyramid } from '@molstar/graphics/gl/compute/cell-cull';
import type { StorageBuffer } from '@molstar/graphics/gl/compute';
import {
  CopyPoseFloats,
  PageRecordWords,
  ParticlePageStore,
  StorePageBits,
  StorePageSize,
  type StorePlacement,
  writeCopyPose,
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

/** a region's room every instance of its member (the tests reading exact regions) */
const FullRooms = { regionRoomShare: 1 };

/**
 * the regions' pose slots and the poses they read: per region its slots written (from its first, its count's within
 * its room; no raster regions)
 */
async function readRegions(culler: ReturnType<typeof createPagedCuller>, counts: StorageBuffer) {
  const regionCount = culler.regionCount;
  const c = new Uint32Array(await counts.readAsync(0, regionCount * 4));
  const bases = new Uint32Array(await culler.bases.readAsync(0, regionCount * 4));
  const indices = new Uint32Array(await culler.indices.readAsync(0, culler.indexRoom * 4));
  const poses = new Uint32Array(await culler.poses.readAsync(0, culler.poses.byteLength));
  const rooms = culler.regionRooms;
  const regions: number[][] = [];
  let total = 0;
  for (let region = 0; region < regionCount; ++region) {
    const n = Math.min(c[region], rooms[region]);
    const slots: number[] = [];
    for (let k = bases[region]; k < bases[region] + n; ++k) slots.push(indices[k]);
    regions.push(slots);
    total += n;
  }
  return { regions, poses, total, counts: c };
}

/** `culler`'s frame version advanced past `v` (a readback of its, not tracked itself, resolved) */
async function advancedPast(culler: { readonly frameVersion: number }, v: number) {
  for (let i = 0; i < 200 && culler.frameVersion === v; ++i) await new Promise((r) => setTimeout(r, 10));
  return culler.frameVersion > v;
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

/**
 * a camera at `eye` looking at `target`, a narrow frustum (`t`, its half width at its near plane 1 away): its planes,
 * its camera plane, two levels by distance
 */
function createView(eye: Vec3, target: Vec3, t = 0.08) {
  const view = Mat4.lookAt(Mat4(), eye, target, Vec3.create(0, 1, 0));
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
      const culler = createPagedCuller(compute, graphics.resources, FullRooms);
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
      const culler = createPagedCuller(compute, graphics.resources, FullRooms);
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

      // occlusion bits for one selected entry at first: the others not tested against the pyramid (drawn in the first
      // phase, none in the second), the count read back, the bits grown, a cull again; then as with bits for all
      const small = createPagedCuller(compute, graphics.resources, { ...FullRooms, occlusionInitialEntries: 1 });
      small.sync(store, { written: [], released: [], runs: [], extras: [], grown: false }, [2, 3, 1.5, 4]);
      small.cull(view, pyramidOf(1));
      const first = await total(small.counts);
      expect(first).toBeGreaterThan(unoccluded);
      expect(first).toBeLessThan(visible);
      expect(small.cullSecond(pyramidOf(0))).toBe(true);
      expect(await total(small.secondCounts)).toBe(visible - first);
      for (let i = 0; i < 200 && small.frameVersion === 0; ++i) await new Promise((r) => setTimeout(r, 10));
      expect(small.frameVersion).toBe(1);
      small.cull(view, pyramidOf(1));
      expect(await total(small.counts)).toBe(unoccluded);
      expect(small.cullSecond(pyramidOf(0))).toBe(true);
      expect(await total(small.secondCounts)).toBe(visible - unoccluded);
      small.destroy();

      // no second phase to follow (Hi-Z's two phases off, a stereo eye): the first's count read where it culled, the
      // bits grown, a cull again; then as with bits for all, its occluded left out
      const single = createPagedCuller(compute, graphics.resources, { ...FullRooms, occlusionInitialEntries: 1 });
      single.sync(store, { written: [], released: [], runs: [], extras: [], grown: false }, [2, 3, 1.5, 4]);
      single.cull({ ...view, secondPhase: false }, pyramidOf(1));
      expect(await total(single.counts)).toBe(first);
      expect(single.cullSecond(pyramidOf(0))).toBe(false);
      expect(await advancedPast(single, 0)).toBe(true);
      single.cull({ ...view, secondPhase: false }, pyramidOf(1));
      expect(await total(single.counts)).toBe(unoccluded);
      single.destroy();

      // a second phase to follow that does not run (no pyramid of the frame): the first's count read at its end
      const dropped = createPagedCuller(compute, graphics.resources, { ...FullRooms, occlusionInitialEntries: 1 });
      dropped.sync(store, { written: [], released: [], runs: [], extras: [], grown: false }, [2, 3, 1.5, 4]);
      dropped.cull(view, pyramidOf(1));
      dropped.endFirst();
      expect(dropped.cullSecond(pyramidOf(0))).toBe(false);
      expect(await advancedPast(dropped, 0)).toBe(true);
      dropped.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
      culler.destroy();
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("culls the selected entries' groups that survive a test whole, their room grown where it held fewer", async () => {
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
      // three pages: three selected entries
      store.add(
        1,
        createList(20_000, (i) => (i % 2 ? 1 : 3), true),
        0,
        placement,
      );
      const view = createView(Vec3.create(30, 40, 250), Vec3.create(50, 50, 40));
      const total = async (culler: ReturnType<typeof createPagedCuller>) => {
        let n = 0;
        for (const c of new Uint32Array(await culler.counts.readAsync(0, culler.regionCount * 4))) n += c;
        return n;
      };
      const full = createPagedCuller(graphics.compute!, graphics.resources, FullRooms);
      full.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      full.cull(view);
      const visible = await total(full);
      expect(visible).toBeGreaterThan(100);
      // room for one entry's groups: the first's alone culled, the count read back, the room grown, a cull again asked
      // for; then every instance
      const small = createPagedCuller(graphics.compute!, graphics.resources, {
        ...FullRooms,
        survivorsInitialEntries: 1,
      });
      small.sync(store, { written: [], released: [], runs: [], extras: [], grown: false }, [2, 3, 1.5, 4]);
      small.cull(view);
      const first = await total(small);
      expect(first).toBeLessThan(visible);
      expect(await advancedPast(small, 0)).toBe(true);
      small.cull(view);
      expect(await total(small)).toBe(visible);
      small.destroy();
      full.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it('leaves out a copy the pyramid hides whole, the second phase adding it where the frame shows it', async () => {
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
        -1,
        createList(300, (i) => (i < 150 ? 1 : 3), true),
        0,
        placement,
        undefined,
        true,
      );
      // three copies one behind the other, their bounds known (a copy tested whole)
      const center = Vec3.create(50, 50, 40);
      const poses = new Float32Array(3 * CopyPoseFloats);
      for (let k = 0; k < 3; ++k) writeCopyPose(poses, k * 8, Vec3.create(0, 0, -60 * k), Quat.identity(), center);
      store.setCopies(-1, { tops: Int32Array.of(5, 6, 7), poses, center, spacing: 5, radius: 100 });
      const culler = createPagedCuller(graphics.compute!, graphics.resources, FullRooms);
      culler.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      // a wide view (each copy's sphere within it: a sphere past its edge not tested against the pyramid)
      const eye = Vec3.create(50, 50, 250);
      const target = Vec3.create(50, 50, 0);
      const t = 1;
      const view = createView(eye, target, t);
      /** a pyramid of one texel of `depth` (reversed: 1 at the near plane) over a viewport of a pixel */
      const pyramidOf = (depth: number): OcclusionPyramid => {
        const texture = graphics.resources.texture('image-float32', 'r', 'float', 'nearest');
        texture.load({ array: Float32Array.of(depth), width: 1, height: 1 });
        return {
          texture,
          levelCount: 1,
          minLevel: 0,
          width: 1,
          height: 1,
          view: Mat4.lookAt(Mat4(), eye, target, Vec3.create(0, 1, 0)),
          projection: Mat4.perspective(Mat4(), -t, t, t, -t, 1, 400),
          near: 1,
          far: 400,
          version: 1,
        } as OcclusionPyramid;
      };
      const total = async (counts: StorageBuffer) => {
        let n = 0;
        for (const c of new Uint32Array(await counts.readAsync(0, culler.regionCount * 4))) n += c;
        return n;
      };
      const copies = async () => new Set((await culler.readSelectedEntries()).map((e) => e.copy));
      culler.cull(view);
      const visible = await total(culler.counts);
      expect(visible).toBeGreaterThan(300);
      expect([...(await copies())].sort()).toEqual([0, 1, 2]);
      // all at the near plane: every copy hidden whole, none of its entries culled (none drawn)
      culler.cull(view, pyramidOf(1));
      expect(await copies()).toEqual(new Set());
      expect(await total(culler.counts)).toBe(0);
      // the second phase against a pyramid at the far plane: the copies added after the first's entries, all drawn
      expect(culler.cullSecond(pyramidOf(0))).toBe(true);
      expect([...(await copies())].sort()).toEqual([0, 1, 2]);
      expect(await total(culler.secondCounts)).toBe(visible);
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
      const culler = createPagedCuller(compute, graphics.resources, FullRooms);
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
      const culler = createPagedCuller(compute, graphics.resources, FullRooms);
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
      const culler = createPagedCuller(compute, graphics.resources, FullRooms);
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
      const merged = createMergedValues('spheres', members[0], members, true, false, {
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
  it('grows its extras in place: the pages synced before read theirs still, the tables not laid out again', async () => {
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
      const reaches = [2, 3, 1.5, 4];
      const culler = createPagedCuller(graphics.compute!, graphics.resources, { ...FullRooms, wantsReadMs: 0 });
      // alternating targets: a member a particle, in the extras (1024 particles' room at first)
      const alternating = (i: number) => [1, 2, 3][i % 3];
      store.add(1, createList(1000, alternating, true), 0, placement);
      culler.sync(store, store.takeChanges(), reaches);
      const view = createView(Vec3.create(50, 50, 300), Vec3.create(50, 50, 40), 1.2);
      culler.cull(view);
      const before = await readRegions(culler, culler.counts);
      expect(before.total).toBe(1000);
      store.add(2, createList(1000, alternating, true), 1000, placement);
      expect(store.extraCapacity).toBeGreaterThan(1024);
      culler.sync(store, store.takeChanges(), reaches);
      culler.cull(view);
      // every instance in its member's region, the first page's as before (its extras copied, not lost)
      const { regions, poses, total } = await readRegions(culler, culler.counts);
      expect(total).toBe(2000);
      const words = new Uint32Array(poses.buffer);
      for (let region = 0; region < regions.length; ++region) {
        for (const slot of regions[region]) expect(words[slot * PagedPoseWords + 7]).toBe(Math.floor(region / 2));
      }
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("culls a template's copies: a copy's instances its own, their poses composed into the frame's slots", async () => {
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
      // leaves drawn as they are (target 2), a template of two targets (1, 3) placed by three copies, one turned
      store.add(
        0,
        createList(10, () => 2, true),
        0,
        placement,
      );
      store.add(
        -1,
        createList(300, (i) => (i < 150 ? 1 : 3), true),
        0,
        placement,
        undefined,
        true,
      );
      const turn = Quat.setAxisAngle(Quat(), Vec3.create(0, 0, 1), Math.PI / 2);
      const tops = Int32Array.of(5, 6, 7);
      const copyPoses = new Float32Array(3 * CopyPoseFloats);
      const center = Vec3.create(50, 50, 40);
      writeCopyPose(copyPoses, 0, Vec3.create(0, 0, 0), Quat.identity(), center);
      writeCopyPose(copyPoses, 8, Vec3.create(200, 0, 0), turn, center);
      writeCopyPose(copyPoses, 16, Vec3.create(0, 200, 0), Quat.identity(), center);
      store.setCopies(-1, { tops, poses: copyPoses });
      expect(store.pool.counts[store.pool.get(1)]).toBe(150 * 3);
      const reaches = [2, 3, 1.5, 4];
      const culler = createPagedCuller(compute, graphics.resources, FullRooms);
      culler.sync(store, store.takeChanges(), reaches);
      expect(culler.instanceTotal).toBe(10 + 300 * 3);

      /** the world point of each pose the regions hold (p = t + s R c), by member, and their count */
      const drawn = async () => {
        const { regions, poses } = await readRegions(culler, culler.counts);
        const words = new Uint32Array(poses.buffer);
        const floats = new Float32Array(poses.buffer);
        const out: Drawn[] = [];
        const levels = 2;
        for (let region = 0; region < store.pool.count * levels; ++region) {
          const m = Math.floor(region / levels);
          for (const slot of regions[region]) {
            const o = slot * PagedPoseWords;
            const q = Quat.create(...unpackSnorm(words[o + 4]), ...unpackSnorm(words[o + 5]));
            Quat.normalize(q, q);
            const c = Vec3.create(
              store.pool.centers[m * 3],
              store.pool.centers[m * 3 + 1],
              store.pool.centers[m * 3 + 2],
            );
            const p = Vec3.scaleAndAdd(
              Vec3(),
              Vec3.create(floats[o], floats[o + 1], floats[o + 2]),
              Vec3.transformQuat(c, c, q),
              floats[o + 3],
            );
            // the segment the member in the low 16 bits, a copy's instance its copy plus one above (none 0)
            expect(words[o + 7] & 0xffff).toBe(m);
            out.push({ m, p, copy: words[o + 7] >>> 16 });
          }
        }
        return out;
      };
      /**
       * the CPU's: the leaves' composed poses and the copies' (`composeCopies`), their points by member, a copy's its
       * copy plus one (tops 5, 6, 7 copies 0, 1, 2)
       */
      const expected = (inView: (p: Vec3) => boolean) => {
        const out: Drawn[] = [];
        const add = (poses: Float32Array, members: Int32Array, tops?: Int32Array) => {
          for (let i = 0; i < members.length; ++i) {
            const m = members[i];
            const q = Quat.create(poses[i * 8 + 4], poses[i * 8 + 5], poses[i * 8 + 6], poses[i * 8 + 7]);
            const c = Vec3.create(
              store.pool.centers[m * 3],
              store.pool.centers[m * 3 + 1],
              store.pool.centers[m * 3 + 2],
            );
            const p = Vec3.scaleAndAdd(
              Vec3(),
              Vec3.create(poses[i * 8], poses[i * 8 + 1], poses[i * 8 + 2]),
              Vec3.transformQuat(c, c, q),
              poses[i * 8 + 3],
            );
            if (inView(p)) out.push({ m, p, copy: tops ? tops[i] - 4 : 0 });
          }
        };
        const leaves = store.compose(store.getBlock(0)!);
        add(leaves.poses, leaves.members);
        const copies = store.composeCopies(-1);
        add(copies.poses, copies.members, copies.tops);
        return out;
      };
      /** each of `a` matched by one of `b` of its member and copy within 0.05 Å (the GPU's quaternions snorm16) */
      const matched = (a: Drawn[], b: Drawn[]) => {
        const used = new Set<number>();
        let n = 0;
        for (const x of a) {
          const k = b.findIndex(
            (y, i) => !used.has(i) && y.m === x.m && y.copy === x.copy && Vec3.distance(x.p, y.p) < 0.05,
          );
          if (k < 0) continue;
          used.add(k);
          n += 1;
        }
        return n;
      };

      // every copy in view
      culler.cull(createView(Vec3.create(130, 130, 300), Vec3.create(130, 130, 0), 1.2));
      const all = await drawn();
      expect(all.length).toBe(10 + 900);
      const ex = expected(() => true);
      expect(matched(ex, all)).toBe(910);
      // the first copy's alone in view (a narrow frustum about it): the others culled whole
      culler.cull(createView(Vec3.create(50, 50, 300), Vec3.create(50, 50, 40), 0.15));
      const one = await drawn();
      expect(one.length).toBeGreaterThan(50);
      expect(one.every(({ p }) => p[0] < 150 && p[1] < 150)).toBe(true);

      // the frame's slots by what a cull needs: a room of 100 at first, the copies' instances past it not drawn, the
      // cursor read back, a cull again and a draw asked for (pending work); not grown where a pyramid is to come (an
      // app's first frames: nothing occluded), grown where none is (Hi-Z off), the next cull's every instance
      const small = createPagedCuller(compute, graphics.resources, { ...FullRooms, frameInitialSlots: 100 });
      small.sync(store, { written: [], released: [], runs: [], extras: [], grown: false }, reaches);
      const wide = createView(Vec3.create(130, 130, 300), Vec3.create(130, 130, 0), 1.2);
      /** the culler's frame version advanced past `v` (its readback, not tracked itself, resolved) */
      const advanced = async (v: number) => {
        for (let i = 0; i < 200 && small.frameVersion === v; ++i) await new Promise((r) => setTimeout(r, 10));
        return small.frameVersion > v;
      };
      const tracked = graphics.resources.pending.trackedCount;
      small.cull({ ...wide, pyramidToCome: true });
      expect((await readRegions(small, small.counts)).total).toBe(10 + 100);
      expect(await advanced(0)).toBe(true);
      expect(graphics.resources.pending.trackedCount).toBeGreaterThan(tracked);
      // a pyramid still to come: the room of 100 again (not grown); then none to come (Hi-Z off): grown
      small.cull({ ...wide, pyramidToCome: true });
      expect((await readRegions(small, small.counts)).total).toBe(10 + 100);
      expect(await advanced(1)).toBe(true);
      small.cull(wide);
      expect((await readRegions(small, small.counts)).total).toBe(10 + 100);
      expect(await advanced(2)).toBe(true);
      small.cull(wide);
      expect((await readRegions(small, small.counts)).total).toBe(10 + 900);
      small.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("selects a template's copies near first", async () => {
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
        -1,
        createList(300, (i) => (i < 150 ? 1 : 3), true),
        0,
        placement,
        undefined,
        true,
      );
      // three copies one behind the other along z
      const center = Vec3.create(50, 50, 40);
      const copyPoses = new Float32Array(3 * CopyPoseFloats);
      for (let k = 0; k < 3; ++k) {
        writeCopyPose(copyPoses, k * 8, Vec3.create(0, 0, -120 * k), Quat.identity(), center);
      }
      store.setCopies(-1, { tops: Int32Array.of(5, 6, 7), poses: copyPoses });
      const culler = createPagedCuller(graphics.compute!, graphics.resources, FullRooms);
      culler.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      /** the copies of the selected entries of a cull of `view`, in their order */
      const orderOf = async (view: ReturnType<typeof createView>) => {
        culler.cull(view);
        // a copy's pages one after another: its runs
        return (await culler.readSelectedEntries()).map((e) => e.copy).filter((c, i, a) => i === 0 || c !== a[i - 1]);
      };
      // from +z the first copy nearest, from -z the last; the cull's workgroups by that order (where more than run at
      // once, as the isg cell's millions: the room the near copies')
      expect(await orderOf(createView(Vec3.create(50, 50, 200), Vec3.create(50, 50, 0), 1.2))).toEqual([0, 1, 2]);
      expect(await orderOf(createView(Vec3.create(50, 50, -300), Vec3.create(50, 50, 0), 1.2))).toEqual([2, 1, 0]);
      culler.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("grows the frame's slots within the device's limit what the persistent slots leave", async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const poseBytes = PagedPoseWords * 4;
      // the persistent slots (the leaves' 10, the template's 300) fit, and some of the copies' 900 instances: not all,
      // nor two copies' 600
      const limit = 700 * poseBytes + 4;
      const compute = Object.create(graphics.compute!, { maxStorageBufferSize: { value: limit } });
      const store = new ParticlePageStore();
      store.add(
        0,
        createList(10, () => 2, true),
        0,
        placement,
      );
      store.add(
        -1,
        createList(300, (i) => (i < 150 ? 1 : 3), true),
        0,
        placement,
        undefined,
        true,
      );
      const center = Vec3.create(50, 50, 40);
      const copyPoses = new Float32Array(3 * CopyPoseFloats);
      writeCopyPose(copyPoses, 0, Vec3.create(0, 0, 0), Quat.identity(), center);
      writeCopyPose(copyPoses, 8, Vec3.create(200, 0, 0), Quat.identity(), center);
      writeCopyPose(copyPoses, 16, Vec3.create(0, 200, 0), Quat.identity(), center);
      store.setCopies(-1, { tops: Int32Array.of(5, 6, 7), poses: copyPoses });
      const culler = createPagedCuller(compute, graphics.resources, { ...FullRooms, frameInitialSlots: 100 });
      culler.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      const wide = createView(Vec3.create(130, 130, 300), Vec3.create(130, 130, 0), 1.2);
      // the room of 100 overflowed, grown to half again what the cull needed: capped by the limit
      culler.cull(wide);
      expect((await readRegions(culler, culler.counts)).total).toBe(10 + 100);
      for (let i = 0; i < 200 && culler.frameVersion === 0; ++i) await new Promise((r) => setTimeout(r, 10));
      expect(culler.frameVersion).toBeGreaterThan(0);
      culler.cull(wide);
      expect(culler.poses.byteLength).toBeLessThanOrEqual(limit);
      const { total } = await readRegions(culler, culler.counts);
      expect(total).toBeGreaterThan(10 + 100);
      expect(total).toBeLessThan(10 + 900);
      // the room at its cap: overflowed again, no cull asked for past the one after the grown (a still view's frames
      // not culled again each)
      culler.cull(wide);
      await readRegions(culler, culler.counts);
      for (let i = 0; i < 20; ++i) await new Promise((r) => setTimeout(r, 10));
      const version = culler.frameVersion;
      culler.cull(wide);
      await readRegions(culler, culler.counts);
      for (let i = 0; i < 20; ++i) await new Promise((r) => setTimeout(r, 10));
      expect(culler.frameVersion).toBe(version);
      // at its cap a cull needing fewer than the last (two copies in view, not three: what a view's own pyramid does
      // after a zoom's first cull of another's) asks for another; the same again none
      const two = createView(Vec3.create(150, 50, 300), Vec3.create(150, 50, 0), 0.6);
      culler.cull(two);
      await readRegions(culler, culler.counts);
      expect(await advancedPast(culler, version)).toBe(true);
      culler.cull(two);
      await readRegions(culler, culler.counts);
      for (let i = 0; i < 20; ++i) await new Promise((r) => setTimeout(r, 10));
      const settled = culler.frameVersion;
      culler.cull(two);
      await readRegions(culler, culler.counts);
      for (let i = 0; i < 20; ++i) await new Promise((r) => setTimeout(r, 10));
      expect(culler.frameVersion).toBe(settled);

      // without a pyramid no entry is tested against the occlusion bits: their room of one entry no bar to growing
      const fewBits = createPagedCuller(compute, graphics.resources, {
        ...FullRooms,
        frameInitialSlots: 100,
        occlusionInitialEntries: 1,
      });
      fewBits.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      fewBits.cull(wide);
      expect((await readRegions(fewBits, fewBits.counts)).total).toBe(10 + 100);
      expect(await advancedPast(fewBits, 0)).toBe(true);
      fewBits.cull(wide);
      expect((await readRegions(fewBits, fewBits.counts)).total).toBeGreaterThan(10 + 100);
      fewBits.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
      culler.destroy();
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("writes a template's leaf pages where the cull wants them (read back), its copies coarse meanwhile, kept", async () => {
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
      // a template of 300 particles and its coarse level of 40, placed by three copies
      store.add(
        -1,
        createList(300, (i) => (i < 150 ? 1 : 3), true),
        0,
        placement,
        undefined,
        true,
      );
      store.add(
        -2,
        createList(40, (i) => (i < 20 ? 1 : 3), true),
        0,
        placement,
        undefined,
        true,
      );
      const center = Vec3.create(50, 50, 40);
      const poses = new Float32Array(3 * CopyPoseFloats);
      writeCopyPose(poses, 0, Vec3.create(0, 0, 0), Quat.identity(), center);
      writeCopyPose(poses, 8, Vec3.create(200, 0, 0), Quat.identity(), center);
      writeCopyPose(poses, 16, Vec3.create(0, 200, 0), Quat.identity(), center);
      const of = { tops: Int32Array.of(5, 6, 7), poses, levels: 1, center, spacing: 5, radius: 100 };
      store.setCopies(-1, { ...of, level: 0 });
      store.setCopies(-2, { ...of, level: 1 });
      const culler = createPagedCuller(graphics.compute!, graphics.resources, { ...FullRooms, wantsReadMs: 0 });
      culler.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      // the leaf block's pages not written before a cull wants them
      const leafPages = store.getBlock(-1)!.pages.length;
      expect(culler.deferredPages).toBe(leafPages);
      // every copy in view; far (a unit's pixels at a distance of 1: the spacing of 5 about 0.02 px at 260): its coarse
      // level drawn, the leaves not wanted; near (a thousand times): the leaves wanted, the copies coarse (lacking
      // them), written by the cull after the wants are read back, drawn there
      const wide = createView(Vec3.create(130, 130, 300), Vec3.create(130, 130, 0), 1.2);
      const at = (scale: number) => ({ ...wide, coarse: { scale, ortho: false, pixels: 4 } });
      culler.cull(at(1));
      expect((await readRegions(culler, culler.counts)).total).toBe(3 * 40);
      expect(culler.deferredPages).toBe(leafPages);
      // the values' blocks the written pages' alone: the coarse level's, not the leaves'
      const levelPages = store.getBlock(-2)!.pages.length;
      expect(culler.valueBlocks).toBe(levelPages);
      const version = culler.frameVersion;
      culler.cull(at(1000));
      expect((await readRegions(culler, culler.counts)).total).toBe(3 * 40);
      expect(await advancedPast(culler, version)).toBe(true);
      culler.cull(at(1000));
      expect(culler.deferredPages).toBe(0);
      expect(culler.valueBlocks).toBe(levelPages + leafPages);
      expect((await readRegions(culler, culler.counts)).total).toBe(3 * 300);
      // far again: the coarse level, the leaves kept; near again drawn
      culler.cull(at(1));
      expect(culler.deferredPages).toBe(0);
      expect((await readRegions(culler, culler.counts)).total).toBe(3 * 40);
      culler.cull(at(1000));
      expect(culler.deferredPages).toBe(0);
      const near = await readRegions(culler, culler.counts);
      expect(near.total).toBe(3 * 300);
      // every slot a pose of the template's pages (its store indices)
      const pages = new Set(store.getBlock(-1)!.pages);
      for (const slots of near.regions) {
        for (const k of slots) expect(pages.has(near.poses[k * PagedPoseWords + 6] >> StorePageBits)).toBe(true);
      }
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
      culler.destroy();
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("draws a template's copies coarser while a frame needs more than its budget, finer again where it fits", async () => {
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
      // a template of 300 particles and its coarse level of 40, placed by three copies
      store.add(
        -1,
        createList(300, (i) => (i < 150 ? 1 : 3), true),
        0,
        placement,
        undefined,
        true,
      );
      store.add(
        -2,
        createList(40, (i) => (i < 20 ? 1 : 3), true),
        0,
        placement,
        undefined,
        true,
      );
      const center = Vec3.create(50, 50, 40);
      const poses = new Float32Array(3 * CopyPoseFloats);
      writeCopyPose(poses, 0, Vec3.create(0, 0, 0), Quat.identity(), center);
      writeCopyPose(poses, 8, Vec3.create(200, 0, 0), Quat.identity(), center);
      writeCopyPose(poses, 16, Vec3.create(0, 200, 0), Quat.identity(), center);
      const of = { tops: Int32Array.of(5, 6, 7), poses, levels: 1, center, spacing: 5, radius: 100 };
      store.setCopies(-1, { ...of, level: 0 });
      store.setCopies(-2, { ...of, level: 1 });
      const culler = createPagedCuller(graphics.compute!, graphics.resources, { ...FullRooms, wantsReadMs: 0 });
      culler.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      // near (the spacing about 17 px): every copy the template whole, its leaves written after a read of the wants
      const wide = createView(Vec3.create(130, 130, 300), Vec3.create(130, 130, 0), 1.2);
      const near = { ...wide, coarse: { scale: 1000, ortho: false, pixels: 4 } };
      const version = culler.frameVersion;
      culler.cull(near);
      expect(await advancedPast(culler, version)).toBe(true);
      culler.cull(near);
      expect((await readRegions(culler, culler.counts)).total).toBe(3 * 300);
      expect(culler.coarseBias).toBe(0);
      /** culled within `max` till its bias settles (a read back changing it culls again), its instances then */
      const settle = async (max: number) => {
        for (let i = 0; i < 2 * PagedMaxCoarseBias; ++i) {
          const v = culler.frameVersion;
          culler.cull({ ...near, maxInstances: max });
          await readRegions(culler, culler.counts);
          if (!(await advancedPast(culler, v))) break;
        }
        return (await readRegions(culler, culler.counts)).total;
      };
      // a budget of 500 (the leaves' 900 past it): a step coarser each read back (the coarse pixels doubled) till the
      // copies draw their coarse level, kept there (its 120 times eight past the budget)
      expect(await settle(500)).toBe(3 * 40);
      expect(culler.coarseBias).toBeGreaterThan(0);
      expect(culler.frameNeed).toBe(3 * 40);
      // a budget of 2,000: finer again (its 120 times eight within it) till the template whole, kept (900 times eight
      // past it); none: no bias
      expect(await settle(2000)).toBe(3 * 300);
      expect(await settle(0)).toBe(3 * 300);
      expect(culler.coarseBias).toBe(0);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
      culler.destroy();
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it('writes a wanted leaf block within a budget a cull, its copies drawn coarse till all its pages are', async () => {
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
      // a template of three pages (20,000 particles) and its coarse level of 40, placed by three copies
      store.add(
        -1,
        createList(20_000, (i) => (i % 2 ? 1 : 3), true),
        0,
        placement,
        undefined,
        true,
      );
      store.add(
        -2,
        createList(40, (i) => (i < 20 ? 1 : 3), true),
        0,
        placement,
        undefined,
        true,
      );
      const center = Vec3.create(50, 50, 40);
      const poses = new Float32Array(3 * CopyPoseFloats);
      writeCopyPose(poses, 0, Vec3.create(0, 0, 0), Quat.identity(), center);
      writeCopyPose(poses, 8, Vec3.create(200, 0, 0), Quat.identity(), center);
      writeCopyPose(poses, 16, Vec3.create(0, 200, 0), Quat.identity(), center);
      const of = { tops: Int32Array.of(5, 6, 7), poses, levels: 1, center, spacing: 5, radius: 100 };
      store.setCopies(-1, { ...of, level: 0 });
      store.setCopies(-2, { ...of, level: 1 });
      const leafPages = store.getBlock(-1)!.pages.length;
      expect(leafPages).toBe(3);
      // a page's positions and rotations a cull (grid positions and snorm16 rotations: 16 bytes a particle); a frame's
      // room of 1,000 (the coarse level's 120 within it, the leaves' 60,000 not)
      const culler = createPagedCuller(graphics.compute!, graphics.resources, {
        ...FullRooms,
        leafUploadBytes: 8192 * 16,
        frameInitialSlots: 1000,
        wantsReadMs: 0,
      });
      culler.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      const wide = createView(Vec3.create(130, 130, 300), Vec3.create(130, 130, 0), 1.2);
      const near = { ...wide, coarse: { scale: 1000, ortho: false, pixels: 4 } };
      // near, every copy at level 0 and every page in view: wanted (the copies coarse), read back; then two culls a page
      // each, the copies at their coarse level meanwhile (a page in view not written), another cull asked for; the
      // third the last page: the leaves drawn
      const wanted = culler.frameVersion;
      culler.cull(near);
      expect((await readRegions(culler, culler.counts)).total).toBe(3 * 40);
      expect(await advancedPast(culler, wanted)).toBe(true);
      for (let i = 0; i < 2; ++i) {
        const version = culler.frameVersion;
        culler.cull(near);
        expect(culler.deferredPages).toBe(leafPages - i - 1);
        expect((await readRegions(culler, culler.counts)).total).toBe(3 * 40);
        expect(culler.frameVersion).toBeGreaterThan(version);
      }
      /** the culler's frame version advanced past `v` (its readback, not tracked itself, resolved) */
      const advanced = async (v: number) => {
        for (let i = 0; i < 200 && culler.frameVersion === v; ++i) await new Promise((r) => setTimeout(r, 10));
        return culler.frameVersion > v;
      };
      // the cull making it resident past the frame's room (the phases' 875, the fallback's 125): a cull asked for, the
      // room not grown by its counts (its pyramid's not of these pages); the copies past it dropped, each drawn at its
      // coarse level into the fallback's (no holes: 875 leaves and the coarse level's 3 x 40); the next's grown, the one
      // after it every leaf
      const levelPages = new Set(store.getBlock(-2)!.pages);
      const fallen = async () => {
        const d = await readRegions(culler, culler.counts);
        let coarse = 0;
        for (const slots of d.regions) {
          for (const k of slots) if (levelPages.has(d.poses[k * PagedPoseWords + 6] >> StorePageBits)) coarse += 1;
        }
        return { total: d.total, coarse };
      };
      let version = culler.frameVersion;
      culler.cull(near);
      expect(culler.deferredPages).toBe(0);
      expect(await fallen()).toEqual({ total: 875 + 3 * 40, coarse: 3 * 40 });
      expect(await advanced(version)).toBe(true);
      version = culler.frameVersion;
      culler.cull(near);
      expect(await fallen()).toEqual({ total: 875 + 3 * 40, coarse: 3 * 40 });
      expect(await advanced(version)).toBe(true);
      culler.cull(near);
      const drawn = await readRegions(culler, culler.counts);
      expect(drawn.total).toBe(3 * 20_000);
      const pages = new Set(store.getBlock(-1)!.pages);
      for (const slots of drawn.regions) {
        for (const k of slots) expect(pages.has(drawn.poses[k * PagedPoseWords + 6] >> StorePageBits)).toBe(true);
      }
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
      culler.destroy();
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it("writes the leaf pages a near copy's view reaches, lets go of those no view wanted for a while", async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    const realNow = performance.now.bind(performance);
    let offset = 0;
    vi.spyOn(performance, 'now').mockImplementation(() => realNow() + offset);
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const store = new ParticlePageStore();
      // a template of three pages apart in space (a page's particles 300 Å from the next's) and its coarse level
      const template = createList(3 * 8192, () => 1, false);
      const xyz = template.coordinates as Float32Array;
      for (let i = 0; i < template.count; ++i) xyz[i * 3] = Math.floor(i / 8192) * 300 + (i % 100);
      store.add(-1, { ...template, coordinates: GridPositions.fromFloat(xyz) }, 0, placement, undefined, true);
      store.add(
        -2,
        createList(40, () => 1, true),
        0,
        placement,
        undefined,
        true,
      );
      const center = Vec3.create(350, 50, 45);
      const poses = new Float32Array(CopyPoseFloats);
      writeCopyPose(poses, 0, Vec3.create(350, 50, 45), Quat.identity(), center);
      const of = { tops: Int32Array.of(5), poses, levels: 1, center, spacing: 5, radius: 400 };
      store.setCopies(-1, { ...of, level: 0 });
      store.setCopies(-2, { ...of, level: 1 });
      const pages = store.getBlock(-1)!.pages;
      expect(pages.length).toBe(3);
      const culler = createPagedCuller(graphics.compute!, graphics.resources, { ...FullRooms, wantsReadMs: 0 });
      culler.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      /** near (a thousand pixels a unit at 1), a narrow view of the template's page `k` */
      const at = (k: number) => ({
        ...createView(Vec3.create(k * 300 + 50, 50, 300), Vec3.create(k * 300 + 50, 50, 45), 0.15),
        coarse: { scale: 1000, ortho: false, pixels: 4 },
      });
      /** the store pages of the poses drawn */
      const drawnPages = async () => {
        const { regions, poses: p } = await readRegions(culler, culler.counts);
        const out = new Set<number>();
        for (const slots of regions) for (const k of slots) out.add(p[k * PagedPoseWords + 6] >> StorePageBits);
        return out;
      };
      const levelPage = store.getBlock(-2)!.pages[0];
      /** a cull of `view`, its wants read back, a cull again: what it then draws */
      const settle = async (view: ReturnType<typeof at>) => {
        const v = culler.frameVersion;
        culler.cull(view);
        await readRegions(culler, culler.counts);
        if (await advancedPast(culler, v)) culler.cull(view);
        return drawnPages();
      };
      // the first page in view: wanted, the copy coarse meanwhile; it alone written, the copy's leaves drawn from it
      culler.cull(at(0));
      expect([...(await drawnPages())]).toEqual([levelPage]);
      expect([...(await settle(at(0)))]).toEqual([pages[0]]);
      expect(culler.deferredPages).toBe(2);
      // the third: written too, drawn
      expect([...(await settle(at(2)))]).toEqual([pages[2]]);
      expect(culler.deferredPages).toBe(1);
      // past the release time the first let go (its value block), the third written again (wanted)
      offset += 20_000;
      expect([...(await settle(at(2)))]).toEqual([pages[2]]);
      expect(culler.deferredPages).toBe(2);
      expect(culler.valueBlocks).toBe(store.getBlock(-2)!.pages.length + 1);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
      culler.destroy();
    } finally {
      vi.restoreAllMocks();
      plugin.dispose();
    }
  }, 120_000);

  it('keeps its value blocks within a budget: the leaf pages the last wants did not want let go for room', async () => {
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
      // a template of three pages apart in space and its coarse level (one page), as above
      const template = createList(3 * 8192, () => 1, false);
      const xyz = template.coordinates as Float32Array;
      for (let i = 0; i < template.count; ++i) xyz[i * 3] = Math.floor(i / 8192) * 300 + (i % 100);
      store.add(-1, { ...template, coordinates: GridPositions.fromFloat(xyz) }, 0, placement, undefined, true);
      store.add(
        -2,
        createList(40, () => 1, true),
        0,
        placement,
        undefined,
        true,
      );
      const center = Vec3.create(350, 50, 45);
      const poses = new Float32Array(CopyPoseFloats);
      writeCopyPose(poses, 0, Vec3.create(350, 50, 45), Quat.identity(), center);
      const of = { tops: Int32Array.of(5), poses, levels: 1, center, spacing: 5, radius: 400 };
      store.setCopies(-1, { ...of, level: 0 });
      store.setCopies(-2, { ...of, level: 1 });
      const pages = store.getBlock(-1)!.pages;
      // room for two value blocks (a page's grid positions and float32 rotations, the template's: 24 bytes a particle):
      // the coarse level's and one leaf page
      const pageBytes = 8192 * 24;
      const culler = createPagedCuller(graphics.compute!, graphics.resources, {
        ...FullRooms,
        wantsReadMs: 0,
        valueBytes: 2 * pageBytes,
      });
      culler.sync(store, store.takeChanges(), [2, 3, 1.5, 4]);
      const at = (k: number) => ({
        ...createView(Vec3.create(k * 300 + 50, 50, 300), Vec3.create(k * 300 + 50, 50, 45), 0.15),
        coarse: { scale: 1000, ortho: false, pixels: 4 },
      });
      const drawnPages = async () => {
        const { regions, poses: p } = await readRegions(culler, culler.counts);
        const out = new Set<number>();
        for (const slots of regions) for (const k of slots) out.add(p[k * PagedPoseWords + 6] >> StorePageBits);
        return out;
      };
      const settle = async (view: ReturnType<typeof at>) => {
        const v = culler.frameVersion;
        culler.cull(view);
        await readRegions(culler, culler.counts);
        if (await advancedPast(culler, v)) culler.cull(view);
        return drawnPages();
      };
      expect([...(await settle(at(0)))]).toEqual([pages[0]]);
      expect(culler.valueBlocks).toBe(2);
      // the third wanted: the first (not wanted by the last read) let go for its room, well before the release time
      expect([...(await settle(at(2)))]).toEqual([pages[2]]);
      expect(culler.valueBlocks).toBe(2);
      expect(culler.deferredPages).toBe(2);
      // the first again: written again, the third let go
      expect([...(await settle(at(0)))]).toEqual([pages[0]]);
      expect(culler.valueBlocks).toBe(2);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
      culler.destroy();
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it('sizes the regions by what a cull needs', async () => {
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
      store.add(
        1,
        createList(20_000, () => 1, true),
        0,
        placement,
      );
      const reaches = [2, 3, 1.5, 4];
      const changes = store.takeChanges();
      const culler = createPagedCuller(compute, graphics.resources);
      culler.sync(store, changes, reaches);
      const full = createPagedCuller(compute, graphics.resources, FullRooms);
      full.sync(store, changes, reaches);
      // every particle in view, about half at each level, all of them the ghost's (no clip objects): three regions
      const view = { ...createView(Vec3.create(50, 50, 280), Vec3.create(50, 50, 45), 1.2), ghost: true };
      full.cull(view);
      const expected = await readRegions(full, full.counts);
      expect(expected.regions.length).toBe(3);
      expect(expected.regions[2].length).toBe(20_000);
      expect(expected.total).toBe(40_000);
      // a region's first room: a sixteenth of its member's 20,000 instances, at least `PagedRegionMinRoom`
      const room = Math.max(PagedRegionMinRoom, Math.ceil(20_000 * PagedRegionRoomShare));
      expect(expected.regions.every((r) => r.length > room)).toBe(true);
      // the first rooms within a budget in all: the shares scaled down (a half a region: 30,000, within 15,000 a
      // quarter, each at least `PagedRegionMinRoom`)
      const budgeted = createPagedCuller(compute, graphics.resources, {
        regionRoomShare: 1 / 2,
        regionInitialSlots: 15_000,
      });
      budgeted.sync(store, changes, reaches);
      budgeted.cull({ ...view, pyramidToCome: true });
      expect(Array.from(budgeted.regionRooms)).toEqual([5000, 5000, 5000]);
      budgeted.destroy();
      const advanced = async (v: number) => {
        for (let i = 0; i < 200 && culler.frameVersion === v; ++i) await new Promise((r) => setTimeout(r, 10));
        return culler.frameVersion > v;
      };
      // a sixteenth of the member's instances (at least `PagedRegionMinRoom`) a region at first: those past them not
      // written, their counts read back, a cull again; not grown where a pyramid is to come, grown where none is
      culler.cull({ ...view, pyramidToCome: true });
      expect(Array.from(culler.regionRooms)).toEqual([room, room, room]);
      expect(culler.indexRoom).toBe(3 * room);
      const first = await readRegions(culler, culler.counts);
      expect(Array.from(first.counts)).toEqual(expected.regions.map((r) => r.length));
      expect(first.total).toBe(3 * room);
      // what was written the instances of the full regions'
      for (let r = 0; r < 3; ++r) {
        const all = new Set(expected.regions[r]);
        expect(first.regions[r].every((s) => all.has(s))).toBe(true);
      }
      expect(await advanced(0)).toBe(true);
      culler.cull({ ...view, pyramidToCome: true });
      expect(Array.from(culler.regionRooms)).toEqual([room, room, room]);
      expect(await advanced(1)).toBe(true);
      culler.cull(view);
      expect(await advanced(2)).toBe(true);
      culler.cull(view);
      const grown = await readRegions(culler, culler.counts);
      expect(grown.total).toBe(40_000);
      // the ghost's grown to its member's every instance
      expect(culler.regionRooms[2]).toBe(20_000);
      for (let r = 0; r < 3; ++r) {
        expect(culler.regionRooms[r]).toBe(
          Math.min(20_000, Math.max(room, Math.ceil(expected.regions[r].length * 1.5))),
        );
        expect([...grown.regions[r]].sort((a, b) => a - b)).toEqual([...expected.regions[r]].sort((a, b) => a - b));
      }
      culler.destroy();
      full.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);
});

/** a drawn instance: its member, its point, its copy plus one (0 none) */
type Drawn = { m: number; p: Vec3; copy: number };

/** two snorm16 of a word (`pack2x16snorm`) */
function unpackSnorm(word: number): [number, number] {
  const lo = (word << 16) >> 16;
  const hi = word >> 16;
  return [Math.max(-1, lo / 32767), Math.max(-1, hi / 32767)];
}
