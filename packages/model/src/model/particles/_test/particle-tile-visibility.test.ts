/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Box3D } from '@molstar/core/math/geometry/primitives/box3d';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Clip } from '@molstar/core/util/clip';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import type { ParticleTarget } from '../particle-list.js';
import { type ParticleTileCamera, ParticleTileVisibility } from '../particle-tile-visibility.js';
import { type ParticleTile, ParticleTileCatalog, ParticleTileGrid, type ParticleTileSet } from '../particle-tiles.js';

type B = [number, number, number, number, number, number];
const box = ([a, b, c, d, e, f]: B) => Box3D.create(Vec3.create(a, b, c), Vec3.create(d, e, f));

/** tiles of one target of `radius` (unknown where null), their boxes */
function createSet(tiles: { box: B }[], radius: number | null = 1): ParticleTileSet {
  const map = new Map<number, ParticleTile>();
  tiles.forEach((t, id) => {
    const b = box(t.box);
    map.set(id, {
      id,
      box: b,
      count: 1,
      base: id,
      targets: Int32Array.of(1),
    });
  });
  return {
    catalog: {
      targets: new Map([[1, { id: 1, kind: 'sites', radius: radius ?? undefined }]]),
      targetInfo: new Map([[1, { name: 'one' }]]),
    },
    grid: ParticleTileGrid.create(Vec3.create(-1, -1, -1), Vec3.create(1, 1, 1), [1, 1, 1]),
    tiles: map,
    count: map.size,
    frameCount: 1,
    provider: {
      getTile: () => {
        throw new Error('not here');
      },
      getTarget: () => {
        throw new Error('not here');
      },
    },
  };
}

const HalfHeight = 256;

/** looking down -z from `eye`, a 90° field of view (perspective) or ±50 Å (orthographic), a viewport 512 high */
function camera(eye: Vec3, orthographic = false, near = 1): ParticleTileCamera {
  const view = Mat4.lookAt(Mat4(), eye, Vec3.create(eye[0], eye[1], eye[2] - 1), Vec3.create(0, 1, 0));
  const projection = orthographic
    ? Mat4.ortho(Mat4(), -50, 50, 50, -50, near, 1e6)
    : Mat4.perspective(Mat4(), -near, near, near, -near, near, 1e6);
  return {
    position: eye,
    direction: Vec3.create(0, 0, -1),
    viewProjection: Mat4.mul(Mat4(), projection, view),
    orthographic,
    // the focal scale (1 at 90°; 1 / 50 Å orthographic) times half the viewport's height
    pixelScale: orthographic ? HalfHeight / 50 : HalfHeight,
    near,
  };
}

const ids = (v: ParticleTileVisibility) => v.tiles.map((t) => t.id);
const eye = Vec3.create(0, 0, 100);

describe('particle tile visibility', () => {
  it('leaves out a tile out of the frustum and one sub-pixel, orders by size on screen', () => {
    const set = createSet([
      { box: [-5, -5, -60, 5, 5, -50] },
      { box: [-5, -5, 0, 5, 5, 10] },
      { box: [500, -5, 0, 510, 5, 10] }, // beside the view
      { box: [-5, -5, -1e5, 5, 5, -1e5 + 10] }, // a radius of 1 Å at 100,000 Å: sub-pixel
    ]);
    const v = ParticleTileVisibility.get(set, camera(eye));
    expect(ids(v)).toEqual([1, 0]);
    expect(v.culled).toEqual({ frustum: 1, clipped: 0, size: 1, occluded: 0 });
    // 1 Å at 89 Å (the reach's nearest point), 256 pixels at a distance of 1
    expect(v.tiles[0].distance).toBeCloseTo(89);
    expect(v.tiles[0].size).toBeCloseTo(256 / 89);
    // a smaller threshold keeps the far one
    expect(ids(ParticleTileVisibility.get(set, camera(eye), { minPixels: 1e-3 }))).toEqual([1, 0, 3]);
  });

  it("keeps a tile whose box is out of view but its targets' reach in, keeps a tile of an unknown radius", () => {
    // the frustum's side at z 0 is x 100: a box from 100.5, its particles' targets reaching 2 Å past it
    const tiles = [{ box: [100.5, -5, -5, 110, 5, 5] as B }];
    expect(ids(ParticleTileVisibility.get(createSet(tiles, 2), camera(eye)))).toEqual([0]);
    expect(ids(ParticleTileVisibility.get(createSet(tiles, 0.1), camera(eye)))).toEqual([]);
    const unknown = ParticleTileVisibility.get(
      createSet([{ box: [-5, -5, -1e5, 5, 5, -1e5 + 10] }], null),
      camera(eye),
    );
    expect(unknown.tiles[0].size).toBe(Infinity);
  });

  it('sizes by the scale alone where orthographic', () => {
    const set = createSet([{ box: [-5, -5, -1e5, 5, 5, -1e5 + 10] }, { box: [-5, -5, 0, 5, 5, 10] }]);
    const v = ParticleTileVisibility.get(set, camera(eye, true));
    expect(v.tiles.map((t) => t.size)).toEqual([256 / 50, 256 / 50]);
    // equal sizes: the nearer first
    expect(ids(v)).toEqual([1, 0]);
  });

  /**
   * a column of cells 40 Å wide (`columns` of them along x, from x -20 per column) and 20 Å deep along z, from z -90:
   * cell z 4 at -10..10 (the front, seen from z 100), 1 at -70..-50; the tiles `densities` gives (by x, z), each its
   * density (1/Å) in one block, none for 0
   */
  function cells(densities: [number, number, number][], columns = 1): ParticleTileSet {
    const set = createSet([]);
    // as many rows along y as columns, the rows past the first empty unless set below
    const grid = ParticleTileGrid.create(
      Vec3.create(-20, -20, -90),
      Vec3.create(-20 + 40 * columns, -20 + 40 * columns, 10),
      [columns, columns, 5],
    );
    const tiles = new Map<number, ParticleTile>();
    for (const [x, z, density] of densities) {
      const id = ParticleTileGrid.getId(grid, x, columns > 1 ? 1 : 0, z);
      tiles.set(id, {
        id,
        box: ParticleTileGrid.getBox(grid, id, Box3D()),
        count: 1,
        base: id,
        targets: Int32Array.of(1),
        opacity: density > 0 ? { dimensions: [1, 1, 1], density: Float32Array.of(density) } : undefined,
      });
    }
    return { ...set, grid, tiles, count: tiles.size };
  }
  const zIds = (set: ParticleTileSet, v: ParticleTileVisibility) =>
    ids(v).map((id) => Math.floor(id / (set.grid.dimensions[0] * set.grid.dimensions[1])));

  it('leaves out a tile behind a dense one, keeps it behind one not drawn opaque or of no opacity', () => {
    // 1/Å over 20 Å: 20
    const set = cells([
      [0, 4, 1],
      [0, 1, 0],
    ]);
    const v = ParticleTileVisibility.get(set, camera(eye));
    expect(zIds(set, v)).toEqual([4]);
    expect(v.culled.occluded).toBe(1);
    expect(ParticleTileVisibility.get(set, camera(eye), { isOpaque: (t) => t.id !== 4 }).tiles.length).toBe(2);
    expect(ParticleTileVisibility.get(set, camera(eye), { occlusion: false }).tiles.length).toBe(2);
    const clear = cells([
      [0, 4, 0],
      [0, 1, 0],
    ]);
    expect(ParticleTileVisibility.get(clear, camera(eye)).tiles.length).toBe(2);
  });

  it('orders the least occluded first: a tile seen through a thin one after a smaller one in clear view', () => {
    // two columns: in the first a thin tile (0.1/Å over 20 Å: 2, about 14% let through) in front of another, nearer
    // and so larger on screen than the second column's farther tile in clear view
    const set = cells(
      [
        [0, 4, 0.1],
        [0, 2, 0.1],
        [1, 1, 0.1],
      ],
      2,
    );
    const v = ParticleTileVisibility.get(set, camera(eye));
    const byId = new Map(v.tiles.map((t) => [t.id, t]));
    // ids: x + 2 (the second row) + 4 z
    expect(byId.get(18)!.exposure).toBe(1);
    // the far one in clear view but for the near thin tile's edge (perspective), the one behind it seen through it
    // but around its edges: more than exp(-2) of it
    expect(byId.get(7)!.exposure).toBeGreaterThan(0.5);
    expect(byId.get(10)!.exposure).toBeLessThan(0.5);
    expect(byId.get(10)!.exposure).toBeGreaterThan(Math.exp(-2));
    expect(byId.get(10)!.size).toBeGreaterThan(byId.get(7)!.size);
    expect(ids(v)).toEqual([18, 7, 10]);
    // without occlusion all exposed: by size
    expect(ids(ParticleTileVisibility.get(set, camera(eye), { occlusion: false }))).toEqual([18, 10, 7]);
  });

  it('gathers the optical depth of the tiles in front: two that each let much through hide together', () => {
    // 0.2/Å over 20 Å: 4, under 7 alone; two: 8
    const one = cells([
      [0, 4, 0.2],
      [0, 1, 0],
    ]);
    expect(ParticleTileVisibility.get(one, camera(eye)).tiles.length).toBe(2);
    const two = cells([
      [0, 4, 0.2],
      [0, 3, 0.2],
      [0, 1, 0],
    ]);
    const v = ParticleTileVisibility.get(two, camera(eye));
    expect(zIds(two, v).sort()).toEqual([3, 4]);
    expect(v.culled.occluded).toBe(1);
    // a higher threshold: not
    expect(ParticleTileVisibility.get(two, camera(eye), { opticalDepth: 9 }).tiles.length).toBe(3);
  });

  it('hides only what is behind the back of what made a ray opaque: what is drawn in it may be anywhere in it', () => {
    // the tile just behind the front one: its reach 1 Å into the front's cell, kept; the next hidden
    const set = cells([
      [0, 4, 1],
      [0, 3, 0],
      [0, 2, 0],
    ]);
    expect(zIds(set, ParticleTileVisibility.get(set, camera(eye))).sort()).toEqual([3, 4]);
    // seen from behind it is the other way round
    const back = Vec3.create(0, 0, -200);
    const fromBehind: ParticleTileCamera = {
      ...camera(back),
      direction: Vec3.create(0, 0, 1),
      viewProjection: Mat4.mul(
        Mat4(),
        Mat4.perspective(Mat4(), -1, 1, 1, -1, 1, 1e6),
        Mat4.lookAt(Mat4(), back, Vec3.create(0, 0, 0), Vec3.create(0, 1, 0)),
      ),
    };
    const flipped = cells([
      [0, 4, 0],
      [0, 1, 1],
    ]);
    expect(zIds(flipped, ParticleTileVisibility.get(flipped, fromBehind))).toEqual([1]);
  });

  it('keeps a tile behind a hole, hides it behind neighbours meeting at a seam', () => {
    // three columns of 40 Å, from x -20: the front's middle one missing, the tile behind it seen through
    const holed = cells(
      [
        [0, 4, 1],
        [2, 4, 1],
        [1, 1, 0],
      ],
      3,
    );
    const at = (x: number) => Vec3.create(x, 40, 100);
    expect(ParticleTileVisibility.get(holed, camera(at(40))).tiles.length).toBe(3);
    // all three in front: the seams between them no gap
    const closed = cells(
      [
        [0, 4, 1],
        [1, 4, 1],
        [2, 4, 1],
        [1, 1, 0],
      ],
      3,
    );
    expect(ParticleTileVisibility.get(closed, camera(at(40))).culled.occluded).toBe(1);
    // the seam off the rays: the camera moved by 0.37 Å
    expect(ParticleTileVisibility.get(closed, camera(at(40.37))).culled.occluded).toBe(1);
    // orthographic: a front of 3 × 3 tiles (the one behind as wide as one, its reach past it)
    const slab = cells([[1, 1, 0]], 3);
    const grid = slab.grid;
    const tiles = new Map(slab.tiles);
    for (let y = 0; y < 3; ++y) {
      for (let x = 0; x < 3; ++x) {
        const id = ParticleTileGrid.getId(grid, x, y, 4);
        const box = ParticleTileGrid.getBox(grid, id, Box3D());
        const opacity = { dimensions: [1, 1, 1] as const, density: Float32Array.of(1) };
        tiles.set(id, { id, box, count: 1, base: id, targets: Int32Array.of(1), opacity });
      }
    }
    const ortho = camera(at(40), true);
    const square = { ...slab, tiles };
    expect(ParticleTileVisibility.get(square, ortho).culled.occluded).toBe(1);
  });

  /** a drawing (opaque) clipped by `objects` (`Clip.Params`' objects, their defaults otherwise) */
  const clipped = (
    objects: Partial<Omit<PD.Values<typeof Clip.Params>['objects'][number], 'rotation'>> &
      {
        rotation?: { axis: Vec3; angle: number };
      }[],
    opaque = true,
  ) => {
    const defaults = PD.getDefaultValues(Clip.Params.objects.element);
    const clip = Clip.getClip({ variant: 'pixel', objects: objects.map((o) => ({ ...defaults, ...o })) });
    return { clip: clip.objects, opaque };
  };
  // the front tile z -10 to 10, the one behind z -70 to -50, across x and y -20 to 20; the eye at z 100
  const front = () =>
    cells([
      [0, 4, 1],
      [0, 1, 0],
    ]);

  it('cuts the blocks by clip planes, leaves out a tile all clipped away, inverted the other side', () => {
    const set = front();
    // a plane clips where its frame's y is at least its position's: turned about z by 90°, x below 0 clipped (the
    // front's left half gone, the tile behind seen past it)
    const half = clipped([{ type: 'plane', rotation: { axis: Vec3.create(0, 0, 1), angle: 90 } }]);
    expect(ParticleTileVisibility.get(set, camera(eye), { drawings: [half] }).tiles.length).toBe(2);
    // turned about x by 90°: z above -30 clipped, the front tile all clipped away
    const cut = {
      type: 'plane' as const,
      position: Vec3.create(0, 0, -30),
      rotation: { axis: Vec3.create(1, 0, 0), angle: 90 },
    };
    const v = ParticleTileVisibility.get(set, camera(eye), { drawings: [clipped([cut])] });
    expect(zIds(set, v)).toEqual([1]);
    expect(v.culled.clipped).toBe(1);
    // inverted: what is below -30 clipped, the tile behind all clipped away (the front kept, hiding nothing more)
    const inverted = ParticleTileVisibility.get(set, camera(eye), { drawings: [clipped([{ ...cut, invert: true }])] });
    expect(zIds(set, inverted)).toEqual([4]);
    expect(inverted.culled.clipped).toBe(1);
  });

  it('cuts the blocks by any clip shape: a hole through the front seen through, a shape clipping a tile whole', () => {
    const set = front();
    expect(zIds(set, ParticleTileVisibility.get(set, camera(eye)))).toEqual([4]);
    // a sphere of radius 15 about the front's center: a hole through it along the view, the tile behind seen through
    const hole = clipped([{ type: 'sphere', scale: Vec3.create(30, 30, 30) }]);
    expect(ParticleTileVisibility.get(set, camera(eye), { drawings: [hole] }).tiles.length).toBe(2);
    // a cylinder along z (its axis the frame's y: turned about x), the same hole
    const tube = clipped([
      { type: 'cylinder', scale: Vec3.create(30, 100, 30), rotation: { axis: Vec3.create(1, 0, 0), angle: 90 } },
    ]);
    expect(ParticleTileVisibility.get(set, camera(eye), { drawings: [tube] }).tiles.length).toBe(2);
    // a small sphere off the view's axis: no hole the rays pass, the tile behind still hidden
    const dent = clipped([{ type: 'sphere', position: Vec3.create(18, 18, 9), scale: Vec3.create(2, 2, 2) }]);
    expect(zIds(set, ParticleTileVisibility.get(set, camera(eye), { drawings: [dent] }))).toEqual([4]);
    // a cube around the front tile: all of it clipped away, left out
    const box = clipped([{ type: 'cube', scale: Vec3.create(50, 50, 30) }]);
    const v = ParticleTileVisibility.get(set, camera(eye), { drawings: [box] });
    expect(zIds(set, v)).toEqual([1]);
    expect(v.culled.clipped).toBe(1);
    // inverted, keeping only its inside: a small sphere about the front's center, the tile behind clipped away whole
    const core = clipped([{ type: 'sphere', invert: true, scale: Vec3.create(6, 6, 6) }]);
    const kept = ParticleTileVisibility.get(set, camera(eye), { drawings: [core] });
    expect(zIds(set, kept)).toEqual([4]);
    expect(kept.culled.clipped).toBe(1);
    // a rod along z through both: the front's rim gone, the tile behind seen around the rod
    const rod = clipped([
      {
        type: 'cylinder',
        invert: true,
        scale: Vec3.create(16, 300, 16),
        rotation: { axis: Vec3.create(1, 0, 0), angle: 90 },
      },
    ]);
    expect(ParticleTileVisibility.get(set, camera(eye), { drawings: [rod] }).tiles.length).toBe(2);
  });

  it('splits a block a clip surface cuts as often as told: unsplit, it hides nothing', () => {
    const set = front();
    // a small sphere off the view's axis: split, the rest of the front hides the tile behind
    const dent = clipped([{ type: 'sphere', position: Vec3.create(18, 18, 9), scale: Vec3.create(2, 2, 2) }]);
    expect(zIds(set, ParticleTileVisibility.get(set, camera(eye), { drawings: [dent] }))).toEqual([4]);
    expect(zIds(set, ParticleTileVisibility.get(set, camera(eye), { drawings: [dent], clipSplitDepth: 1 }))).toEqual([
      4,
    ]);
    // not split: the blocks it cuts add nothing, the tile behind kept
    const unsplit = ParticleTileVisibility.get(set, camera(eye), { drawings: [dent], clipSplitDepth: 0 });
    expect(unsplit.tiles.length).toBe(2);
  });

  it('hides by what every opaque drawing keeps, shows what any keeps', () => {
    const set = front();
    const hole = clipped([{ type: 'sphere', scale: Vec3.create(30, 30, 30) }]);
    const whole = { opaque: true };
    // another drawing of the front whole hides what is behind (the hole's only one of them)
    expect(ParticleTileVisibility.get(set, camera(eye), { drawings: [whole] }).tiles.length).toBe(1);
    expect(ParticleTileVisibility.get(set, camera(eye), { drawings: [hole, whole] }).tiles.length).toBe(2);
    // drawn not opaque: nothing hides
    expect(ParticleTileVisibility.get(set, camera(eye), { drawings: [{ opaque: false }] }).tiles.length).toBe(2);
    // a tile clipped away by one drawing, kept by the other: shown
    const box = clipped([{ type: 'cube', scale: Vec3.create(50, 50, 30) }]);
    expect(zIds(set, ParticleTileVisibility.get(set, camera(eye), { drawings: [box, { opaque: false }] }))).toEqual([
      4, 1,
    ]);
  });

  it('clips a box to its largest axis-aligned part on the kept side', () => {
    const out = Box3D();
    const b = box([0, 0, 0, 10, 10, 10]);
    expect(ParticleTileVisibility.clipBox(b, [Plane3D.create(Vec3.create(0, -1, 0), 4)], out)).toBe(true);
    expect([...out.min, ...out.max]).toEqual([0, 0, 0, 10, 4, 10]);
    // oblique, x + y <= 15 (normalized): shrunk along x to where every y keeps it, x <= 5
    const n = Vec3.normalize(Vec3(), Vec3.create(-1, -1, 0));
    expect(ParticleTileVisibility.clipBox(b, [Plane3D.create(n, 15 / Math.SQRT2)], out)).toBe(true);
    expect(out.max[0]).toBeCloseTo(5);
    expect(out.max[1]).toBe(10);
    expect(ParticleTileVisibility.clipBox(b, [Plane3D.create(n, 5 / Math.SQRT2)], out)).toBe(false);
  });

  it('gathers from the near plane where the camera is inside a tile, sizes a tile around it as Infinity', () => {
    const inside = Vec3.create(0, 0, 0);
    // from the near plane (1 Å) to the cell's back at -10: 9 Å, 1/Å: 9
    const dense = cells([
      [0, 4, 1],
      [0, 1, 0],
    ]);
    const v = ParticleTileVisibility.get(dense, camera(inside));
    expect(zIds(dense, v)).toEqual([4]);
    expect(v.tiles[0].size).toBe(Infinity);
    expect(v.tiles[0].distance).toBe(0);
    // 0.5/Å: 4.5
    const thin = cells([
      [0, 4, 0.5],
      [0, 1, 0],
    ]);
    expect(ParticleTileVisibility.get(thin, camera(inside)).tiles.length).toBe(2);
  });

  it('gives the catalog a sites target radius around its center', () => {
    const target: ParticleTarget = {
      kind: 'sites',
      sites: { count: 2, positions: Float32Array.of(-3, 0, 0, 3, 0, 0), radii: Float32Array.of(1, 2) },
    };
    expect(ParticleTileCatalog.getTargetRadius(target)).toBe(5);
    expect(ParticleTileCatalog.getTargetInfo(1, target).radius).toBe(5);
  });
});
