/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Box3D } from '@molstar/core/math/geometry/primitives/box3d';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { type ParticleList, type ParticleTarget, ParticleValues } from '../particle-list.js';
import { ParticleDensityGrid, type ParticleTileMoleculeSink, ParticleTileOpaque } from '../particle-tile-opaque.js';
import { ParticleTileVisibility } from '../particle-tile-visibility.js';
import { ParticleTileCatalog, ParticleTiles } from '../particle-tiles.js';

const boxOf = (b: Box3D | undefined) => (b ? [...b.min, ...b.max].map((v) => Math.round(v * 1000) / 1000) : b);
const cube = (size: number) => Box3D.create(Vec3.create(0, 0, 0), Vec3.create(size, size, size));

/**
 * molecules of radius 10 every 20 Å (a voxel's center each) in a cube of `size` where `keep`, each a cross-section of
 * `crossSection` (π 10², about 314, by default): an optical depth of crossSection / 8000 an Å
 */
function molecules(size: number, keep: (x: number, y: number, z: number) => boolean, crossSection = Math.PI * 100) {
  return (sink: ParticleTileMoleculeSink) => {
    const m = size / 20;
    for (let k = 0; k < m; ++k) {
      for (let j = 0; j < m; ++j) {
        for (let i = 0; i < m; ++i) {
          const x = i * 20 + 10,
            y = j * 20 + 10,
            z = k * 20 + 10;
          if (keep(x, y, z)) sink(x, y, z, 10, crossSection);
        }
      }
    }
  };
}

describe('particle tile opacity', () => {
  const blocks = (o: ReturnType<typeof ParticleTileOpaque.computeOpacity>) =>
    o ? [...o.density].map((d) => Math.round(d * 1e5) / 1e5) : o;
  // 314 Å² every 8000 Å³
  const uniform = Math.round((Math.PI * 100 * 1e5) / 8000) / 1e5;

  it('gives a block the optical depth a length of its content: 0.0393 an Å for molecules of 314 Å² every 8000 Å³', () => {
    const o = ParticleTileOpaque.computeOpacity(
      cube(160),
      molecules(160, () => true),
      { blocks: 2 },
    )!;
    expect(o.dimensions).toEqual([2, 2, 2]);
    expect(blocks(o)).toEqual(new Array(8).fill(uniform));
    // 4 blocks an axis by default
    expect(
      ParticleTileOpaque.computeOpacity(
        cube(160),
        molecules(160, () => true),
      )!.dimensions,
    ).toEqual([4, 4, 4]);
    expect(ParticleTileOpaque.computeOpacity(cube(160), () => {})).toBeUndefined();
  });

  it("counts a molecule in its center's block: a half-filled cell's empty half none, the far face the last block", () => {
    // the cell's x < 80 filled: blocks 0, 2, 4, 6 (x first) full, the others empty
    const half = ParticleTileOpaque.computeOpacity(
      cube(160),
      molecules(160, (x) => x < 80),
      { blocks: 2 },
    )!;
    expect(blocks(half)).toEqual([uniform, 0, uniform, 0, uniform, 0, uniform, 0]);
    // a center on the far face: the last block (80³: 1000 / 512,000); outside the cell nothing
    const face = ParticleTileOpaque.computeOpacity(cube(160), (sink) => sink(160, 160, 160, 1, 1000), { blocks: 2 })!;
    expect(blocks(face)).toEqual([0, 0, 0, 0, 0, 0, 0, Math.round((1000 / 512000) * 1e5) / 1e5]);
    expect(ParticleTileOpaque.computeOpacity(cube(160), (sink) => sink(161, 0, 0, 1, 1000))).toBeUndefined();
  });

  it("gives a tile's molecules: a sphere of its target's volume at each particle", () => {
    const list = {
      count: 3,
      targets: Int32Array.of(7, 7, 8),
      coordinates: Float32Array.of(100, 0, 0, 0, 50, 0, 1, 2, 3),
    } as unknown as ParticleList;
    const out: number[][] = [];
    const sink: ParticleTileMoleculeSink = (x, y, z, r, c) =>
      out.push([x, y, z, r, c].map((v) => Math.round(v * 1e4) / 1e4));
    ParticleTileOpaque.forEachMolecule(list, (id) => (id === 7 ? (4 / 3) * Math.PI * 1000 : undefined), sink);
    const c = Math.round(Math.PI * 100 * 1e4) / 1e4;
    expect(out).toEqual([
      [100, 0, 0, 10, c],
      [0, 50, 0, 10, c],
    ]);
  });

  it("estimates a target's volume by its mass, else its spheres", () => {
    const carbons = {
      kind: 'sites',
      sites: {
        count: 2,
        positions: new Float32Array(6),
        radii: Float32Array.of(1.7, 1.7),
        elements: Uint8Array.of(6, 6),
      },
    } as unknown as ParticleTarget;
    // two carbons and a hydrogen each, 1.21 Å³ a dalton
    expect(ParticleTileCatalog.getTargetVolume(carbons)).toBeCloseTo(2 * (12.011 + 1.008) * 1.21, 6);
    expect(ParticleTileCatalog.getTargetInfo(1, carbons).volume).toBeCloseTo(31.506, 3);
    const beads: ParticleTarget = {
      kind: 'sites',
      sites: { count: 1, positions: new Float32Array(3), radii: Float32Array.of(5) },
    };
    expect(ParticleTileCatalog.getTargetVolume(beads)).toBeCloseTo((4 / 3) * Math.PI * 125 * 0.7, 6);
  });

  it("computes a tile set's opacities, which hide the tiles past their neighbours", async () => {
    // a lattice of 8³ particles 10 Å apart, each of a target of 600 carbons (about 9,450 Å³: a sphere of 13 Å, a
    // cross-section of 540 Å² every 1000 Å³: 0.54 an Å)
    const count = 512;
    const coordinates = new Float32Array(count * 3);
    for (let i = 0; i < count; ++i) {
      coordinates.set([(i % 8) * 10 + 5, (Math.floor(i / 8) % 8) * 10 + 5, Math.floor(i / 64) * 10 + 5], i * 3);
    }
    const sites = {
      count: 600,
      positions: new Float32Array(1800),
      radii: new Float32Array(600).fill(9),
      elements: new Uint8Array(600).fill(6),
    };
    const list = {
      label: 'lattice',
      count,
      targets: Int32Array.from({ length: count }, () => 1),
      targetInfo: new Map([[1, { name: 'one' }]]),
      targetMapping: new Map([[1, { kind: 'sites', sites }]]),
      coordinates,
    } as unknown as ParticleList;
    // four layers along z, 17.5 Å each, 2 × 2 × 1 particles a block
    const set = await ParticleTileOpaque.computeOpacities(ParticleTiles.fromList(list, [2, 2, 4]), {
      blocks: 2,
    }).run();
    expect(set.tiles.size).toBe(16);
    // every block some (the molecules spread over the blocks they cover)
    for (const tile of set.tiles.values()) expect(Math.min(...tile.opacity!.density)).toBeGreaterThan(0);

    // seen from +z close up, the view within the lattice's outline (its tiles' reach past it uncovered): the front layer
    // (deep enough alone: 0.54 an Å over 17.5 Å) and the next (its reach, 9 Å past its box, within the front's cells:
    // what is drawn in those may be behind it), the two past them hidden
    const eye = Vec3.create(40, 40, 200);
    const view = Mat4.lookAt(Mat4(), eye, Vec3.create(40, 40, 0), Vec3.create(0, 1, 0));
    const projection = Mat4.perspective(Mat4(), -0.1, 0.1, 0.1, -0.1, 1, 1e5);
    const v = ParticleTileVisibility.get(set, {
      position: eye,
      direction: Vec3.create(0, 0, -1),
      viewProjection: Mat4.mul(Mat4(), projection, view),
      pixelScale: 256 / 0.1,
      near: 1,
    });
    expect(v.tiles.length).toBe(8);
    expect(v.culled.occluded).toBe(8);
  });

  it("puts the shown tiles' opacities side by side in one density grid over the set's grid", () => {
    // a set of 2 x 1 x 1 cells, 2 blocks a tile: the first tile's blocks 1..8, the second's 11..18
    const set = {
      grid: { origin: Vec3.create(0, 0, 0), cellSize: Vec3.create(10, 20, 30), dimensions: [2, 1, 1] as const },
      tiles: new Map([
        [
          0,
          {
            id: 0,
            opacity: { dimensions: [2, 2, 2] as const, density: Float32Array.from({ length: 8 }, (_, i) => 1 + i) },
          },
        ],
        [
          1,
          {
            id: 1,
            opacity: { dimensions: [2, 2, 2] as const, density: Float32Array.from({ length: 8 }, (_, i) => 11 + i) },
          },
        ],
      ]),
    } as unknown as Parameters<typeof ParticleTileOpaque.getDensityGrid>[0];
    const grid = ParticleTileOpaque.getDensityGrid(set, [0, 1])!;
    expect(grid.dimensions).toEqual([4, 2, 2]);
    expect(Array.from(grid.min)).toEqual([0, 0, 0]);
    expect(Array.from(grid.max)).toEqual([20, 20, 30]);
    // x fastest: a row the first tile's two blocks then the second's
    const at = (x: number, y: number, z: number) => grid.density[x + y * 4 + z * 8];
    expect([at(0, 0, 0), at(1, 0, 0), at(2, 0, 0), at(3, 0, 0)]).toEqual([1, 2, 11, 12]);
    expect([at(0, 1, 1), at(3, 1, 1)]).toEqual([7, 18]);
    // the tiles shown alone: the others' cells none
    const first = ParticleTileOpaque.getDensityGrid(set, [0])!;
    expect(first.density[2]).toBe(0);
    expect(first.density[0]).toBe(1);
    // none with an opacity: none
    expect(ParticleTileOpaque.getDensityGrid(set, [])).toBeUndefined();

    // updated in place from the last grid (the tiles let go cleared, those added written): as made again, a new object
    const one = ParticleTileOpaque.updateDensityGrid(set, [0])!;
    expect(Array.from(one.density)).toEqual(Array.from(first.density));
    const both = ParticleTileOpaque.updateDensityGrid(set, [0, 1], { grid: one, tiles: [0] })!;
    expect(both).not.toBe(one);
    expect(both.density).toBe(one.density);
    expect(Array.from(both.density)).toEqual(Array.from(grid.density));
    // each update a version of its own, from the last's, with the cells it wrote: the second tile added (x 2-4)
    expect(both.update).toEqual({ from: one.version, boxes: Int32Array.of(2, 0, 0, 4, 2, 2) });
    expect(both.version).toBeGreaterThan(one.version!);
    // the tiles with opacities it holds, kept by the changes
    expect([one.tileCount, both.tileCount]).toEqual([1, 2]);
    const second = ParticleTileOpaque.updateDensityGrid(set, [1], { grid: both, tiles: [0, 1] })!;
    expect(Array.from(second.density)).toEqual(Array.from(ParticleTileOpaque.getDensityGrid(set, [1])!.density));
    // the first let go: its cells cleared (x 0-2)
    expect(second.update).toEqual({ from: both.version, boxes: Int32Array.of(0, 0, 0, 2, 2, 2) });
    expect(second.tileCount).toBe(1);
    expect(ParticleTileOpaque.updateDensityGrid(set, [], { grid: second, tiles: [1] })).toBeUndefined();
  });

  it("gives a box near-cubic blocks: the longest side's count, the others' of about that size", () => {
    const box = Box3D.create(Vec3.create(0, 0, 0), Vec3.create(100, 50, 26));
    expect(ParticleTileOpaque.getGridDimensions(box, 4)).toEqual([4, 2, 1]);
    expect(ParticleTileOpaque.getGridDimensions(box)).toEqual([64, 32, 17]);
    // flat: one block across
    expect(ParticleTileOpaque.getGridDimensions(Box3D.create(Vec3.create(0, 0, 0), Vec3.create(10, 10, 0)), 8)).toEqual(
      [8, 8, 1],
    );
    // a point: one block
    expect(ParticleTileOpaque.getGridDimensions(cube(0))).toEqual([1, 1, 1]);
  });

  it("computes a list's density grid over its particles' box, its targets' volumes as the catalog's", () => {
    // 8 particles at a cube's corners, a target of a known volume; one more of a target without one
    const corners = [0, 1].flatMap((z) => [0, 1].flatMap((y) => [0, 1].map((x) => [x * 40, y * 40, z * 40])));
    const coordinates = Float32Array.from([...corners.flat(), 20, 20, 20]);
    const targets = Int32Array.from([...corners.map(() => 1), 2]);
    const target: ParticleTarget = { kind: 'shape', shape: undefined as never, center: Vec3() } as never;
    const list = {
      count: 9,
      coordinates,
      targets,
      targetInfo: new Map([
        [1, { name: 'one' }],
        [2, { name: 'two' }],
      ]),
      targetMapping: new Map([
        [1, target],
        [2, target],
      ]),
      customProperties: { add: () => {} },
      _propertyData: Object.create(null),
    } as unknown as ParticleList;
    const volume = (id: number) => (id === 1 ? 1000 : undefined);
    const grid = ParticleDensityGrid.compute(list, { longest: 2, volume })!;
    expect(grid.dimensions).toEqual([2, 2, 2]);
    expect(Array.from(grid.min)).toEqual([0, 0, 0]);
    expect(Array.from(grid.max)).toEqual([40, 40, 40]);
    // a corner's molecule (radius 6.2 Å) in its corner's block alone, its cross-section over the block's 8000 Å³
    const r = Math.cbrt((3 * 1000) / (4 * Math.PI));
    for (const d of grid.density) expect(d).toBeCloseTo((Math.PI * r * r) / 8000, 6);
    // known: kept, not computed again
    ParticleDensityGrid.set(list, grid);
    expect(ParticleDensityGrid.get(list)).toBe(grid);
    // none of a volume: none
    expect(ParticleDensityGrid.compute(list, { volume: () => undefined })).toBeUndefined();
  });

  it('computes a list grid as a task in chunks, the same as at once', async () => {
    // past a chunk (2^20 particles): a particle every 8 Å of a line, folded into a cube
    const n = (1 << 20) + 1000;
    const coordinates = new Float32Array(n * 3);
    for (let i = 0; i < n; ++i) coordinates.set([(i * 7) % 401, (i * 13) % 397, (i * 29) % 389], i * 3);
    const list = {
      count: n,
      coordinates,
      targets: ParticleValues.ofRuns(Int32Array.from([0, n]), Int32Array.from([1])),
      targetInfo: new Map([[1, { name: 'one' }]]),
      customProperties: { add: () => {} },
      _propertyData: Object.create(null),
    } as unknown as ParticleList;
    const volume = () => 500;
    const once = ParticleDensityGrid.compute(list, { longest: 8, volume })!;
    const chunked = (await ParticleDensityGrid.computeTask(list, { longest: 8, volume }).run())!;
    expect(chunked.dimensions).toEqual(once.dimensions);
    expect(Array.from(chunked.min)).toEqual(Array.from(once.min));
    expect(Array.from(chunked.density)).toEqual(Array.from(once.density));
  });
});
