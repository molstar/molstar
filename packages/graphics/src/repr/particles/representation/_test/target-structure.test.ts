/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, StructureElement, Unit } from '@molstar/model/model/structure';
import type { Location } from '@molstar/model/model/location';
import { createGroupSize } from '@molstar/graphics/geo/geometry/size-data';
import { getPhysicalRadius } from '@molstar/graphics/theme/size/physical';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import type { ParticleSites } from '@molstar/model/model/particles/particle-list';
import type { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import type { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { ElementIterator } from '../../../structure/visual/util/element.js';
import { DefaultBlobDensityProps } from '../../../structure/visual/util/blob-surface.js';
import { createStructureTargetGeometry, createStructureTargetSizes } from '../target/structure.js';
import { createBlobSurfaceTargetGeometry } from '../target/blob-surface.js';
import { createSitesBlobSurfaceGeometry, createSitesTargetGeometry, createSitesTargetSizes } from '../target/sites.js';

let _1cbs: Promise<Structure> | undefined;
/** the fixture parsed once for the file */
function load1cbs() {
  _1cbs ??= _load1cbs();
  return _1cbs;
}

async function _load1cbs() {
  const data = fs.readFileSync(
    path.resolve(import.meta.dirname, '../../../../../../../data/examples/1cbs_updated.cif'),
    'utf8',
  );
  const parsed = await CIF.parse(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

describe('createStructureTargetSizes', () => {
  it("writes each element's physical radius in the structure's serial order", async () => {
    const structure = await load1cbs();
    expect(structure.units.length).toBeGreaterThan(1);
    for (const scale of [1, 1.5]) {
      const theme = { color: {} as any, size: { granularity: 'group', props: { scale }, size: () => 1 } as any };
      const sizes = createStructureTargetSizes(structure, theme).tSize.ref.value.array;
      // as through a location per element, in the serial mapping's order
      const expected = createGroupSize(
        ElementIterator.fromStructure(structure),
        (l: Location) =>
          getPhysicalRadius((l as StructureElement.Location).unit, (l as StructureElement.Location).element) * scale,
      ).tSize.ref.value.array;
      expect(Array.from(sizes)).toEqual(Array.from(expected));
    }
  });
});

/** The structure's elements as sites: positions and physical radii in its serial order, a unit a chain */
function sitesOf(structure: Structure): ParticleSites {
  const count = structure.elementCount;
  const positions = new Float32Array(count * 3);
  const radii = new Float32Array(count);
  const offsets = new Int32Array(structure.units.length + 1);
  const v = [0, 0, 0] as any;
  let g = 0;
  for (let u = 0; u < structure.units.length; ++u) {
    const unit = structure.units[u];
    offsets[u + 1] = offsets[u] + unit.elements.length;
    for (let i = 0, il = unit.elements.length; i < il; ++i, ++g) {
      unit.conformation.position(unit.elements[i], v);
      positions.set(v, g * 3);
      radii[g] = getPhysicalRadius(unit, unit.elements[i]);
    }
  }
  return { count, positions, radii, chains: { offsets, labels: structure.units.map((u) => `${u.id}`) } };
}

/** The structure's units (all, or those given), with fast boundaries (their boxes' spheres): as a target structure */
function targetOf(structure: Structure, units = structure.units) {
  const builder = new Structure.StructureBuilder({ representativeModel: structure.representativeModel });
  for (const u of units) {
    builder.addUnit(
      u.kind,
      u.model,
      u.conformation.operator,
      u.elements,
      u.traits | Unit.Trait.FastBoundary,
      u.invariantId,
    );
  }
  return builder.getStructure();
}

describe('sites targets', () => {
  const ctx = { runtime: SyncRuntimeContext } as any;
  const theme = (scale: number) =>
    ({ color: {} as any, size: { granularity: 'group', props: { scale }, size: () => 1 } as any }) as any;
  const graphics = { capabilities: { fragDepth: true } } as any;
  const props = (tryUseImpostor: boolean) => ({ type: 'spacefill', sizeFactor: 1.2, detail: 1, tryUseImpostor });

  it("are sized as a structure's elements, to a step of the packed size", async () => {
    const structure = await load1cbs();
    const sites = sitesOf(structure);
    // the radii as float32: a size at a half of the packing's step (0.01) can round to its neighbor
    const unpack = (a: ArrayLike<number>, i: number) => a[i * 4] * 65536 + a[i * 4 + 1] * 256 + a[i * 4 + 2];
    for (const scale of [1, 1.5]) {
      const a = createSitesTargetSizes(sites, theme(scale)).tSize.ref.value.array;
      const b = createStructureTargetSizes(structure, theme(scale)).tSize.ref.value.array;
      expect(a.length).toBe(b.length);
      let maxStep = 0;
      for (let i = 0; i < sites.count; ++i) maxStep = Math.max(maxStep, Math.abs(unpack(a, i) - unpack(b, i)));
      expect(maxStep).toBeLessThanOrEqual(1);
    }
  });

  it("are the spheres of a structure's elements: impostor centers and groups", async () => {
    const structure = await load1cbs();
    const a = createSitesTargetGeometry(sitesOf(structure), theme(1), props(true), graphics) as Spheres;
    const b = createStructureTargetGeometry(ctx, structure, theme(1), props(true), graphics) as Spheres;
    expect(a.kind).toBe('spheres');
    expect(b.kind).toBe('spheres');
    expect(a.sphereCount).toBe(b.sphereCount);
    expect(Array.from(a.positionGroup.ref.value.array)).toEqual(Array.from(b.positionGroup.ref.value.array));
  });

  it('are bounded as a target structure, a chain a unit', async () => {
    const structure = await load1cbs();
    expect(structure.units.length).toBeGreaterThan(1);
    for (const target of [targetOf(structure, structure.units.slice(0, 1)), targetOf(structure)]) {
      const a = createSitesTargetGeometry(sitesOf(target), theme(1), props(true), graphics) as Spheres;
      const b = createStructureTargetGeometry(ctx, target, theme(1), props(true), graphics) as Spheres;
      for (let k = 0; k < 3; ++k) expect(a.boundingSphere.center[k]).toBeCloseTo(b.boundingSphere.center[k], 4);
      expect(a.boundingSphere.radius).toBeCloseTo(b.boundingSphere.radius, 4);
      expect(a.boundingSphere.extrema?.map((e) => Array.from(e))).toEqual(
        b.boundingSphere.extrema?.map((e) => Array.from(e)),
      );
    }
  });

  it("are the sphere mesh of a structure's elements", async () => {
    const structure = await load1cbs();
    const a = createSitesTargetGeometry(sitesOf(structure), theme(1.5), props(false)) as Mesh;
    const b = createStructureTargetGeometry(ctx, structure, theme(1.5), props(false)) as Mesh;
    expect(a.kind).toBe('mesh');
    expect(a.vertexCount).toBe(b.vertexCount);
    expect(a.triangleCount).toBe(b.triangleCount);
    const va = a.vertexBuffer.ref.value;
    const vb = b.vertexBuffer.ref.value;
    for (let i = 0; i < a.vertexCount * 3; ++i) expect(va[i]).toBeCloseTo(vb[i], 4);
    expect(Array.from(a.groupBuffer.ref.value.subarray(0, a.vertexCount))).toEqual(
      Array.from(b.groupBuffer.ref.value.subarray(0, b.vertexCount)),
    );
  });

  it('are the blob surface of a target structure', async () => {
    const chain = targetOf(await load1cbs());
    const { blobSize, blobMethod, resolution, adjustResolution, blobShape, radiusOffset, smoothness } =
      DefaultBlobDensityProps;
    const blob = { blobSize, blobMethod, resolution, adjustResolution, blobShape, radiusOffset, smoothness };
    const a = await createSitesBlobSurfaceGeometry(ctx, sitesOf(chain), theme(1), blob);
    const b = await createBlobSurfaceTargetGeometry(ctx, chain, theme(1), blob);
    expect(a.vertexCount).toBeGreaterThan(0);
    expect(a.vertexCount).toBe(b.vertexCount);
    expect(a.triangleCount).toBe(b.triangleCount);
    const va = a.vertexBuffer.ref.value;
    const vb = b.vertexBuffer.ref.value;
    for (let i = 0; i < a.vertexCount * 3; ++i) expect(va[i]).toBeCloseTo(vb[i], 3);
  });
});
