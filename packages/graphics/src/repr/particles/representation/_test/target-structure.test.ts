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
import { ParticleSites } from '@molstar/model/model/particles/particle-list';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Sphere3D } from '@molstar/core/math/geometry';
import { getFastBoundary } from '@molstar/core/math/geometry/boundary';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { OrderedSet } from '@molstar/core/data/int';
import type { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import type { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { ElementIterator } from '../../../structure/visual/util/element.js';
import { DefaultBlobDensityProps } from '../../../structure/visual/util/blob-surface.js';
import { createStructureTargetGeometry, createStructureTargetSizes } from '../target/structure.js';
import { createBlobSurfaceTargetGeometry } from '../target/blob-surface.js';
import {
  createSitesBlobSurfaceGeometry,
  createSitesTargetGeometry,
  createSitesTargetSizes,
  getSitesSignatureDefines,
  setSitesTargetSiteColor,
} from '../target/sites.js';
import { createSiteColorValues } from '@molstar/graphics/gl/renderable/spheres';
import { getStructureSignatureDefines, setStructureTargetSiteColor } from '../target/site-color.js';
import { ElementSymbolColors, elementSymbolColor } from '@molstar/graphics/theme/color/element-symbol';
import { Color } from '@molstar/core/util/color';

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

describe('structure targets colored near the camera', () => {
  const on = (carbonColor: 'particle' | 'element') => ({
    siteColor: { name: 'element' as const, params: { carbonColor, radius: [3, 9] as [number, number] } },
  });

  it('by element, a texel an element in the serial order, a carbon alpha 0 where it keeps the particle color', async () => {
    const structure = await load1cbs();
    const values = createSiteColorValues();
    setStructureTargetSiteColor(structure, on('particle'), values);
    expect(values.dSiteColor.ref.value).toBe('group');
    expect(values.tSiteColor.ref.value.count).toBe(structure.elementCount);
    expect(Array.from(values.uSiteColorRange.ref.value)).toEqual([3, 9]);
    const array = values.tSiteColor.ref.value.array;
    let g = 0;
    let carbons = 0;
    for (const unit of structure.units) {
      const symbols = (unit as Unit.Atomic).model.atomicHierarchy.atoms.type_symbol;
      for (let i = 0; i < unit.elements.length; ++i, ++g) {
        const symbol = symbols.value(unit.elements[i]);
        const texel = Array.from(array.subarray(g * 4, g * 4 + 4));
        if (symbol === 'C') {
          carbons += 1;
          expect(texel[3]).toBe(0);
        } else {
          expect(texel).toEqual([...Color.toRgb(elementSymbolColor(ElementSymbolColors, symbol)), 255]);
        }
      }
    }
    expect(carbons).toBeGreaterThan(0);
    // carbons their element's, or off
    setStructureTargetSiteColor(structure, on('element'), values);
    expect(Array.from(values.tSiteColor.ref.value.array).filter((_, k) => k % 4 === 3 && _ === 0)).toEqual([]);
    setStructureTargetSiteColor(structure, { siteColor: { name: 'off', params: {} } }, values);
    expect(values.dSiteColor.ref.value).toBe('off');
    expect(getStructureSignatureDefines(on('particle'))).toEqual({ dSiteColor: ['group'] });
    expect(getStructureSignatureDefines({ siteColor: { name: 'off', params: {} } })).toEqual({ dSiteColor: ['off'] });
  });
});

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

  it("is bound as the fast boundaries of its chains' decoded positions are, from the quantized integers", () => {
    // two chains, a third empty, odd scales and offsets: the float32 positions' boxes, decoded from the extremes
    const n = 40;
    const data = new Int16Array(n * 4);
    for (let i = 0; i < n; ++i) data.set([((i * 37) % 201) - 100, ((i * 53) % 301) - 150, (i * 11) % 97, 0], i * 4);
    const quantized = { data, origin: Vec3.create(103.17, -7.3, 2.2), scale: 0.0137, typeRadii: Float32Array.of(1.5) };
    const chains = { offsets: Int32Array.of(0, 15, 15, n), labels: ['A', 'B', 'C'] };
    const sites = ParticleSites.ofQuantized(n, quantized, undefined, chains);
    const sphere = createSitesTargetGeometry(sites, theme(1), props(true), graphics).boundingSphere;
    // the reference: the fast boundaries (`getFastBoundary`) of the decoded float32 positions a chain, merged as the target does
    const positions = sites.positions;
    const helper = new BoundaryHelper('98');
    const spheres = [
      [0, 15],
      [15, n],
    ].map(([a, b]) => {
      const x = new Float32Array(n);
      const y = new Float32Array(n);
      const z = new Float32Array(n);
      for (let i = 0; i < n; ++i) [x[i], y[i], z[i]] = [positions[i * 3], positions[i * 3 + 1], positions[i * 3 + 2]];
      return getFastBoundary({ x, y, z, indices: OrderedSet.ofBounds(a, b) }).sphere;
    });
    helper.reset();
    for (const s of spheres) helper.includeSphere(s);
    helper.finishedIncludeStep();
    for (const s of spheres) helper.radiusSphere(s);
    const expected = Sphere3D.expand(Sphere3D(), helper.getSphere(), 1.5 * 1.2 + 0.05);
    expect(Array.from(sphere.center)).toEqual(Array.from(expected.center));
    expect(sphere.radius).toBe(expected.radius);
  });

  it('declares the defines its impostors take from its data, for the precompile', () => {
    const data = Int16Array.of(0, 0, 0, 0, 1, 0, 0, 0);
    const quantized = { data, origin: Vec3(), scale: 1, typeRadii: Float32Array.of(1.5) };
    const withElements = ParticleSites.ofQuantized(2, { ...quantized, typeElements: Uint8Array.of(8) });
    const without = ParticleSites.ofQuantized(2, quantized);
    const on = {
      siteColor: {
        name: 'element' as const,
        params: { carbonColor: 'particle' as const, radius: [4, 12] as [number, number] },
      },
    };
    const off = { siteColor: { name: 'off' as const, params: {} } };
    const webgpu = { capabilities: { backend: 'webgpu' } } as any;
    const webgl = { capabilities: { backend: 'webgl' } } as any;
    expect(getSitesSignatureDefines(withElements, on, webgpu)).toEqual({ dSites: [true], dSiteColor: ['siteType'] });
    // where quantized sites are not drawn as they are: by group
    expect(getSitesSignatureDefines(withElements, on, webgl)).toEqual({ dSites: [false], dSiteColor: ['group'] });
    expect(getSitesSignatureDefines(withElements, off, webgpu)).toEqual({ dSites: [true], dSiteColor: ['off'] });
    expect(getSitesSignatureDefines(without, on, webgpu)).toEqual({ dSites: [true], dSiteColor: ['off'] });
  });

  it("is bound alike from a pack's extremes of its chains, with no pass over the sites", () => {
    const n = 40;
    const data = new Int16Array(n * 4);
    for (let i = 0; i < n; ++i) data.set([((i * 37) % 201) - 100, ((i * 53) % 301) - 150, (i * 11) % 97, 0], i * 4);
    const quantized = { data, origin: Vec3.create(103.17, -7.3, 2.2), scale: 0.0137, typeRadii: Float32Array.of(1.5) };
    for (const offsets of [Int32Array.of(0, 15, 15, n), undefined]) {
      const chains = offsets ? { offsets, labels: ['A', 'B', 'C'] } : undefined;
      const runs = offsets ?? Int32Array.of(0, n);
      const extremes = new Int16Array((runs.length - 1) * 6);
      for (let c = 0; c + 1 < runs.length; ++c)
        ParticleSites.getQuantizedExtremes(data, runs[c], runs[c + 1], extremes, c * 6);
      const scanned = ParticleSites.ofQuantized(n, quantized, undefined, chains);
      const given = ParticleSites.ofQuantized(n, { ...quantized, extremes }, undefined, chains);
      const a = createSitesTargetGeometry(scanned, theme(1), props(true), graphics).boundingSphere;
      const b = createSitesTargetGeometry(given, theme(1), props(true), graphics).boundingSphere;
      expect(Array.from(b.center)).toEqual(Array.from(a.center));
      expect(b.radius).toBe(a.radius);
    }
  });

  it("are colored near the camera by element: a table by type (a carbon alpha 0, its particle's), or a texel a site", () => {
    // C, O, S: three types, the sites' types 0, 1, 2, 0
    const data = Int16Array.of(0, 0, 0, 0, 1, 0, 0, 1, 2, 0, 0, 2, 3, 0, 0, 0);
    const quantized = { data, origin: Vec3(), scale: 1, typeRadii: Float32Array.of(1.7, 1.52, 1.8) };
    const withElements = ParticleSites.ofQuantized(4, { ...quantized, typeElements: Uint8Array.of(6, 8, 16) });
    const rgb = (c: Color) => Color.toRgb(c);
    const range = {
      siteColor: {
        name: 'element' as const,
        params: { carbonColor: 'particle' as const, radius: [3, 9] as [number, number] },
      },
    };
    const byType = createSiteColorValues();
    setSitesTargetSiteColor(withElements, range, byType, true);
    expect(byType.dSiteColor.ref.value).toBe('siteType');
    expect(byType.tSiteColor.ref.value.count).toBe(3);
    const table = Array.from(byType.tSiteColor.ref.value.array.subarray(0, 12));
    expect(table).toEqual([
      ...rgb(ElementSymbolColors.C),
      0,
      ...rgb(ElementSymbolColors.O),
      255,
      ...rgb(ElementSymbolColors.S),
      255,
    ]);
    expect(Array.from(byType.uSiteColorRange.ref.value)).toEqual([3, 9]);
    // a site each (spheres not drawn from the types), carbons their element's
    const bySite = createSiteColorValues();
    const element = {
      siteColor: { ...range.siteColor, params: { ...range.siteColor.params, carbonColor: 'element' as const } },
    };
    setSitesTargetSiteColor(withElements, element, bySite, false);
    expect(bySite.dSiteColor.ref.value).toBe('group');
    expect(bySite.tSiteColor.ref.value.count).toBe(4);
    expect(Array.from(bySite.tSiteColor.ref.value.array.subarray(12, 16))).toEqual([
      ...rgb(ElementSymbolColors.C),
      255,
    ]);
    // off where asked, or where the elements are not known
    setSitesTargetSiteColor(withElements, { siteColor: { name: 'off', params: {} } }, bySite, false);
    expect(bySite.dSiteColor.ref.value).toBe('off');
    const unknown = createSiteColorValues();
    setSitesTargetSiteColor(ParticleSites.ofQuantized(4, quantized), range, unknown, true);
    expect(unknown.dSiteColor.ref.value).toBe('off');
  });

  it('drawn from quantized sites: the sites as they are, sized by a table of their types', () => {
    const data = Int16Array.of(0, 0, 0, 1, 4, 0, 0, 0, 0, 8, 0, 1);
    const sites = ParticleSites.ofQuantized(3, {
      data,
      origin: Vec3.create(1, 2, 3),
      scale: 0.5,
      typeRadii: Float32Array.of(1.5, 1.75),
    });
    const spheres = createSitesTargetGeometry(sites, theme(1), props(true), graphics) as Spheres;
    expect(spheres.sites!.data).toBe(data);
    const unpack = (a: ArrayLike<number>, i: number) => (a[i * 4] * 65536 + a[i * 4 + 1] * 256 + a[i * 4 + 2]) / 100;
    // the table: a size a type, the theme's scale applied
    const byType = createSitesTargetSizes(sites, theme(2), undefined, true);
    expect(byType.sizeType.ref.value).toBe('siteType');
    expect(byType.uSizeIndexType.ref.value).toBe(5);
    expect(byType.tSize.ref.value.count).toBe(2);
    // to a step of the packing (0.01)
    const close = (a: number[], b: number[]) =>
      a.forEach((v, i) => expect(Math.abs(v - b[i])).toBeLessThanOrEqual(0.011));
    close(
      [0, 1].map((i) => unpack(byType.tSize.ref.value.array, i)),
      [3, 3.5],
    );
    // a site each otherwise (a float copy of the centers)
    const bySite = createSitesTargetSizes(sites, theme(2));
    expect(bySite.sizeType.ref.value).toBe('group');
    close(
      [0, 1, 2].map((i) => unpack(bySite.tSize.ref.value.array, i)),
      [3.5, 3, 3.5],
    );
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
