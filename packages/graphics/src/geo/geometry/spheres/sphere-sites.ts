/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ValueCell } from '@molstar/core/util';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { InstanceTransforms } from '@molstar/core/math/geometry/instance-transforms';
import { type DataArray, calculateInvariantBoundingSphere } from '@molstar/graphics/gl/renderable/util';
import type { GroupMapping } from '../../util.js';
import { int32ToFloatBits } from '../center-group.js';
import { LocationIterator, PositionLocation } from '../../util/location-iterator.js';
import { type TransformData, getInstanceRows, getInstanceTransforms } from '../transform-data.js';

/**
 * Quantized sphere centers: sphere `i`'s center `origin + q * scale`, `q` int16 x, y, z of `data` (four int16 a sphere,
 * the fourth not read: a particle pack's site type), its group `i`. 8 bytes a sphere, drawn from the array as it is (a
 * view of a particle pack's sites, no copy; the shaders' `tSite`, `dSites`); WebGL2 and WebGPU
 */
export interface SphereSites {
  readonly data: Int16Array;
  readonly origin: Vec3;
  readonly scale: number;
  /** the spheres every `stride`-th of a whole's (a coarse read): their levels of detail from that stride (`getLodLevelsValue`'s `base`); 1 without */
  readonly stride?: number;
}

export namespace SphereSites {
  /** Sphere `i`'s center into `out` */
  export function getCenter(sites: SphereSites, i: number, out: Vec3): Vec3 {
    const { data, origin, scale } = sites;
    const o = i * 4;
    out[0] = origin[0] + data[o] * scale;
    out[1] = origin[1] + data[o + 1] * scale;
    out[2] = origin[2] + data[o + 2] * scale;
    return out;
  }

  /** The `count` centers, xyz each (a transient copy, where a reader takes them whole) */
  export function getCenters(sites: SphereSites, count: number) {
    const out = new Float32Array(count * 3);
    const p = Vec3();
    for (let i = 0; i < count; ++i) Vec3.toArray(getCenter(sites, i, p), out, i * 3);
    return out;
  }
}

/** The data array of `sites` (`tSite`): the view as it is, `count` elements of four int16 */
export function getSitesDataArray(sites: SphereSites, count: number): DataArray<Int16Array> {
  return { array: sites.data, count };
}

/**
 * The bounding sphere (from the decoded centers, where not set) and the group mapping (sphere `i` of group `i`: offsets
 * `i * 4` as `positionGroup`'s) of the sites in `cell`, as `CenterGroupData`'s
 */
export function createSitesData(cell: ValueCell<SphereSites>, getCount: () => number) {
  const boundingSphere = Sphere3D();
  let groupMapping: GroupMapping | undefined;
  let mappingKey = '';
  let sphereKey = '';
  const key = () => `${cell.ref.version}:${getCount()}`;

  return {
    getBoundingSphere() {
      const k = key();
      if (k !== sphereKey) {
        const count = getCount();
        Sphere3D.copy(
          boundingSphere,
          calculateInvariantBoundingSphere(SphereSites.getCenters(cell.ref.value, count), count, 1),
        );
        sphereKey = k;
      }
      return boundingSphere;
    },
    getGroupMapping() {
      const k = key();
      if (!groupMapping || k !== mappingKey) {
        const count = getCount();
        const indices = new Int32Array(count);
        const offsets = new Int32Array(count + 1);
        for (let i = 0; i < count; ++i) {
          indices[i] = i * 4;
          offsets[i + 1] = i + 1;
        }
        groupMapping = { indices, offsets };
        mappingKey = k;
      }
      return groupMapping;
    },
    setBoundingSphere(sphere: Sphere3D) {
      Sphere3D.copy(boundingSphere, sphere);
      sphereKey = key();
    },
    hasBoundingSphere() {
      return sphereKey === key();
    },
  };
}

/** The centers of the `count` spheres of `sites` per instance of `transform` (the sphere as the iterator's group) */
export function createSitesIterator(sites: SphereSites, count: number, transform: TransformData): LocationIterator {
  const instanceCount = transform.instanceCount.ref.value;
  const location = PositionLocation();
  const p = location.position;
  const t = getInstanceTransforms(transform);
  // the instance's row: the grid stores them in its order (`getInstanceRows`)
  const rows = getInstanceRows(transform);
  const local = Vec3();
  const getLocation = (groupIndex: number, instanceIndex: number) => {
    if (instanceIndex < 0) {
      SphereSites.getCenter(sites, groupIndex, p);
    } else {
      const row = rows ? rows[instanceIndex] : instanceIndex;
      InstanceTransforms.getPoint(t, row, SphereSites.getCenter(sites, groupIndex, local), p);
    }
    return location;
  };
  return LocationIterator(count, instanceCount, 1, getLocation);
}

/**
 * The centers and groups of the `count` spheres of `data` (four int16 a sphere) by `origin` and `scale`, as
 * `positionGroup` holds them (`PositionGroupStride` words, the group the sphere as its bits): a transient copy, for
 * readers of that layout (loci bounds, export)
 */
export function decodeSitesPositionGroup(data: Int16Array, count: number, origin: Vec3, scale: number) {
  const out = new Float32Array(count * 4);
  for (let i = 0; i < count; ++i) {
    const o = i * 4;
    out[o] = origin[0] + data[o] * scale;
    out[o + 1] = origin[1] + data[o + 1] * scale;
    out[o + 2] = origin[2] + data[o + 2] * scale;
    out[o + 3] = int32ToFloatBits(i);
  }
  return out;
}
