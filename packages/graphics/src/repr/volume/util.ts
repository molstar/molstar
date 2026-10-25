/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Grid, Volume } from '@molstar/model/model/volume';
import type { Loci } from '@molstar/model/model/loci';
import { equalEps } from '@molstar/core/math/linear-algebra/3d/common';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { packIntToRGBArray } from '@molstar/core/util/number-packing';
import { SetUtils } from '@molstar/core/util/set';
import { Box3D } from '@molstar/core/math/geometry';
import { toHalfFloat } from '@molstar/core/util/number-conversion';
import { clamp } from '@molstar/core/math/interpolate';
import { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import { Tensor } from '@molstar/core/math/linear-algebra/tensor';
import { Interval } from '@molstar/core/data/int/interval';
import type { SortedArray } from '@molstar/core/data/int/sorted-array';
import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import type { Image } from '@molstar/graphics/geo/geometry/image/image';
import { type BrickRange, BrickSize } from '@molstar/graphics/geo/geometry/direct-volume/direct-volume';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { GraphicsFeature } from '@molstar/graphics/gl/capabilities';
import { isComputeIsosurfaceSupported } from '@molstar/graphics/gl/compute/marching-cubes/compute';
import { getAtlasIsosurfaceExtraction } from '@molstar/graphics/gl/compute/marching-cubes/atlas';
import { missingCapabilities } from '../util.js';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3set = Vec3.set;
const v3normalize = Vec3.normalize;
const v3sub = Vec3.sub;
const v3addScalar = Vec3.addScalar;
const v3scale = Vec3.scale;
const v3toArray = Vec3.toArray;

export function eachVolumeLoci(
  loci: Loci,
  volume: Volume,
  props: { isoValue?: Volume.IsoValue; segments?: SortedArray } | undefined,
  apply: (interval: Interval) => boolean,
) {
  let changed = false;
  const cellCount = volume.grid.cells.data.length;

  if (Volume.isLoci(loci)) {
    if (!Volume.areEquivalent(loci.volume, volume)) return false;
    if (Interval.is(loci.instances)) {
      const start = Interval.start(loci.instances) * cellCount;
      const end = Interval.end(loci.instances) * cellCount;
      if (apply(Interval.ofBounds(start, end))) changed = true;
    } else {
      for (let i = 0, il = loci.instances.length; i < il; ++i) {
        const offset = loci.instances[i] * cellCount;
        if (apply(Interval.ofBounds(offset, offset + cellCount))) changed = true;
      }
    }
  } else if (Volume.Isosurface.isLoci(loci)) {
    if (!Volume.areEquivalent(loci.volume, volume)) return false;
    if (props?.isoValue) {
      if (!Volume.IsoValue.areSame(loci.isoValue, props.isoValue, volume.grid.stats)) return false;
      if (Interval.is(loci.instances)) {
        const start = Interval.start(loci.instances) * cellCount;
        const end = Interval.end(loci.instances) * cellCount;
        if (apply(Interval.ofBounds(start, end))) changed = true;
      } else {
        for (let i = 0, il = loci.instances.length; i < il; ++i) {
          const offset = loci.instances[i] * cellCount;
          if (apply(Interval.ofBounds(offset, offset + cellCount))) changed = true;
        }
      }
    } else {
      const {
        stats,
        cells: { data },
      } = volume.grid;
      const eps = stats.sigma;
      const v = Volume.IsoValue.toAbsolute(loci.isoValue, stats).absoluteValue;
      for (let i = 0, il = data.length; i < il; ++i) {
        if (equalEps(v, data[i], eps)) {
          OrderedSet.forEach(loci.instances, (j) => {
            const offset = j * cellCount;
            if (apply(Interval.ofSingleton(offset + i))) changed = true;
          });
        }
      }
    }
  } else if (Volume.Cell.isLoci(loci)) {
    if (!Volume.areEquivalent(loci.volume, volume)) return false;
    for (const { indices, instances } of loci.elements) {
      if (Interval.is(indices)) {
        OrderedSet.forEach(instances, (j) => {
          const offset = j * cellCount;
          if (apply(Interval.offset(indices, offset))) changed = true;
        });
      } else {
        OrderedSet.forEach(indices, (v) => {
          OrderedSet.forEach(instances, (j) => {
            const offset = j * cellCount;
            if (apply(Interval.ofSingleton(offset + v))) changed = true;
          });
        });
      }
    }
  } else if (Volume.Segment.isLoci(loci)) {
    if (!Volume.areEquivalent(loci.volume, volume)) return false;
    if (props?.segments) {
      for (const { segments, instances } of loci.elements) {
        if (OrderedSet.areIntersecting(segments, props.segments)) {
          OrderedSet.forEach(instances, (j) => {
            const offset = j * cellCount;
            if (apply(Interval.ofBounds(offset, offset + cellCount))) changed = true;
          });
        }
      }
    } else {
      const segmentation = Volume.Segmentation.get(volume);
      if (segmentation) {
        for (const { segments, instances } of loci.elements) {
          // this element's sets, and the instances it names (perf-audit-2 07/N3)
          const set = new Set<number>();
          for (let i = 0, _i = OrderedSet.size(segments); i < _i; i++) {
            const o = OrderedSet.getAt(segments, i);
            SetUtils.add(set, segmentation.segments.get(o)!);
          }
          const d = volume.grid.cells.data;
          for (let i = 0, il = d.length; i < il; ++i) {
            if (set.has(d[i])) {
              for (let j = 0, _j = OrderedSet.size(instances); j < _j; j++) {
                const offset = OrderedSet.getAt(instances, j) * cellCount;
                if (apply(Interval.ofSingleton(i + offset))) changed = true;
              }
            }
          }
        }
      }
    }
  }
  return changed;
}

export function createVolumeCellLocationIterator(volume: Volume, key: number): LocationIterator {
  const [xn, yn, zn] = volume.grid.cells.space.dimensions;
  const groupCount = xn * yn * zn;
  const instanceCount = volume.instances.length;
  const location = Volume.Cell.Location(volume);
  const getLocation = (groupIndex: number, instanceIndex: number) => {
    location.cell = groupIndex as Volume.CellIndex;
    location.instance = instanceIndex as Volume.InstanceIndex;
    return location;
  };
  return LocationIterator(groupCount, instanceCount, 1, getLocation);
}

export function createVolumeSliceCellLocationIterator(volume: Volume, key: number, image: Image): LocationIterator {
  const mapping = image.meta.mapping as any;
  if (!mapping) return createVolumeCellLocationIterator(volume, key);

  const groupCount = mapping.cell.length;
  const instanceCount = Volume.isPeriodic(volume) ? 1 : volume.instances.length;
  const location = Volume.Cell.Location(volume);
  const getLocation = (groupIndex: number, instanceIndex: number) => {
    location.cell = groupIndex as Volume.CellIndex;
    location.instance = instanceIndex as Volume.InstanceIndex;
    return location;
  };
  const nonInstanceable = Volume.isPeriodic(volume);
  return LocationIterator(groupCount, instanceCount, 1, getLocation, nonInstanceable);
}

//

export function getVolumeTexture2dLayout(dim: Vec3, padding = 0) {
  const area = dim[0] * dim[1] * dim[2];
  const squareDim = Math.sqrt(area);
  const powerOfTwoSize = Math.pow(2, Math.ceil(Math.log(squareDim) / Math.log(2)));

  let width = dim[0] + padding;
  let height = dim[1] + padding;
  let rows = 1;
  let columns = width;
  if (powerOfTwoSize < width * dim[2]) {
    columns = Math.floor(powerOfTwoSize / width);
    rows = Math.ceil(dim[2] / columns);
    width *= columns;
    height *= rows;
  } else {
    width *= dim[2];
  }
  return {
    width,
    height,
    columns,
    rows,
    powerOfTwoSize: height < powerOfTwoSize ? powerOfTwoSize : powerOfTwoSize * 2,
  };
}

export function createVolumeTexture2d(
  volume: Volume,
  variant: 'normals' | 'groups' | 'data',
  padding = 0,
  type: 'byte' | 'float' | 'halfFloat' = 'byte',
) {
  const {
    cells: { space, data },
    stats: { max, min },
  } = volume.grid;
  const dim = space.dimensions as Vec3;
  const { dataOffset: o } = space;
  const { width, height } = getVolumeTexture2dLayout(dim, padding);

  const itemSize = variant === 'data' ? 1 : 4;
  const array =
    type === 'byte'
      ? new Uint8Array(width * height * itemSize)
      : type === 'halfFloat'
        ? new Uint16Array(width * height * itemSize)
        : new Float32Array(width * height * itemSize);
  const textureImage = { array, width, height };

  const diff = max - min;
  const [xn, yn, zn] = dim;
  const xnp = xn + padding;
  const ynp = yn + padding;

  const n0 = Vec3();
  const n1 = Vec3();

  const xn1 = xn - 1;
  const yn1 = yn - 1;
  const zn1 = zn - 1;

  for (let z = 0; z < zn; ++z) {
    for (let y = 0; y < yn; ++y) {
      for (let x = 0; x < xn; ++x) {
        const column = Math.floor(((z * xnp) % width) / xnp);
        const row = Math.floor((z * xnp) / width);
        const px = column * xnp + x;
        const index = itemSize * (row * ynp * width + y * width + px);
        const offset = o(x, y, z);

        let value: number;
        if (type === 'byte') {
          value = Math.round(((data[offset] - min) / diff) * 255);
        } else if (type === 'halfFloat') {
          value = toHalfFloat((data[offset] - min) / diff);
        } else {
          value = (data[offset] - min) / diff;
        }

        if (variant === 'data') {
          array[index] = value;
        } else {
          if (variant === 'groups') {
            if (type === 'halfFloat') {
              let group = clamp(Math.round(offset), 0, 16777216 - 1) + 1;
              array[index + 2] = toHalfFloat(group % 256);
              group = Math.floor(group / 256);
              array[index + 1] = toHalfFloat(group % 256);
              group = Math.floor(group / 256);
              array[index] = toHalfFloat(group % 256);
            } else {
              packIntToRGBArray(offset, array, index);
            }
          } else {
            v3set(
              n0,
              data[o(Math.max(0, x - 1), y, z)],
              data[o(x, Math.max(0, y - 1), z)],
              data[o(x, y, Math.max(0, z - 1))],
            );
            v3set(
              n1,
              data[o(Math.min(xn1, x + 1), y, z)],
              data[o(x, Math.min(yn1, y + 1), z)],
              data[o(x, y, Math.min(zn1, z + 1))],
            );
            v3normalize(n0, v3sub(n0, n0, n1));
            v3addScalar(n0, v3scale(n0, n0, 0.5), 0.5);

            if (type === 'byte') {
              v3toArray(v3scale(n0, n0, 255), array, index);
            } else if (type === 'halfFloat') {
              array[index] = toHalfFloat(n0[0]);
              array[index + 1] = toHalfFloat(n0[1]);
              array[index + 2] = toHalfFloat(n0[2]);
            } else {
              v3toArray(n0, array, index);
            }
          }

          array[index + 3] = value;
        }
      }
    }
  }

  return textureImage;
}

/**
 * The values of a 3d grid and the strides to read them by (`data[x * sx + y * sy + z * sz]`): the data itself where the
 * space is linear (`Tensor.strides`), else a copy in x-fastest order read once through the space's accessors (spaces
 * that override them: periodic and wrapped volumes, which `dataOffset` maps modulo the data)
 */
function getStridedGrid(cells: Tensor): { data: ArrayLike<number>; sx: number; sy: number; sz: number } {
  const { space, data } = cells;
  const strides = Tensor.strides(space);
  if (strides) return { data, sx: strides[0], sy: strides[1], sz: strides[2] };
  const [nx, ny, nz] = space.dimensions;
  const copy = new Float32Array(nx * ny * nz);
  for (let z = 0, i = 0; z < nz; ++z) {
    for (let y = 0; y < ny; ++y) {
      for (let x = 0; x < nx; ++x, ++i) copy[i] = data[space.dataOffset(x, y, z)];
    }
  }
  return { data: copy, sx: 1, sy: nx, sz: nx * ny };
}

/**
 * The value range per brick of `BrickSize`³ cells (`BrickRange`): the grid's values normalized by its min and max (as in
 * `createVolumeTexture3d`), of a brick's cells and the next cell along each axis, times 255 and rounded outwards
 */
export function createVolumeBrickRange(volume: Volume): BrickRange {
  const {
    cells,
    stats: { max, min },
  } = volume.grid;
  const [nx, ny, nz] = cells.space.dimensions as Vec3;
  // offsets from strides, not `dataOffset` per cell (polymorphic, slow)
  const { data, sx, sy, sz } = getStridedGrid(cells);

  const bx = Math.ceil(nx / BrickSize);
  const by = Math.ceil(ny / BrickSize);
  const bz = Math.ceil(nz / BrickSize);
  const bxy = bx * by;
  const count = bxy * bz;
  // per brick the lowest and highest value, normalized, times 255; rounded once at the end
  const lows = new Float32Array(count).fill(Infinity);
  const highs = new Float32Array(count).fill(-Infinity);
  const scale = max > min ? 255 / (max - min) : 0;

  // per x its brick column, and whether it is on a brick's lower face (the next cell of the brick before)
  const brickX = new Int32Array(nx);
  const faceX = new Uint8Array(nx);
  for (let x = 0; x < nx; ++x) {
    brickX[x] = Math.floor(x / BrickSize);
    faceX[x] = x > 0 && x % BrickSize === 0 ? 1 : 0;
  }

  let v = 0;
  const add = (b: number) => {
    if (v < lows[b]) lows[b] = v;
    if (v > highs[b]) highs[b] = v;
  };

  for (let z = 0; z < nz; ++z) {
    const bz0 = Math.floor(z / BrickSize);
    const ez = z > 0 && z % BrickSize === 0;
    for (let y = 0; y < ny; ++y) {
      const by0 = Math.floor(y / BrickSize);
      const ey = y > 0 && y % BrickSize === 0;
      const row = bz0 * bxy + by0 * bx;
      let i = z * sz + y * sy;
      for (let x = 0; x < nx; ++x, i += sx) {
        v = (data[i] - min) * scale;
        const b = row + brickX[x];
        if (v < lows[b]) lows[b] = v;
        if (v > highs[b]) highs[b] = v;
        if (faceX[x] === 0 && !ey && !ez) continue;
        const ex = faceX[x] === 1;
        if (ex) add(b - 1);
        if (ey) {
          add(b - bx);
          if (ex) add(b - bx - 1);
        }
        if (ez) {
          add(b - bxy);
          if (ex) add(b - bxy - 1);
          if (ey) {
            add(b - bxy - bx);
            if (ex) add(b - bxy - bx - 1);
          }
        }
      }
    }
  }

  const array = new Uint8Array(count * 2);
  for (let b = 0; b < count; ++b) {
    // a little wider than the values: half-float data rounds
    array[2 * b] = clamp(Math.floor(lows[b] - 0.01), 0, 255);
    array[2 * b + 1] = clamp(Math.ceil(highs[b] + 0.01), 0, 255);
  }
  return { array, dim: Vec3.create(bx, by, bz) };
}

export function createVolumeTexture3d(volume: Volume, type: 'byte' | 'float' | 'halfFloat' = 'byte') {
  const {
    cells,
    stats: { max, min },
  } = volume.grid;
  const [width, height, depth] = cells.space.dimensions as Vec3;
  // offsets from strides, not `dataOffset` per cell (polymorphic, seven calls a cell)
  const { data, sx, sy, sz } = getStridedGrid(cells);

  const array =
    type === 'byte'
      ? new Uint8Array(width * height * depth * 4)
      : type === 'halfFloat'
        ? new Uint16Array(width * height * depth * 4)
        : new Float32Array(width * height * depth * 4);
  const textureVolume = { array, width, height, depth };
  const diff = max - min;
  const byte = type === 'byte';
  const half = type === 'halfFloat';

  let i = 0;
  for (let z = 0; z < depth; ++z) {
    // the neighbors' offsets, clamped at the grid's faces
    const zm = z > 0 ? -sz : 0;
    const zp = z < depth - 1 ? sz : 0;
    for (let y = 0; y < height; ++y) {
      const ym = y > 0 ? -sy : 0;
      const yp = y < height - 1 ? sy : 0;
      let offset = z * sz + y * sy;
      for (let x = 0; x < width; ++x, offset += sx, i += 4) {
        const xm = x > 0 ? -sx : 0;
        const xp = x < width - 1 ? sx : 0;

        // the normal: the negated central differences, normalized (zero stays zero), mapped to [0, 1]
        let nx = data[offset + xm] - data[offset + xp];
        let ny = data[offset + ym] - data[offset + yp];
        let nz = data[offset + zm] - data[offset + zp];
        const len = nx * nx + ny * ny + nz * nz;
        if (len > 0) {
          const inv = 1 / Math.sqrt(len);
          nx *= inv;
          ny *= inv;
          nz *= inv;
        }
        nx = nx * 0.5 + 0.5;
        ny = ny * 0.5 + 0.5;
        nz = nz * 0.5 + 0.5;
        const value = (data[offset] - min) / diff;

        if (byte) {
          array[i] = nx * 255;
          array[i + 1] = ny * 255;
          array[i + 2] = nz * 255;
          array[i + 3] = Math.round(value * 255);
        } else if (half) {
          array[i] = toHalfFloat(nx);
          array[i + 1] = toHalfFloat(ny);
          array[i + 2] = toHalfFloat(nz);
          array[i + 3] = toHalfFloat(value);
        } else {
          array[i] = nx;
          array[i + 1] = ny;
          array[i + 2] = nz;
          array[i + 3] = value;
        }
      }
    }
  }

  return textureVolume;
}

/**
 * The smoothed segment mask at grid cell (i, j, k): the cell twice and its 6 neighbours, 255 each when in the segment
 * `set`, over 8. Cells outside the grid count as not in the segment: a segment's expanded box reaches past the grid's
 * faces, where a read wrapped into the next row or gave undefined (perf-audit-2 07/N1).
 */
/**
 * The smoothed segment mask (0..255) of a voxel, from its label and its 6 neighbours'. With `bbox` (the voxels that
 * will be asked, `[min, max)`) the membership is computed once per voxel of the box and its 1-voxel border, not 7 times
 * per voxel (perf-audit-2 07/P8); the same values.
 */
export function getSegmentMaskValue(cells: Tensor, set: number[], bbox?: Box3D) {
  const data = cells.data;
  const o = cells.space.dataOffset;
  const [gx, gy, gz] = cells.space.dimensions;
  const mask = (x: number, y: number, z: number) =>
    x < 0 || y < 0 || z < 0 || x >= gx || y >= gy || z >= gz ? 0 : set.includes(data[o(x, y, z)]) ? 255 : 0;

  if (bbox) {
    const x0 = bbox.min[0] - 1,
      y0 = bbox.min[1] - 1,
      z0 = bbox.min[2] - 1;
    const mx = bbox.max[0] - bbox.min[0] + 2,
      my = bbox.max[1] - bbox.min[1] + 2,
      mz = bbox.max[2] - bbox.min[2] + 2;
    const values = new Uint8Array(mx * my * mz);
    for (let z = 0; z < mz; ++z) {
      for (let y = 0; y < my; ++y) {
        for (let x = 0, r = (z * my + y) * mx; x < mx; ++x) values[r + x] = mask(x + x0, y + y0, z + z0);
      }
    }
    const sy = mx,
      sz = mx * my;
    return (i: number, j: number, k: number) => {
      const lx = i - x0,
        ly = j - y0,
        lz = k - z0;
      // a voxel of the box: it and its neighbours are in `values`
      if (lx >= 1 && ly >= 1 && lz >= 1 && lx < mx - 1 && ly < my - 1 && lz < mz - 1) {
        const l = lz * sz + ly * sy + lx;
        const v0 = values[l];
        const n = values[l + 1] + values[l - 1] + values[l + sy] + values[l - sy] + values[l + sz] + values[l - sz];
        return Math.round((v0 + v0 + n) / 8);
      }
      return stencil(i, j, k);
    };
  }
  return stencil;

  function stencil(i: number, j: number, k: number) {
    const v0 = mask(i, j, k);
    const n =
      mask(i + 1, j, k) +
      mask(i - 1, j, k) +
      mask(i, j + 1, k) +
      mask(i, j - 1, k) +
      mask(i, j, k + 1) +
      mask(i, j, k - 1);
    return Math.round((v0 + v0 + n) / 8);
  }
}

export function createSegmentTexture2d(volume: Volume, set: number[], bbox: Box3D, padding = 0) {
  const dim = Box3D.size(Vec3(), bbox);
  const { width, height } = getVolumeTexture2dLayout(dim, padding);
  const value = getSegmentMaskValue(volume.grid.cells, set, bbox);

  const itemSize = 1;
  const array = new Uint8Array(width * height * itemSize);
  const textureImage = { array, width, height };

  const [xn, yn, zn] = dim;

  const xnp = xn + padding;
  const ynp = yn + padding;

  const [minx, miny, minz] = bbox.min;

  for (let z = 0; z < zn; ++z) {
    for (let y = 0; y < yn; ++y) {
      for (let x = 0; x < xn; ++x) {
        const column = Math.floor(((z * xnp) % width) / xnp);
        const row = Math.floor((z * xnp) / width);
        const px = column * xnp + x;
        const index = itemSize * (row * ynp * width + y * width + px);

        array[index] = value(x + minx, y + miny, z + minz);
      }
    }
  }

  return textureImage;
}

/**
 * Create a new tensor that is wrapped by one cell in all dimensions.
 * Reuses the original tensor data with new data accessors.
 * Only intended for isosurface construction.
 */
export function createWrappedTensor(tensor: Tensor): Tensor {
  const { space, data } = tensor;
  const { get, set, add, dataOffset } = space;
  const [xn, yn, zn] = space.dimensions as Vec3;

  const _dimensions = Vec3.create(xn + 1, yn + 1, zn + 1);

  const _get = (data: Tensor.Data, x: number, y: number, z: number) => get(data, x % xn, y % yn, z % zn);
  const _set = (data: Tensor.Data, x: number, y: number, z: number, d: number) => set(data, x % xn, y % yn, z % zn, d);
  const _add = (data: Tensor.Data, x: number, y: number, z: number, d: number) => add(data, x % xn, y % yn, z % zn, d);
  const _dataOffset = (x: number, y: number, z: number) => dataOffset(x % xn, y % yn, z % zn);

  const _space: Tensor.Space = {
    ...space,
    dimensions: _dimensions,
    get: _get,
    set: _set,
    add: _add,
    dataOffset: _dataOffset,
  };

  return {
    data,
    space: _space,
  };
}

/**
 * Create a new volume that is wrapped by one cell in all dimensions.
 * Reuses the original volume grid data with new data accessors.
 * Only intended for isosurface construction.
 */
export function createWrappedVolume(volume: Volume): Volume {
  const { grid } = volume;
  const matrix = Grid.getGridToCartesianTransform(grid);

  const _grid: Grid = {
    ...grid,
    transform: { kind: 'matrix', matrix },
    cells: createWrappedTensor(grid.cells),
  };

  return {
    ...volume,
    grid: _grid,
    parent: volume,
  };
}

export function createPeriodicVolume(volume: Volume, min: Vec3, max: Vec3): Volume {
  const { grid } = volume;
  const { space } = grid.cells;
  const { get, set, add, dataOffset } = space;
  const [xn, yn, zn] = space.dimensions as Vec3;

  const r = Vec3.sub(Vec3(), max, min);
  const _dimensions = Vec3.create(xn * r[0], yn * r[1], zn * r[2]);

  const _get = (data: Tensor.Data, x: number, y: number, z: number) => get(data, x % xn, y % yn, z % zn);
  const _set = (data: Tensor.Data, x: number, y: number, z: number, d: number) => set(data, x % xn, y % yn, z % zn, d);
  const _add = (data: Tensor.Data, x: number, y: number, z: number, d: number) => add(data, x % xn, y % yn, z % zn, d);
  const _dataOffset = (x: number, y: number, z: number) => dataOffset(x % xn, y % yn, z % zn);

  const _space: Tensor.Space = {
    ...space,
    dimensions: _dimensions,
    get: _get,
    set: _set,
    add: _add,
    dataOffset: _dataOffset,
  };

  const matrix = Grid.getGridToCartesianTransform(volume.grid);

  const t = Mat4.getScaling(Vec3(), matrix);
  Vec3.mul(t, t, space.dimensions as Vec3);
  const s = Vec3.mul(Vec3(), min, t);
  Mat4.setTranslation(matrix, s);

  const _transform: Grid.Transform = { kind: 'matrix', matrix };

  const _grid: Grid = {
    ...grid,
    transform: _transform,
    cells: {
      ...grid.cells,
      space: _space,
    },
  };

  return {
    ...volume,
    grid: _grid,
    parent: volume,
    instances: [{ transform: Mat4.identity() }],
  };
}

/**
 * Why the GPU marching cubes of a volume (isosurface, segment) are not used, if they are not: with compute kernels
 * (WebGPU) the values in one storage buffer, otherwise (WebGL) `capabilities` and the 2d atlas of the volume (with
 * `padding`)
 */
export function getVolumeGpuFallbackReason(
  volume: Volume,
  graphics: GraphicsContext,
  padding: number,
  capabilities: readonly GraphicsFeature[],
) {
  const cellCount = volume.grid.cells.data.length;
  // small volumes are about as fast or faster on CPU vs integrated GPU
  if (cellCount < Math.pow(10, 3)) return `the volume has ${cellCount} cells, fewer than 1000 (as fast on the CPU)`;
  if (isComputeIsosurfaceSupported(graphics)) {
    const max = graphics.compute.maxStorageBufferSize;
    if (cellCount * 4 > max) return `the volume values need ${cellCount * 4} bytes, a storage buffer at most ${max}`;
    return;
  }
  if (!getAtlasIsosurfaceExtraction(graphics)) return 'no compute kernels or fragment-shader marching cubes';
  const missing = missingCapabilities(graphics, capabilities);
  if (missing.length) return `missing ${missing.join(', ')}`;
  // the GPU is much more memory contraint, especially true for integrated GPUs,
  // fallback to CPU for large volumes
  const { powerOfTwoSize } = getVolumeTexture2dLayout(volume.grid.cells.space.dimensions as Vec3, padding);
  const max = graphics.resources.limits.maxTextureSize / 2;
  if (powerOfTwoSize > max) return `the 2d texture of the volume needs ${powerOfTwoSize}², at most ${max}²`;
}
