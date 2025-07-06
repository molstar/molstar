/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Volume } from '../../mol-model/volume';
import { Loci } from '../../mol-model/loci';
import { Interval, OrderedSet, SortedArray } from '../../mol-data/int';
import { equalEps } from '../../mol-math/linear-algebra/3d/common';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { packIntToRGBArray } from '../../mol-util/number-packing';
import { SetUtils } from '../../mol-util/set';
import { Box3D } from '../../mol-math/geometry';
import { toHalfFloat } from '../../mol-util/number-conversion';
import { clamp } from '../../mol-math/interpolate';
import { LocationIterator } from '../../mol-geo/util/location-iterator';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3set = Vec3.set;
const v3normalize = Vec3.normalize;
const v3sub = Vec3.sub;
const v3addScalar = Vec3.addScalar;
const v3scale = Vec3.scale;
const v3toArray = Vec3.toArray;

export function eachVolumeLoci(loci: Loci, volume: Volume, props: { isoValue?: Volume.IsoValue, segments?: SortedArray } | undefined, apply: (interval: Interval) => boolean) {
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
            const { stats, cells: { data } } = volume.grid;
            const eps = stats.sigma;
            const v = Volume.IsoValue.toAbsolute(loci.isoValue, stats).absoluteValue;
            for (let i = 0, il = data.length; i < il; ++i) {
                if (equalEps(v, data[i], eps)) {
                    OrderedSet.forEach(loci.instances, j => {
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
                OrderedSet.forEach(instances, j => {
                    const offset = j * cellCount;
                    if (apply(Interval.offset(indices, offset))) changed = true;
                });
            } else {
                OrderedSet.forEach(indices, v => {
                    OrderedSet.forEach(instances, j => {
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
                    OrderedSet.forEach(instances, j => {
                        const offset = j * cellCount;
                        if (apply(Interval.ofBounds(offset, offset + cellCount))) changed = true;
                    });
                }
            }
        } else {
            const segmentation = Volume.Segmentation.get(volume);
            if (segmentation) {
                const set = new Set<number>();
                for (const { segments, instances } of loci.elements) {
                    for (let i = 0, _i = OrderedSet.size(segments); i < _i; i++) {
                        const o = OrderedSet.getAt(segments, i);
                        SetUtils.add(set, segmentation.segments.get(o)!);
                    }
                    const s = Array.from(set.values());
                    const d = volume.grid.cells.data;
                    for (let i = 0, il = d.length; i < il; ++i) {
                        if (s.includes(d[i])) {
                            for (let j = 0, _j = OrderedSet.size(instances); j < _j; j++) {
                                const offset = j * cellCount;
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

export function createVolumeCellLocationIterator(volume: Volume): LocationIterator {
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
    return { width, height, columns, rows, powerOfTwoSize: height < powerOfTwoSize ? powerOfTwoSize : powerOfTwoSize * 2 };
}

export function createVolumeTexture2d(volume: Volume, variant: 'normals' | 'groups' | 'data', padding = 0, type: 'byte' | 'float' | 'halfFloat' = 'byte') {
    const { cells: { space, data }, stats: { max, min } } = volume.grid;
    const dim = space.dimensions as Vec3;
    const { dataOffset: o } = space;
    const { width, height } = getVolumeTexture2dLayout(dim, padding);

    const itemSize = variant === 'data' ? 1 : 4;
    const array = type === 'byte'
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
                const index = itemSize * ((row * ynp * width) + (y * width) + px);
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
                        v3set(n0,
                            data[o(Math.max(0, x - 1), y, z)],
                            data[o(x, Math.max(0, y - 1), z)],
                            data[o(x, y, Math.max(0, z - 1))]
                        );
                        v3set(n1,
                            data[o(Math.min(xn1, x + 1), y, z)],
                            data[o(x, Math.min(yn1, y + 1), z)],
                            data[o(x, y, Math.min(zn1, z + 1))]
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

export function createVolumeTexture3d(volume: Volume, type: 'byte' | 'float' | 'halfFloat' = 'byte') {
    const { cells: { space, data }, stats: { max, min } } = volume.grid;
    const [width, height, depth] = space.dimensions as Vec3;
    const { dataOffset: o } = space;

    const array = type === 'byte'
        ? new Uint8Array(width * height * depth * 4)
        : type === 'halfFloat'
            ? new Uint16Array(width * height * depth * 4)
            : new Float32Array(width * height * depth * 4);
    const textureVolume = { array, width, height, depth };
    const diff = max - min;

    const n0 = Vec3();
    const n1 = Vec3();

    const width1 = width - 1;
    const height1 = height - 1;
    const depth1 = depth - 1;

    let i = 0;
    for (let z = 0; z < depth; ++z) {
        for (let y = 0; y < height; ++y) {
            for (let x = 0; x < width; ++x) {
                const offset = o(x, y, z);

                v3set(n0,
                    data[o(Math.max(0, x - 1), y, z)],
                    data[o(x, Math.max(0, y - 1), z)],
                    data[o(x, y, Math.max(0, z - 1))]
                );
                v3set(n1,
                    data[o(Math.min(width1, x + 1), y, z)],
                    data[o(x, Math.min(height1, y + 1), z)],
                    data[o(x, y, Math.min(depth1, z + 1))]
                );
                v3normalize(n0, v3sub(n0, n0, n1));
                v3addScalar(n0, v3scale(n0, n0, 0.5), 0.5);

                if (type === 'byte') {
                    v3toArray(v3scale(n0, n0, 255), array, i);
                    array[i + 3] = Math.round(((data[offset] - min) / diff) * 255);
                } else if (type === 'halfFloat') {
                    array[i] = toHalfFloat(n0[0]);
                    array[i + 1] = toHalfFloat(n0[1]);
                    array[i + 2] = toHalfFloat(n0[2]);
                    array[i + 3] = toHalfFloat((data[offset] - min) / diff);
                } else {
                    v3toArray(n0, array, i);
                    array[i + 3] = (data[offset] - min) / diff;
                }
                i += 4;
            }
        }
    }

    return textureVolume;
}

export function createSegmentTexture2d(volume: Volume, set: number[], bbox: Box3D, padding = 0) {
    const data = volume.grid.cells.data;
    const dim = Box3D.size(Vec3(), bbox);
    const o = volume.grid.cells.space.dataOffset;
    const { width, height } = getVolumeTexture2dLayout(dim, padding);

    const itemSize = 1;
    const array = new Uint8Array(width * height * itemSize);
    const textureImage = { array, width, height };

    const [xn, yn, zn] = dim;
    const xn1 = xn - 1;
    const yn1 = yn - 1;
    const zn1 = zn - 1;

    const xnp = xn + padding;
    const ynp = yn + padding;

    const [minx, miny, minz] = bbox.min;
    const [maxx, maxy, maxz] = bbox.max;

    for (let z = 0; z < zn; ++z) {
        for (let y = 0; y < yn; ++y) {
            for (let x = 0; x < xn; ++x) {
                const column = Math.floor(((z * xnp) % width) / xnp);
                const row = Math.floor((z * xnp) / width);
                const px = column * xnp + x;
                const index = itemSize * ((row * ynp * width) + (y * width) + px);

                const v0 = set.includes(data[o(x + minx, y + miny, z + minz)]) ? 255 : 0;
                const xp = set.includes(data[o(Math.min(xn1 + maxx, x + 1 + minx), y + miny, z + minz)]) ? 255 : 0;
                const xn = set.includes(data[o(Math.max(0, x - 1 + minx), y + miny, z + minz)]) ? 255 : 0;
                const yp = set.includes(data[o(x + minx, Math.min(yn1 + maxy, y + 1 + miny), z + minz)]) ? 255 : 0;
                const yn = set.includes(data[o(x + minx, Math.max(0, y - 1 + miny), z + minz)]) ? 255 : 0;
                const zp = set.includes(data[o(x + minx, y + miny, Math.min(zn1 + maxz, z + 1 + minz))]) ? 255 : 0;
                const zn = set.includes(data[o(x + minx, y + miny, Math.max(0, z - 1 + minz))]) ? 255 : 0;

                array[index] = Math.round((v0 + v0 + xp + xn + yp + yn + zp + zn) / 8);
            }
        }
    }

    return textureImage;
}
