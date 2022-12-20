/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3set = Vec3.set;
const v3normalize = Vec3.normalize;
const v3sub = Vec3.sub;
const v3addScalar = Vec3.addScalar;
const v3scale = Vec3.scale;
const v3toArray = Vec3.toArray;

export function eachVolumeLoci(loci: Loci, volume: Volume, props: { isoValue?: Volume.IsoValue, segments?: SortedArray } | undefined, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (Volume.isLoci(loci)) {
        if (!Volume.areEquivalent(loci.volume, volume)) return false;
        if (apply(Interval.ofLength(volume.grid.cells.data.length))) changed = true;
    } else if (Volume.Isosurface.isLoci(loci)) {
        if (!Volume.areEquivalent(loci.volume, volume)) return false;
        if (props?.isoValue) {
            if (!Volume.IsoValue.areSame(loci.isoValue, props.isoValue, volume.grid.stats)) return false;
            if (apply(Interval.ofLength(volume.grid.cells.data.length))) changed = true;
        } else {
            const { stats, cells: { data } } = volume.grid;
            const eps = stats.sigma;
            const v = Volume.IsoValue.toAbsolute(loci.isoValue, stats).absoluteValue;
            for (let i = 0, il = data.length; i < il; ++i) {
                if (equalEps(v, data[i], eps)) {
                    if (apply(Interval.ofSingleton(i))) changed = true;
                }
            }
        }
    } else if (Volume.Cell.isLoci(loci)) {
        if (!Volume.areEquivalent(loci.volume, volume)) return false;
        if (Interval.is(loci.indices)) {
            if (apply(loci.indices)) changed = true;
        } else {
            OrderedSet.forEach(loci.indices, v => {
                if (apply(Interval.ofSingleton(v))) changed = true;
            });
        }
    } else if (Volume.Segment.isLoci(loci)) {
        if (!Volume.areEquivalent(loci.volume, volume)) return false;
        if (props?.segments) {
            if (!SortedArray.areIntersecting(loci.segments, props.segments)) return false;
            if (apply(Interval.ofLength(volume.grid.cells.data.length))) changed = true;
        } else {
            const segmentation = Volume.Segmentation.get(volume);
            if (segmentation) {
                const set = new Set<number>();
                for (let i = 0, il = loci.segments.length; i < il; ++i) {
                    SetUtils.add(set, segmentation.segments.get(loci.segments[i])!);
                }
                const s = Array.from(set.values());
                const d = volume.grid.cells.data;
                for (let i = 0, il = d.length; i < il; ++i) {
                    if (s.includes(d[i])) {
                        if (apply(Interval.ofSingleton(i))) changed = true;
                    }
                }
            }
        }
    }
    return changed;
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

export function createVolumeTexture2d(volume: Volume, variant: 'normals' | 'groups' | 'data', padding = 0) {
    const { cells: { space, data }, stats: { max, min } } = volume.grid;
    const dim = space.dimensions as Vec3;
    const { dataOffset: o } = space;
    const { width, height } = getVolumeTexture2dLayout(dim, padding);

    const itemSize = variant === 'data' ? 1 : 4;
    const array = new Uint8Array(width * height * itemSize);
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

                if (variant === 'data') {
                    array[index] = Math.round(((data[offset] - min) / diff) * 255);
                } else {
                    if (variant === 'groups') {
                        packIntToRGBArray(offset, array, index);
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
                        v3toArray(v3scale(n0, n0, 255), array, index);
                    }

                    array[index + 3] = Math.round(((data[offset] - min) / diff) * 255);
                }
            }
        }
    }

    return textureImage;
}

export function createVolumeTexture3d(volume: Volume) {
    const { cells: { space, data }, stats: { max, min } } = volume.grid;
    const [width, height, depth] = space.dimensions as Vec3;
    const { dataOffset: o } = space;

    const array = new Uint8Array(width * height * depth * 4);
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
                v3toArray(v3scale(n0, n0, 255), array, i);

                array[i + 3] = Math.round(((data[offset] - min) / diff) * 255);
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
