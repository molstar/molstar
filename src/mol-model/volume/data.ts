/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SpacegroupCell, Box3D } from 'mol-math/geometry'
import { Tensor, Mat4, Vec3 } from 'mol-math/linear-algebra'

interface VolumeData {
    // The basic unit cell that contains the data.
    readonly cell: SpacegroupCell,
    readonly fractionalBox: Box3D,
    readonly data: Tensor,
    readonly dataStats: Readonly<{
        min: number,
        max: number,
        mean: number,
        sigma: number
    }>
}

namespace VolumeData {
    const _scale = Mat4.zero(), _translate = Mat4.zero();
    export function getGridToCartesianTransform(volume: VolumeData) {
        const { data: { space } } = volume;
        const scale = Mat4.fromScaling(_scale, Vec3.div(Vec3.zero(), Box3D.size(Vec3.zero(), volume.fractionalBox), Vec3.ofArray(space.dimensions)));
        const translate = Mat4.fromTranslation(_translate, volume.fractionalBox.min);
        return Mat4.mul3(Mat4.zero(), volume.cell.fromFractional, translate, scale);
    }
}

type VolumeIsoValue = VolumeIsoValue.Absolute | VolumeIsoValue.Relative

namespace VolumeIsoValue {
    export type Relative = Readonly<{ kind: 'relative', stats: VolumeData['dataStats'], relativeValue: number }>
    export type Absolute = Readonly<{ kind: 'absolute', stats: VolumeData['dataStats'], absoluteValue: number }>

    export function absolute(stats: VolumeData['dataStats'], value: number): Absolute { return { kind: 'absolute', stats, absoluteValue: value }; }
    export function relative(stats: VolumeData['dataStats'], value: number): Relative { return { kind: 'relative', stats, relativeValue: value }; }

    export function calcAbsolute(stats: VolumeData['dataStats'], relativeValue: number): number {
        return relativeValue * stats.sigma + stats.mean
    }

    export function calcRelative(stats: VolumeData['dataStats'], absoluteValue: number): number {
        return (stats.mean - absoluteValue) / stats.sigma
    }

    export function toAbsolute(value: VolumeIsoValue): Absolute {
        if (value.kind === 'absolute') return value;
        return {
            kind: 'absolute',
            stats: value.stats,
            absoluteValue: calcAbsolute(value.stats, value.relativeValue)
        }
    }

    export function toRelative(value: VolumeIsoValue): Relative {
        if (value.kind === 'relative') return value;
        return {
            kind: 'relative',
            stats: value.stats,
            relativeValue: calcRelative(value.stats, value.absoluteValue)
        }
    }
}

export { VolumeData, VolumeIsoValue }