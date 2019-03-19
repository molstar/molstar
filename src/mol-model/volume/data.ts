/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SpacegroupCell, Box3D } from 'mol-math/geometry'
import { Tensor, Mat4, Vec3 } from 'mol-math/linear-algebra'
import { equalEps } from 'mol-math/linear-algebra/3d/common';
import { Histogram } from 'mol-math/misc';

/** The basic unit cell that contains the data. */
interface VolumeData {
    readonly cell: SpacegroupCell,
    readonly fractionalBox: Box3D,
    readonly data: Tensor,
    readonly dataStats: Readonly<{
        min: number,
        max: number,
        mean: number,
        sigma: number,
        histogram?: Histogram
    }>
}

namespace VolumeData {
    export const One: VolumeData = {
        cell: SpacegroupCell.Zero,
        fractionalBox: Box3D.empty(),
        data: Tensor.create(Tensor.Space([1, 1, 1], [0, 1, 2]), Tensor.Data1([0])),
        dataStats: { min: 0, max: 0, mean: 0, sigma: 0, histogram: { bins: [0, 0], counts: [1] } }
    }

    const _scale = Mat4.zero(), _translate = Mat4.zero();
    export function getGridToCartesianTransform(volume: VolumeData) {
        const { data: { space } } = volume;
        const scale = Mat4.fromScaling(_scale, Vec3.div(Vec3.zero(), Box3D.size(Vec3.zero(), volume.fractionalBox), Vec3.ofArray(space.dimensions)));
        const translate = Mat4.fromTranslation(_translate, volume.fractionalBox.min);
        return Mat4.mul3(Mat4.zero(), volume.cell.fromFractional, translate, scale);
    }

    export function areEquivalent(volA: VolumeData, volB: VolumeData) {
        return volA === volB
    }
}

type VolumeIsoValue = VolumeIsoValue.Absolute | VolumeIsoValue.Relative

namespace VolumeIsoValue {
    export type Relative = Readonly<{ kind: 'relative', relativeValue: number }>
    export type Absolute = Readonly<{ kind: 'absolute', absoluteValue: number }>

    export function areSame(a: VolumeIsoValue, b: VolumeIsoValue, stats: VolumeData['dataStats']) {
        return equalEps(toAbsolute(a, stats).absoluteValue, toAbsolute(b, stats).absoluteValue, stats.sigma / 100)
    }

    export function absolute(value: number): Absolute { return { kind: 'absolute', absoluteValue: value }; }
    export function relative(value: number): Relative { return { kind: 'relative', relativeValue: value }; }

    export function calcAbsolute(stats: VolumeData['dataStats'], relativeValue: number): number {
        return relativeValue * stats.sigma + stats.mean
    }

    export function calcRelative(stats: VolumeData['dataStats'], absoluteValue: number): number {
        return stats.sigma === 0 ? 0 : ((absoluteValue - stats.mean) / stats.sigma)
    }

    export function toAbsolute(value: VolumeIsoValue, stats: VolumeData['dataStats']): Absolute {
        return value.kind === 'absolute' ? value : { kind: 'absolute', absoluteValue: VolumeIsoValue.calcAbsolute(stats, value.relativeValue) }
    }

    export function toRelative(value: VolumeIsoValue, stats: VolumeData['dataStats']): Relative {
        return value.kind === 'relative' ? value : { kind: 'relative', relativeValue: VolumeIsoValue.calcRelative(stats, value.absoluteValue) }
    }

    export function toString(value: VolumeIsoValue) {
        return value.kind === 'relative'
            ? `${value.relativeValue} σ`
            : `${value.absoluteValue}`
    }
}

export { VolumeData, VolumeIsoValue }