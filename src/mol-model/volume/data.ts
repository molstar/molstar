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
        const scale = Mat4.fromScaling(_scale, Vec3.div(Vec3.zero(), Box3D.size(volume.fractionalBox), Vec3.ofArray(space.dimensions)));
        const translate = Mat4.fromTranslation(_translate, volume.fractionalBox.min);
        return Mat4.mul3(Mat4.zero(), volume.cell.fromFractional, translate, scale);
    }
}

export { VolumeData }