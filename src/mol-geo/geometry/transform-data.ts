/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util';
import { Mat4 } from '../../mol-math/linear-algebra';
import { fillSerial } from '../../mol-util/array';

export type TransformData = {
    /**
     * final per-instance transform calculated for instance `i` as
     * `aTransform[i] = matrix * transform[i] * extraTransform[i]`
     */
    aTransform: ValueCell<Float32Array>,
    /** global transform, see aTransform */
    matrix: ValueCell<Mat4>,
    /** base per-instance transform, see aTransform */
    transform: ValueCell<Float32Array>,
    /** additional per-instance transform, see aTransform */
    extraTransform: ValueCell<Float32Array>,

    uInstanceCount: ValueCell<number>,
    instanceCount: ValueCell<number>,
    aInstance: ValueCell<Float32Array>,
}

export function createTransform(transformArray: Float32Array, instanceCount: number, transformData?: TransformData): TransformData {
    if (transformData) {
        ValueCell.update(transformData.matrix, transformData.matrix.ref.value);
        ValueCell.update(transformData.transform, transformArray);
        ValueCell.update(transformData.uInstanceCount, instanceCount);
        ValueCell.update(transformData.instanceCount, instanceCount);

        const aTransform = transformData.aTransform.ref.value.length >= instanceCount * 16 ? transformData.aTransform.ref.value : new Float32Array(instanceCount * 16);
        aTransform.set(transformArray);
        ValueCell.update(transformData.aTransform, aTransform);

        // Note that this sets `extraTransform` to identity transforms
        const extraTransform = transformData.extraTransform.ref.value.length >= instanceCount * 16 ? transformData.extraTransform.ref.value : new Float32Array(instanceCount * 16);
        ValueCell.update(transformData.extraTransform, fillIdentityTransform(extraTransform, instanceCount));

        const aInstance = transformData.aInstance.ref.value.length >= instanceCount ? transformData.aInstance.ref.value : new Float32Array(instanceCount);
        ValueCell.update(transformData.aInstance, fillSerial(aInstance, instanceCount));

        updateTransformData(transformData);
        return transformData;
    } else {
        return {
            aTransform: ValueCell.create(new Float32Array(transformArray)),
            matrix: ValueCell.create(Mat4.identity()),
            transform: ValueCell.create(transformArray),
            extraTransform: ValueCell.create(fillIdentityTransform(new Float32Array(instanceCount * 16), instanceCount)),
            uInstanceCount: ValueCell.create(instanceCount),
            instanceCount: ValueCell.create(instanceCount),
            aInstance: ValueCell.create(fillSerial(new Float32Array(instanceCount)))
        };
    }
}

const identityTransform = new Float32Array(16);
Mat4.toArray(Mat4.identity(), identityTransform, 0);

export function createIdentityTransform(transformData?: TransformData): TransformData {
    return createTransform(new Float32Array(identityTransform), 1, transformData);
}

export function fillIdentityTransform(transform: Float32Array, count: number) {
    for (let i = 0; i < count; i++) {
        transform.set(identityTransform, i * 16);
    }
    return transform;
}

/**
 * updates per-instance transform calculated for instance `i` as
 * `aTransform[i] = matrix * transform[i] * extraTransform[i]`
 */
export function updateTransformData(transformData: TransformData) {
    const aTransform = transformData.aTransform.ref.value;
    const instanceCount = transformData.instanceCount.ref.value;
    const matrix = transformData.matrix.ref.value;
    const transform = transformData.transform.ref.value;
    const extraTransform = transformData.extraTransform.ref.value;
    for (let i = 0; i < instanceCount; i++) {
        const i16 = i * 16;
        Mat4.mulOffset(aTransform, extraTransform, transform, i16, i16, i16);
        Mat4.mulOffset(aTransform, matrix, aTransform, i16, 0, i16);
    }
    ValueCell.update(transformData.aTransform, aTransform);
}