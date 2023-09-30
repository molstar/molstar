/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util';
import { Mat4, Mat3 } from '../../mol-math/linear-algebra';
import { fillSerial } from '../../mol-util/array';
import { Sphere3D } from '../../mol-math/geometry';
import { calcInstanceGrid, createEmptyInstanceGrid, InstanceGrid } from '../../mol-math/geometry/instance-grid';

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

    hasReflection: ValueCell<boolean>,

    instanceGrid: ValueCell<InstanceGrid>,
}

const _m3 = Mat3();
const _m4 = Mat4();
function checkReflection(transformArray: Float32Array, instanceCount: number) {
    for (let i = 0; i < instanceCount; i++) {
        Mat3.fromMat4(_m3, Mat4.fromArray(_m4, transformArray, i * 16));
        if (Mat3.determinant(_m3) < 0) return true;
    }
    return false;
}

export function createTransform(transformArray: Float32Array, instanceCount: number, invariantBoundingSphere: Sphere3D | undefined, cellSize: number, batchSize: number, transformData?: TransformData): TransformData {
    const hasReflection = checkReflection(transformArray, instanceCount);

    if (transformData) {
        ValueCell.update(transformData.matrix, transformData.matrix.ref.value);
        const transform = transformData.transform.ref.value.length >= instanceCount * 16 ? transformData.transform.ref.value : new Float32Array(instanceCount * 16);
        transform.set(transformArray);
        ValueCell.update(transformData.transform, transform);
        ValueCell.updateIfChanged(transformData.uInstanceCount, instanceCount);
        ValueCell.updateIfChanged(transformData.instanceCount, instanceCount);

        const aTransform = transformData.aTransform.ref.value.length >= instanceCount * 16 ? transformData.aTransform.ref.value : new Float32Array(instanceCount * 16);
        ValueCell.update(transformData.aTransform, aTransform);

        // Note that this sets `extraTransform` to identity transforms
        const extraTransform = transformData.extraTransform.ref.value.length >= instanceCount * 16 ? transformData.extraTransform.ref.value : new Float32Array(instanceCount * 16);
        ValueCell.update(transformData.extraTransform, fillIdentityTransform(extraTransform, instanceCount));

        const aInstance = transformData.aInstance.ref.value.length >= instanceCount ? transformData.aInstance.ref.value : new Float32Array(instanceCount);
        ValueCell.update(transformData.aInstance, fillSerial(aInstance, instanceCount));

        ValueCell.update(transformData.hasReflection, hasReflection);
    } else {
        transformData = {
            aTransform: ValueCell.create(new Float32Array(instanceCount * 16)),
            matrix: ValueCell.create(Mat4.identity()),
            transform: ValueCell.create(new Float32Array(transformArray)),
            extraTransform: ValueCell.create(fillIdentityTransform(new Float32Array(instanceCount * 16), instanceCount)),
            uInstanceCount: ValueCell.create(instanceCount),
            instanceCount: ValueCell.create(instanceCount),
            aInstance: ValueCell.create(fillSerial(new Float32Array(instanceCount))),
            hasReflection: ValueCell.create(hasReflection),
            instanceGrid: ValueCell.create(createEmptyInstanceGrid()),
        };
    }

    updateTransformData(transformData, invariantBoundingSphere, cellSize, batchSize);
    return transformData;
}

const identityTransform = new Float32Array(16);
Mat4.toArray(Mat4.identity(), identityTransform, 0);

export function createIdentityTransform(transformData?: TransformData): TransformData {
    return createTransform(new Float32Array(identityTransform), 1, undefined, 0, 0, transformData);
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
export function updateTransformData(transformData: TransformData, invariantBoundingSphere: Sphere3D | undefined, cellSize: number, batchSize: number) {
    const aTransform = transformData.aTransform.ref.value;
    const aInstance = transformData.aInstance.ref.value;
    const instanceCount = transformData.instanceCount.ref.value;
    const matrix = transformData.matrix.ref.value;
    const transform = transformData.transform.ref.value;
    const extraTransform = transformData.extraTransform.ref.value;
    for (let i = 0; i < instanceCount; i++) {
        const i16 = i * 16;
        Mat4.mulOffset(aTransform, extraTransform, transform, i16, i16, i16);
        Mat4.mulOffset(aTransform, matrix, aTransform, i16, 0, i16);
        aInstance[i] = i;
    }

    if (invariantBoundingSphere && instanceCount > 0) {
        const instanceGrid = calcInstanceGrid({
            instanceCount,
            instance: aInstance,
            transform: aTransform,
            invariantBoundingSphere
        }, cellSize, batchSize);

        ValueCell.update(transformData.instanceGrid, instanceGrid);
        ValueCell.update(transformData.aInstance, instanceGrid.cellInstance);
        ValueCell.update(transformData.aTransform, instanceGrid.cellTransform);
    } else {
        ValueCell.update(transformData.aInstance, aInstance);
        ValueCell.update(transformData.aTransform, aTransform);
    }
}
