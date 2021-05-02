/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { calcMeshColorSmoothing } from '../../../../mol-geo/geometry/mesh/color-smoothing';
import { calcTextureMeshColorSmoothing } from '../../../../mol-geo/geometry/texture-mesh/color-smoothing';
import { MeshValues } from '../../../../mol-gl/renderable/mesh';
import { TextureMeshValues } from '../../../../mol-gl/renderable/texture-mesh';
import { WebGLContext } from '../../../../mol-gl/webgl/context';
import { Texture } from '../../../../mol-gl/webgl/texture';
import { ValueCell } from '../../../../mol-util';

function isSupportedColorType(x: string): x is 'group' | 'groupInstance' | 'vertex' | 'vertexInstance' {
    return x === 'group' || x === 'groupInstance' || x === 'vertex' || x === 'vertexInstance';
}

export function applyMeshColorSmoothing(webgl: WebGLContext, values: MeshValues, resolution: number) {
    if (!isSupportedColorType(values.dColorType.ref.value)) return;

    const smoothingData = calcMeshColorSmoothing(webgl, {
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.uInstanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionBuffer: values.aPosition.ref.value,
        groupBuffer: values.aGroup.ref.value,
        colorData: values.tColor.ref.value,
        colorType: values.dColorType.ref.value,
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
    }, resolution * 2, 6);

    ValueCell.updateIfChanged(values.dColorType, smoothingData.type);
    ValueCell.update(values.tColorGrid, smoothingData.texture);
    ValueCell.update(values.uColorTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uColorGridDim, smoothingData.gridDim);
    ValueCell.update(values.uColorGridTransform, smoothingData.transform);
    ValueCell.updateIfChanged(values.dColorGridType, '2d');
}

export function applyTextureMeshColorSmoothing(webgl: WebGLContext, values: TextureMeshValues, resolution: number, colorTexture?: Texture) {
    if (!isSupportedColorType(values.dColorType.ref.value)) return;

    const smoothingData = calcTextureMeshColorSmoothing(webgl, {
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.uInstanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionTexture: values.tPosition.ref.value,
        groupTexture: values.tGroup.ref.value,
        colorData: values.tColor.ref.value,
        colorType: values.dColorType.ref.value,
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
    }, resolution * 2, 12, colorTexture);

    ValueCell.updateIfChanged(values.dColorType, smoothingData.type);
    ValueCell.update(values.tColorGrid, smoothingData.texture);
    ValueCell.update(values.uColorTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uColorGridDim, smoothingData.gridDim);
    ValueCell.update(values.uColorGridTransform, smoothingData.transform);
    ValueCell.updateIfChanged(values.dColorGridType, '2d');
}