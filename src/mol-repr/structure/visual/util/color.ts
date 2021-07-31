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
import { smoothstep } from '../../../../mol-math/interpolate';
import { Theme } from '../../../../mol-theme/theme';
import { ValueCell } from '../../../../mol-util';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';

export const ColorSmoothingParams = {
    smoothColors: PD.MappedStatic('auto', {
        auto: PD.Group({}),
        on: PD.Group({
            resolutionFactor: PD.Numeric(2, { min: 0.5, max: 6, step: 0.1 }),
            sampleStride: PD.Numeric(3, { min: 1, max: 12, step: 1 }),
        }),
        off: PD.Group({})
    }),
};
export type ColorSmoothingParams = typeof ColorSmoothingParams

export function getColorSmoothingProps(props: PD.Values<ColorSmoothingParams>, theme: Theme, resolution?: number) {
    if ((props.smoothColors.name === 'on' || (props.smoothColors.name === 'auto' && theme.color.preferSmoothing)) && resolution && resolution < 3) {
        let stride = 3;
        if (props.smoothColors.name === 'on') {
            resolution *= props.smoothColors.params.resolutionFactor;
            stride = props.smoothColors.params.sampleStride;
        } else {
            // https://graphtoy.com/?f1(x,t)=(2-smoothstep(0,1.1,x))*x&coords=0.7,0.6,1.8
            resolution *= 2 - smoothstep(0, 1.1, resolution);
            resolution = Math.max(0.5, resolution);
            if (resolution > 1.2) stride = 2;
        }
        return { resolution, stride };
    };
}

function isSupportedColorType(x: string): x is 'group' | 'groupInstance' {
    return x === 'group' || x === 'groupInstance';
}

export function applyMeshColorSmoothing(values: MeshValues, resolution: number, stride: number, webgl?: WebGLContext, colorTexture?: Texture) {
    if (!isSupportedColorType(values.dColorType.ref.value)) return;

    const smoothingData = calcMeshColorSmoothing({
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
    }, resolution, stride, webgl, colorTexture);

    if (smoothingData.kind === 'volume') {
        ValueCell.updateIfChanged(values.dColorType, smoothingData.type);
        ValueCell.update(values.tColorGrid, smoothingData.texture);
        ValueCell.update(values.uColorTexDim, smoothingData.gridTexDim);
        ValueCell.update(values.uColorGridDim, smoothingData.gridDim);
        ValueCell.update(values.uColorGridTransform, smoothingData.gridTransform);
    } else if (smoothingData.kind === 'vertex') {
        ValueCell.updateIfChanged(values.dColorType, smoothingData.type);
        ValueCell.update(values.tColor, smoothingData.texture);
        ValueCell.update(values.uColorTexDim, smoothingData.texDim);
    }
}

export function applyTextureMeshColorSmoothing(values: TextureMeshValues, resolution: number, stride: number, webgl: WebGLContext, colorTexture?: Texture) {
    if (!isSupportedColorType(values.dColorType.ref.value)) return;

    stride *= 3; // triple because TextureMesh is never indexed (no elements buffer)

    const smoothingData = calcTextureMeshColorSmoothing({
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
    }, resolution, stride, webgl, colorTexture);

    ValueCell.updateIfChanged(values.dColorType, smoothingData.type);
    ValueCell.update(values.tColorGrid, smoothingData.texture);
    ValueCell.update(values.uColorTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uColorGridDim, smoothingData.gridDim);
    ValueCell.update(values.uColorGridTransform, smoothingData.gridTransform);
}