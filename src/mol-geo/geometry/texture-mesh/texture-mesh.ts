/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../../mol-util';
import { Sphere3D } from '../../../mol-math/geometry';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { TransformData } from '../transform-data';
import { createColors } from '../color-data';
import { createMarkers } from '../marker-data';
import { GeometryUtils } from '../geometry';
import { Theme } from '../../../mol-theme/theme';
import { Color } from '../../../mol-util/color';
import { BaseGeometry } from '../base';
import { createEmptyOverpaint } from '../overpaint-data';
import { createEmptyTransparency } from '../transparency-data';
import { TextureMeshValues } from '../../../mol-gl/renderable/texture-mesh';
import { calculateTransformBoundingSphere } from '../../../mol-gl/renderable/util';
import { Texture } from '../../../mol-gl/webgl/texture';
import { Vec2, Vec4 } from '../../../mol-math/linear-algebra';
import { fillSerial } from '../../../mol-util/array';
import { createEmptyClipping } from '../clipping-data';

export interface TextureMesh {
    readonly kind: 'texture-mesh',

    /** Number of vertices in the texture-mesh */
    vertexCount: number,
    /** Number of groups in the texture-mesh */
    groupCount: number,

    readonly geoTextureDim: ValueCell<Vec2>,
    /** texture has vertex positions in XYZ and group id in W */
    readonly vertexGroupTexture: ValueCell<Texture>,
    readonly normalTexture: ValueCell<Texture>,

    readonly boundingSphere: Sphere3D
}

export namespace TextureMesh {
    export function create(vertexCount: number, groupCount: number, vertexGroupTexture: Texture, normalTexture: Texture, boundingSphere: Sphere3D, textureMesh?: TextureMesh): TextureMesh {
        const width = vertexGroupTexture.getWidth();
        const height = vertexGroupTexture.getHeight();
        if (textureMesh) {
            textureMesh.vertexCount = vertexCount;
            textureMesh.groupCount = groupCount;
            ValueCell.update(textureMesh.geoTextureDim, Vec2.set(textureMesh.geoTextureDim.ref.value, width, height));
            ValueCell.update(textureMesh.vertexGroupTexture, vertexGroupTexture);
            ValueCell.update(textureMesh.normalTexture, normalTexture);
            Sphere3D.copy(textureMesh.boundingSphere, boundingSphere);
            return textureMesh;
        } else {
            return {
                kind: 'texture-mesh',
                vertexCount,
                groupCount,
                geoTextureDim: ValueCell.create(Vec2.create(width, height)),
                vertexGroupTexture: ValueCell.create(vertexGroupTexture),
                normalTexture: ValueCell.create(normalTexture),
                boundingSphere: Sphere3D.clone(boundingSphere),
            };
        }
    }

    export function createEmpty(textureMesh?: TextureMesh): TextureMesh {
        return {} as TextureMesh; // TODO
    }

    export const Params = {
        ...BaseGeometry.Params,
        doubleSided: PD.Boolean(false, BaseGeometry.CustomQualityParamInfo),
        flipSided: PD.Boolean(false, BaseGeometry.ShadingCategory),
        flatShaded: PD.Boolean(false, BaseGeometry.ShadingCategory),
    };
    export type Params = typeof Params

    export const Utils: GeometryUtils<TextureMesh, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState: BaseGeometry.createRenderableState,
        updateRenderableState: BaseGeometry.updateRenderableState
    };

    function createValues(textureMesh: TextureMesh, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): TextureMeshValues {
        const { instanceCount, groupCount } = locationIt;
        const color = createColors(locationIt, theme.color);
        const marker = createMarkers(instanceCount * groupCount);
        const overpaint = createEmptyOverpaint();
        const transparency = createEmptyTransparency();
        const clipping = createEmptyClipping();

        const counts = { drawCount: textureMesh.vertexCount, groupCount, instanceCount };

        const transformBoundingSphere = calculateTransformBoundingSphere(textureMesh.boundingSphere, transform.aTransform.ref.value, transform.instanceCount.ref.value);

        return {
            uGeoTexDim: textureMesh.geoTextureDim,
            tPositionGroup: textureMesh.vertexGroupTexture,
            tNormal: textureMesh.normalTexture,

            // aGroup is used as a vertex index here and the group id is retirieved from tPositionGroup
            aGroup: ValueCell.create(fillSerial(new Float32Array(textureMesh.vertexCount))),
            boundingSphere: ValueCell.create(transformBoundingSphere),
            invariantBoundingSphere: ValueCell.create(Sphere3D.clone(textureMesh.boundingSphere)),
            uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(textureMesh.boundingSphere)),

            ...color,
            ...marker,
            ...overpaint,
            ...transparency,
            ...clipping,
            ...transform,

            ...BaseGeometry.createValues(props, counts),
            dDoubleSided: ValueCell.create(props.doubleSided),
            dFlatShaded: ValueCell.create(props.flatShaded),
            dFlipSided: ValueCell.create(props.flipSided),
            dGeoTexture: ValueCell.create(true),
        };
    }

    function createValuesSimple(textureMesh: TextureMesh, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
        const p = { ...PD.getDefaultValues(Params), ...props };
        return createValues(textureMesh, s.transform, s.locationIterator, s.theme, p);
    }

    function updateValues(values: TextureMeshValues, props: PD.Values<Params>) {
        ValueCell.updateIfChanged(values.alpha, props.alpha); // `uAlpha` is set in renderable.render

        ValueCell.updateIfChanged(values.dDoubleSided, props.doubleSided);
        ValueCell.updateIfChanged(values.dFlatShaded, props.flatShaded);
        ValueCell.updateIfChanged(values.dFlipSided, props.flipSided);

        if (values.drawCount.ref.value > values.aGroup.ref.value.length) {
            // console.log('updating vertex ids in aGroup to handle larger drawCount')
            ValueCell.update(values.aGroup, fillSerial(new Float32Array(values.drawCount.ref.value)));
        }
    }

    function updateBoundingSphere(values: TextureMeshValues, textureMesh: TextureMesh) {
        const invariantBoundingSphere = textureMesh.boundingSphere;
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, values.aTransform.ref.value, values.instanceCount.ref.value);
        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere);
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
            ValueCell.update(values.uInvariantBoundingSphere, Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere));
        }
    }
}