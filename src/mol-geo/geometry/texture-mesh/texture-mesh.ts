/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Cai Huiyu <szmun.caihy@gmail.com>
 */

import { ValueCell } from '../../../mol-util';
import { Sphere3D } from '../../../mol-math/geometry';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { LocationIterator, PositionLocation } from '../../../mol-geo/util/location-iterator';
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
import { createNullTexture, Texture } from '../../../mol-gl/webgl/texture';
import { Vec2, Vec3, Vec4 } from '../../../mol-math/linear-algebra';
import { createEmptyClipping } from '../clipping-data';
import { NullLocation } from '../../../mol-model/location';
import { createEmptySubstance } from '../substance-data';
import { RenderableState } from '../../../mol-gl/renderable';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { createEmptyEmissive } from '../emissive-data';

export interface TextureMesh {
    readonly kind: 'texture-mesh',

    /** Number of vertices in the texture-mesh */
    vertexCount: number,
    /** Number of groups in the texture-mesh */
    groupCount: number,

    readonly geoTextureDim: ValueCell<Vec2>,
    readonly vertexTexture: ValueCell<Texture>,
    readonly groupTexture: ValueCell<Texture>,
    readonly normalTexture: ValueCell<Texture>,
    readonly varyingGroup: ValueCell<boolean>,
    readonly doubleBuffer: TextureMesh.DoubleBuffer

    readonly boundingSphere: Sphere3D

    readonly meta: {
        webgl?: WebGLContext
        [k: string]: unknown
    }
}

export namespace TextureMesh {
    export class DoubleBuffer {
        private index = 0;
        private textures: ({ vertex: Texture, group: Texture, normal: Texture } | undefined)[] = [];

        get() {
            return this.textures[this.index];
        }

        set(vertex: Texture, group: Texture, normal: Texture) {
            this.textures[this.index] = Object.assign(this.textures[this.index] || {}, {
                vertex, group, normal
            });
            this.index = (this.index + 1) % 2;
        }

        destroy() {
            for (const buffer of this.textures) {
                buffer!.vertex.destroy();
                buffer!.group.destroy();
                buffer!.normal.destroy();
            }
        }
    }

    export function create(vertexCount: number, groupCount: number, vertexTexture: Texture, groupTexture: Texture, normalTexture: Texture, boundingSphere: Sphere3D, textureMesh?: TextureMesh): TextureMesh {
        const width = vertexTexture.getWidth();
        const height = vertexTexture.getHeight();
        if (textureMesh) {
            textureMesh.vertexCount = vertexCount;
            textureMesh.groupCount = groupCount;
            ValueCell.update(textureMesh.geoTextureDim, Vec2.set(textureMesh.geoTextureDim.ref.value, width, height));
            ValueCell.update(textureMesh.vertexTexture, vertexTexture);
            ValueCell.update(textureMesh.groupTexture, groupTexture);
            ValueCell.update(textureMesh.normalTexture, normalTexture);
            textureMesh.doubleBuffer.set(vertexTexture, groupTexture, normalTexture);
            Sphere3D.copy(textureMesh.boundingSphere, boundingSphere);
            return textureMesh;
        } else {
            return {
                kind: 'texture-mesh',
                vertexCount,
                groupCount,
                geoTextureDim: ValueCell.create(Vec2.create(width, height)),
                vertexTexture: ValueCell.create(vertexTexture),
                groupTexture: ValueCell.create(groupTexture),
                normalTexture: ValueCell.create(normalTexture),
                varyingGroup: ValueCell.create(false),
                doubleBuffer: new DoubleBuffer(),
                boundingSphere: Sphere3D.clone(boundingSphere),
                meta: {}
            };
        }
    }

    export function createEmpty(textureMesh?: TextureMesh): TextureMesh {
        const vt = textureMesh ? textureMesh.vertexTexture.ref.value : createNullTexture();
        const gt = textureMesh ? textureMesh.groupTexture.ref.value : createNullTexture();
        const nt = textureMesh ? textureMesh.normalTexture.ref.value : createNullTexture();
        const bs = textureMesh ? textureMesh.boundingSphere : Sphere3D();
        return create(0, 0, vt, gt, nt, bs, textureMesh);
    }

    export const Params = {
        ...BaseGeometry.Params,
        doubleSided: PD.Boolean(false, BaseGeometry.CustomQualityParamInfo),
        flipSided: PD.Boolean(false, BaseGeometry.ShadingCategory),
        flatShaded: PD.Boolean(false, BaseGeometry.ShadingCategory),
        ignoreLight: PD.Boolean(false, BaseGeometry.ShadingCategory),
        celShaded: PD.Boolean(false, BaseGeometry.ShadingCategory),
        xrayShaded: PD.Select<boolean | 'inverted'>(false, [[false, 'Off'], [true, 'On'], ['inverted', 'Inverted']], BaseGeometry.ShadingCategory),
        transparentBackfaces: PD.Select('off', PD.arrayToOptions(['off', 'on', 'opaque'] as const), BaseGeometry.ShadingCategory),
        bumpFrequency: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
        bumpAmplitude: PD.Numeric(1, { min: 0, max: 5, step: 0.1 }, BaseGeometry.ShadingCategory),
    };
    export type Params = typeof Params

    export const Utils: GeometryUtils<TextureMesh, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState,
        updateRenderableState,
        createPositionIterator,
    };

    const TextureMeshName = 'texture-mesh';

    function createPositionIterator(textureMesh: TextureMesh, transform: TransformData): LocationIterator {
        const webgl = textureMesh.meta.webgl;
        if (!webgl) return LocationIterator(1, 1, 1, () => NullLocation);

        if (!webgl.namedFramebuffers[TextureMeshName]) {
            webgl.namedFramebuffers[TextureMeshName] = webgl.resources.framebuffer();
        }
        const framebuffer = webgl.namedFramebuffers[TextureMeshName];
        const [width, height] = textureMesh.geoTextureDim.ref.value;
        const vertices = new Float32Array(width * height * 4);
        framebuffer.bind();
        textureMesh.vertexTexture.ref.value.attachFramebuffer(framebuffer, 0);
        webgl.readPixels(0, 0, width, height, vertices);

        const normals = new Float32Array(width * height * 4);
        framebuffer.bind();
        textureMesh.normalTexture.ref.value.attachFramebuffer(framebuffer, 0);
        webgl.readPixels(0, 0, width, height, normals);

        const groupCount = textureMesh.vertexCount;
        const instanceCount = transform.instanceCount.ref.value;
        const location = PositionLocation();
        const p = location.position;
        const n = location.normal;
        const m = transform.aTransform.ref.value;
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            if (instanceIndex < 0) {
                Vec3.fromArray(p, vertices, groupIndex * 4);
                Vec3.fromArray(n, normals, groupIndex * 4);
            } else {
                Vec3.transformMat4Offset(p, vertices, m, 0, groupIndex * 4, instanceIndex * 16);
                Vec3.transformDirectionOffset(n, normals, m, 0, groupIndex * 4, instanceIndex * 16);
            }
            return location;
        };
        return LocationIterator(groupCount, instanceCount, 1, getLocation);
    }

    function createValues(textureMesh: TextureMesh, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): TextureMeshValues {
        const { instanceCount, groupCount } = locationIt;
        const positionIt = Utils.createPositionIterator(textureMesh, transform);

        const color = createColors(locationIt, positionIt, theme.color);
        const marker = props.instanceGranularity
            ? createMarkers(instanceCount, 'instance')
            : createMarkers(instanceCount * groupCount, 'groupInstance');
        const overpaint = createEmptyOverpaint();
        const transparency = createEmptyTransparency();
        const emissive = createEmptyEmissive();
        const substance = createEmptySubstance();
        const clipping = createEmptyClipping();

        const counts = { drawCount: textureMesh.vertexCount, vertexCount: textureMesh.vertexCount, groupCount, instanceCount };

        const invariantBoundingSphere = Sphere3D.clone(textureMesh.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform.aTransform.ref.value, instanceCount, 0);

        return {
            dGeometryType: ValueCell.create('textureMesh'),

            uGeoTexDim: textureMesh.geoTextureDim,
            tPosition: textureMesh.vertexTexture,
            tGroup: textureMesh.groupTexture,
            tNormal: textureMesh.normalTexture,
            dVaryingGroup: textureMesh.varyingGroup,

            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),

            ...color,
            ...marker,
            ...overpaint,
            ...transparency,
            ...emissive,
            ...substance,
            ...clipping,
            ...transform,

            ...BaseGeometry.createValues(props, counts),
            uDoubleSided: ValueCell.create(props.doubleSided),
            dFlatShaded: ValueCell.create(props.flatShaded),
            dFlipSided: ValueCell.create(props.flipSided),
            dIgnoreLight: ValueCell.create(props.ignoreLight),
            dCelShaded: ValueCell.create(props.celShaded),
            dXrayShaded: ValueCell.create(props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off'),
            dTransparentBackfaces: ValueCell.create(props.transparentBackfaces),
            uBumpFrequency: ValueCell.create(props.bumpFrequency),
            uBumpAmplitude: ValueCell.create(props.bumpAmplitude),

            meta: ValueCell.create(textureMesh.meta),
        };
    }

    function createValuesSimple(textureMesh: TextureMesh, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
        const p = { ...PD.getDefaultValues(Params), ...props };
        return createValues(textureMesh, s.transform, s.locationIterator, s.theme, p);
    }

    function updateValues(values: TextureMeshValues, props: PD.Values<Params>) {
        BaseGeometry.updateValues(values, props);
        ValueCell.updateIfChanged(values.uDoubleSided, props.doubleSided);
        ValueCell.updateIfChanged(values.dFlatShaded, props.flatShaded);
        ValueCell.updateIfChanged(values.dFlipSided, props.flipSided);
        ValueCell.updateIfChanged(values.dIgnoreLight, props.ignoreLight);
        ValueCell.updateIfChanged(values.dCelShaded, props.celShaded);
        ValueCell.updateIfChanged(values.dXrayShaded, props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off');
        ValueCell.updateIfChanged(values.dTransparentBackfaces, props.transparentBackfaces);
        ValueCell.updateIfChanged(values.uBumpFrequency, props.bumpFrequency);
        ValueCell.updateIfChanged(values.uBumpAmplitude, props.bumpAmplitude);
    }

    function updateBoundingSphere(values: TextureMeshValues, textureMesh: TextureMesh) {
        const invariantBoundingSphere = Sphere3D.clone(textureMesh.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, values.aTransform.ref.value, values.instanceCount.ref.value, 0);

        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere);
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
            ValueCell.update(values.uInvariantBoundingSphere, Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere));
        }
    }

    function createRenderableState(props: PD.Values<Params>): RenderableState {
        const state = BaseGeometry.createRenderableState(props);
        updateRenderableState(state, props);
        return state;
    }

    function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
        BaseGeometry.updateRenderableState(state, props);
        state.opaque = state.opaque && !props.xrayShaded;
        state.writeDepth = state.opaque;
    }
}