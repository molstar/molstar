/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../../mol-util';
import { GeometryUtils } from '../geometry';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { TransformData } from '../transform-data';
import { LocationIterator, PositionLocation } from '../../../mol-geo/util/location-iterator';
import { Theme } from '../../../mol-theme/theme';
import { SpheresValues } from '../../../mol-gl/renderable/spheres';
import { createColors } from '../color-data';
import { createMarkers } from '../marker-data';
import { TextureImage, calculateInvariantBoundingSphere, calculateTransformBoundingSphere, createTextureImage } from '../../../mol-gl/renderable/util';
import { Sphere3D } from '../../../mol-math/geometry';
import { createSizes, getMaxSize } from '../size-data';
import { Color } from '../../../mol-util/color';
import { BaseGeometry } from '../base';
import { createEmptyOverpaint } from '../overpaint-data';
import { createEmptyTransparency } from '../transparency-data';
import { hashFnv32a } from '../../../mol-data/util';
import { GroupMapping, createGroupMapping } from '../../util';
import { createEmptyClipping } from '../clipping-data';
import { Vec2, Vec3, Vec4 } from '../../../mol-math/linear-algebra';
import { RenderableState } from '../../../mol-gl/renderable';
import { createEmptySubstance } from '../substance-data';
import { createEmptyEmissive } from '../emissive-data';

export interface Spheres {
    readonly kind: 'spheres',

    /** Number of spheres */
    sphereCount: number,

    /** Center buffer as array of xyz values wrapped in a value cell */
    readonly centerBuffer: ValueCell<Float32Array>,
    /** Group buffer as array of group ids for each vertex wrapped in a value cell */
    readonly groupBuffer: ValueCell<Float32Array>,

    /** Bounding sphere of the spheres */
    readonly boundingSphere: Sphere3D
    /** Maps group ids to sphere indices */
    readonly groupMapping: GroupMapping

    setBoundingSphere(boundingSphere: Sphere3D): void

    shaderData: Spheres.ShaderData
}

export namespace Spheres {
    export interface ShaderData {
        readonly positionGroup: ValueCell<TextureImage<Float32Array>>
        readonly texDim: ValueCell<Vec2>
        readonly lodLevels: ValueCell<LodLevelsValue>
        readonly sizeFactor: ValueCell<number>
        update(props?: { lodLevels: LodLevels, sizeFactor: number }): void
    }

    export function create(centers: Float32Array, groups: Float32Array, sphereCount: number, spheres?: Spheres): Spheres {
        return spheres ?
            update(centers, groups, sphereCount, spheres) :
            fromArrays(centers, groups, sphereCount);
    }

    export function createEmpty(spheres?: Spheres): Spheres {
        const cb = spheres ? spheres.centerBuffer.ref.value : new Float32Array(0);
        const gb = spheres ? spheres.groupBuffer.ref.value : new Float32Array(0);
        return create(cb, gb, 0, spheres);
    }

    function hashCode(spheres: Spheres) {
        return hashFnv32a([
            spheres.sphereCount,
            spheres.centerBuffer.ref.version,
            spheres.groupBuffer.ref.version
        ]);
    }

    function fromArrays(centers: Float32Array, groups: Float32Array, sphereCount: number): Spheres {
        const boundingSphere = Sphere3D();
        let groupMapping: GroupMapping;

        let currentHash = -1;
        let currentGroup = -1;

        const positionGroup = ValueCell.create(createTextureImage(1, 4, Float32Array));
        const texDim = ValueCell.create(Vec2.create(0, 0));
        const lodLevels = ValueCell.create([] as LodLevelsValue);
        const sizeFactor = ValueCell.create(0);

        const spheres = {
            kind: 'spheres' as const,
            sphereCount,
            centerBuffer: ValueCell.create(centers),
            groupBuffer: ValueCell.create(groups),
            get boundingSphere() {
                const newHash = hashCode(spheres);
                if (newHash !== currentHash) {
                    const b = calculateInvariantBoundingSphere(spheres.centerBuffer.ref.value, spheres.sphereCount, 1);
                    Sphere3D.copy(boundingSphere, b);
                    currentHash = newHash;
                }
                return boundingSphere;
            },
            get groupMapping() {
                if (spheres.groupBuffer.ref.version !== currentGroup) {
                    groupMapping = createGroupMapping(spheres.groupBuffer.ref.value, spheres.sphereCount);
                    currentGroup = spheres.groupBuffer.ref.version;
                }
                return groupMapping;
            },
            setBoundingSphere(sphere: Sphere3D) {
                Sphere3D.copy(boundingSphere, sphere);
                currentHash = hashCode(spheres);
            },
            shaderData: {
                positionGroup,
                texDim,
                lodLevels,
                sizeFactor,
                update(props?: { lodLevels: LodLevels, sizeFactor: number }) {
                    const lodLevelsProp = props?.lodLevels ?? getLodLevels(lodLevels.ref.value);
                    const sizeFactorProp = props?.sizeFactor ?? sizeFactor.ref.value;

                    const strides = getStrides(lodLevelsProp, sizeFactorProp);
                    const pgt = createTextureImage(spheres.sphereCount, 4, Float32Array, positionGroup.ref.value.array);
                    const offsets = getStrideOffsetsAndSetPositionGroup(pgt, spheres.centerBuffer.ref.value, spheres.groupBuffer.ref.value, spheres.sphereCount, strides);
                    const newLodLevels = offsets ? getLodLevelsValue(lodLevelsProp, sizeFactorProp, offsets, spheres.sphereCount) : [];

                    ValueCell.update(positionGroup, pgt);
                    ValueCell.update(texDim, Vec2.set(texDim.ref.value, pgt.width, pgt.height));
                    ValueCell.update(lodLevels, newLodLevels);
                    ValueCell.update(sizeFactor, sizeFactorProp);
                }
            },
        };
        spheres.shaderData.update();
        return spheres;
    }

    function update(centers: Float32Array, groups: Float32Array, sphereCount: number, spheres: Spheres) {
        spheres.sphereCount = sphereCount;
        ValueCell.update(spheres.centerBuffer, centers);
        ValueCell.update(spheres.groupBuffer, groups);
        spheres.shaderData.update();
        return spheres;
    }

    function getStrideOffsetsAndSetPositionGroup(out: TextureImage<Float32Array>, centers: Float32Array, groups: Float32Array, count: number, strides: number[]) {
        const { array } = out;
        if (strides.length === 0) {
            for (let i = 0; i < count; ++i) {
                array[i * 4 + 0] = centers[i * 3 + 0];
                array[i * 4 + 1] = centers[i * 3 + 1];
                array[i * 4 + 2] = centers[i * 3 + 2];
                array[i * 4 + 3] = groups[i];
            }
            return;
        }

        const offsets = [0];

        let o = 0;
        for (let i = 0, il = strides.length; i < il; ++i) {
            const s = strides[i];
            for (let j = 0; j < count; ++j) {
                let handled = false;
                for (let k = 0; k < i; ++k) {
                    if (j % strides[k] === 0) {
                        handled = true;
                        break;
                    }
                }
                if (!handled && j % s === 0) {
                    array[o * 4 + 0] = centers[j * 3 + 0];
                    array[o * 4 + 1] = centers[j * 3 + 1];
                    array[o * 4 + 2] = centers[j * 3 + 2];
                    array[o * 4 + 3] = groups[j];
                    o += 1;
                }
            }
            offsets.push(o * 6);
        }

        return offsets;
    }

    type LodLevels = {
        minDistance: number
        maxDistance: number
        overlap: number
        stride: number
        scaleBias: number
    }[]

    function areLodLevelsEqual(a: LodLevels, b: LodLevels) {
        if (a.length !== b.length) return false;
        for (let i = 0, il = a.length; i < il; ++i) {
            if (a[i].maxDistance !== b[i].maxDistance) return false;
            if (a[i].minDistance !== b[i].minDistance) return false;
            if (a[i].overlap !== b[i].overlap) return false;
            if (a[i].stride !== b[i].stride) return false;
            if (a[i].scaleBias !== b[i].scaleBias) return false;
        }
        return true;
    }

    type LodLevelsValue = [minDistance: number, maxDistance: number, overlap: number, count: number, scale: number, stride: number, scaleBias: number][];

    function getLodLevelsValue(prop: LodLevels, sizeFactor: number, offsets: number[], count: number): LodLevelsValue {
        return prop.map((l, i) => {
            const stride = getAdjustedStride(l, sizeFactor);
            return [
                l.minDistance,
                l.maxDistance,
                l.overlap,
                offsets[offsets.length - 1 - i],
                Math.pow(Math.min(count, stride), 1 / l.scaleBias),
                l.stride,
                l.scaleBias,
            ];
        });
    }

    function getLodLevels(lodLevelsValue: LodLevelsValue): LodLevels {
        return lodLevelsValue.map(l => ({
            minDistance: l[0],
            maxDistance: l[1],
            overlap: l[2],
            stride: l[5],
            scaleBias: l[6],
        }));
    }

    function getAdjustedStride(lodLevel: LodLevels[0], sizeFactor: number) {
        return Math.max(1, Math.round(lodLevel.stride / Math.pow(sizeFactor, lodLevel.scaleBias)));
    }

    function getStrides(lodLevels: LodLevels, sizeFactor: number) {
        return lodLevels.map(l => getAdjustedStride(l, sizeFactor)).reverse();
    }

    export const Params = {
        ...BaseGeometry.Params,
        sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
        doubleSided: PD.Boolean(false, BaseGeometry.CustomQualityParamInfo),
        ignoreLight: PD.Boolean(false, BaseGeometry.ShadingCategory),
        celShaded: PD.Boolean(false, BaseGeometry.ShadingCategory),
        xrayShaded: PD.Select<boolean | 'inverted'>(false, [[false, 'Off'], [true, 'On'], ['inverted', 'Inverted']], BaseGeometry.ShadingCategory),
        transparentBackfaces: PD.Select('off', PD.arrayToOptions(['off', 'on', 'opaque'] as const), BaseGeometry.ShadingCategory),
        solidInterior: PD.Boolean(true, BaseGeometry.ShadingCategory),
        clipPrimitive: PD.Boolean(false, { ...BaseGeometry.ShadingCategory, description: 'Clip whole sphere instead of cutting it.' }),
        approximate: PD.Boolean(false, { ...BaseGeometry.ShadingCategory, description: 'Faster rendering, but has artifacts.' }),
        alphaThickness: PD.Numeric(0, { min: 0, max: 20, step: 1 }, { ...BaseGeometry.ShadingCategory, description: 'If not zero, adjusts alpha for radius.' }),
        bumpFrequency: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
        bumpAmplitude: PD.Numeric(1, { min: 0, max: 5, step: 0.1 }, BaseGeometry.ShadingCategory),
        lodLevels: PD.ObjectList({
            minDistance: PD.Numeric(0),
            maxDistance: PD.Numeric(0),
            overlap: PD.Numeric(0),
            stride: PD.Numeric(0),
            scaleBias: PD.Numeric(3, { min: 0.1, max: 10, step: 0.1 }),
        }, o => `${o.stride}`, {
            ...BaseGeometry.CullingLodCategory,
            defaultValue: [] as LodLevels
        })
    };
    export type Params = typeof Params

    export const Utils: GeometryUtils<Spheres, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState,
        updateRenderableState,
        createPositionIterator
    };

    function createPositionIterator(spheres: Spheres, transform: TransformData): LocationIterator {
        const groupCount = spheres.sphereCount;
        const instanceCount = transform.instanceCount.ref.value;
        const location = PositionLocation();
        const p = location.position;
        const v = spheres.centerBuffer.ref.value;
        const m = transform.aTransform.ref.value;
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            if (instanceIndex < 0) {
                Vec3.fromArray(p, v, groupIndex * 3);
            } else {
                Vec3.transformMat4Offset(p, v, m, 0, groupIndex * 3, instanceIndex * 16);
            }
            return location;
        };
        return LocationIterator(groupCount, instanceCount, 1, getLocation);
    }

    function createValues(spheres: Spheres, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): SpheresValues {
        const { instanceCount, groupCount } = locationIt;
        const positionIt = createPositionIterator(spheres, transform);

        const color = createColors(locationIt, positionIt, theme.color);
        const size = createSizes(locationIt, theme.size);
        const marker = props.instanceGranularity
            ? createMarkers(instanceCount, 'instance')
            : createMarkers(instanceCount * groupCount, 'groupInstance');
        const overpaint = createEmptyOverpaint();
        const transparency = createEmptyTransparency();
        const emissive = createEmptyEmissive();
        const material = createEmptySubstance();
        const clipping = createEmptyClipping();

        const counts = { drawCount: spheres.sphereCount * 2 * 3, vertexCount: spheres.sphereCount * 6, groupCount, instanceCount };

        const padding = spheres.boundingSphere.radius ? getMaxSize(size) * props.sizeFactor : 0;
        const invariantBoundingSphere = Sphere3D.expand(Sphere3D(), spheres.boundingSphere, padding);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform.aTransform.ref.value, instanceCount, 0);

        spheres.shaderData.update({ lodLevels: props.lodLevels, sizeFactor: props.sizeFactor });

        return {
            dGeometryType: ValueCell.create('spheres'),

            uTexDim: spheres.shaderData.texDim,
            tPositionGroup: spheres.shaderData.positionGroup,

            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),
            ...color,
            ...size,
            ...marker,
            ...overpaint,
            ...transparency,
            ...emissive,
            ...material,
            ...clipping,
            ...transform,

            padding: ValueCell.create(padding),

            ...BaseGeometry.createValues(props, counts),
            uSizeFactor: spheres.shaderData.sizeFactor,
            uDoubleSided: ValueCell.create(props.doubleSided),
            dIgnoreLight: ValueCell.create(props.ignoreLight),
            dCelShaded: ValueCell.create(props.celShaded),
            dXrayShaded: ValueCell.create(props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off'),
            dTransparentBackfaces: ValueCell.create(props.transparentBackfaces),
            dSolidInterior: ValueCell.create(props.solidInterior),
            dClipPrimitive: ValueCell.create(props.clipPrimitive),
            dApproximate: ValueCell.create(props.approximate),
            uAlphaThickness: ValueCell.create(props.alphaThickness),
            uBumpFrequency: ValueCell.create(props.bumpFrequency),
            uBumpAmplitude: ValueCell.create(props.bumpAmplitude),

            lodLevels: spheres.shaderData.lodLevels,
            centerBuffer: spheres.centerBuffer,
            groupBuffer: spheres.groupBuffer,
        };
    }

    function createValuesSimple(spheres: Spheres, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
        const p = { ...PD.getDefaultValues(Params), ...props };
        return createValues(spheres, s.transform, s.locationIterator, s.theme, p);
    }

    function updateValues(values: SpheresValues, props: PD.Values<Params>) {
        BaseGeometry.updateValues(values, props);
        ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor);
        ValueCell.updateIfChanged(values.uDoubleSided, props.doubleSided);
        ValueCell.updateIfChanged(values.dIgnoreLight, props.ignoreLight);
        ValueCell.updateIfChanged(values.dCelShaded, props.celShaded);
        ValueCell.updateIfChanged(values.dXrayShaded, props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off');
        ValueCell.updateIfChanged(values.dTransparentBackfaces, props.transparentBackfaces);
        ValueCell.updateIfChanged(values.dSolidInterior, props.solidInterior);
        ValueCell.updateIfChanged(values.dClipPrimitive, props.clipPrimitive);
        ValueCell.updateIfChanged(values.dApproximate, props.approximate);
        ValueCell.updateIfChanged(values.uAlphaThickness, props.alphaThickness);
        ValueCell.updateIfChanged(values.uBumpFrequency, props.bumpFrequency);
        ValueCell.updateIfChanged(values.uBumpAmplitude, props.bumpAmplitude);

        const lodLevels = getLodLevels(values.lodLevels.ref.value as LodLevelsValue);
        if (!areLodLevelsEqual(props.lodLevels, lodLevels)) {
            const count = values.uVertexCount.ref.value / 6;
            const strides = getStrides(props.lodLevels, props.sizeFactor);
            const offsets = getStrideOffsetsAndSetPositionGroup(values.tPositionGroup.ref.value, values.centerBuffer.ref.value, values.groupBuffer.ref.value, count, strides);
            const lodLevels = offsets ? getLodLevelsValue(props.lodLevels, props.sizeFactor, offsets, count) : [];
            ValueCell.update(values.tPositionGroup, values.tPositionGroup.ref.value);
            ValueCell.update(values.lodLevels, lodLevels);
        }
    }

    function updateBoundingSphere(values: SpheresValues, spheres: Spheres) {
        const padding = spheres.boundingSphere.radius
            ? getMaxSize(values) * values.uSizeFactor.ref.value
            : 0;
        const invariantBoundingSphere = Sphere3D.expand(Sphere3D(), spheres.boundingSphere, padding);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, values.aTransform.ref.value, values.instanceCount.ref.value, 0);

        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere);
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
            ValueCell.update(values.uInvariantBoundingSphere, Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere));
        }
        ValueCell.update(values.padding, padding);
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