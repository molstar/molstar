/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../../mol-util';
import { GeometryUtils } from '../geometry';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { TransformData } from '../transform-data';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { Theme } from '../../../mol-theme/theme';
import { SpheresValues } from '../../../mol-gl/renderable/spheres';
import { createColors } from '../color-data';
import { createMarkers } from '../marker-data';
import { calculateInvariantBoundingSphere, calculateTransformBoundingSphere } from '../../../mol-gl/renderable/util';
import { Sphere3D } from '../../../mol-math/geometry';
import { createSizes, getMaxSize } from '../size-data';
import { Color } from '../../../mol-util/color';
import { BaseGeometry } from '../base';
import { createEmptyOverpaint } from '../overpaint-data';
import { createEmptyTransparency } from '../transparency-data';
import { hashFnv32a } from '../../../mol-data/util';
import { GroupMapping, createGroupMapping } from '../../util';
import { createEmptyClipping } from '../clipping-data';
import { Vec4 } from '../../../mol-math/linear-algebra';

export interface Spheres {
    readonly kind: 'spheres',

    /** Number of spheres */
    sphereCount: number,

    /** Center buffer as array of xyz values wrapped in a value cell */
    readonly centerBuffer: ValueCell<Float32Array>,
    /** Mapping buffer as array of xy values wrapped in a value cell */
    readonly mappingBuffer: ValueCell<Float32Array>,
    /** Index buffer as array of center index triplets wrapped in a value cell */
    readonly indexBuffer: ValueCell<Uint32Array>,
    /** Group buffer as array of group ids for each vertex wrapped in a value cell */
    readonly groupBuffer: ValueCell<Float32Array>,

    /** Bounding sphere of the spheres */
    readonly boundingSphere: Sphere3D
    /** Maps group ids to sphere indices */
    readonly groupMapping: GroupMapping

    setBoundingSphere(boundingSphere: Sphere3D): void
}

export namespace Spheres {
    export function create(centers: Float32Array, mappings: Float32Array, indices: Uint32Array, groups: Float32Array, sphereCount: number, spheres?: Spheres): Spheres {
        return spheres ?
            update(centers, mappings, indices, groups, sphereCount, spheres) :
            fromArrays(centers, mappings, indices, groups, sphereCount);
    }

    export function createEmpty(spheres?: Spheres): Spheres {
        const cb = spheres ? spheres.centerBuffer.ref.value : new Float32Array(0);
        const mb = spheres ? spheres.mappingBuffer.ref.value : new Float32Array(0);
        const ib = spheres ? spheres.indexBuffer.ref.value : new Uint32Array(0);
        const gb = spheres ? spheres.groupBuffer.ref.value : new Float32Array(0);
        return create(cb, mb, ib, gb, 0, spheres);
    }

    function hashCode(spheres: Spheres) {
        return hashFnv32a([
            spheres.sphereCount,
            spheres.centerBuffer.ref.version, spheres.mappingBuffer.ref.version,
            spheres.indexBuffer.ref.version, spheres.groupBuffer.ref.version
        ]);
    }

    function fromArrays(centers: Float32Array, mappings: Float32Array, indices: Uint32Array, groups: Float32Array, sphereCount: number): Spheres {

        const boundingSphere = Sphere3D();
        let groupMapping: GroupMapping;

        let currentHash = -1;
        let currentGroup = -1;

        const spheres = {
            kind: 'spheres' as const,
            sphereCount,
            centerBuffer: ValueCell.create(centers),
            mappingBuffer: ValueCell.create(mappings),
            indexBuffer: ValueCell.create(indices),
            groupBuffer: ValueCell.create(groups),
            get boundingSphere() {
                const newHash = hashCode(spheres);
                if (newHash !== currentHash) {
                    const b = calculateInvariantBoundingSphere(spheres.centerBuffer.ref.value, spheres.sphereCount * 4, 4);
                    Sphere3D.copy(boundingSphere, b);
                    currentHash = newHash;
                }
                return boundingSphere;
            },
            get groupMapping() {
                if (spheres.groupBuffer.ref.version !== currentGroup) {
                    groupMapping = createGroupMapping(spheres.groupBuffer.ref.value, spheres.sphereCount, 4);
                    currentGroup = spheres.groupBuffer.ref.version;
                }
                return groupMapping;
            },
            setBoundingSphere(sphere: Sphere3D) {
                Sphere3D.copy(boundingSphere, sphere);
                currentHash = hashCode(spheres);
            }
        };
        return spheres;
    }

    function update(centers: Float32Array, mappings: Float32Array, indices: Uint32Array, groups: Float32Array, sphereCount: number, spheres: Spheres) {
        spheres.sphereCount = sphereCount;
        ValueCell.update(spheres.centerBuffer, centers);
        ValueCell.update(spheres.mappingBuffer, mappings);
        ValueCell.update(spheres.indexBuffer, indices);
        ValueCell.update(spheres.groupBuffer, groups);
        return spheres;
    }

    export const Params = {
        ...BaseGeometry.Params,
        sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
        doubleSided: PD.Boolean(false, BaseGeometry.CustomQualityParamInfo),
        ignoreLight: PD.Boolean(false, BaseGeometry.ShadingCategory),
    };
    export type Params = typeof Params

    export const Utils: GeometryUtils<Spheres, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState: BaseGeometry.createRenderableState,
        updateRenderableState: BaseGeometry.updateRenderableState
    };

    function createValues(spheres: Spheres, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): SpheresValues {
        const { instanceCount, groupCount } = locationIt;
        if (instanceCount !== transform.instanceCount.ref.value) {
            throw new Error('instanceCount values in TransformData and LocationIterator differ');
        }

        const color = createColors(locationIt, theme.color);
        const size = createSizes(locationIt, theme.size);
        const marker = createMarkers(instanceCount * groupCount);
        const overpaint = createEmptyOverpaint();
        const transparency = createEmptyTransparency();
        const clipping = createEmptyClipping();

        const counts = { drawCount: spheres.sphereCount * 2 * 3, groupCount, instanceCount };

        const padding = getMaxSize(size);
        const invariantBoundingSphere = Sphere3D.expand(Sphere3D(), spheres.boundingSphere, padding);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform.aTransform.ref.value, instanceCount);

        return {
            aPosition: spheres.centerBuffer,
            aMapping: spheres.mappingBuffer,
            aGroup: spheres.groupBuffer,
            elements: spheres.indexBuffer,
            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),
            ...color,
            ...size,
            ...marker,
            ...overpaint,
            ...transparency,
            ...clipping,
            ...transform,

            padding: ValueCell.create(padding),

            ...BaseGeometry.createValues(props, counts),
            uSizeFactor: ValueCell.create(props.sizeFactor),
            dDoubleSided: ValueCell.create(props.doubleSided),
            dIgnoreLight: ValueCell.create(props.ignoreLight),
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
        ValueCell.updateIfChanged(values.dDoubleSided, props.doubleSided);
        ValueCell.updateIfChanged(values.dIgnoreLight, props.ignoreLight);
    }

    function updateBoundingSphere(values: SpheresValues, spheres: Spheres) {
        const padding = getMaxSize(values);
        const invariantBoundingSphere = Sphere3D.expand(Sphere3D(), spheres.boundingSphere, padding);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, values.aTransform.ref.value, values.instanceCount.ref.value);

        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere);
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
            ValueCell.update(values.uInvariantBoundingSphere, Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere));
        }
        ValueCell.update(values.padding, padding);
    }
}