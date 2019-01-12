/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import { Geometry } from '../geometry';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { TransformData, createIdentityTransform } from '../transform-data';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { Theme } from 'mol-theme/theme';
import { SpheresValues } from 'mol-gl/renderable/spheres';
import { createColors, createValueColor } from '../color-data';
import { createMarkers } from '../marker-data';
import { calculateBoundingSphere } from 'mol-gl/renderable/util';
import { ColorNames } from 'mol-util/color/tables';
import { Sphere3D } from 'mol-math/geometry';
import { createSizes, createValueSize } from '../size-data';

/** Spheres */
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
}

export namespace Spheres {
    export function createEmpty(spheres?: Spheres): Spheres {
        const cb = spheres ? spheres.centerBuffer.ref.value : new Float32Array(0)
        const mb = spheres ? spheres.mappingBuffer.ref.value : new Float32Array(0)
        const ib = spheres ? spheres.indexBuffer.ref.value : new Uint32Array(0)
        const gb = spheres ? spheres.groupBuffer.ref.value : new Float32Array(0)
        return {
            kind: 'spheres',
            sphereCount: 0,
            centerBuffer: spheres ? ValueCell.update(spheres.centerBuffer, cb) : ValueCell.create(cb),
            mappingBuffer: spheres ? ValueCell.update(spheres.mappingBuffer, mb) : ValueCell.create(mb),
            indexBuffer: spheres ? ValueCell.update(spheres.indexBuffer, ib) : ValueCell.create(ib),
            groupBuffer: spheres ? ValueCell.update(spheres.groupBuffer, gb) : ValueCell.create(gb)
        }
    }

    export const Params = {
        ...Geometry.Params,
        sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
        doubleSided: PD.Boolean(false),
    }
    export type Params = typeof Params

    export function createValues(spheres: Spheres, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): SpheresValues {
        const { instanceCount, groupCount } = locationIt
        if (instanceCount !== transform.instanceCount.ref.value) {
            throw new Error('instanceCount values in TransformData and LocationIterator differ')
        }

        const color = createColors(locationIt, theme.color)
        const size = createSizes(locationIt, theme.size)
        const marker = createMarkers(instanceCount * groupCount)

        const counts = { drawCount: spheres.sphereCount * 2 * 3, groupCount, instanceCount }

        const padding = 5 // TODO get max sphere size
        const { boundingSphere, invariantBoundingSphere } = calculateBoundingSphere(
            spheres.centerBuffer.ref.value, spheres.sphereCount * 4,
            transform.aTransform.ref.value, instanceCount, padding
        )

        return {
            aPosition: spheres.centerBuffer,
            aMapping: spheres.mappingBuffer,
            aGroup: spheres.groupBuffer,
            elements: spheres.indexBuffer,
            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            ...color,
            ...size,
            ...marker,
            ...transform,

            padding: ValueCell.create(padding),

            ...Geometry.createValues(props, counts),
            uSizeFactor: ValueCell.create(props.sizeFactor),
            dDoubleSided: ValueCell.create(props.doubleSided),
        }
    }

    export function createValuesSimple(spheres: Spheres, props: Partial<PD.Values<Params>>, colorValue = ColorNames.grey, sizeValue = 1, transform?: TransformData): SpheresValues {
        const p = { ...PD.getDefaultValues(Params), ...props }
        if (!transform) transform = createIdentityTransform()
        const instanceCount = transform.instanceCount.ref.value
        const groupCount = 1
        const color = createValueColor(colorValue)
        const size = createValueSize(sizeValue)
        const marker = createMarkers(instanceCount * groupCount)

        const counts = { drawCount: spheres.sphereCount * 2 * 3, groupCount, instanceCount }

        const { boundingSphere, invariantBoundingSphere } = calculateBoundingSphere(
            spheres.centerBuffer.ref.value, spheres.sphereCount * 4,
            transform.aTransform.ref.value, instanceCount, sizeValue
        )

        return {
            aPosition: spheres.centerBuffer,
            aMapping: spheres.mappingBuffer,
            aGroup: spheres.groupBuffer,
            elements: spheres.indexBuffer,
            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            ...color,
            ...size,
            ...marker,
            ...transform,

            padding: ValueCell.create(sizeValue),

            ...Geometry.createValues(p, counts),
            uSizeFactor: ValueCell.create(p.sizeFactor),
            dDoubleSided: ValueCell.create(p.doubleSided),
        }
    }

    export function updateValues(values: SpheresValues, props: PD.Values<Params>) {
        Geometry.updateValues(values, props)
        ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor)
        ValueCell.updateIfChanged(values.dDoubleSided, props.doubleSided)
    }

    export function updateBoundingSphere(values: SpheresValues, spheres: Spheres) {
        const { boundingSphere, invariantBoundingSphere } = calculateBoundingSphere(
            values.aPosition.ref.value, spheres.sphereCount * 4,
            values.aTransform.ref.value, values.instanceCount.ref.value, values.padding.ref.value
        )
        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere)
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere)
        }
    }
}