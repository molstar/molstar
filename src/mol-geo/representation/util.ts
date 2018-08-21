/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { BaseValues } from 'mol-gl/renderable/schema';
import { MeshValues, RenderableState } from 'mol-gl/renderable';
import { defaults } from 'mol-util';
import { Structure } from 'mol-model/structure';
import { fillSerial } from 'mol-util/array';
import { Mat4 } from 'mol-math/linear-algebra';

export const DefaultBaseProps = {
    alpha: 1,
    visible: true,
    depthMask: true,
    useFog: true,
    quality: 'auto' as VisualQuality
}
export type BaseProps = typeof DefaultBaseProps

export const DefaultMeshProps = {
    ...DefaultBaseProps,
    doubleSided: false,
    flipSided: false,
    flatShaded: false,
}
export type MeshProps = typeof DefaultMeshProps

const identityTransform = new Float32Array(16)
Mat4.toArray(Mat4.identity(), identityTransform, 0)
export function createIdentityTransform(transforms?: ValueCell<Float32Array>) {
    return transforms ? ValueCell.update(transforms, identityTransform) : ValueCell.create(identityTransform)
}

type Counts = { drawCount: number, groupCount: number, instanceCount: number }

export function createBaseValues(props: Required<BaseProps>, counts: Counts) {
    return {
        uAlpha: ValueCell.create(props.alpha),
        uInstanceCount: ValueCell.create(counts.instanceCount),
        uGroupCount: ValueCell.create(counts.groupCount),
        aInstance: ValueCell.create(fillSerial(new Float32Array(counts.instanceCount))),
        drawCount: ValueCell.create(counts.drawCount),
        instanceCount: ValueCell.create(counts.instanceCount),
    }
}

export function createMeshValues(props: Required<MeshProps>, counts: Counts) {
    return {
        ...createBaseValues(props, counts),
        dDoubleSided: ValueCell.create(props.doubleSided),
        dFlatShaded: ValueCell.create(props.flatShaded),
        dFlipSided: ValueCell.create(props.flipSided),
        dUseFog: ValueCell.create(props.useFog),
    }
}

export function createRenderableState(props: Required<BaseProps>): RenderableState {
    return {
        visible: props.visible,
        depthMask: props.depthMask
    }
}

export function updateBaseValues(values: BaseValues, props: Required<BaseProps>) {
    ValueCell.updateIfChanged(values.uAlpha, props.alpha)
}

export function updateMeshValues(values: MeshValues, props: Required<MeshProps>) {
    updateBaseValues(values, props)
    ValueCell.updateIfChanged(values.dDoubleSided, props.doubleSided)
    ValueCell.updateIfChanged(values.dFlatShaded, props.flatShaded)
    ValueCell.updateIfChanged(values.dFlipSided, props.flipSided)
    ValueCell.updateIfChanged(values.dUseFog, props.useFog)
}

export function updateRenderableState(state: RenderableState, props: Required<BaseProps>) {
    state.visible = props.visible
    state.depthMask = props.depthMask
}

export type VisualQuality = 'custom' | 'auto' | 'highest' | 'high' | 'medium' | 'low' | 'lowest'

interface QualityProps {
    quality: VisualQuality
    detail: number
    radialSegments: number
}

export function getQualityProps(props: Partial<QualityProps>, structure: Structure) {
    let quality = defaults(props.quality, 'auto' as VisualQuality)
    let detail = 1
    let radialSegments = 12

    if (quality === 'auto') {
        const score = structure.elementCount
        if (score > 500_000) {
            quality = 'lowest'
        } else if (score > 100_000) {
            quality = 'low'
        } else if (score > 30_000) {
            quality = 'medium'
        } else {
            quality = 'high'
        }
    }

    switch (quality) {
        case 'highest':
            detail = 2
            radialSegments = 36
            break
        case 'high':
            detail = 1
            radialSegments = 24
            break
        case 'medium':
            detail = 0
            radialSegments = 12
            break
        case 'low':
            detail = 0
            radialSegments = 5
            break
        case 'lowest':
            detail = 0
            radialSegments = 3
            break
        case 'custom':
            detail = defaults(props.detail, 1)
            radialSegments = defaults(props.radialSegments, 12)
            break
    }

    return {
        detail,
        radialSegments
    }
}