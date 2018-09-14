/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mesh } from './mesh/mesh';
import { Point } from './point/point';
import { PointValues, MeshValues, RenderableState } from 'mol-gl/renderable';
import { MeshRenderObject, PointRenderObject } from 'mol-gl/render-object';
import { ValueCell } from 'mol-util';
import { BaseValues } from 'mol-gl/renderable/schema';

export type GeometryKindType = { 'mesh': Mesh, 'point': Point }
export type GeometryKind = keyof GeometryKindType
export type Geometry = Helpers.ValueOf<GeometryKindType>

export type GeometryDef = {
    'mesh': { 'props': MeshProps, 'values': MeshValues, 'renderObject': MeshRenderObject }
    'point': { 'props': PointProps, 'values': PointValues, 'renderObject': PointRenderObject }
}

export namespace Geometry {
    export function getDrawCount(geometry: Geometry) {
        switch (geometry.kind) {
            case 'mesh': return geometry.triangleCount * 3
            case 'point': return geometry.vertexCount
        }
    }
}

//

export const VisualQualityInfo = {
    'custom': {},
    'auto': {},
    'highest': {},
    'high': {},
    'medium': {},
    'low': {},
    'lowest': {},
}
export type VisualQuality = keyof typeof VisualQualityInfo
export const VisualQualityNames = Object.keys(VisualQualityInfo)

//

export const DefaultBaseProps = {
    alpha: 1,
    visible: true,
    depthMask: true,
    useFog: false,
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

export const DefaultPointProps = {
    ...DefaultBaseProps,
    pointSizeAttenuation: true
}
export type PointProps = typeof DefaultPointProps

type Counts = { drawCount: number, groupCount: number, instanceCount: number }

export function createBaseValues(props: BaseProps, counts: Counts) {
    return {
        uAlpha: ValueCell.create(props.alpha),
        uGroupCount: ValueCell.create(counts.groupCount),
        drawCount: ValueCell.create(counts.drawCount),
        dUseFog: ValueCell.create(props.useFog),
    }
}

export function createMeshValues(props: MeshProps, counts: Counts) {
    return {
        ...createBaseValues(props, counts),
        dDoubleSided: ValueCell.create(props.doubleSided),
        dFlatShaded: ValueCell.create(props.flatShaded),
        dFlipSided: ValueCell.create(props.flipSided),
    }
}

export function createPointValues(props: PointProps, counts: Counts) {
    return {
        ...createBaseValues(props, counts),
        dPointSizeAttenuation: ValueCell.create(props.pointSizeAttenuation),
    }
}

export function createRenderableState(props: BaseProps): RenderableState {
    return {
        visible: props.visible,
        depthMask: props.depthMask
    }
}

export function updateBaseValues(values: BaseValues, props: BaseProps) {
    ValueCell.updateIfChanged(values.uAlpha, props.alpha)
    ValueCell.updateIfChanged(values.dUseFog, props.useFog)
}

export function updateMeshValues(values: MeshValues, props: MeshProps) {
    updateBaseValues(values, props)
    ValueCell.updateIfChanged(values.dDoubleSided, props.doubleSided)
    ValueCell.updateIfChanged(values.dFlatShaded, props.flatShaded)
    ValueCell.updateIfChanged(values.dFlipSided, props.flipSided)
}

export function updatePointValues(values: PointValues, props: PointProps) {
    updateBaseValues(values, props)
    ValueCell.updateIfChanged(values.dPointSizeAttenuation, props.pointSizeAttenuation)
}

export function updateRenderableState(state: RenderableState, props: BaseProps) {
    state.visible = props.visible
    state.depthMask = props.depthMask
}