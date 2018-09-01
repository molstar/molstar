/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit, Structure } from 'mol-model/structure';
import { Mat4 } from 'mol-math/linear-algebra'

import { createUniformColor, ColorData, createGroupColor, createGroupInstanceColor, createInstanceColor } from '../../../../util/color-data';
import { createUniformSize, SizeData, createGroupSize, createGroupInstanceSize, createInstanceSize } from '../../../../util/size-data';
import { ValueCell } from 'mol-util';
import { LocationIterator } from '../../../../util/location-iterator';
import { Mesh } from '../../../../mesh/mesh';
import { MeshValues, PointValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../../util/mesh-data';
import { MeshProps, createMeshValues, createRenderableState, createIdentityTransform, TransformData, createPointValues } from '../../../util';
import { StructureProps } from '../..';
import { createMarkers } from '../../../../util/marker-data';
import { createMeshRenderObject, createPointRenderObject } from 'mol-gl/render-object';
import { ColorThemeProps, ColorTheme } from 'mol-view/theme/color';
import { SizeThemeProps, SizeTheme } from 'mol-view/theme/size';
import { RuntimeContext } from 'mol-task';
import { PointProps } from 'mol-geo/representation/structure/representation/point';
import { fillSerial } from 'mol-util/array';

export function createTransforms({ units }: Unit.SymmetryGroup, transformData?: TransformData) {
    const unitCount = units.length
    const n = unitCount * 16
    const array = transformData && transformData.aTransform && transformData.aTransform.ref.value.length >= n ? transformData.aTransform.ref.value : new Float32Array(n)
    for (let i = 0; i < unitCount; i++) {
        Mat4.toArray(units[i].conformation.operator.matrix, array, i * 16)
    }
    if (transformData) {
        ValueCell.update(transformData.aTransform, array)
        return transformData
    } else {
        return { aTransform: ValueCell.create(array) }
    }
}

export function createColors(ctx: RuntimeContext, locationIt: LocationIterator, props: ColorThemeProps, colorData?: ColorData): Promise<ColorData> {
    const colorTheme = ColorTheme(props)
    switch (colorTheme.kind) {
        case 'uniform': return createUniformColor(ctx, locationIt, colorTheme.color, colorData)
        case 'group': return createGroupColor(ctx, locationIt, colorTheme.color, colorData)
        case 'groupInstance': return createGroupInstanceColor(ctx, locationIt, colorTheme.color, colorData)
        case 'instance': return createInstanceColor(ctx, locationIt, colorTheme.color, colorData)
    }
}

export async function createSizes(ctx: RuntimeContext, locationIt: LocationIterator, props: SizeThemeProps, sizeData?: SizeData): Promise<SizeData> {
    const sizeTheme = SizeTheme(props)
    switch (sizeTheme.kind) {
        case 'uniform': return createUniformSize(ctx, locationIt, sizeTheme.size, sizeData)
        case 'group': return createGroupSize(ctx, locationIt, sizeTheme.size, sizeData)
        case 'groupInstance': return createGroupInstanceSize(ctx, locationIt, sizeTheme.size, sizeData)
        case 'instance': return createInstanceSize(ctx, locationIt, sizeTheme.size, sizeData)
    }
}

// mesh

type StructureMeshProps = Required<MeshProps & StructureProps>

async function _createMeshValues(ctx: RuntimeContext, transforms: TransformData, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps): Promise<MeshValues> {
    const { instanceCount, groupCount } = locationIt
    const color = await createColors(ctx, locationIt, props.colorTheme)
    const marker = createMarkers(instanceCount * groupCount)

    const counts = { drawCount: mesh.triangleCount * 3, groupCount, instanceCount }

    return {
        ...getMeshData(mesh),
        ...color,
        ...marker,
        ...transforms,
        elements: mesh.indexBuffer,
        ...createMeshValues(props, counts)
    }
}

export async function createComplexMeshValues(ctx: RuntimeContext, structure: Structure, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps): Promise<MeshValues> {
    const transforms = createIdentityTransform()
    return _createMeshValues(ctx, transforms, mesh, locationIt, props)
}

export async function createUnitsMeshValues(ctx: RuntimeContext, group: Unit.SymmetryGroup, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps): Promise<MeshValues> {
    const transforms = createTransforms(group)
    return _createMeshValues(ctx, transforms, mesh, locationIt, props)
}

export async function createComplexMeshRenderObject(ctx: RuntimeContext, structure: Structure, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps) {
    const values = await createComplexMeshValues(ctx, structure, mesh, locationIt, props)
    const state = createRenderableState(props)
    return createMeshRenderObject(values, state)
}

export async function createUnitsMeshRenderObject(ctx: RuntimeContext, group: Unit.SymmetryGroup, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps) {
    const values = await createUnitsMeshValues(ctx, group, mesh, locationIt, props)
    const state = createRenderableState(props)
    return createMeshRenderObject(values, state)
}

export async function updateComplexMeshRenderObject(ctx: RuntimeContext, structure: Structure, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps): Promise<MeshValues> {
    const transforms = createIdentityTransform()
    return _createMeshValues(ctx, transforms, mesh, locationIt, props)
}

// point

type StructurePointProps = Required<PointProps & StructureProps>

async function _createPointValues(ctx: RuntimeContext, transforms: TransformData, vertices: ValueCell<Float32Array>, locationIt: LocationIterator, props: StructurePointProps): Promise<PointValues> {
    const { instanceCount, groupCount } = locationIt
    const color = await createColors(ctx, locationIt, props.colorTheme)
    const size = await createSizes(ctx, locationIt, props.sizeTheme)
    const marker = createMarkers(instanceCount * groupCount)

    const counts = { drawCount: groupCount, groupCount, instanceCount }

    return {
        aPosition: vertices,
        aGroup: ValueCell.create(fillSerial(new Float32Array(groupCount))),
        ...color,
        ...size,
        ...marker,
        ...transforms,
        ...createPointValues(props, counts)
    }
}

export async function createUnitsPointValues(ctx: RuntimeContext, group: Unit.SymmetryGroup, vertices: ValueCell<Float32Array>, locationIt: LocationIterator, props: StructurePointProps): Promise<PointValues> {
    const transforms = createTransforms(group)
    return _createPointValues(ctx, transforms, vertices, locationIt, props)
}

export async function createUnitsPointRenderObject(ctx: RuntimeContext, group: Unit.SymmetryGroup, vertices: ValueCell<Float32Array>, locationIt: LocationIterator, props: StructurePointProps) {
    const values = await createUnitsPointValues(ctx, group, vertices, locationIt, props)
    const state = createRenderableState(props)
    return createPointRenderObject(values, state)
}