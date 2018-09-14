/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit, Structure } from 'mol-model/structure';
import { createUniformColor, ColorData, createGroupColor, createGroupInstanceColor, createInstanceColor, ColorType } from '../../../../util/color-data';
import { createUniformSize, SizeData, createGroupSize, createGroupInstanceSize, createInstanceSize, SizeType } from '../../../../util/size-data';
import { LocationIterator } from '../../../../util/location-iterator';
import { Mesh } from '../../../../geometry/mesh/mesh';
import { MeshValues, PointValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../../util/mesh-data';
import { StructureProps } from '../..';
import { createMarkers } from '../../../../util/marker-data';
import { createMeshRenderObject, createPointRenderObject } from 'mol-gl/render-object';
import { ColorThemeProps, ColorTheme } from 'mol-view/theme/color';
import { SizeThemeProps, SizeTheme } from 'mol-view/theme/size';
import { RuntimeContext } from 'mol-task';
import { PointProps } from 'mol-geo/representation/structure/representation/point';
import { TransformData, createIdentityTransform, createTransforms } from '../../../../util/transform-data';
import { Point } from '../../../../geometry/point/point';
import { getPointData } from '../../../../util/point-data';
import { MeshProps, createMeshValues, createRenderableState, createPointValues } from '../../../../geometry/geometry';

function getGranularity(locationIt: LocationIterator, granularity: ColorType | SizeType) {
    // Always use 'group' granularity for 'complex' location iterators,
    // i.e. for which an instance may include multiple units
    return granularity === 'instance' && locationIt.isComplex ? 'group' : granularity
}

export function createColors(ctx: RuntimeContext, locationIt: LocationIterator, props: ColorThemeProps, colorData?: ColorData): Promise<ColorData> {
    const colorTheme = ColorTheme(props)
    switch (getGranularity(locationIt, colorTheme.granularity)) {
        case 'uniform': return createUniformColor(ctx, locationIt, colorTheme.color, colorData)
        case 'group': return createGroupColor(ctx, locationIt, colorTheme.color, colorData)
        case 'groupInstance': return createGroupInstanceColor(ctx, locationIt, colorTheme.color, colorData)
        case 'instance': return createInstanceColor(ctx, locationIt, colorTheme.color, colorData)
    }
}

export async function createSizes(ctx: RuntimeContext, locationIt: LocationIterator, props: SizeThemeProps, sizeData?: SizeData): Promise<SizeData> {
    const sizeTheme = SizeTheme(props)
    switch (getGranularity(locationIt, sizeTheme.granularity)) {
        case 'uniform': return createUniformSize(ctx, locationIt, sizeTheme.size, sizeData)
        case 'group': return createGroupSize(ctx, locationIt, sizeTheme.size, sizeData)
        case 'groupInstance': return createGroupInstanceSize(ctx, locationIt, sizeTheme.size, sizeData)
        case 'instance': return createInstanceSize(ctx, locationIt, sizeTheme.size, sizeData)
    }
}

// mesh

type StructureMeshProps = MeshProps & StructureProps

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

type StructurePointProps = PointProps & StructureProps

async function _createPointValues(ctx: RuntimeContext, transforms: TransformData, point: Point, locationIt: LocationIterator, props: StructurePointProps): Promise<PointValues> {
    const { instanceCount, groupCount } = locationIt
    const color = await createColors(ctx, locationIt, props.colorTheme)
    const size = await createSizes(ctx, locationIt, props.sizeTheme)
    const marker = createMarkers(instanceCount * groupCount)

    const counts = { drawCount: point.vertexCount, groupCount, instanceCount }

    return {
        ...getPointData(point),
        ...color,
        ...size,
        ...marker,
        ...transforms,
        ...createPointValues(props, counts)
    }
}

export async function createUnitsPointValues(ctx: RuntimeContext, group: Unit.SymmetryGroup, point: Point, locationIt: LocationIterator, props: StructurePointProps): Promise<PointValues> {
    const transforms = createTransforms(group)
    return _createPointValues(ctx, transforms, point, locationIt, props)
}

export async function createUnitsPointRenderObject(ctx: RuntimeContext, group: Unit.SymmetryGroup, point: Point, locationIt: LocationIterator, props: StructurePointProps) {
    const values = await createUnitsPointValues(ctx, group, point, locationIt, props)
    const state = createRenderableState(props)
    return createPointRenderObject(values, state)
}