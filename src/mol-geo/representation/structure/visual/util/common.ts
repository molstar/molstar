/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { LocationIterator } from '../../../../util/location-iterator';
import { Mesh } from '../../../../geometry/mesh/mesh';
import { StructureProps } from '../..';
import { createMeshRenderObject, createPointsRenderObject, createLinesRenderObject, createDirectVolumeRenderObject } from 'mol-gl/render-object';
import { RuntimeContext } from 'mol-task';
import { TransformData, createIdentityTransform, createTransform } from '../../../../geometry/transform-data';
import { Points } from '../../../../geometry/points/points';
import { createRenderableState } from '../../../../geometry/geometry';
import { Mat4 } from 'mol-math/linear-algebra';
import { Lines } from '../../../../geometry/lines/lines';
import { DirectVolume } from '../../../../geometry/direct-volume/direct-volume';

export function createUnitsTransform({ units }: Unit.SymmetryGroup, transformData?: TransformData) {
    const unitCount = units.length
    const n = unitCount * 16
    const array = transformData && transformData.aTransform.ref.value.length >= n ? transformData.aTransform.ref.value : new Float32Array(n)
    for (let i = 0; i < unitCount; i++) {
        Mat4.toArray(units[i].conformation.operator.matrix, array, i * 16)
    }
    return createTransform(array, unitCount, transformData)
}

// mesh

type StructureMeshProps = Mesh.Props & StructureProps

export async function createComplexMeshRenderObject(ctx: RuntimeContext, structure: Structure, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps) {
    const transform = createIdentityTransform()
    const values = await Mesh.createValues(ctx, mesh, transform, locationIt, props)
    const state = createRenderableState(props)
    return createMeshRenderObject(values, state)
}

export async function createUnitsMeshRenderObject(ctx: RuntimeContext, group: Unit.SymmetryGroup, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps) {
    const transform = createUnitsTransform(group)
    const values = await Mesh.createValues(ctx, mesh, transform, locationIt, props)
    const state = createRenderableState(props)
    return createMeshRenderObject(values, state)
}

// points

type StructurePointsProps = Points.Props & StructureProps

export async function createUnitsPointsRenderObject(ctx: RuntimeContext, group: Unit.SymmetryGroup, points: Points, locationIt: LocationIterator, props: StructurePointsProps) {
    const transform = createUnitsTransform(group)
    const values = await Points.createValues(ctx, points, transform, locationIt, props)
    const state = createRenderableState(props)
    return createPointsRenderObject(values, state)
}

// lines

type StructureLinesProps = Lines.Props & StructureProps

export async function createUnitsLinesRenderObject(ctx: RuntimeContext, group: Unit.SymmetryGroup, lines: Lines, locationIt: LocationIterator, props: StructureLinesProps) {
    const transform = createUnitsTransform(group)
    const values = await Lines.createValues(ctx, lines, transform, locationIt, props)
    console.log('values', values)
    const state = createRenderableState(props)
    return createLinesRenderObject(values, state)
}

// direct-volume

type StructureDirectVolumeProps = DirectVolume.Props & StructureProps

export async function createUnitsDirectVolumeRenderObject(ctx: RuntimeContext, group: Unit.SymmetryGroup, directVolume: DirectVolume, locationIt: LocationIterator, props: StructureDirectVolumeProps) {
    // TODO transform support
    // const transform = createUnitsTransform(group)
    const values = await DirectVolume.createValues(ctx, directVolume, props)
    const state = createRenderableState(props)
    return createDirectVolumeRenderObject(values, state)
}