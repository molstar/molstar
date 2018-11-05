/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { StructureProps } from '../../index';
import { createMeshRenderObject, createPointsRenderObject, createLinesRenderObject, createDirectVolumeRenderObject } from 'mol-gl/render-object';
import { Mat4 } from 'mol-math/linear-algebra';
import { TransformData, createTransform, createIdentityTransform } from 'mol-geo/geometry/transform-data';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { createRenderableState } from 'mol-geo/geometry/geometry';
import { Points } from 'mol-geo/geometry/points/points';
import { Lines } from 'mol-geo/geometry/lines/lines';
import { DirectVolume } from 'mol-geo/geometry/direct-volume/direct-volume';
import { VisualContext } from 'mol-repr';

export function createUnitsTransform({ units }: Unit.SymmetryGroup, transformData?: TransformData) {
    const unitCount = units.length
    const n = unitCount * 16
    const array = transformData && transformData.aTransform.ref.value.length >= n ? transformData.aTransform.ref.value : new Float32Array(n)
    for (let i = 0; i < unitCount; i++) {
        Mat4.toArray(units[i].conformation.operator.matrix, array, i * 16)
    }
    return createTransform(array, unitCount, transformData)
}

export const UnitKindInfo = {
    'atomic': {},
    'spheres': {},
    'gaussians': {},
}
export type UnitKind = keyof typeof UnitKindInfo
export const UnitKindNames = Object.keys(UnitKindInfo)
export const UnitKindOptions = UnitKindNames.map(n => [n, n] as [UnitKind, string])

export function includesUnitKind(unitKinds: UnitKind[], unit: Unit) {
    for (let i = 0, il = unitKinds.length; i < il; ++i) {
        if (Unit.isAtomic(unit) && unitKinds[i] === 'atomic') return true
        if (Unit.isSpheres(unit) && unitKinds[i] === 'spheres') return true
        if (Unit.isGaussians(unit) && unitKinds[i] === 'gaussians') return true
    }
    return false
}

// mesh

type StructureMeshProps = Mesh.Props & StructureProps

export async function createComplexMeshRenderObject(ctx: VisualContext, structure: Structure, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps) {
    const transform = createIdentityTransform()
    const values = await Mesh.createValues(ctx.runtime, mesh, transform, locationIt, props)
    const state = createRenderableState(props)
    return createMeshRenderObject(values, state)
}

export async function createUnitsMeshRenderObject(ctx: VisualContext, group: Unit.SymmetryGroup, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps) {
    const transform = createUnitsTransform(group)
    const values = await Mesh.createValues(ctx.runtime, mesh, transform, locationIt, props)
    const state = createRenderableState(props)
    return createMeshRenderObject(values, state)
}

// points

type StructurePointsProps = Points.Props & StructureProps

export async function createUnitsPointsRenderObject(ctx: VisualContext, group: Unit.SymmetryGroup, points: Points, locationIt: LocationIterator, props: StructurePointsProps) {
    const transform = createUnitsTransform(group)
    const values = await Points.createValues(ctx.runtime, points, transform, locationIt, props)
    const state = createRenderableState(props)
    return createPointsRenderObject(values, state)
}

// lines

type StructureLinesProps = Lines.Props & StructureProps

export async function createUnitsLinesRenderObject(ctx: VisualContext, group: Unit.SymmetryGroup, lines: Lines, locationIt: LocationIterator, props: StructureLinesProps) {
    const transform = createUnitsTransform(group)
    const values = await Lines.createValues(ctx.runtime, lines, transform, locationIt, props)
    const state = createRenderableState(props)
    return createLinesRenderObject(values, state)
}

// direct-volume

type StructureDirectVolumeProps = DirectVolume.Props & StructureProps

export async function createUnitsDirectVolumeRenderObject(ctx: VisualContext, group: Unit.SymmetryGroup, directVolume: DirectVolume, locationIt: LocationIterator, props: StructureDirectVolumeProps) {
    const transform = createUnitsTransform(group)
    const values = await DirectVolume.createValues(ctx.runtime, directVolume, transform, locationIt, props)
    const state = createRenderableState(props)
    return createDirectVolumeRenderObject(values, state)
}