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
import { SizeProps } from 'mol-geo/geometry/size-data';
import { ColorProps } from 'mol-geo/geometry/color-data';

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

export function sizeChanged(oldProps: SizeProps, newProps: SizeProps) {
    return (
        oldProps.sizeTheme !== newProps.sizeTheme ||
        oldProps.sizeValue !== newProps.sizeValue ||
        oldProps.sizeFactor !== newProps.sizeFactor
    )
}

export function colorChanged(oldProps: ColorProps, newProps: ColorProps) {
    return (
        oldProps.colorTheme !== newProps.colorTheme ||
        oldProps.colorValue !== newProps.colorValue
    )
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
    const state = createRenderableState(props)
    return createLinesRenderObject(values, state)
}

// direct-volume

type StructureDirectVolumeProps = DirectVolume.Props & StructureProps

export async function createUnitsDirectVolumeRenderObject(ctx: RuntimeContext, group: Unit.SymmetryGroup, directVolume: DirectVolume, locationIt: LocationIterator, props: StructureDirectVolumeProps) {
    const transform = createUnitsTransform(group)
    const values = await DirectVolume.createValues(ctx, directVolume, transform, locationIt, props)
    const state = createRenderableState(props)
    return createDirectVolumeRenderObject(values, state)
}