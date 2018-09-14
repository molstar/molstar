/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { LocationIterator } from '../../../../util/location-iterator';
import { Mesh } from '../../../../geometry/mesh/mesh';
import { StructureProps } from '../..';
import { createMeshRenderObject, createPointRenderObject } from 'mol-gl/render-object';
import { RuntimeContext } from 'mol-task';
import { PointProps } from 'mol-geo/representation/structure/representation/point';
import { TransformData, createIdentityTransform, createTransform } from '../../../../geometry/transform-data';
import { Point } from '../../../../geometry/point/point';
import { createRenderableState } from '../../../../geometry/geometry';
import { Mat4 } from 'mol-math/linear-algebra';

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

// point

type StructurePointProps = PointProps & StructureProps

export async function createUnitsPointRenderObject(ctx: RuntimeContext, group: Unit.SymmetryGroup, point: Point, locationIt: LocationIterator, props: StructurePointProps) {
    const transform = createUnitsTransform(group)
    const values = await Point.createValues(ctx, point, transform, locationIt, props)
    const state = createRenderableState(props)
    return createPointRenderObject(values, state)
}