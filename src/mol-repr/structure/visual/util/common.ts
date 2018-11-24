/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure, ElementIndex, StructureElement } from 'mol-model/structure';
import { createMeshRenderObject, createPointsRenderObject, createLinesRenderObject, createDirectVolumeRenderObject } from 'mol-gl/render-object';
import { Mat4 } from 'mol-math/linear-algebra';
import { TransformData, createTransform, createIdentityTransform } from 'mol-geo/geometry/transform-data';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { createRenderableState } from 'mol-geo/geometry/geometry';
import { Points } from 'mol-geo/geometry/points/points';
import { Lines } from 'mol-geo/geometry/lines/lines';
import { DirectVolume } from 'mol-geo/geometry/direct-volume/direct-volume';
import { Theme } from 'mol-theme/theme';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { StructureMeshParams, StructurePointsParams, StructureLinesParams, StructureDirectVolumeParams } from 'mol-repr/structure/representation';
import { OrderedSet, SortedArray } from 'mol-data/int';
import { EmptyLoci, Loci } from 'mol-model/loci';

/** Return a Loci for the elements of a whole residue the elementIndex belongs to. */
export function getResidueLoci(structure: Structure, unit: Unit, elementIndex: ElementIndex): Loci {
    const { elements, model } = unit
    if (OrderedSet.indexOf(elements, elementIndex) !== -1) {
        const { index, offsets } = model.atomicHierarchy.residueAtomSegments
        const rI = index[elementIndex]
        const _indices: number[] = []
        for (let i = offsets[rI], il = offsets[rI + 1]; i < il; ++i) {
            const unitIndex = OrderedSet.indexOf(elements, i)
            if (unitIndex !== -1) _indices.push(unitIndex)
        }
        const indices = OrderedSet.ofSortedArray<StructureElement.UnitIndex>(SortedArray.ofSortedArray(_indices))
        return StructureElement.Loci(structure, [{ unit, indices }])
    }
    return EmptyLoci
}

//

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

export function createComplexMeshRenderObject(structure: Structure, mesh: Mesh, locationIt: LocationIterator, theme: Theme, props: PD.Values<StructureMeshParams>) {
    const transform = createIdentityTransform()
    const values = Mesh.createValues(mesh, transform, locationIt, theme, props)
    const state = createRenderableState(props)
    return createMeshRenderObject(values, state)
}

export function createUnitsMeshRenderObject(group: Unit.SymmetryGroup, mesh: Mesh, locationIt: LocationIterator, theme: Theme, props: PD.Values<StructureMeshParams>) {
    const transform = createUnitsTransform(group)
    const values = Mesh.createValues(mesh, transform, locationIt, theme, props)
    const state = createRenderableState(props)
    return createMeshRenderObject(values, state)
}

// points

export function createUnitsPointsRenderObject(group: Unit.SymmetryGroup, points: Points, locationIt: LocationIterator, theme: Theme, props: PD.Values<StructurePointsParams>) {
    const transform = createUnitsTransform(group)
    const values = Points.createValues(points, transform, locationIt, theme, props)
    const state = createRenderableState(props)
    return createPointsRenderObject(values, state)
}

// lines

export function createUnitsLinesRenderObject(group: Unit.SymmetryGroup, lines: Lines, locationIt: LocationIterator, theme: Theme, props: PD.Values<StructureLinesParams>) {
    const transform = createUnitsTransform(group)
    const values = Lines.createValues(lines, transform, locationIt, theme, props)
    const state = createRenderableState(props)
    return createLinesRenderObject(values, state)
}

// direct-volume

export function createUnitsDirectVolumeRenderObject(group: Unit.SymmetryGroup, directVolume: DirectVolume, locationIt: LocationIterator, theme: Theme, props: PD.Values<StructureDirectVolumeParams>) {
    const transform = createUnitsTransform(group)
    const values = DirectVolume.createValues(directVolume, transform, locationIt, theme, props)
    const state = createRenderableState(props)
    return createDirectVolumeRenderObject(values, state)
}