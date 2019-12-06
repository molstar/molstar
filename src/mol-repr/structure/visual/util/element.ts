/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { Unit, StructureElement, Structure } from '../../../../mol-model/structure';
import { Loci, EmptyLoci } from '../../../../mol-model/loci';
import { Interval, OrderedSet } from '../../../../mol-data/int';
import { Mesh } from '../../../../mol-geo/geometry/mesh/mesh';
import { sphereVertexCount } from '../../../../mol-geo/primitive/sphere';
import { MeshBuilder } from '../../../../mol-geo/geometry/mesh/mesh-builder';
import { addSphere } from '../../../../mol-geo/geometry/mesh/builder/sphere';
import { PickingId } from '../../../../mol-geo/geometry/picking';
import { LocationIterator } from '../../../../mol-geo/util/location-iterator';
import { VisualContext } from '../../../../mol-repr/visual';
import { Theme } from '../../../../mol-theme/theme';
import { StructureGroup } from '../../../../mol-repr/structure/units-visual';
import { Spheres } from '../../../../mol-geo/geometry/spheres/spheres';
import { SpheresBuilder } from '../../../../mol-geo/geometry/spheres/spheres-builder';
import { isHydrogen } from './common';

export interface ElementSphereMeshProps {
    detail: number,
    sizeFactor: number,
    ignoreHydrogens: boolean,
}

export function createElementSphereMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: ElementSphereMeshProps, mesh?: Mesh): Mesh {
    const { detail, sizeFactor } = props

    const { elements } = unit;
    const elementCount = elements.length;
    const vertexCount = elementCount * sphereVertexCount(detail)
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 2, mesh)

    const v = Vec3.zero()
    const pos = unit.conformation.invariantPosition
    const l = StructureElement.Location.create()
    l.unit = unit

    for (let i = 0; i < elementCount; i++) {
        if (props.ignoreHydrogens && isHydrogen(unit, elements[i])) continue

        l.element = elements[i]
        pos(elements[i], v)

        builderState.currentGroup = i
        addSphere(builderState, v, theme.size.size(l) * sizeFactor, detail)
    }

    return MeshBuilder.getMesh(builderState)
}

export interface ElementSphereImpostorProps {
    ignoreHydrogens: boolean,
}

export function createElementSphereImpostor(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: ElementSphereImpostorProps, spheres?: Spheres): Spheres {

    const { elements } = unit;
    const elementCount = elements.length;
    const builder = SpheresBuilder.create(elementCount, elementCount / 2, spheres)

    const v = Vec3.zero()
    const pos = unit.conformation.invariantPosition

    for (let i = 0; i < elementCount; i++) {
        if (props.ignoreHydrogens && isHydrogen(unit, elements[i])) continue

        pos(elements[i], v)
        builder.add(v[0], v[1], v[2], i)
    }

    return builder.getSpheres()
}

export function eachElement(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false
    if (!StructureElement.Loci.is(loci)) return false
    const { structure, group } = structureGroup
    if (!Structure.areEquivalent(loci.structure, structure)) return false
    const elementCount = group.elements.length
    for (const e of loci.elements) {
        const unitIdx = group.unitIndexMap.get(e.unit.id)
        if (unitIdx !== undefined) {
            if (Interval.is(e.indices)) {
                const start = unitIdx * elementCount + Interval.start(e.indices)
                const end = unitIdx * elementCount + Interval.end(e.indices)
                if (apply(Interval.ofBounds(start, end))) changed = true
            } else {
                for (let i = 0, _i = e.indices.length; i < _i; i++) {
                    const start = e.indices[i];
                    let endI = i + 1;
                    while (endI < _i && e.indices[i] === start) endI++;
                    i = endI - 1;
                    const end = e.indices[i];
                    changed = apply(Interval.ofRange(start, end)) || changed;
                }
            }
        }
    }
    return changed
}

export function getElementLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId
    if (id === objectId) {
        const { structure, group } = structureGroup
        const unit = group.units[instanceId]
        const indices = OrderedSet.ofSingleton(groupId as StructureElement.UnitIndex)
        return StructureElement.Loci(structure, [{ unit, indices }])
    }
    return EmptyLoci
}

//

export function eachSerialElement(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false
    if (!StructureElement.Loci.is(loci)) return false
    if (!Structure.areEquivalent(loci.structure, structure)) return false
    const { unitElementCount } = structure.serialMapping
    for (const e of loci.elements) {
        const unitIdx = structure.unitIndexMap.get(e.unit.id)
        if (unitIdx !== undefined) {
            if (Interval.is(e.indices)) {
                const start = unitElementCount[unitIdx] + Interval.start(e.indices)
                const end = unitElementCount[unitIdx] + Interval.end(e.indices)
                if (apply(Interval.ofBounds(start, end))) changed = true
            } else {
                for (let i = 0, _i = e.indices.length; i < _i; i++) {
                    const idx = unitElementCount[unitIdx] + e.indices[i]
                    if (apply(Interval.ofSingleton(idx))) changed = true
                }
            }
        }
    }
    return changed
}

export function getSerialElementLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        const { unitIndices, unitElementCount } = structure.serialMapping
        const unitIdx = unitIndices[groupId]
        const unit = structure.units[unitIdx]
        const idx = groupId - unitElementCount[unitIdx]
        const indices = OrderedSet.ofSingleton(idx as StructureElement.UnitIndex)
        return StructureElement.Loci(structure, [{ unit, indices }])
    }
    return EmptyLoci
}

//

export namespace ElementIterator {
    export function fromGroup(group: Unit.SymmetryGroup): LocationIterator {
        const groupCount = group.elements.length
        const instanceCount = group.units.length
        const location = StructureElement.Location.create()
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex]
            location.unit = unit
            location.element = unit.elements[groupIndex]
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation)
    }

    export function fromStructure(structure: Structure): LocationIterator {
        const { units, elementCount } = structure
        const groupCount = elementCount
        const instanceCount = 1
        const { unitIndices, elementIndices } = structure.serialMapping
        const location = StructureElement.Location.create()
        const getLocation = (groupIndex: number) => {
            location.unit = units[unitIndices[groupIndex]]
            location.element = elementIndices[groupIndex]
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation, true)
    }
}