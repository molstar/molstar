/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { Unit, StructureElement, Structure, ElementIndex } from '../../../../mol-model/structure';
import { Loci, EmptyLoci } from '../../../../mol-model/loci';
import { Interval, OrderedSet, SortedArray } from '../../../../mol-data/int';
import { Mesh } from '../../../../mol-geo/geometry/mesh/mesh';
import { sphereVertexCount } from '../../../../mol-geo/primitive/sphere';
import { MeshBuilder } from '../../../../mol-geo/geometry/mesh/mesh-builder';
import { addSphere } from '../../../../mol-geo/geometry/mesh/builder/sphere';
import { PickingId } from '../../../../mol-geo/geometry/picking';
import { LocationIterator } from '../../../../mol-geo/util/location-iterator';
import { VisualContext } from '../../../../mol-repr/visual';
import { Theme } from '../../../../mol-theme/theme';
import { Spheres } from '../../../../mol-geo/geometry/spheres/spheres';
import { SpheresBuilder } from '../../../../mol-geo/geometry/spheres/spheres-builder';
import { isTrace, isH, StructureGroup } from './common';
import { Sphere3D } from '../../../../mol-math/geometry';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3add = Vec3.add;

type ElementProps = {
    ignoreHydrogens: boolean,
    traceOnly: boolean,
}

export type ElementSphereMeshProps = {
    detail: number,
    sizeFactor: number,
} & ElementProps

export function makeElementIgnoreTest(structure: Structure, unit: Unit, props: ElementProps): undefined | ((i: ElementIndex) => boolean) {
    const { ignoreHydrogens, traceOnly } = props;

    const { atomicNumber } = unit.model.atomicHierarchy.derived.atom;
    const isCoarse = Unit.isCoarse(unit);

    const { child } = structure;
    const childUnit = child?.unitMap.get(unit.id);
    if (child && !childUnit) throw new Error('expected childUnit to exist if child exists');

    if (!child && !ignoreHydrogens && !traceOnly) return;

    return (element: ElementIndex) => {
        return (
            (!!childUnit && !SortedArray.has(childUnit.elements, element)) ||
            (!isCoarse && ignoreHydrogens && isH(atomicNumber, element)) ||
            (traceOnly && !isTrace(unit, element))
        );
    };
}

export function createElementSphereMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: ElementSphereMeshProps, mesh?: Mesh): Mesh {
    const { child } = structure;
    const childUnit = child?.unitMap.get(unit.id);
    if (child && !childUnit) return Mesh.createEmpty(mesh);

    const { detail, sizeFactor } = props;

    const { elements } = unit;
    const elementCount = elements.length;
    const vertexCount = elementCount * sphereVertexCount(detail);
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 2, mesh);

    const v = Vec3();
    const pos = unit.conformation.invariantPosition;
    const ignore = makeElementIgnoreTest(structure, unit, props);
    const l = StructureElement.Location.create(structure, unit);
    const themeSize = theme.size.size;
    const center = Vec3();
    let maxSize = 0;
    let count = 0;

    for (let i = 0; i < elementCount; i++) {
        if (ignore && ignore(elements[i])) continue;

        l.element = elements[i];
        pos(elements[i], v);
        v3add(center, center, v);
        count += 1;

        builderState.currentGroup = i;
        const size = themeSize(l);
        if (size > maxSize) maxSize = size;

        addSphere(builderState, v, size * sizeFactor, detail);
    }

    // re-use boundingSphere if it has not changed much
    let boundingSphere: Sphere3D;
    Vec3.scale(center, center, 1 / count);
    if (mesh && Vec3.distance(center, mesh.boundingSphere.center) / mesh.boundingSphere.radius < 1.0) {
        boundingSphere = Sphere3D.clone(mesh.boundingSphere);
    } else {
        boundingSphere = Sphere3D.expand(Sphere3D(), (childUnit ?? unit).boundary.sphere, maxSize * sizeFactor + 0.05);
    }

    const m = MeshBuilder.getMesh(builderState);
    m.setBoundingSphere(boundingSphere);

    return m;
}

export type ElementSphereImpostorProps = {
    sizeFactor: number,
} & ElementProps

export function createElementSphereImpostor(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: ElementSphereImpostorProps, spheres?: Spheres): Spheres {
    const { child } = structure;
    const childUnit = child?.unitMap.get(unit.id);
    if (child && !childUnit) return Spheres.createEmpty(spheres);

    const { elements } = unit;
    const elementCount = elements.length;
    const builder = SpheresBuilder.create(elementCount, elementCount / 2, spheres);

    const v = Vec3();
    const pos = unit.conformation.invariantPosition;
    const ignore = makeElementIgnoreTest(structure, unit, props);

    const l = StructureElement.Location.create(structure, unit);
    const themeSize = theme.size.size;
    const center = Vec3();
    let maxSize = 0;
    let count = 0;

    for (let i = 0; i < elementCount; i++) {
        if (ignore?.(elements[i])) continue;

        pos(elements[i], v);
        builder.add(v[0], v[1], v[2], i);
        v3add(center, center, v);
        count += 1;

        l.element = elements[i];
        const size = themeSize(l);
        if (size > maxSize) maxSize = size;
    }

    // re-use boundingSphere if it has not changed much
    let boundingSphere: Sphere3D;
    Vec3.scale(center, center, 1 / count);
    if (spheres && Vec3.distance(center, spheres.boundingSphere.center) / spheres.boundingSphere.radius < 1.0) {
        boundingSphere = Sphere3D.clone(spheres.boundingSphere);
    } else {
        boundingSphere = Sphere3D.expand(Sphere3D(), (childUnit ?? unit).boundary.sphere, maxSize * props.sizeFactor + 0.05);
    }

    const s = builder.getSpheres();
    s.setBoundingSphere(boundingSphere);

    return s;
}

export function eachElement(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (!StructureElement.Loci.is(loci)) return false;
    const { structure, group } = structureGroup;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;
    const elementCount = group.elements.length;
    const { unitIndexMap } = group;
    for (const e of loci.elements) {
        const unitIdx = unitIndexMap.get(e.unit.id);
        if (unitIdx !== undefined) {
            const offset = unitIdx * elementCount; // to target unit instance
            if (Interval.is(e.indices)) {
                const start = offset + Interval.start(e.indices);
                const end = offset + Interval.end(e.indices);
                if (apply(Interval.ofBounds(start, end))) changed = true;
            } else {
                for (let i = 0, _i = e.indices.length; i < _i; i++) {
                    const start = e.indices[i];
                    let endI = i + 1;
                    while (endI < _i && e.indices[endI] === start) endI++;
                    i = endI - 1;
                    const end = e.indices[i];
                    changed = apply(Interval.ofRange(offset + start, offset + end)) || changed;
                }
            }
        }
    }
    return changed;
}

export function getElementLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId;
    if (id === objectId) {
        const { structure, group } = structureGroup;
        const unit = group.units[instanceId];
        const indices = OrderedSet.ofSingleton(groupId as StructureElement.UnitIndex);
        return StructureElement.Loci(structure.target, [{ unit, indices }]);
    }
    return EmptyLoci;
}

//

export function eachSerialElement(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (!StructureElement.Loci.is(loci)) return false;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;
    const { cumulativeUnitElementCount } = structure.serialMapping;
    for (const e of loci.elements) {
        const unitIdx = structure.unitIndexMap.get(e.unit.id);
        if (unitIdx !== undefined) {
            if (Interval.is(e.indices)) {
                const start = cumulativeUnitElementCount[unitIdx] + Interval.start(e.indices);
                const end = cumulativeUnitElementCount[unitIdx] + Interval.end(e.indices);
                if (apply(Interval.ofBounds(start, end))) changed = true;
            } else {
                for (let i = 0, _i = e.indices.length; i < _i; i++) {
                    const idx = cumulativeUnitElementCount[unitIdx] + e.indices[i];
                    if (apply(Interval.ofSingleton(idx))) changed = true;
                }
            }
        }
    }
    return changed;
}

export function getSerialElementLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        const { unitIndices, cumulativeUnitElementCount } = structure.serialMapping;
        const unitIdx = unitIndices[groupId];
        const unit = structure.units[unitIdx];
        const idx = groupId - cumulativeUnitElementCount[unitIdx];
        const indices = OrderedSet.ofSingleton(idx as StructureElement.UnitIndex);
        return StructureElement.Loci(structure, [{ unit, indices }]);
    }
    return EmptyLoci;
}

//

export namespace ElementIterator {
    export function fromGroup(structureGroup: StructureGroup): LocationIterator {
        const { group, structure } = structureGroup;
        const groupCount = group.elements.length;
        const instanceCount = group.units.length;
        const location = StructureElement.Location.create(structure);
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex];
            location.unit = unit;
            location.element = unit.elements[groupIndex];
            return location;
        };
        return LocationIterator(groupCount, instanceCount, 1, getLocation);
    }

    export function fromStructure(structure: Structure): LocationIterator {
        const { units, elementCount } = structure;
        const groupCount = elementCount;
        const instanceCount = 1;
        const { unitIndices, elementIndices } = structure.serialMapping;
        const location = StructureElement.Location.create(structure);
        const getLocation = (groupIndex: number) => {
            location.unit = units[unitIndices[groupIndex]];
            location.element = elementIndices[groupIndex];
            return location;
        };
        return LocationIterator(groupCount, instanceCount, 1, getLocation, true);
    }
}