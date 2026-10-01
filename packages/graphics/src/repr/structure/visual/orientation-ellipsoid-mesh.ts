/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { UnitsMeshParams, type UnitsVisual, UnitsMeshVisual } from '../units-visual.js';
import type { VisualUpdateState } from '../../util.js';
import type { VisualContext } from '../../visual.js';
import { Unit, Structure, StructureElement } from '@molstar/model/model/structure';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { addEllipsoid } from '@molstar/graphics/geo/geometry/mesh/builder/ellipsoid';
import { Axes3D, Sphere3D } from '@molstar/core/math/geometry';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { OrderedSet, Interval } from '@molstar/core/data/int';
import { EmptyLoci, type Loci } from '@molstar/model/model/loci';
import type { UnitIndex } from '@molstar/model/model/structure/structure/element/element';
import { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import { MoleculeType } from '@molstar/model/model/structure/model/types';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import type { StructureGroup } from './util/common.js';

export const OrientationEllipsoidMeshParams = {
    ...UnitsMeshParams,
    sizeFactor: PD.Numeric(1, { min: 0, max: 2, step: 0.1 }),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
};
export type OrientationEllipsoidMeshParams = typeof OrientationEllipsoidMeshParams

export function OrientationEllipsoidMeshVisual(materialId: number): UnitsVisual<OrientationEllipsoidMeshParams> {
    return UnitsMeshVisual<OrientationEllipsoidMeshParams>({
        defaultProps: PD.getDefaultValues(OrientationEllipsoidMeshParams),
        createGeometry: createOrientationEllipsoidMesh,
        createLocationIterator: UnitIterator,
        getLoci: getUnitLoci,
        eachLocation: eachUnit,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<OrientationEllipsoidMeshParams>, currentProps: PD.Values<OrientationEllipsoidMeshParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail
            );
        }
    }, materialId);
}

//

export interface OrientationEllipsoidMeshProps {
    detail: number,
    sizeFactor: number,
}

function isUnitApplicable(unit: Unit) {
    if (Unit.Traits.is(unit.traits, Unit.Trait.MultiChain)) return false;
    if (Unit.Traits.is(unit.traits, Unit.Trait.Partitioned)) return false;
    if (Unit.isCoarse(unit)) return true;
    if (unit.elements.length === 0) return false;
    unit.model.atomicHierarchy.derived.residue.moleculeType;
    const rI = unit.residueIndex[unit.elements[0]];
    const mt = unit.model.atomicHierarchy.derived.residue.moleculeType[rI];
    if (mt === MoleculeType.Ion) return false;
    if (mt === MoleculeType.Water) return false;
    return true;
}

export function createOrientationEllipsoidMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: OrientationEllipsoidMeshProps, mesh?: Mesh): Mesh {
    if (!isUnitApplicable(unit)) return Mesh.createEmpty(mesh);

    const { detail, sizeFactor } = props;

    const vertexCount = 256;
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 2, mesh);
    const axes = unit.principalAxes.boxAxes;
    const { origin, dirA, dirB } = axes;

    const size = Axes3D.size(Vec3(), axes);
    Vec3.scale(size, size, sizeFactor / 2);
    const radiusScale = Vec3.create(size[2], size[1], size[0]);

    builderState.currentGroup = 0;
    addEllipsoid(builderState, origin, dirA, dirB, radiusScale, detail + 1);

    const m = MeshBuilder.getMesh(builderState);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    m.setBoundingSphere(sphere);

    return m;
}

//

function UnitIterator(structureGroup: StructureGroup): LocationIterator {
    const { group, structure } = structureGroup;
    const groupCount = 1;
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

function getUnitLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId } = pickingId;
    if (id === objectId) {
        const { structure, group } = structureGroup;
        const unit = group.units[instanceId];
        const indices = OrderedSet.ofBounds(0, unit.elements.length) as OrderedSet<UnitIndex>;
        return StructureElement.Loci(structure, [{ unit, indices }]);
    }
    return EmptyLoci;
}

function eachUnit(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (!StructureElement.Loci.is(loci)) return false;
    const { structure, group } = structureGroup;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;
    const elementCount = group.elements.length;
    for (const e of loci.elements) {
        const unitIdx = group.unitIndexMap.get(e.unit.id);
        if (unitIdx !== undefined) {
            if (OrderedSet.size(e.indices) === elementCount) {
                if (apply(Interval.ofSingleton(unitIdx))) changed = true;
            }
        }
    }
    return changed;
}