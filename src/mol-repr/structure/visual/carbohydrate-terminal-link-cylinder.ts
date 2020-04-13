/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Structure, StructureElement, Unit } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { createLinkCylinderMesh, LinkCylinderParams, LinkCylinderStyle } from './util/link';
import { UnitsMeshParams } from '../units-visual';
import { ComplexVisual, ComplexMeshVisual } from '../complex-visual';
import { VisualUpdateState } from '../../util';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { OrderedSet, Interval } from '../../../mol-data/int';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { getElementIdx, MetalsSet } from '../../../mol-model/structure/structure/unit/bonds/common';
import { getAltResidueLociFromId, getAltResidueLoci } from './util/common';

function createCarbohydrateTerminalLinkCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<CarbohydrateTerminalLinkParams>, mesh?: Mesh) {
    const { terminalLinks, elements } = structure.carbohydrates;
    const { terminalLinkSizeFactor } = props;

    const location = StructureElement.Location.create(structure);

    const builderProps = {
        linkCount: terminalLinks.length,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const l = terminalLinks[edgeIndex];
            if (l.fromCarbohydrate) {
                Vec3.copy(posA, elements[l.carbohydrateIndex].geometry.center);
                l.elementUnit.conformation.position(l.elementUnit.elements[l.elementIndex], posB);
            } else {
                l.elementUnit.conformation.position(l.elementUnit.elements[l.elementIndex], posA);
                Vec3.copy(posB, elements[l.carbohydrateIndex].geometry.center);
            }
        },
        radius: (edgeIndex: number) => {
            const l = terminalLinks[edgeIndex];
            if (l.fromCarbohydrate) {
                const carb = elements[l.carbohydrateIndex];
                const ring = carb.unit.rings.all[carb.ringIndex];
                location.unit = carb.unit;
                location.element = carb.unit.elements[ring[0]];
            } else {
                location.unit = l.elementUnit;
                location.element = l.elementUnit.elements[l.elementIndex];
            }
            return theme.size.size(location) * terminalLinkSizeFactor;
        },
        style: (edgeIndex: number) => {
            const l = terminalLinks[edgeIndex];
            const eI = l.elementUnit.elements[l.elementIndex];
            const beI = getElementIdx(l.elementUnit.model.atomicHierarchy.atoms.type_symbol.value(eI));
            return MetalsSet.has(beI) ? LinkCylinderStyle.Dashed : LinkCylinderStyle.Solid;
        }
    };

    return createLinkCylinderMesh(ctx, builderProps, props, mesh);
}

export const CarbohydrateTerminalLinkParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
    terminalLinkSizeFactor: PD.Numeric(0.2, { min: 0, max: 3, step: 0.01 }),
};
export type CarbohydrateTerminalLinkParams = typeof CarbohydrateTerminalLinkParams

export function CarbohydrateTerminalLinkVisual(materialId: number): ComplexVisual<CarbohydrateTerminalLinkParams> {
    return ComplexMeshVisual<CarbohydrateTerminalLinkParams>({
        defaultProps: PD.getDefaultValues(CarbohydrateTerminalLinkParams),
        createGeometry: createCarbohydrateTerminalLinkCylinderMesh,
        createLocationIterator: CarbohydrateTerminalLinkIterator,
        getLoci: getTerminalLinkLoci,
        eachLocation: eachTerminalLink,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<CarbohydrateTerminalLinkParams>, currentProps: PD.Values<CarbohydrateTerminalLinkParams>) => {
            state.createGeometry = (
                newProps.terminalLinkSizeFactor !== currentProps.terminalLinkSizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.linkCap !== currentProps.linkCap
            );
        }
    }, materialId);
}

function CarbohydrateTerminalLinkIterator(structure: Structure): LocationIterator {
    const { elements, terminalLinks } = structure.carbohydrates;
    const groupCount = terminalLinks.length;
    const instanceCount = 1;
    const location = StructureElement.Location.create(structure);
    const getLocation = (groupIndex: number) => {
        const terminalLink = terminalLinks[groupIndex];
        if (terminalLink.fromCarbohydrate) {
            const carb = elements[terminalLink.carbohydrateIndex];
            const ring = carb.unit.rings.all[carb.ringIndex];
            location.unit = carb.unit;
            location.element = carb.unit.elements[ring[0]];
        } else {
            location.unit = terminalLink.elementUnit;
            location.element = terminalLink.elementUnit.elements[terminalLink.elementIndex];
        }
        return location;
    };
    return LocationIterator(groupCount, instanceCount, getLocation, true);
}

function getTerminalLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        const { terminalLinks, elements } = structure.carbohydrates;
        const l = terminalLinks[groupId];
        const carb = elements[l.carbohydrateIndex];

        return StructureElement.Loci.union(
            getAltResidueLociFromId(structure, carb.unit, carb.residueIndex, carb.altId),
            getAltResidueLoci(structure, l.elementUnit, l.elementUnit.elements[l.elementIndex])
        );
    }
    return EmptyLoci;
}

function eachTerminalLink(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (!StructureElement.Loci.is(loci)) return false;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;

    const { getTerminalLinkIndices } = structure.carbohydrates;
    for (const { unit, indices } of loci.elements) {
        if (!Unit.isAtomic(unit)) continue;

        OrderedSet.forEach(indices, v => {
            // TODO avoid duplicate calls to apply
            const linkIndices = getTerminalLinkIndices(unit, unit.elements[v]);
            for (let i = 0, il = linkIndices.length; i < il; ++i) {
                if (apply(Interval.ofSingleton(linkIndices[i]))) changed = true;
            }
        });
    }
    return changed;
}