/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, StructureElement, Unit } from '../../../mol-model/structure';
import { Loci, EmptyLoci } from '../../../mol-model/loci';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { createLinkCylinderMesh, LinkCylinderParams } from './util/link';
import { OrderedSet, Interval } from '../../../mol-data/int';
import { ComplexMeshVisual, ComplexVisual } from '../complex-visual';
import { UnitsMeshParams } from '../units-visual';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { VisualUpdateState } from '../../util';
import { VisualContext } from '../../../mol-repr/visual';
import { Theme } from '../../../mol-theme/theme';
import { getAltResidueLociFromId } from './util/common';

function createCarbohydrateLinkCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<CarbohydrateLinkParams>, mesh?: Mesh) {
    const { links, elements } = structure.carbohydrates;
    const { linkSizeFactor } = props;

    const location = StructureElement.Location.create(structure);

    const builderProps = {
        linkCount: links.length,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const l = links[edgeIndex];
            Vec3.copy(posA, elements[l.carbohydrateIndexA].geometry.center);
            Vec3.copy(posB, elements[l.carbohydrateIndexB].geometry.center);
        },
        radius: (edgeIndex: number) => {
            const l = links[edgeIndex];
            const carbA = elements[l.carbohydrateIndexA];
            const ringA = carbA.unit.rings.all[carbA.ringIndex];
            location.unit = carbA.unit;
            location.element = carbA.unit.elements[ringA[0]];
            return theme.size.size(location) * linkSizeFactor;
        },
    };

    return createLinkCylinderMesh(ctx, builderProps, props, mesh);
}

export const CarbohydrateLinkParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
    linkSizeFactor: PD.Numeric(0.3, { min: 0, max: 3, step: 0.01 }),
};
export type CarbohydrateLinkParams = typeof CarbohydrateLinkParams

export function CarbohydrateLinkVisual(materialId: number): ComplexVisual<CarbohydrateLinkParams> {
    return ComplexMeshVisual<CarbohydrateLinkParams>({
        defaultProps: PD.getDefaultValues(CarbohydrateLinkParams),
        createGeometry: createCarbohydrateLinkCylinderMesh,
        createLocationIterator: CarbohydrateLinkIterator,
        getLoci: getLinkLoci,
        eachLocation: eachCarbohydrateLink,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<CarbohydrateLinkParams>, currentProps: PD.Values<CarbohydrateLinkParams>) => {
            state.createGeometry = (
                newProps.linkSizeFactor !== currentProps.linkSizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.linkCap !== currentProps.linkCap
            );
        }
    }, materialId);
}

function CarbohydrateLinkIterator(structure: Structure): LocationIterator {
    const { elements, links } = structure.carbohydrates;
    const groupCount = links.length;
    const instanceCount = 1;
    const location = StructureElement.Location.create(structure);
    const getLocation = (groupIndex: number) => {
        const link = links[groupIndex];
        const carbA = elements[link.carbohydrateIndexA];
        const ringA = carbA.unit.rings.all[carbA.ringIndex];
        location.unit = carbA.unit;
        location.element = carbA.unit.elements[ringA[0]];
        return location;
    };
    return LocationIterator(groupCount, instanceCount, getLocation, true);
}

function getLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        const { links, elements } = structure.carbohydrates;
        const l = links[groupId];
        const carbA = elements[l.carbohydrateIndexA];
        const carbB = elements[l.carbohydrateIndexB];
        return StructureElement.Loci.union(
            getAltResidueLociFromId(structure, carbA.unit, carbA.residueIndex, carbA.altId),
            getAltResidueLociFromId(structure, carbB.unit, carbB.residueIndex, carbB.altId)
        );
    }
    return EmptyLoci;
}

function eachCarbohydrateLink(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (!StructureElement.Loci.is(loci)) return false;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;

    const { getLinkIndices } = structure.carbohydrates;
    for (const { unit, indices } of loci.elements) {
        if (!Unit.isAtomic(unit)) continue;

        OrderedSet.forEach(indices, v => {
            // TODO avoid duplicate calls to apply
            const linkIndices = getLinkIndices(unit, unit.elements[v]);
            for (let i = 0, il = linkIndices.length; i < il; ++i) {
                if (apply(Interval.ofSingleton(linkIndices[i]))) changed = true;
            }
        });
    }
    return changed;
}