/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Structure, StructureElement, Link } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { BitFlags } from '../../../mol-util';
import { createLinkCylinderMesh, LinkCylinderParams, LinkIterator } from './util/link';
import { ComplexMeshParams, ComplexVisual, ComplexMeshVisual } from '../complex-visual';
import { VisualUpdateState } from '../../util';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Interval, OrderedSet } from '../../../mol-data/int';
import { isHydrogen } from './util/common';

function createInterUnitLinkCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InterUnitLinkParams>, mesh?: Mesh) {
    const links = structure.links
    const { bondCount, bonds } = links
    const { sizeFactor, sizeAspectRatio, ignoreHydrogens } = props

    if (!bondCount) return Mesh.createEmpty(mesh)

    const location = StructureElement.Location.create()

    const builderProps = {
        linkCount: bondCount,
        referencePosition: (edgeIndex: number) => null, // TODO
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const b = bonds[edgeIndex]
            const uA = b.unitA, uB = b.unitB
            uA.conformation.position(uA.elements[b.indexA], posA)
            uB.conformation.position(uB.elements[b.indexB], posB)
        },
        order: (edgeIndex: number) => bonds[edgeIndex].order,
        flags: (edgeIndex: number) => BitFlags.create(bonds[edgeIndex].flag),
        radius: (edgeIndex: number) => {
            const b = bonds[edgeIndex]
            location.unit = b.unitA
            location.element = b.unitA.elements[b.indexA]
            return theme.size.size(location) * sizeFactor * sizeAspectRatio
        },
        ignore: ignoreHydrogens ? (edgeIndex: number) => {
            const b = bonds[edgeIndex]
            const uA = b.unitA, uB = b.unitB
            return isHydrogen(uA, uA.elements[b.indexA]) || isHydrogen(uB, uB.elements[b.indexB])
        } : () => false
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const InterUnitLinkParams = {
    ...ComplexMeshParams,
    ...LinkCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2/3, { min: 0, max: 3, step: 0.01 }),
    ignoreHydrogens: PD.Boolean(false),
}
export type InterUnitLinkParams = typeof InterUnitLinkParams

export function InterUnitLinkVisual(materialId: number): ComplexVisual<InterUnitLinkParams> {
    return ComplexMeshVisual<InterUnitLinkParams>({
        defaultProps: PD.getDefaultValues(InterUnitLinkParams),
        createGeometry: createInterUnitLinkCylinderMesh,
        createLocationIterator: LinkIterator.fromStructure,
        getLoci: getLinkLoci,
        eachLocation: eachLink,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<InterUnitLinkParams>, currentProps: PD.Values<InterUnitLinkParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.sizeAspectRatio !== currentProps.sizeAspectRatio ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.linkScale !== currentProps.linkScale ||
                newProps.linkSpacing !== currentProps.linkSpacing ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens
            )
        }
    }, materialId)
}

function getLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        const bond = structure.links.bonds[groupId]
        return Link.Loci(structure, [
            Link.Location(
                bond.unitA, bond.indexA as StructureElement.UnitIndex,
                bond.unitB, bond.indexB as StructureElement.UnitIndex
            ),
            Link.Location(
                bond.unitB, bond.indexB as StructureElement.UnitIndex,
                bond.unitA, bond.indexA as StructureElement.UnitIndex
            )
        ])
    }
    return EmptyLoci
}

function eachLink(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false
    if (Link.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false
        for (const b of loci.links) {
            const idx = structure.links.getBondIndex(b.aIndex, b.aUnit, b.bIndex, b.bUnit)
            if (idx !== -1) {
                if (apply(Interval.ofSingleton(idx))) changed = true
            }
        }
    } else if (StructureElement.Loci.is(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false
        // TODO mark link only when both of the link elements are in a StructureElement.Loci
        for (const e of loci.elements) {
            OrderedSet.forEach(e.indices, v => {
                const indices = structure.links.getBondIndices(v, e.unit)
                for (let i = 0, il = indices.length; i < il; ++i) {
                    if (apply(Interval.ofSingleton(indices[i]))) changed = true
                }
            })
        }
    }
    return changed
}