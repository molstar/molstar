/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Unit, Structure, StructureElement, Link } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { BitFlags } from '../../../mol-util';
import { createLinkCylinderMesh, LinkCylinderParams, LinkIterator } from './util/link';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual, StructureGroup } from '../units-visual';
import { VisualUpdateState } from '../../util';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Interval, OrderedSet } from '../../../mol-data/int';
import { isHydrogen } from './util/common';

function createIntraUnitLinkCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<IntraUnitLinkParams>, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    const location = StructureElement.Location.create(unit)

    const elements = unit.elements;
    const links = unit.links
    const { edgeCount, a, b, edgeProps, offset } = links
    const { order: _order, flags: _flags } = edgeProps
    const { sizeFactor, sizeAspectRatio, ignoreHydrogens } = props

    if (!edgeCount) return Mesh.createEmpty(mesh)

    const vRef = Vec3.zero()
    const pos = unit.conformation.invariantPosition

    const builderProps = {
        linkCount: edgeCount * 2,
        referencePosition: (edgeIndex: number) => {
            let aI = a[edgeIndex], bI = b[edgeIndex];

            if (aI > bI) [aI, bI] = [bI, aI]
            if (offset[aI + 1] - offset[aI] === 1) [aI, bI] = [bI, aI]
            // TODO prefer reference atoms in rings

            for (let i = offset[aI], il = offset[aI + 1]; i < il; ++i) {
                const _bI = b[i]
                if (_bI !== bI && _bI !== aI) return pos(elements[_bI], vRef)
            }
            for (let i = offset[bI], il = offset[bI + 1]; i < il; ++i) {
                const _aI = a[i]
                if (_aI !== aI && _aI !== bI) return pos(elements[_aI], vRef)
            }
            return null
        },
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            pos(elements[a[edgeIndex]], posA)
            pos(elements[b[edgeIndex]], posB)
        },
        order: (edgeIndex: number) => _order[edgeIndex],
        flags: (edgeIndex: number) => BitFlags.create(_flags[edgeIndex]),
        radius: (edgeIndex: number) => {
            location.element = elements[a[edgeIndex]]
            return theme.size.size(location) * sizeFactor * sizeAspectRatio
        },
        ignore: ignoreHydrogens ? (edgeIndex: number) => {
            return isHydrogen(unit, elements[a[edgeIndex]]) || isHydrogen(unit, elements[b[edgeIndex]])
        } : () => false
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const IntraUnitLinkParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2/3, { min: 0, max: 3, step: 0.01 }),
    ignoreHydrogens: PD.Boolean(false),
}
export type IntraUnitLinkParams = typeof IntraUnitLinkParams

export function IntraUnitLinkVisual(materialId: number): UnitsVisual<IntraUnitLinkParams> {
    return UnitsMeshVisual<IntraUnitLinkParams>({
        defaultProps: PD.getDefaultValues(IntraUnitLinkParams),
        createGeometry: createIntraUnitLinkCylinderMesh,
        createLocationIterator: LinkIterator.fromGroup,
        getLoci: getLinkLoci,
        eachLocation: eachLink,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<IntraUnitLinkParams>, currentProps: PD.Values<IntraUnitLinkParams>) => {
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

function getLinkLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId
    if (id === objectId) {
        const { structure, group } = structureGroup
        const unit = group.units[instanceId]
        if (Unit.isAtomic(unit)) {
            return Link.Loci(structure, [
                Link.Location(
                    unit, unit.links.a[groupId] as StructureElement.UnitIndex,
                    unit, unit.links.b[groupId] as StructureElement.UnitIndex
                ),
                Link.Location(
                    unit, unit.links.b[groupId] as StructureElement.UnitIndex,
                    unit, unit.links.a[groupId] as StructureElement.UnitIndex
                )
            ])
        }
    }
    return EmptyLoci
}

function eachLink(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false
    if (Link.isLoci(loci)) {
        const { structure, group } = structureGroup
        if (!Structure.areEquivalent(loci.structure, structure)) return false
        const unit = group.units[0]
        if (!Unit.isAtomic(unit)) return false
        const groupCount = unit.links.edgeCount * 2
        for (const b of loci.links) {
            const unitIdx = group.unitIndexMap.get(b.aUnit.id)
            if (unitIdx !== undefined) {
                const idx = unit.links.getDirectedEdgeIndex(b.aIndex, b.bIndex)
                if (idx !== -1) {
                    if (apply(Interval.ofSingleton(unitIdx * groupCount + idx))) changed = true
                }
            }
        }
    } else if (StructureElement.Loci.is(loci)) {
        const { structure, group } = structureGroup
        if (!Structure.areEquivalent(loci.structure, structure)) return false
        const unit = group.units[0]
        if (!Unit.isAtomic(unit)) return false
        const groupCount = unit.links.edgeCount * 2
        for (const e of loci.elements) {
            const unitIdx = group.unitIndexMap.get(e.unit.id)
            if (unitIdx !== undefined) {
                const { offset, b } = unit.links
                OrderedSet.forEach(e.indices, v => {
                    for (let t = offset[v], _t = offset[v + 1]; t < _t; t++) {
                        if (OrderedSet.has(e.indices, b[t])) {
                            if (apply(Interval.ofSingleton(unitIdx * groupCount + t))) changed = true
                        }
                    }
                })
            }
        }
    }
    return changed
}