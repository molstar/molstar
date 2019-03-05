/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit, Link, StructureElement, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { createLinkCylinderMesh, LinkIterator, LinkCylinderParams } from './util/link';
import { Vec3 } from 'mol-math/linear-algebra';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { UnitsMeshVisual, UnitsMeshParams, StructureGroup } from '../units-visual';
import { Interval, OrderedSet } from 'mol-data/int';
import { BitFlags } from 'mol-util';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { PickingId } from 'mol-geo/geometry/picking';
import { VisualContext } from 'mol-repr/visual';
import { Theme } from 'mol-theme/theme';

function createIntraUnitLinkCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<IntraUnitLinkParams>, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    const location = StructureElement.create(unit)

    const elements = unit.elements;
    const links = unit.links
    const { edgeCount, a, b, edgeProps, offset } = links
    const { order: _order, flags: _flags } = edgeProps
    const { sizeFactor, sizeAspectRatio } = props

    if (!edgeCount) return Mesh.createEmpty(mesh)

    const vRef = Vec3.zero()
    const pos = unit.conformation.invariantPosition

    const builderProps = {
        linkCount: edgeCount * 2,
        referencePosition: (edgeIndex: number) => {
            let aI = a[edgeIndex], bI = b[edgeIndex];
            if (aI > bI) [aI, bI] = [bI, aI]
            for (let i = offset[aI], il = offset[aI + 1]; i < il; ++i) {
                if (b[i] !== bI) return pos(elements[b[i]], vRef)
            }
            for (let i = offset[bI], il = offset[bI + 1]; i < il; ++i) {
                if (a[i] !== aI) return pos(elements[a[i]], vRef)
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
        }
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const IntraUnitLinkParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2/3, { min: 0, max: 3, step: 0.01 }),
}
export type IntraUnitLinkParams = typeof IntraUnitLinkParams

export function IntraUnitLinkVisual(): UnitsVisual<IntraUnitLinkParams> {
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
                newProps.linkSpacing !== currentProps.linkSpacing
            )
        }
    })
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
        if (loci.structure !== structure) return false
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
    } else if (StructureElement.isLoci(loci)) {
        const { structure, group } = structureGroup
        if (loci.structure !== structure) return false
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