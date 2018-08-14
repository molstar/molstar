/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit, Link, StructureElement } from 'mol-model/structure';
import { UnitsVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { DefaultLinkCylinderProps, LinkCylinderProps, createLinkCylinderMesh } from './util/link';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { Vec3 } from 'mol-math/linear-algebra';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { LinkIterator } from './util/location-iterator';
import { UnitsMeshVisual, DefaultUnitsMeshProps } from '../units-visual';
import { Interval } from 'mol-data/int';
import { SizeThemeProps } from 'mol-view/theme/size';

async function createIntraUnitLinkCylinderMesh(ctx: RuntimeContext, unit: Unit, props: LinkCylinderProps, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    const elements = unit.elements;
    const links = unit.links
    const { edgeCount, a, b, edgeProps, offset } = links
    const { order: _order, flags: _flags } = edgeProps

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
        flags: (edgeIndex: number) => _flags[edgeIndex]
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const DefaultIntraUnitLinkProps = {
    ...DefaultUnitsMeshProps,
    ...DefaultLinkCylinderProps,
    sizeTheme: { name: 'physical', factor: 0.3 } as SizeThemeProps,
}
export type IntraUnitLinkProps = typeof DefaultIntraUnitLinkProps

export function IntraUnitLinkVisual(): UnitsVisual<IntraUnitLinkProps> {
    return UnitsMeshVisual<IntraUnitLinkProps>({
        defaultProps: DefaultIntraUnitLinkProps,
        createMesh: createIntraUnitLinkCylinderMesh,
        createLocationIterator: LinkIterator.fromGroup,
        getLoci: getLinkLoci,
        mark: markLink,
        setUpdateState: () => {}
    })
}

function getLinkLoci(pickingId: PickingId, group: Unit.SymmetryGroup, id: number) {
    const { objectId, instanceId, elementId } = pickingId
    const unit = group.units[instanceId]
    if (id === objectId && Unit.isAtomic(unit)) {
        return Link.Loci([
            Link.Location(
                unit, unit.links.a[elementId] as StructureElement.UnitIndex,
                unit, unit.links.b[elementId] as StructureElement.UnitIndex
            )
        ])
    }
    return EmptyLoci
}

function markLink(loci: Loci, group: Unit.SymmetryGroup, apply: (interval: Interval) => boolean) {
    const unit = group.units[0]

    let changed = false
    if (Unit.isAtomic(unit) && Link.isLoci(loci)) {
        for (const b of loci.links) {
            const unitIdx = Unit.findUnitById(b.aUnit.id, group.units)
            if (unitIdx !== -1) {
                const idx = unit.links.getDirectedEdgeIndex(b.aIndex, b.bIndex)
                if (idx !== -1) {
                    if (apply(Interval.ofSingleton(idx))) changed = true
                }
            }
        }
    }
    return changed
}