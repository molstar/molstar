/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit, Link, StructureElement, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../index';
import { VisualUpdateState } from '../../util';
import { LinkCylinderProps, createLinkCylinderMesh, LinkIterator, LinkCylinderParams } from './util/link';
import { Vec3 } from 'mol-math/linear-algebra';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { Interval } from 'mol-data/int';
import { SizeThemeName, SizeThemeOptions } from 'mol-theme/size';
import { BitFlags } from 'mol-util';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { PickingId } from 'mol-geo/geometry/picking';
import { VisualContext } from 'mol-repr';
import { Theme } from 'mol-geo/geometry/geometry';

async function createIntraUnitLinkCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: LinkCylinderProps, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    const location = StructureElement.create(unit)

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
        flags: (edgeIndex: number) => BitFlags.create(_flags[edgeIndex]),
        radius: (edgeIndex: number) => {
            location.element = elements[a[edgeIndex]]
            return theme.size.size(location)
        }
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const IntraUnitLinkParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
    sizeTheme: PD.Select<SizeThemeName>('Size Theme', '', 'physical', SizeThemeOptions),
    sizeValue: PD.Numeric('Size Value', '', 0.2, 0, 10, 0.1),
    sizeFactor: PD.Numeric('Size Factor', '', 1, 0, 10, 0.1),
}
export const DefaultIntraUnitLinkProps = PD.getDefaultValues(IntraUnitLinkParams)
export type IntraUnitLinkProps = typeof DefaultIntraUnitLinkProps

export function IntraUnitLinkVisual(): UnitsVisual<IntraUnitLinkProps> {
    return UnitsMeshVisual<IntraUnitLinkProps>({
        defaultProps: DefaultIntraUnitLinkProps,
        createGeometry: createIntraUnitLinkCylinderMesh,
        createLocationIterator: LinkIterator.fromGroup,
        getLoci: getLinkLoci,
        mark: markLink,
        setUpdateState: (state: VisualUpdateState, newProps: LinkCylinderProps, currentProps: LinkCylinderProps) => {
            state.createGeometry = newProps.radialSegments !== currentProps.radialSegments
        }
    })
}

function getLinkLoci(pickingId: PickingId, group: Unit.SymmetryGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId
    const unit = group.units[instanceId]
    if (id === objectId && Unit.isAtomic(unit)) {
        return Link.Loci([
            Link.Location(
                unit, unit.links.a[groupId] as StructureElement.UnitIndex,
                unit, unit.links.b[groupId] as StructureElement.UnitIndex
            )
        ])
    }
    return EmptyLoci
}

function markLink(loci: Loci, group: Unit.SymmetryGroup, apply: (interval: Interval) => boolean) {
    const unit = group.units[0]

    let changed = false
    if (Unit.isAtomic(unit) && Link.isLoci(loci)) {
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
    }
    return changed
}