/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Structure, StructureElement, Bond, Unit } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { BitFlags } from '../../../mol-util';
import { createBondCylinderMesh, BondCylinderParams, BondIterator } from './util/bond';
import { ComplexMeshParams, ComplexVisual, ComplexMeshVisual } from '../complex-visual';
import { VisualUpdateState } from '../../util';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Interval, OrderedSet } from '../../../mol-data/int';
import { isHydrogen } from './util/common';

const tmpRefPosBondIt = new Bond.ElementBondIterator()
function setRefPosition(pos: Vec3, structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    tmpRefPosBondIt.setElement(structure, unit, index)
    while (tmpRefPosBondIt.hasNext) {
        const bA = tmpRefPosBondIt.move()
        bA.otherUnit.conformation.position(bA.otherUnit.elements[bA.otherIndex], pos)
        return pos
    }
    return null
}

const tmpRef = Vec3()
const tmpLoc = StructureElement.Location.create()

function createInterUnitBondCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InterUnitBondParams>, mesh?: Mesh) {
    const bonds = structure.interUnitBonds
    const { edgeCount, edges } = bonds
    const { sizeFactor, sizeAspectRatio, ignoreHydrogens } = props

    if (!edgeCount) return Mesh.createEmpty(mesh)

    const builderProps = {
        bondCount: edgeCount,
        referencePosition: (edgeIndex: number) => {
            const b = edges[edgeIndex]
            let unitA: Unit, unitB: Unit
            let indexA: StructureElement.UnitIndex, indexB: StructureElement.UnitIndex
            if (b.unitA.id < b.unitB.id) {
                unitA = b.unitA, unitB = b.unitB
                indexA = b.indexA, indexB = b.indexB
            } else if (b.unitA.id > b.unitB.id) {
                unitA = b.unitB, unitB = b.unitA
                indexA = b.indexB, indexB = b.indexA
            } else {
                throw new Error('same units in createInterUnitBondCylinderMesh')
            }
            return setRefPosition(tmpRef, structure, unitA, indexA) || setRefPosition(tmpRef, structure, unitB, indexB)
        },
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const b = edges[edgeIndex]
            const uA = b.unitA, uB = b.unitB
            uA.conformation.position(uA.elements[b.indexA], posA)
            uB.conformation.position(uB.elements[b.indexB], posB)
        },
        order: (edgeIndex: number) => edges[edgeIndex].props.order,
        flags: (edgeIndex: number) => BitFlags.create(edges[edgeIndex].props.flag),
        radius: (edgeIndex: number) => {
            const b = edges[edgeIndex]
            tmpLoc.unit = b.unitA
            tmpLoc.element = b.unitA.elements[b.indexA]
            const sizeA = theme.size.size(tmpLoc)
            tmpLoc.unit = b.unitB
            tmpLoc.element = b.unitB.elements[b.indexB]
            const sizeB = theme.size.size(tmpLoc)
            return Math.min(sizeA, sizeB) * sizeFactor * sizeAspectRatio
        },
        ignore: ignoreHydrogens ? (edgeIndex: number) => {
            const b = edges[edgeIndex]
            const uA = b.unitA, uB = b.unitB
            return isHydrogen(uA, uA.elements[b.indexA]) || isHydrogen(uB, uB.elements[b.indexB])
        } : () => false
    }

    return createBondCylinderMesh(ctx, builderProps, props, mesh)
}

export const InterUnitBondParams = {
    ...ComplexMeshParams,
    ...BondCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2/3, { min: 0, max: 3, step: 0.01 }),
    ignoreHydrogens: PD.Boolean(false),
}
export type InterUnitBondParams = typeof InterUnitBondParams

export function InterUnitBondVisual(materialId: number): ComplexVisual<InterUnitBondParams> {
    return ComplexMeshVisual<InterUnitBondParams>({
        defaultProps: PD.getDefaultValues(InterUnitBondParams),
        createGeometry: createInterUnitBondCylinderMesh,
        createLocationIterator: BondIterator.fromStructure,
        getLoci: getBondLoci,
        eachLocation: eachBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<InterUnitBondParams>, currentProps: PD.Values<InterUnitBondParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.sizeAspectRatio !== currentProps.sizeAspectRatio ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.bondScale !== currentProps.bondScale ||
                newProps.bondSpacing !== currentProps.bondSpacing ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.bondCap !== currentProps.bondCap
            )
        }
    }, materialId)
}

function getBondLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        const bond = structure.interUnitBonds.edges[groupId]
        return Bond.Loci(structure, [
            Bond.Location(
                bond.unitA, bond.indexA as StructureElement.UnitIndex,
                bond.unitB, bond.indexB as StructureElement.UnitIndex
            ),
            Bond.Location(
                bond.unitB, bond.indexB as StructureElement.UnitIndex,
                bond.unitA, bond.indexA as StructureElement.UnitIndex
            )
        ])
    }
    return EmptyLoci
}

function eachBond(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false
    if (Bond.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false
        for (const b of loci.bonds) {
            const idx = structure.interUnitBonds.getBondIndexFromLocation(b)
            if (idx !== -1) {
                if (apply(Interval.ofSingleton(idx))) changed = true
            }
        }
    } else if (StructureElement.Loci.is(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false
        if (loci.elements.length === 1) return false // only a single unit

        const map = new Map<number, OrderedSet<StructureElement.UnitIndex>>()
        for (const e of loci.elements) map.set(e.unit.id, e.indices)

        for (const e of loci.elements) {
            const { unit } = e
            if (!Unit.isAtomic(unit)) continue
            structure.interUnitBonds.getConnectedUnits(unit).forEach(b => {
                const otherLociIndices = map.get(b.unitB.id)
                if (otherLociIndices) {
                    OrderedSet.forEach(e.indices, v => {
                        if (!b.connectedIndices.includes(v)) return
                        b.getEdges(v).forEach(bi => {
                            if (OrderedSet.has(otherLociIndices, bi.indexB)) {
                                const idx = structure.interUnitBonds.getEdgeIndex(v, unit, bi.indexB, b.unitB)
                                if (apply(Interval.ofSingleton(idx))) changed = true
                            }
                        })
                    })
                }
            })
        }
    }
    return changed
}