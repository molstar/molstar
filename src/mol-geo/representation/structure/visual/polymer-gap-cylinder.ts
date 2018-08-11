/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit } from 'mol-model/structure';
import { UnitsVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../shape/mesh';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { getPolymerGapCount, PolymerGapIterator } from './util/polymer';
import { getElementLoci, markElement } from './util/element';
import { Vec3 } from 'mol-math/linear-algebra';
import { StructureElementIterator } from './util/location-iterator';
import { UnitsMeshVisual, DefaultUnitsMeshProps } from '../units-visual';

async function createPolymerGapCylinderMesh(ctx: RuntimeContext, unit: Unit, props: {}, mesh?: Mesh) {
    const polymerGapCount = getPolymerGapCount(unit)
    if (!polymerGapCount) return Mesh.createEmpty(mesh)
    console.log('polymerGapCount', polymerGapCount)

    // TODO better vertex count estimates
    const builder = MeshBuilder.create(polymerGapCount * 30, polymerGapCount * 30 / 2, mesh)

    const { elements } = unit
    const pos = unit.conformation.invariantPosition
    const pA = Vec3.zero()
    const pB = Vec3.zero()

    let i = 0
    const polymerGapIt = PolymerGapIterator(unit)
    while (polymerGapIt.hasNext) {
        // TODO size theme
        const { centerA, centerB } = polymerGapIt.move()
        if (centerA.element === centerB.element) {
            builder.setId(centerA.element)
            pos(elements[centerA.element], pA)
            builder.addSphere(pA, 0.6, 0)
        } else {
            pos(elements[centerA.element], pA)
            pos(elements[centerB.element], pB)
            builder.setId(centerA.element)
            builder.addFixedCountDashedCylinder(pA, pB, 0.5, 10, { radiusTop: 0.2, radiusBottom: 0.2 })
            builder.setId(centerB.element)
            builder.addFixedCountDashedCylinder(pB, pA, 0.5, 10, { radiusTop: 0.2, radiusBottom: 0.2 })
        }

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Gap mesh', current: i, max: polymerGapCount });
        }
        ++i
    }

    return builder.getMesh()
}

export const DefaultPolymerGapProps = {
    ...DefaultUnitsMeshProps
}
export type PolymerGapProps = typeof DefaultPolymerGapProps

export function PolymerGapVisual(): UnitsVisual<PolymerGapProps> {
    return UnitsMeshVisual<PolymerGapProps>({
        defaultProps: DefaultPolymerGapProps,
        createMesh: createPolymerGapCylinderMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement
    })
}