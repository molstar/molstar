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
import { getPolymerElementCount, PolymerBackboneIterator } from './util/polymer';
import { getElementLoci, markElement } from './util/element';
import { Vec3 } from 'mol-math/linear-algebra';
import { StructureElementIterator } from './util/location-iterator';
import { DefaultUnitsMeshProps, UnitsMeshVisual } from '../units-visual';

async function createPolymerBackboneCylinderMesh(ctx: RuntimeContext, unit: Unit, props: {}, mesh?: Mesh) {
    const polymerElementCount = getPolymerElementCount(unit)
    if (!polymerElementCount) return Mesh.createEmpty(mesh)
    console.log('polymerElementCount backbone', polymerElementCount)

    // TODO better vertex count estimates
    const builder = MeshBuilder.create(polymerElementCount * 30, polymerElementCount * 30 / 2, mesh)

    const { elements } = unit
    const pos = unit.conformation.invariantPosition
    const pA = Vec3.zero()
    const pB = Vec3.zero()

    let i = 0
    const polymerBackboneIt = PolymerBackboneIterator(unit)
    while (polymerBackboneIt.hasNext) {
        // TODO size theme
        const { centerA, centerB } = polymerBackboneIt.move()
        pos(elements[centerA.element], pA)
        pos(elements[centerB.element], pB)
        builder.setId(centerA.element)
        builder.addCylinder(pA, pB, 0.5, { radiusTop: 0.2, radiusBottom: 0.2 })
        builder.setId(centerB.element)
        builder.addCylinder(pB, pA, 0.5, { radiusTop: 0.2, radiusBottom: 0.2 })

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Backbone mesh', current: i, max: polymerElementCount });
        }
        ++i
    }

    return builder.getMesh()
}

export const DefaultPolymerBackboneProps = {
    ...DefaultUnitsMeshProps
}
export type PolymerBackboneProps = typeof DefaultPolymerBackboneProps

export function PolymerBackboneVisual(): UnitsVisual<PolymerBackboneProps> {
    return UnitsMeshVisual<PolymerBackboneProps>({
        defaultProps: DefaultPolymerBackboneProps,
        createMesh: createPolymerBackboneCylinderMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: () => {}
    })
}