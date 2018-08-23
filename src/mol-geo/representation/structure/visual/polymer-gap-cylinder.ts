/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement } from 'mol-model/structure';
import { UnitsVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../mesh/mesh';
import { MeshBuilder } from '../../../mesh/mesh-builder';
import { getPolymerGapCount, PolymerGapIterator } from './util/polymer';
import { getElementLoci, markElement, StructureElementIterator } from './util/element';
import { Vec3 } from 'mol-math/linear-algebra';
import { UnitsMeshVisual, DefaultUnitsMeshProps } from '../units-visual';
import { SizeThemeProps, SizeTheme } from 'mol-view/theme/size';
import { CylinderProps } from '../../../primitive/cylinder';
import { addSphere } from '../../../mesh/builder/sphere';
import { addFixedCountDashedCylinder } from '../../../mesh/builder/cylinder';

const segmentCount = 10

export interface PolymerGapCylinderProps {
    sizeTheme: SizeThemeProps
    radialSegments: number
}

async function createPolymerGapCylinderMesh(ctx: RuntimeContext, unit: Unit, props: PolymerGapCylinderProps, mesh?: Mesh) {
    const polymerGapCount = getPolymerGapCount(unit)
    if (!polymerGapCount) return Mesh.createEmpty(mesh)

    const sizeTheme = SizeTheme(props.sizeTheme)
    const { radialSegments } = props

    const vertexCountEstimate = segmentCount * radialSegments * 2 * polymerGapCount * 2
    const builder = MeshBuilder.create(vertexCountEstimate, vertexCountEstimate / 10, mesh)

    const { elements } = unit
    const pos = unit.conformation.invariantPosition
    const pA = Vec3.zero()
    const pB = Vec3.zero()
    const l = StructureElement.create(unit)
    const cylinderProps: CylinderProps = {
        radiusTop: 1, radiusBottom: 1, topCap: true, bottomCap: true, radialSegments
    }

    let i = 0
    const polymerGapIt = PolymerGapIterator(unit)
    while (polymerGapIt.hasNext) {
        const { centerA, centerB } = polymerGapIt.move()
        if (centerA.element === centerB.element) {
            builder.setGroup(centerA.element)
            pos(elements[centerA.element], pA)
            addSphere(builder, pA, 0.6, 0)
        } else {
            const elmA = elements[centerA.element]
            const elmB = elements[centerB.element]
            pos(elmA, pA)
            pos(elmB, pB)

            l.element = elmA
            cylinderProps.radiusTop = cylinderProps.radiusBottom = sizeTheme.size(l)
            builder.setGroup(centerA.element)
            addFixedCountDashedCylinder(builder, pA, pB, 0.5, segmentCount, cylinderProps)

            l.element = elmB
            cylinderProps.radiusTop = cylinderProps.radiusBottom = sizeTheme.size(l)
            builder.setGroup(centerB.element)
            addFixedCountDashedCylinder(builder, pB, pA, 0.5, segmentCount, cylinderProps)
        }

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Gap mesh', current: i, max: polymerGapCount });
        }
        ++i
    }

    return builder.getMesh()
}

export const DefaultPolymerGapProps = {
    ...DefaultUnitsMeshProps,
    radialSegments: 16
}
export type PolymerGapProps = typeof DefaultPolymerGapProps

export function PolymerGapVisual(): UnitsVisual<PolymerGapProps> {
    return UnitsMeshVisual<PolymerGapProps>({
        defaultProps: DefaultPolymerGapProps,
        createMesh: createPolymerGapCylinderMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: () => {}
    })
}