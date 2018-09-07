/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit } from 'mol-model/structure';
import { UnitsVisual, MeshUpdateState } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../mesh/mesh';
import { MeshBuilder } from '../../../mesh/mesh-builder';
import { PolymerGapIterator, PolymerGapLocationIterator } from './util/polymer';
import { getElementLoci, markElement } from './util/element';
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
    const polymerGapCount = unit.gapElements.length
    if (!polymerGapCount) return Mesh.createEmpty(mesh)

    const sizeTheme = SizeTheme(props.sizeTheme)
    const { radialSegments } = props

    const vertexCountEstimate = segmentCount * radialSegments * 2 * polymerGapCount * 2
    const builder = MeshBuilder.create(vertexCountEstimate, vertexCountEstimate / 10, mesh)

    const pos = unit.conformation.invariantPosition
    const pA = Vec3.zero()
    const pB = Vec3.zero()
    const cylinderProps: CylinderProps = {
        radiusTop: 1, radiusBottom: 1, topCap: true, bottomCap: true, radialSegments
    }

    let i = 0
    const polymerGapIt = PolymerGapIterator(unit)
    while (polymerGapIt.hasNext) {
        const { centerA, centerB } = polymerGapIt.move()
        if (centerA.element === centerB.element) {
            builder.setGroup(i)
            pos(centerA.element, pA)
            addSphere(builder, pA, 0.6, 0)
        } else {
            pos(centerA.element, pA)
            pos(centerB.element, pB)

            cylinderProps.radiusTop = cylinderProps.radiusBottom = sizeTheme.size(centerA)
            builder.setGroup(i)
            addFixedCountDashedCylinder(builder, pA, pB, 0.5, segmentCount, cylinderProps)

            cylinderProps.radiusTop = cylinderProps.radiusBottom = sizeTheme.size(centerB)
            builder.setGroup(i + 1)
            addFixedCountDashedCylinder(builder, pB, pA, 0.5, segmentCount, cylinderProps)
        }

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Gap mesh', current: i, max: polymerGapCount });
        }
        i += 2
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
        createLocationIterator: PolymerGapLocationIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: MeshUpdateState, newProps: PolymerGapProps, currentProps: PolymerGapProps) => {
            state.createMesh = newProps.radialSegments !== currentProps.radialSegments
        }
    })
}