/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { PolymerGapIterator, PolymerGapLocationIterator, markPolymerGapElement, getPolymerGapElementLoci } from './util/polymer';
import { Vec3 } from 'mol-math/linear-algebra';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { LinkCylinderParams } from './util/link';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { CylinderProps } from 'mol-geo/primitive/cylinder';
import { addSphere } from 'mol-geo/geometry/mesh/builder/sphere';
import { addFixedCountDashedCylinder } from 'mol-geo/geometry/mesh/builder/cylinder';
import { VisualContext } from 'mol-repr/visual';
import { Theme } from 'mol-theme/theme';

const segmentCount = 10

export const PolymerGapCylinderParams = {
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    radialSegments: PD.Numeric(16, { min: 3, max: 56, step: 1 }),
}
export const DefaultPolymerGapCylinderProps = PD.getDefaultValues(PolymerGapCylinderParams)
export type PolymerGapCylinderProps = typeof DefaultPolymerGapCylinderProps

function createPolymerGapCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PolymerGapCylinderProps, mesh?: Mesh) {
    const polymerGapCount = unit.gapElements.length
    if (!polymerGapCount) return Mesh.createEmpty(mesh)

    const { sizeFactor, radialSegments } = props

    const vertexCountEstimate = segmentCount * radialSegments * 2 * polymerGapCount * 2
    const builderState = MeshBuilder.createState(vertexCountEstimate, vertexCountEstimate / 10, mesh)

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
            builderState.currentGroup = i
            pos(centerA.element, pA)
            addSphere(builderState, pA, 0.6, 0)
        } else {
            pos(centerA.element, pA)
            pos(centerB.element, pB)

            cylinderProps.radiusTop = cylinderProps.radiusBottom = theme.size.size(centerA) * sizeFactor
            builderState.currentGroup = i
            addFixedCountDashedCylinder(builderState, pA, pB, 0.5, segmentCount, cylinderProps)

            cylinderProps.radiusTop = cylinderProps.radiusBottom = theme.size.size(centerB) * sizeFactor
            builderState.currentGroup = i + 1
            addFixedCountDashedCylinder(builderState, pB, pA, 0.5, segmentCount, cylinderProps)
        }

        i += 2
    }

    return MeshBuilder.getMesh(builderState)
}

export const InterUnitLinkParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
}
export const DefaultIntraUnitLinkProps = PD.getDefaultValues(InterUnitLinkParams)
export type IntraUnitLinkProps = typeof DefaultIntraUnitLinkProps

export const PolymerGapParams = {
    ...UnitsMeshParams,
    ...PolymerGapCylinderParams
}
export type PolymerGapParams = typeof PolymerGapParams

export function PolymerGapVisual(): UnitsVisual<PolymerGapParams> {
    return UnitsMeshVisual<PolymerGapParams>({
        defaultProps: PD.getDefaultValues(PolymerGapParams),
        createGeometry: createPolymerGapCylinderMesh,
        createLocationIterator: PolymerGapLocationIterator.fromGroup,
        getLoci: getPolymerGapElementLoci,
        mark: markPolymerGapElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<PolymerGapParams>, currentProps: PD.Values<PolymerGapParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments
            )
        }
    })
}