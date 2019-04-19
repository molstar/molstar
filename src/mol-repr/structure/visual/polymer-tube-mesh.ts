/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { PolymerTraceIterator, createCurveSegmentState, interpolateCurveSegment, PolymerLocationIterator, getPolymerElementLoci, eachPolymerElement, interpolateSizes } from './util/polymer';
import { isNucleic } from 'mol-model/structure/model/types';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { addTube } from 'mol-geo/geometry/mesh/builder/tube';
import { VisualContext } from 'mol-repr/visual';
import { Theme } from 'mol-theme/theme';

export const PolymerTubeMeshParams = {
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    linearSegments: PD.Numeric(8, { min: 1, max: 48, step: 1 }),
    radialSegments: PD.Numeric(16, { min: 3, max: 56, step: 1 }),
}
export const DefaultPolymerTubeMeshProps = PD.getDefaultValues(PolymerTubeMeshParams)
export type PolymerTubeMeshProps = typeof DefaultPolymerTubeMeshProps

// TODO handle polymer ends properly

function createPolymerTubeMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PolymerTubeMeshProps, mesh?: Mesh) {
    const polymerElementCount = unit.polymerElements.length

    if (!polymerElementCount) return Mesh.createEmpty(mesh)
    const { sizeFactor, linearSegments, radialSegments } = props

    const vertexCount = linearSegments * radialSegments * polymerElementCount + (radialSegments + 1) * polymerElementCount * 2
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 10, mesh)

    const state = createCurveSegmentState(linearSegments)
    const { curvePoints, normalVectors, binormalVectors, widthValues, heightValues } = state

    let i = 0
    const polymerTraceIt = PolymerTraceIterator(unit, structure)
    while (polymerTraceIt.hasNext) {
        const v = polymerTraceIt.move()
        builderState.currentGroup = i

        const isNucleicType = isNucleic(v.moleculeType)
        const tension = isNucleicType ? 0.5 : 0.9
        const shift = isNucleicType ? 0.3 : 0.5

        interpolateCurveSegment(state, v, tension, shift)

        let s0 = theme.size.size(v.centerPrev) * sizeFactor
        let s1 = theme.size.size(v.center) * sizeFactor
        let s2 = theme.size.size(v.centerNext) * sizeFactor

        interpolateSizes(state, s0, s1, s2, s0, s1, s2, shift)
        addTube(builderState, curvePoints, normalVectors, binormalVectors, linearSegments, radialSegments, widthValues, heightValues, 1, v.first, v.last)

        ++i
    }

    return MeshBuilder.getMesh(builderState)
}

export const PolymerTubeParams = {
    ...UnitsMeshParams,
    ...PolymerTubeMeshParams
}
export type PolymerTubeParams = typeof PolymerTubeParams

export function PolymerTubeVisual(materialId: number): UnitsVisual<PolymerTubeParams> {
    return UnitsMeshVisual<PolymerTubeParams>({
        defaultProps: PD.getDefaultValues(PolymerTubeParams),
        createGeometry: createPolymerTubeMesh,
        createLocationIterator: PolymerLocationIterator.fromGroup,
        getLoci: getPolymerElementLoci,
        eachLocation: eachPolymerElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<PolymerTubeParams>, currentProps: PD.Values<PolymerTubeParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.linearSegments !== currentProps.linearSegments ||
                newProps.radialSegments !== currentProps.radialSegments
            )
        }
    }, materialId)
}