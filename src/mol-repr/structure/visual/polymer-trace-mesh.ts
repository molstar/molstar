/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { PolymerTraceIterator, createCurveSegmentState, interpolateCurveSegment, PolymerLocationIterator, getPolymerElementLoci, markPolymerElement, interpolateSizes } from './util/polymer';
import { SecondaryStructureType, isNucleic } from 'mol-model/structure/model/types';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { addSheet } from 'mol-geo/geometry/mesh/builder/sheet';
import { addTube } from 'mol-geo/geometry/mesh/builder/tube';
import { VisualContext } from 'mol-repr/visual';
import { Theme } from 'mol-theme/theme';

export const PolymerTraceMeshParams = {
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    linearSegments: PD.Numeric(8, { min: 1, max: 48, step: 1 }),
    radialSegments: PD.Numeric(16, { min: 3, max: 56, step: 1 }),
    aspectRatio: PD.Numeric(5, { min: 0.1, max: 5, step: 0.1 }),
    arrowFactor: PD.Numeric(1.5, { min: 0.1, max: 5, step: 0.1 }),
}
export const DefaultPolymerTraceMeshProps = PD.getDefaultValues(PolymerTraceMeshParams)
export type PolymerTraceMeshProps = typeof DefaultPolymerTraceMeshProps

// TODO handle polymer ends properly

function createPolymerTraceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PolymerTraceMeshProps, mesh?: Mesh) {
    const polymerElementCount = unit.polymerElements.length

    if (!polymerElementCount) return Mesh.createEmpty(mesh)
    const { sizeFactor, linearSegments, radialSegments, aspectRatio, arrowFactor } = props

    const vertexCount = linearSegments * radialSegments * polymerElementCount + (radialSegments + 1) * polymerElementCount * 2
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 10, mesh)

    const isCoarse = Unit.isCoarse(unit)
    const state = createCurveSegmentState(linearSegments)
    const { curvePoints, normalVectors, binormalVectors, widthValues, heightValues } = state

    let i = 0
    const polymerTraceIt = PolymerTraceIterator(unit)
    while (polymerTraceIt.hasNext) {
        const v = polymerTraceIt.move()
        builderState.currentGroup = i

        const isNucleicType = isNucleic(v.moleculeType)
        const isSheet = SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Beta)
        const isHelix = SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Helix)
        const tension = isNucleicType ? 0.5 : 0.9
        const shift = isNucleicType ? 0.3 : 0.5

        interpolateCurveSegment(state, v, tension, shift)

        let w0 = theme.size.size(v.centerPrev) * sizeFactor
        let w1 = theme.size.size(v.center) * sizeFactor
        let w2 = theme.size.size(v.centerNext) * sizeFactor
        if (isCoarse) {
            w0 *= aspectRatio / 2
            w1 *= aspectRatio / 2
            w2 *= aspectRatio / 2
        }

        if (isSheet) {
            const h1 = w1 * aspectRatio
            const arrowHeight = v.secStrucLast ? h1 * arrowFactor : 0
            addSheet(builderState, curvePoints, normalVectors, binormalVectors, linearSegments, w1, h1, arrowHeight, v.secStrucFirst, v.secStrucLast)
        } else {
            let h0: number, h1: number, h2: number
            if (isHelix && !v.isCoarseBackbone) {
                h0 = w0 * aspectRatio
                h1 = w1 * aspectRatio
                h2 = w2 * aspectRatio
            } else if (isNucleicType && !v.isCoarseBackbone) {
                h0 = w0 * aspectRatio;
                [w0, h0] = [h0, w0]
                h1 = w1 * aspectRatio;
                [w1, h1] = [h1, w1]
                h2 = w2 * aspectRatio;
                [w2, h2] = [h2, w2]
            } else {
                h0 = w0
                h1 = w1
                h2 = w2
            }
            interpolateSizes(state, w0, w1, w2, h0, h1, h2, shift)
            addTube(builderState, curvePoints, normalVectors, binormalVectors, linearSegments, radialSegments, widthValues, heightValues, 1, v.secStrucFirst || v.coarseBackboneFirst, v.secStrucLast || v.coarseBackboneLast)
        }

        ++i
    }

    return MeshBuilder.getMesh(builderState)
}

export const PolymerTraceParams = {
    ...UnitsMeshParams,
    ...PolymerTraceMeshParams
}
export type PolymerTraceParams = typeof PolymerTraceParams

export function PolymerTraceVisual(): UnitsVisual<PolymerTraceParams> {
    return UnitsMeshVisual<PolymerTraceParams>({
        defaultProps: PD.getDefaultValues(PolymerTraceParams),
        createGeometry: createPolymerTraceMesh,
        createLocationIterator: PolymerLocationIterator.fromGroup,
        getLoci: getPolymerElementLoci,
        mark: markPolymerElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<PolymerTraceParams>, currentProps: PD.Values<PolymerTraceParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.linearSegments !== currentProps.linearSegments ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.aspectRatio !== currentProps.aspectRatio ||
                newProps.arrowFactor !== currentProps.arrowFactor
            )
        }
    })
}