/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit } from 'mol-model/structure';
import { UnitsVisual, MeshUpdateState } from '..';
import { RuntimeContext } from 'mol-task'
import { markElement, getElementLoci } from './util/element';
import { Mesh } from '../../../shape/mesh';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { getPolymerElementCount, PolymerTraceIterator, createCurveSegmentState, interpolateCurveSegment } from './util/polymer';
import { SecondaryStructureType, MoleculeType } from 'mol-model/structure/model/types';
import { StructureElementIterator } from './util/location-iterator';
import { UnitsMeshVisual, DefaultUnitsMeshProps } from '../units-visual';
import { SizeThemeProps, SizeTheme } from 'mol-view/theme/size';

export interface PolymerTraceMeshProps {
    sizeTheme: SizeThemeProps
    linearSegments: number
    radialSegments: number
    aspectRatio: number
    arrowFactor: number
}

// TODO handle polymer ends properly

async function createPolymerTraceMesh(ctx: RuntimeContext, unit: Unit, props: PolymerTraceMeshProps, mesh?: Mesh) {
    const polymerElementCount = getPolymerElementCount(unit)
    if (!polymerElementCount) return Mesh.createEmpty(mesh)

    const sizeTheme = SizeTheme(props.sizeTheme)
    const { linearSegments, radialSegments, aspectRatio, arrowFactor } = props

    const vertexCount = linearSegments * radialSegments * polymerElementCount + (radialSegments + 1) * polymerElementCount * 2
    const builder = MeshBuilder.create(vertexCount, vertexCount / 10, mesh)

    const state = createCurveSegmentState(linearSegments)
    const { curvePoints, normalVectors, binormalVectors } = state

    let i = 0
    const polymerTraceIt = PolymerTraceIterator(unit)
    while (polymerTraceIt.hasNext) {
        const v = polymerTraceIt.move()
        builder.setId(v.center.element)

        const isNucleic = v.moleculeType === MoleculeType.DNA || v.moleculeType === MoleculeType.RNA
        const isSheet = SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Beta)
        const isHelix = SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Helix)
        const tension = (isNucleic || isSheet) ? 0.5 : 0.9

        interpolateCurveSegment(state, v, tension)

        let width = sizeTheme.size(v.center)

        if (isSheet) {
            const height = width * aspectRatio
            const arrowHeight = v.secStrucChange ? height * arrowFactor : 0
            builder.addSheet(curvePoints, normalVectors, binormalVectors, linearSegments, width, height, arrowHeight, true, true)
        } else {
            let height: number
            if (isHelix) {
                height = width * aspectRatio
            } else if (isNucleic) {
                height = width * aspectRatio;
                [width, height] = [height, width]
            } else {
                height = width
            }
            builder.addTube(curvePoints, normalVectors, binormalVectors, linearSegments, radialSegments, width, height, 1, true, true)
        }

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Polymer trace mesh', current: i, max: polymerElementCount });
        }
        ++i
    }

    return builder.getMesh()
}

export const DefaultPolymerTraceProps = {
    ...DefaultUnitsMeshProps,
    linearSegments: 8,
    radialSegments: 12,
    aspectRatio: 5,
    arrowFactor: 1.5
}
export type PolymerTraceProps = typeof DefaultPolymerTraceProps

export function PolymerTraceVisual(): UnitsVisual<PolymerTraceProps> {
    return UnitsMeshVisual<PolymerTraceProps>({
        defaultProps: DefaultPolymerTraceProps,
        createMesh: createPolymerTraceMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: MeshUpdateState, newProps: PolymerTraceProps, currentProps: PolymerTraceProps) => {
            state.createMesh = (
                newProps.linearSegments !== currentProps.linearSegments ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.aspectRatio !== currentProps.aspectRatio ||
                newProps.arrowFactor !== currentProps.arrowFactor
            )
        }
    })
}