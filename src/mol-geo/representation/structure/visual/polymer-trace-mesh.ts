/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit } from 'mol-model/structure';
import { UnitsVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { markElement, getElementLoci } from './util/element';
import { Mesh } from '../../../shape/mesh';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { getPolymerElementCount, PolymerTraceIterator, createCurveSegmentState, interpolateCurveSegment } from './util/polymer';
import { SecondaryStructureType, MoleculeType } from 'mol-model/structure/model/types';
import { StructureElementIterator } from './util/location-iterator';
import { UnitsMeshVisual, DefaultUnitsMeshProps } from '../units-visual';

// TODO handle polymer ends properly

async function createPolymerTraceMesh(ctx: RuntimeContext, unit: Unit, props: {}, mesh?: Mesh) {
    const polymerElementCount = getPolymerElementCount(unit)
    console.log('polymerElementCount trace', polymerElementCount)
    if (!polymerElementCount) return Mesh.createEmpty(mesh)

    // TODO better vertex count estimates
    const builder = MeshBuilder.create(polymerElementCount * 30, polymerElementCount * 30 / 2, mesh)
    const linearSegments = 8
    const radialSegments = 12

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

        // console.log('ELEMENT', i)
        interpolateCurveSegment(state, v, tension)

        let width = 0.2, height = 0.2

        // TODO size theme
        if (isSheet) {
            width = 0.15; height = 1.0
            const arrowHeight = v.secStrucChange ? 1.7 : 0
            builder.addSheet(curvePoints, normalVectors, binormalVectors, linearSegments, width, height, arrowHeight, true, true)
        } else {
            if (isHelix) {
                width = 0.2; height = 1.0
            } else if (isNucleic) {
                width = 1.5; height = 0.3
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
    ...DefaultUnitsMeshProps
}
export type PolymerTraceProps = typeof DefaultPolymerTraceProps

export function PolymerTraceVisual(): UnitsVisual<PolymerTraceProps> {
    return UnitsMeshVisual<PolymerTraceProps>({
        defaultProps: DefaultPolymerTraceProps,
        createMesh: createPolymerTraceMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: () => {}
    })
}