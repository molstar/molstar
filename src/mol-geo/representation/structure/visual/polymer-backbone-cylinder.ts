/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement } from 'mol-model/structure';
import { UnitsVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../shape/mesh';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { getPolymerElementCount, PolymerBackboneIterator } from './util/polymer';
import { getElementLoci, markElement } from './util/element';
import { Vec3 } from 'mol-math/linear-algebra';
import { StructureElementIterator } from './util/location-iterator';
import { DefaultUnitsMeshProps, UnitsMeshVisual } from '../units-visual';
import { SizeThemeProps, SizeTheme } from 'mol-view/theme/size';
import { CylinderProps } from '../../../primitive/cylinder';

export interface PolymerBackboneCylinderProps {
    sizeTheme: SizeThemeProps
    radialSegments: number
}

async function createPolymerBackboneCylinderMesh(ctx: RuntimeContext, unit: Unit, props: PolymerBackboneCylinderProps, mesh?: Mesh) {
    const polymerElementCount = getPolymerElementCount(unit)
    if (!polymerElementCount) return Mesh.createEmpty(mesh)

    const sizeTheme = SizeTheme(props.sizeTheme)
    const { radialSegments } = props

    const vertexCountEstimate = radialSegments * 2 * polymerElementCount * 2
    const builder = MeshBuilder.create(vertexCountEstimate, vertexCountEstimate / 10, mesh)

    const { elements } = unit
    const pos = unit.conformation.invariantPosition
    const pA = Vec3.zero()
    const pB = Vec3.zero()
    const l = StructureElement.create(unit)
    const cylinderProps: CylinderProps = { radiusTop: 1, radiusBottom: 1, radialSegments }

    let i = 0
    const polymerBackboneIt = PolymerBackboneIterator(unit)
    while (polymerBackboneIt.hasNext) {
        const { centerA, centerB } = polymerBackboneIt.move()
        const elmA = elements[centerA.element]
        const elmB = elements[centerB.element]
        pos(elmA, pA)
        pos(elmB, pB)

        l.element = elmA
        cylinderProps.radiusTop = cylinderProps.radiusBottom = sizeTheme.size(l)
        builder.setId(centerA.element)
        builder.addCylinder(pA, pB, 0.5, cylinderProps)

        l.element = elmB
        cylinderProps.radiusTop = cylinderProps.radiusBottom = sizeTheme.size(l)
        builder.setId(centerB.element)
        builder.addCylinder(pB, pA, 0.5, cylinderProps)

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Backbone mesh', current: i, max: polymerElementCount });
        }
        ++i
    }

    return builder.getMesh()
}

export const DefaultPolymerBackboneProps = {
    ...DefaultUnitsMeshProps,
    radialSegments: 16
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