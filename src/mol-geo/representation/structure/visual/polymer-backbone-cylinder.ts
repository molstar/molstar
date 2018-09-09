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
import { PolymerBackboneIterator } from './util/polymer';
import { getElementLoci, markElement, StructureElementIterator } from './util/element';
import { Vec3 } from 'mol-math/linear-algebra';
import { DefaultUnitsMeshProps, UnitsMeshVisual } from '../units-visual';
import { SizeThemeProps, SizeTheme } from 'mol-view/theme/size';
import { CylinderProps } from '../../../primitive/cylinder';
import { OrderedSet } from 'mol-data/int';
import { addCylinder } from '../../../mesh/builder/cylinder';

export interface PolymerBackboneCylinderProps {
    sizeTheme: SizeThemeProps
    radialSegments: number
}

async function createPolymerBackboneCylinderMesh(ctx: RuntimeContext, unit: Unit, props: PolymerBackboneCylinderProps, mesh?: Mesh) {
    const polymerElementCount = unit.polymerElements.length
    if (!polymerElementCount) return Mesh.createEmpty(mesh)

    const sizeTheme = SizeTheme(props.sizeTheme)
    const { radialSegments } = props

    const vertexCountEstimate = radialSegments * 2 * polymerElementCount * 2
    const builder = MeshBuilder.create(vertexCountEstimate, vertexCountEstimate / 10, mesh)

    const { elements } = unit
    const pos = unit.conformation.invariantPosition
    const pA = Vec3.zero()
    const pB = Vec3.zero()
    const cylinderProps: CylinderProps = { radiusTop: 1, radiusBottom: 1, radialSegments }

    let i = 0
    const polymerBackboneIt = PolymerBackboneIterator(unit)
    while (polymerBackboneIt.hasNext) {
        const { centerA, centerB } = polymerBackboneIt.move()
        pos(centerA.element, pA)
        pos(centerB.element, pB)

        cylinderProps.radiusTop = cylinderProps.radiusBottom = sizeTheme.size(centerA)
        builder.setGroup(OrderedSet.indexOf(elements, centerA.element))
        addCylinder(builder, pA, pB, 0.5, cylinderProps)

        cylinderProps.radiusTop = cylinderProps.radiusBottom = sizeTheme.size(centerB)
        builder.setGroup(OrderedSet.indexOf(elements, centerB.element))
        addCylinder(builder, pB, pA, 0.5, cylinderProps)

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
        // TODO create a specialized location iterator
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: MeshUpdateState, newProps: PolymerBackboneProps, currentProps: PolymerBackboneProps) => {
            state.createMesh = newProps.radialSegments !== currentProps.radialSegments
        }
    })
}