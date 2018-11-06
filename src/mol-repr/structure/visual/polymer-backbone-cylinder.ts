/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../index';
import { VisualUpdateState } from '../../util';
import { PolymerBackboneIterator } from './util/polymer';
import { getElementLoci, markElement, StructureElementIterator } from './util/element';
import { Vec3 } from 'mol-math/linear-algebra';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { SizeTheme, SizeThemeOptions, SizeThemeName } from 'mol-theme/size';
import { OrderedSet } from 'mol-data/int';
import { paramDefaultValues, NumberParam, SelectParam } from 'mol-util/parameter';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { CylinderProps } from 'mol-geo/primitive/cylinder';
import { addCylinder } from 'mol-geo/geometry/mesh/builder/cylinder';
import { VisualContext } from 'mol-repr';

export const PolymerBackboneCylinderParams = {
    sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
    sizeValue: NumberParam('Size Value', '', 1, 0, 10, 0.1),
    radialSegments: NumberParam('Radial Segments', '', 16, 3, 56, 1),
}
export const DefaultPolymerBackboneCylinderProps = paramDefaultValues(PolymerBackboneCylinderParams)
export type PolymerBackboneCylinderProps = typeof DefaultPolymerBackboneCylinderProps

async function createPolymerBackboneCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, props: PolymerBackboneCylinderProps, mesh?: Mesh) {
    const polymerElementCount = unit.polymerElements.length
    if (!polymerElementCount) return Mesh.createEmpty(mesh)

    const sizeTheme = SizeTheme({ name: props.sizeTheme, value: props.sizeValue })
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

        if (i % 10000 === 0 && ctx.runtime.shouldUpdate) {
            await ctx.runtime.update({ message: 'Backbone mesh', current: i, max: polymerElementCount });
        }
        ++i
    }

    return builder.getMesh()
}

export const PolymerBackboneParams = {
    ...UnitsMeshParams,
    ...PolymerBackboneCylinderParams,
}
export const DefaultPolymerBackboneProps = paramDefaultValues(PolymerBackboneParams)
export type PolymerBackboneProps = typeof DefaultPolymerBackboneProps

export function PolymerBackboneVisual(): UnitsVisual<PolymerBackboneProps> {
    return UnitsMeshVisual<PolymerBackboneProps>({
        defaultProps: DefaultPolymerBackboneProps,
        createGeometry: createPolymerBackboneCylinderMesh,
        // TODO create a specialized location iterator
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: PolymerBackboneProps, currentProps: PolymerBackboneProps) => {
            state.createGeometry = newProps.radialSegments !== currentProps.radialSegments
        }
    })
}