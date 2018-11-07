/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../index';
import { VisualUpdateState } from '../../util';
import { PolymerGapIterator, PolymerGapLocationIterator, markPolymerGapElement, getPolymerGapElementLoci } from './util/polymer';
import { Vec3 } from 'mol-math/linear-algebra';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { SizeThemeOptions, SizeThemeName } from 'mol-theme/size';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { LinkCylinderParams } from './util/link';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { CylinderProps } from 'mol-geo/primitive/cylinder';
import { addSphere } from 'mol-geo/geometry/mesh/builder/sphere';
import { addFixedCountDashedCylinder } from 'mol-geo/geometry/mesh/builder/cylinder';
import { VisualContext } from 'mol-repr';
import { Theme } from 'mol-geo/geometry/geometry';

const segmentCount = 10

export const PolymerGapCylinderParams = {
    radialSegments: PD.Numeric('Radial Segments', '', 16, 3, 56, 1),
}
export const DefaultPolymerGapCylinderProps = PD.getDefaultValues(PolymerGapCylinderParams)
export type PolymerGapCylinderProps = typeof DefaultPolymerGapCylinderProps

async function createPolymerGapCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PolymerGapCylinderProps, mesh?: Mesh) {
    const polymerGapCount = unit.gapElements.length
    if (!polymerGapCount) return Mesh.createEmpty(mesh)

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

            cylinderProps.radiusTop = cylinderProps.radiusBottom = theme.size.size(centerA)
            builder.setGroup(i)
            addFixedCountDashedCylinder(builder, pA, pB, 0.5, segmentCount, cylinderProps)

            cylinderProps.radiusTop = cylinderProps.radiusBottom = theme.size.size(centerB)
            builder.setGroup(i + 1)
            addFixedCountDashedCylinder(builder, pB, pA, 0.5, segmentCount, cylinderProps)
        }

        if (i % 10000 === 0 && ctx.runtime.shouldUpdate) {
            await ctx.runtime.update({ message: 'Gap mesh', current: i, max: polymerGapCount });
        }
        i += 2
    }

    return builder.getMesh()
}

export const InterUnitLinkParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
    sizeTheme: PD.Select<SizeThemeName>('Size Theme', '', 'physical', SizeThemeOptions),
    sizeValue: PD.Numeric('Size Value', '', 1, 0, 20, 0.1),
    sizeFactor: PD.Numeric('Size Factor', '', 0.3, 0, 10, 0.1),
}
export const DefaultIntraUnitLinkProps = PD.getDefaultValues(InterUnitLinkParams)
export type IntraUnitLinkProps = typeof DefaultIntraUnitLinkProps

export const PolymerGapParams = {
    ...UnitsMeshParams,
    ...PolymerGapCylinderParams
}
export const DefaultPolymerGapProps = PD.getDefaultValues(PolymerGapParams)
export type PolymerGapProps = typeof DefaultPolymerGapProps

export function PolymerGapVisual(): UnitsVisual<PolymerGapProps> {
    return UnitsMeshVisual<PolymerGapProps>({
        defaultProps: DefaultPolymerGapProps,
        createGeometry: createPolymerGapCylinderMesh,
        createLocationIterator: PolymerGapLocationIterator.fromGroup,
        getLoci: getPolymerGapElementLoci,
        mark: markPolymerGapElement,
        setUpdateState: (state: VisualUpdateState, newProps: PolymerGapProps, currentProps: PolymerGapProps) => {
            state.createGeometry = newProps.radialSegments !== currentProps.radialSegments
        }
    })
}