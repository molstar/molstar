/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual, VisualUpdateState } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../geometry/mesh/mesh';
import { UnitsMeshVisual, DefaultUnitsMeshProps } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { computeMarchingCubes } from '../../../util/marching-cubes/algorithm';
import { computeGaussianDensity, DefaultGaussianDensityProps } from './util/gaussian';

async function createGaussianSurfaceMesh(ctx: RuntimeContext, unit: Unit, structure: Structure, props: GaussianSurfaceProps, mesh?: Mesh): Promise<Mesh> {
    const { smoothness } = props
    const { transform, field } = await computeGaussianDensity(unit, structure, props).runAsChild(ctx)

    const surface = await computeMarchingCubes({
        isoLevel: Math.exp(-smoothness),
        scalarField: field,
        oldSurface: mesh
    }).runAsChild(ctx)

    Mesh.transformImmediate(surface, transform)
    Mesh.computeNormalsImmediate(surface)

    return surface;
}

export const DefaultGaussianSurfaceProps = {
    ...DefaultUnitsMeshProps,
    ...DefaultGaussianDensityProps,

    flipSided: true, // TODO should not be required
}
export type GaussianSurfaceProps = typeof DefaultGaussianSurfaceProps

export function GaussianSurfaceVisual(): UnitsVisual<GaussianSurfaceProps> {
    return UnitsMeshVisual<GaussianSurfaceProps>({
        defaultProps: DefaultGaussianSurfaceProps,
        createMesh: createGaussianSurfaceMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: GaussianSurfaceProps, currentProps: GaussianSurfaceProps) => {
            if (newProps.resolutionFactor !== currentProps.resolutionFactor) state.createGeometry = true
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true
        }
    })
}