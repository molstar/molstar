/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from 'mol-geo/util/marching-cubes/algorithm';
import { VisualContext } from 'mol-repr/visual';
import { Theme } from 'mol-theme/theme';
import { GaussianDensityProps, computeUnitGaussianDensity, GaussianDensityParams } from './util/gaussian';

async function createGaussianSurfaceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityProps, mesh?: Mesh): Promise<Mesh> {
    const { smoothness } = props
    const { transform, field, idField } = await computeUnitGaussianDensity(unit, props, ctx.webgl).runInContext(ctx.runtime)

    const params = {
        isoLevel: Math.exp(-smoothness),
        scalarField: field,
        idField
    }
    const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime)

    Mesh.transformImmediate(surface, transform)
    Mesh.computeNormalsImmediate(surface)
    Mesh.uniformTriangleGroup(surface)

    return surface;
}

export const GaussianSurfaceParams = {
    ...UnitsMeshParams,
    ...GaussianDensityParams,
}
export type GaussianSurfaceParams = typeof GaussianSurfaceParams

export function GaussianSurfaceVisual(): UnitsVisual<GaussianSurfaceParams> {
    return UnitsMeshVisual<GaussianSurfaceParams>({
        defaultProps: PD.getDefaultValues(GaussianSurfaceParams),
        createGeometry: createGaussianSurfaceMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<GaussianSurfaceParams>, currentProps: PD.Values<GaussianSurfaceParams>) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true
            if (newProps.useGpu !== currentProps.useGpu) state.createGeometry = true
        }
    })
}