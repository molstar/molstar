/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../index';
import { VisualUpdateState } from '../../util';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { GaussianDensityProps, GaussianDensityParams } from 'mol-model/structure/structure/unit/gaussian-density';
import { paramDefaultValues } from 'mol-util/parameter';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from 'mol-geo/util/marching-cubes/algorithm';
import { VisualContext } from 'mol-repr';
import { Theme } from 'mol-geo/geometry/geometry';

async function createGaussianSurfaceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityProps, mesh?: Mesh): Promise<Mesh> {
    const { smoothness } = props
    const { transform, field, idField } = await unit.computeGaussianDensity(props, ctx.runtime, ctx.webgl)

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
export const DefaultGaussianSurfaceProps = paramDefaultValues(GaussianSurfaceParams)
export type GaussianSurfaceProps = typeof DefaultGaussianSurfaceProps

export function GaussianSurfaceVisual(): UnitsVisual<GaussianSurfaceProps> {
    return UnitsMeshVisual<GaussianSurfaceProps>({
        defaultProps: DefaultGaussianSurfaceProps,
        createGeometry: createGaussianSurfaceMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: GaussianSurfaceProps, currentProps: GaussianSurfaceProps) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true
            if (newProps.useGpu !== currentProps.useGpu) state.createGeometry = true
            if (newProps.ignoreCache !== currentProps.ignoreCache) state.createGeometry = true
        }
    })
}