/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../units-visual';
import { MolecularSurfaceCalculationParams } from '../../../mol-math/geometry/molecular-surface';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { computeUnitMolecularSurface, MolecularSurfaceProps } from './util/molecular-surface';
import { computeMarchingCubesMesh } from '../../../mol-geo/util/marching-cubes/algorithm';
import { ElementIterator, getElementLoci, eachElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { CommonSurfaceParams } from './util/common';
import { Sphere3D } from '../../../mol-math/geometry';

export const MolecularSurfaceMeshParams = {
    ...UnitsMeshParams,
    ...MolecularSurfaceCalculationParams,
    ...CommonSurfaceParams
};
export type MolecularSurfaceMeshParams = typeof MolecularSurfaceMeshParams

//

async function createMolecularSurfaceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: MolecularSurfaceProps, mesh?: Mesh): Promise<Mesh> {
    const { transform, field, idField } = await computeUnitMolecularSurface(structure, unit, props).runInContext(ctx.runtime);
    const params = {
        isoLevel: props.probeRadius,
        scalarField: field,
        idField
    };
    const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);

    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) Mesh.uniformTriangleGroup(surface);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, props.probeRadius);
    surface.setBoundingSphere(sphere);

    return surface;
}

export function MolecularSurfaceMeshVisual(materialId: number): UnitsVisual<MolecularSurfaceMeshParams> {
    return UnitsMeshVisual<MolecularSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(MolecularSurfaceMeshParams),
        createGeometry: createMolecularSurfaceMesh,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<MolecularSurfaceMeshParams>, currentProps: PD.Values<MolecularSurfaceMeshParams>) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true;
            if (newProps.probeRadius !== currentProps.probeRadius) state.createGeometry = true;
            if (newProps.probePositions !== currentProps.probePositions) state.createGeometry = true;
            if (newProps.ignoreHydrogens !== currentProps.ignoreHydrogens) state.createGeometry = true;
            if (newProps.includeParent !== currentProps.includeParent) state.createGeometry = true;
        }
    }, materialId);
}