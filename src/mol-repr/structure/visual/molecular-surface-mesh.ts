/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { UnitsMeshVisual, UnitsMeshParams, UnitsTextureMeshParams } from '../units-visual';
import { StructureElementIterator, getElementLoci, eachElement } from './util/element';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from 'mol-geo/util/marching-cubes/algorithm';
import { VisualContext } from 'mol-repr/visual';
import { Theme } from 'mol-theme/theme';
import { computeUnitMolecularSurface } from './util/molecular-surface';
import { MolecularSurfaceCalculationParams, MolecularSurfaceCalculationProps } from 'mol-math/geometry/molecular-surface';

export const MolecularSurfaceMeshParams = {
    ...UnitsMeshParams,
    ...UnitsTextureMeshParams,
    ...MolecularSurfaceCalculationParams,
}
export type MolecularSurfaceMeshParams = typeof MolecularSurfaceMeshParams

//

async function createMolecularSurfaceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: MolecularSurfaceCalculationProps, mesh?: Mesh): Promise<Mesh> {

    const { transform, field, idField } = await computeUnitMolecularSurface(unit, props).runInContext(ctx.runtime)
    const params = {
        isoLevel: props.probeRadius,
        scalarField: field,
        idField
    }
    const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime)

    Mesh.transformImmediate(surface, transform)
    Mesh.computeNormalsImmediate(surface)
    Mesh.uniformTriangleGroup(surface)

    return surface
}

export function MolecularSurfaceMeshVisual(materialId: number): UnitsVisual<MolecularSurfaceMeshParams> {
    return UnitsMeshVisual<MolecularSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(MolecularSurfaceMeshParams),
        createGeometry: createMolecularSurfaceMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<MolecularSurfaceMeshParams>, currentProps: PD.Values<MolecularSurfaceMeshParams>) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true
            if (newProps.probeRadius !== currentProps.probeRadius) state.createGeometry = true
            if (newProps.probePositions !== currentProps.probePositions) state.createGeometry = true
        }
    }, materialId)
}