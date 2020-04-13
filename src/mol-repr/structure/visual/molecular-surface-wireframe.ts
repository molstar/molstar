/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsVisual, UnitsLinesVisual, UnitsLinesParams } from '../units-visual';
import { MolecularSurfaceCalculationParams } from '../../../mol-math/geometry/molecular-surface';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { computeUnitMolecularSurface, MolecularSurfaceProps } from './util/molecular-surface';
import { computeMarchingCubesLines } from '../../../mol-geo/util/marching-cubes/algorithm';
import { ElementIterator, getElementLoci, eachElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { CommonSurfaceParams } from './util/common';
import { Sphere3D } from '../../../mol-math/geometry';

export const MolecularSurfaceWireframeParams = {
    ...UnitsLinesParams,
    ...MolecularSurfaceCalculationParams,
    ...CommonSurfaceParams,
    sizeFactor: PD.Numeric(1.5, { min: 0, max: 10, step: 0.1 }),
};
export type MolecularSurfaceWireframeParams = typeof MolecularSurfaceWireframeParams

//

async function createMolecularSurfaceWireframe(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: MolecularSurfaceProps, lines?: Lines): Promise<Lines> {
    const { transform, field, idField } = await computeUnitMolecularSurface(structure, unit, props).runInContext(ctx.runtime);
    const params = {
        isoLevel: props.probeRadius,
        scalarField: field,
        idField
    };
    const wireframe = await computeMarchingCubesLines(params, lines).runAsChild(ctx.runtime);

    Lines.transform(wireframe, transform);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, props.probeRadius);
    wireframe.setBoundingSphere(sphere);

    return wireframe;
}

export function MolecularSurfaceWireframeVisual(materialId: number): UnitsVisual<MolecularSurfaceWireframeParams> {
    return UnitsLinesVisual<MolecularSurfaceWireframeParams>({
        defaultProps: PD.getDefaultValues(MolecularSurfaceWireframeParams),
        createGeometry: createMolecularSurfaceWireframe,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<MolecularSurfaceWireframeParams>, currentProps: PD.Values<MolecularSurfaceWireframeParams>) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true;
            if (newProps.probeRadius !== currentProps.probeRadius) state.createGeometry = true;
            if (newProps.probePositions !== currentProps.probePositions) state.createGeometry = true;
            if (newProps.ignoreHydrogens !== currentProps.ignoreHydrogens) state.createGeometry = true;
            if (newProps.includeParent !== currentProps.includeParent) state.createGeometry = true;
        }
    }, materialId);
}