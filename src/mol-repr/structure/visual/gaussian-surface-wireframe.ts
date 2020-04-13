/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { computeUnitGaussianDensity, GaussianDensityParams, GaussianDensityProps } from './util/gaussian';
import { computeMarchingCubesLines } from '../../../mol-geo/util/marching-cubes/algorithm';
import { UnitsLinesParams, UnitsVisual, UnitsLinesVisual } from '../units-visual';
import { ElementIterator, getElementLoci, eachElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { Sphere3D } from '../../../mol-math/geometry';

async function createGaussianWireframe(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityProps, lines?: Lines): Promise<Lines> {
    const { smoothness } = props;
    const { transform, field, idField } = await computeUnitGaussianDensity(structure, unit, props, ctx.webgl).runInContext(ctx.runtime);

    const params = {
        isoLevel: Math.exp(-smoothness),
        scalarField: field,
        idField
    };
    const wireframe = await computeMarchingCubesLines(params, lines).runAsChild(ctx.runtime);

    Lines.transform(wireframe, transform);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, props.radiusOffset);
    wireframe.setBoundingSphere(sphere);

    return wireframe;
}

export const GaussianWireframeParams = {
    ...UnitsLinesParams,
    ...GaussianDensityParams,
    sizeFactor: PD.Numeric(1.5, { min: 0, max: 10, step: 0.1 }),
    lineSizeAttenuation: PD.Boolean(false),
    ignoreHydrogens: PD.Boolean(false),
};
export type GaussianWireframeParams = typeof GaussianWireframeParams

export function GaussianWireframeVisual(materialId: number): UnitsVisual<GaussianWireframeParams> {
    return UnitsLinesVisual<GaussianWireframeParams>({
        defaultProps: PD.getDefaultValues(GaussianWireframeParams),
        createGeometry: createGaussianWireframe,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<GaussianWireframeParams>, currentProps: PD.Values<GaussianWireframeParams>) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true;
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true;
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true;
            if (newProps.useGpu !== currentProps.useGpu) state.createGeometry = true;
            if (newProps.ignoreHydrogens !== currentProps.ignoreHydrogens) state.createGeometry = true;
            if (newProps.traceOnly !== currentProps.traceOnly) state.createGeometry = true;
            if (newProps.includeParent !== currentProps.includeParent) state.createGeometry = true;
        }
    }, materialId);
}