/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { computeStructureGaussianDensity, computeUnitGaussianDensity, GaussianDensityParams, GaussianDensityProps } from './util/gaussian';
import { computeMarchingCubesLines } from '../../../mol-geo/util/marching-cubes/algorithm';
import { UnitsLinesParams, UnitsVisual, UnitsLinesVisual } from '../units-visual';
import { ElementIterator, getElementLoci, eachElement, getSerialElementLoci, eachSerialElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { Sphere3D } from '../../../mol-math/geometry';
import { ComplexLinesParams, ComplexLinesVisual, ComplexVisual } from '../complex-visual';

const SharedParams = {
    ...GaussianDensityParams,
    sizeFactor: PD.Numeric(3, { min: 0, max: 10, step: 0.1 }),
    lineSizeAttenuation: PD.Boolean(false),
};
type SharedParams = typeof SharedParams

export const GaussianWireframeParams = {
    ...UnitsLinesParams,
    ...SharedParams,
};
export type GaussianWireframeParams = typeof GaussianWireframeParams

export const StructureGaussianWireframeParams = {
    ...ComplexLinesParams,
    ...SharedParams,
};
export type StructureGaussianWireframeParams = typeof StructureGaussianWireframeParams

async function createGaussianWireframe(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityProps, lines?: Lines): Promise<Lines> {
    const { smoothness } = props;
    const { transform, field, idField, maxRadius } = await computeUnitGaussianDensity(structure, unit, theme.size, props).runInContext(ctx.runtime);

    const params = {
        isoLevel: Math.exp(-smoothness),
        scalarField: field,
        idField
    };
    const wireframe = await computeMarchingCubesLines(params, lines).runAsChild(ctx.runtime);

    Lines.transform(wireframe, transform);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, maxRadius);
    wireframe.setBoundingSphere(sphere);

    return wireframe;
}


export function GaussianWireframeVisual(materialId: number): UnitsVisual<GaussianWireframeParams> {
    return UnitsLinesVisual<GaussianWireframeParams>({
        defaultProps: PD.getDefaultValues(GaussianWireframeParams),
        createGeometry: createGaussianWireframe,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<GaussianWireframeParams>, currentProps: PD.Values<GaussianWireframeParams>) => {
            state.createGeometry = (
                newProps.resolution !== currentProps.resolution ||
                newProps.radiusOffset !== currentProps.radiusOffset ||
                newProps.smoothness !== currentProps.smoothness ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.includeParent !== currentProps.includeParent ||
                newProps.floodfill !== currentProps.floodfill
            );
        }
    }, materialId);
}

//

async function createStructureGaussianWireframe(ctx: VisualContext, structure: Structure, theme: Theme, props: GaussianDensityProps, lines?: Lines): Promise<Lines> {
    const { smoothness } = props;
    const { transform, field, idField, maxRadius } = await computeStructureGaussianDensity(structure, theme.size, props).runInContext(ctx.runtime);

    const params = {
        isoLevel: Math.exp(-smoothness),
        scalarField: field,
        idField
    };
    const wireframe = await computeMarchingCubesLines(params, lines).runAsChild(ctx.runtime);

    Lines.transform(wireframe, transform);

    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, maxRadius);
    wireframe.setBoundingSphere(sphere);

    return wireframe;
}

export function StructureGaussianWireframeVisual(materialId: number): ComplexVisual<StructureGaussianWireframeParams> {
    return ComplexLinesVisual<StructureGaussianWireframeParams>({
        defaultProps: PD.getDefaultValues(StructureGaussianWireframeParams),
        createGeometry: createStructureGaussianWireframe,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<StructureGaussianWireframeParams>, currentProps: PD.Values<StructureGaussianWireframeParams>) => {
            state.createGeometry = (
                newProps.resolution !== currentProps.resolution ||
                newProps.radiusOffset !== currentProps.radiusOffset ||
                newProps.smoothness !== currentProps.smoothness ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.includeParent !== currentProps.includeParent ||
                newProps.floodfill !== currentProps.floodfill
            );
        }
    }, materialId);
}
