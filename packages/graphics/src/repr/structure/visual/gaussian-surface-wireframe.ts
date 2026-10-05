/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { VisualContext } from '../../visual.js';
import type { Unit, Structure } from '@molstar/model/model/structure';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Lines } from '@molstar/graphics/geo/geometry/lines/lines';
import { computeStructureGaussianDensity, computeUnitGaussianDensity, GaussianDensityParams, type GaussianDensityProps } from './util/gaussian.js';
import { computeMarchingCubesLines } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import { UnitsLinesParams, type UnitsVisual, UnitsLinesVisual } from '../units-visual.js';
import { ElementIterator, getElementLoci, eachElement, getSerialElementLoci, eachSerialElement } from './util/element.js';
import type { VisualUpdateState } from '../../util.js';
import { Sphere3D } from '@molstar/core/math/geometry';
import { ComplexLinesParams, ComplexLinesVisual, type ComplexVisual } from '../complex-visual.js';
import { Tensor } from '@molstar/core/math/linear-algebra/tensor';

const SharedParams = {
    ...GaussianDensityParams,
    sizeFactor: PD.Numeric(3, { min: 0, max: 10, step: 0.1 }),
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
    const { smoothness, floodfill, radiusOffset } = props;
    const { transform, field, idField, maxRadius, radiusFactor } = await computeUnitGaussianDensity(structure, unit, theme.size, props).runInContext(ctx.runtime);

    const isoLevel = Math.exp(-smoothness) / radiusFactor;
    const params = {
        isoLevel,
        scalarField: floodfill !== 'off' ? Tensor.createFloodfilled(field, isoLevel, floodfill) : field,
        idField
    };
    const wireframe = await computeMarchingCubesLines(params, lines).runAsChild(ctx.runtime);

    Lines.transform(wireframe, transform);

    const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, maxRadius + extraRadius);
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
    const { smoothness, floodfill, radiusOffset } = props;
    const { transform, field, idField, maxRadius, radiusFactor } = await computeStructureGaussianDensity(structure, theme.size, props).runInContext(ctx.runtime);

    const isoLevel = Math.exp(-smoothness) / radiusFactor;
    const params = {
        isoLevel,
        scalarField: floodfill !== 'off' ? Tensor.createFloodfilled(field, isoLevel, floodfill) : field,
        idField
    };
    const wireframe = await computeMarchingCubesLines(params, lines).runAsChild(ctx.runtime);

    Lines.transform(wireframe, transform);

    const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, maxRadius + extraRadius);
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
