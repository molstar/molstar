/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { computeMarchingCubesLines } from '../../../mol-geo/util/marching-cubes/algorithm';
import { UnitsLinesParams, UnitsVisual, UnitsLinesVisual } from '../units-visual';
import { ComplexLinesParams, ComplexLinesVisual, ComplexVisual } from '../complex-visual';
import { ElementIterator, getElementLoci, eachElement, getSerialElementLoci, eachSerialElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { BlobDensityParams, computeUnitBlobSurface, computeStructureBlobSurface, shouldUpdateBlobGeometry } from './util/blob-surface';
import { Sphere3D } from '../../../mol-math/geometry';

const SharedParams = {
    ...BlobDensityParams,
    sizeFactor: PD.Numeric(3, { min: 0, max: 10, step: 0.1 }),
};
type SharedParams = typeof SharedParams

export const BlobSurfaceWireframeParams = {
    ...UnitsLinesParams,
    ...SharedParams,
};
export type BlobSurfaceWireframeParams = typeof BlobSurfaceWireframeParams

export const StructureBlobSurfaceWireframeParams = {
    ...ComplexLinesParams,
    ...SharedParams,
};
export type StructureBlobSurfaceWireframeParams = typeof StructureBlobSurfaceWireframeParams

//

async function createBlobSurfaceWireframe(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<BlobSurfaceWireframeParams>, lines?: Lines): Promise<Lines> {
    const { smoothness, radiusOffset } = props;
    const { transform, field, idField, radiusFactor, maxRadius } = await computeUnitBlobSurface(structure, unit, theme.size, props).runInContext(ctx.runtime);

    const isoLevel = Math.exp(-smoothness) / radiusFactor;
    const wireframe = await computeMarchingCubesLines({ isoLevel, scalarField: field, idField }, lines).runAsChild(ctx.runtime);

    Lines.transform(wireframe, transform);

    const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, maxRadius + extraRadius);
    wireframe.setBoundingSphere(sphere);

    return wireframe;
}

export function BlobSurfaceWireframeVisual(materialId: number): UnitsVisual<BlobSurfaceWireframeParams> {
    return UnitsLinesVisual<BlobSurfaceWireframeParams>({
        defaultProps: PD.getDefaultValues(BlobSurfaceWireframeParams),
        createGeometry: createBlobSurfaceWireframe,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<BlobSurfaceWireframeParams>, currentProps: PD.Values<BlobSurfaceWireframeParams>) => {
            state.createGeometry = shouldUpdateBlobGeometry(newProps, currentProps);
        }
    }, materialId);
}

//

async function createStructureBlobSurfaceWireframe(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<StructureBlobSurfaceWireframeParams>, lines?: Lines): Promise<Lines> {
    const { smoothness, radiusOffset } = props;
    const { transform, field, idField, radiusFactor, maxRadius } = await computeStructureBlobSurface(structure, theme.size, props).runInContext(ctx.runtime);

    const isoLevel = Math.exp(-smoothness) / radiusFactor;
    const wireframe = await computeMarchingCubesLines({ isoLevel, scalarField: field, idField }, lines).runAsChild(ctx.runtime);

    Lines.transform(wireframe, transform);

    const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, maxRadius + extraRadius);
    wireframe.setBoundingSphere(sphere);

    return wireframe;
}

export function StructureBlobSurfaceWireframeVisual(materialId: number): ComplexVisual<StructureBlobSurfaceWireframeParams> {
    return ComplexLinesVisual<StructureBlobSurfaceWireframeParams>({
        defaultProps: PD.getDefaultValues(StructureBlobSurfaceWireframeParams),
        createGeometry: createStructureBlobSurfaceWireframe,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<StructureBlobSurfaceWireframeParams>, currentProps: PD.Values<StructureBlobSurfaceWireframeParams>) => {
            state.createGeometry = shouldUpdateBlobGeometry(newProps, currentProps);
        }
    }, materialId);
}
