/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../units-visual';
import { ComplexVisual, ComplexMeshParams, ComplexMeshVisual } from '../complex-visual';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from '../../../mol-geo/util/marching-cubes/algorithm';
import { ElementIterator, getElementLoci, eachElement, getSerialElementLoci, eachSerialElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { BlobDensityParams, BlobDensityProps, computeUnitBlobSurface, computeStructureBlobSurface } from './util/blob-surface';
import { Sphere3D } from '../../../mol-math/geometry';
import { ValueCell } from '../../../mol-util/value-cell';

export const BlobSurfaceMeshParams = {
    ...UnitsMeshParams,
    ...BlobDensityParams,
};
export type BlobSurfaceMeshParams = typeof BlobSurfaceMeshParams

export const StructureBlobSurfaceMeshParams = {
    ...ComplexMeshParams,
    ...BlobDensityParams,
};
export type StructureBlobSurfaceMeshParams = typeof StructureBlobSurfaceMeshParams

function shouldUpdateGeometry(newProps: BlobDensityProps, currentProps: BlobDensityProps) {
    return (
        newProps.blobSize !== currentProps.blobSize ||
        newProps.blobMethod.name !== currentProps.blobMethod.name ||
        (newProps.blobMethod.name === 'clustering' && currentProps.blobMethod.name === 'clustering' &&
            newProps.blobMethod.params.iterations !== currentProps.blobMethod.params.iterations) ||
        newProps.resolution !== currentProps.resolution ||
        newProps.resolutionFactor !== currentProps.resolutionFactor ||
        newProps.radiusOffset !== currentProps.radiusOffset ||
        newProps.smoothness !== currentProps.smoothness ||
        newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
        newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
        newProps.traceOnly !== currentProps.traceOnly ||
        newProps.includeParent !== currentProps.includeParent
    );
}

//

async function createBlobSurfaceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<BlobSurfaceMeshParams>, mesh?: Mesh): Promise<Mesh> {
    const { smoothness, radiusOffset } = props;
    const { transform, field, idField, radiusFactor, maxRadius } = await computeUnitBlobSurface(structure, unit, theme.size, props).runInContext(ctx.runtime);

    const isoLevel = Math.exp(-smoothness) / radiusFactor;
    const surface = await computeMarchingCubesMesh({ isoLevel, scalarField: field, idField }, mesh).runAsChild(ctx.runtime);

    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) {
        Mesh.uniformTriangleGroup(surface);
        ValueCell.updateIfChanged(surface.varyingGroup, false);
    } else {
        ValueCell.updateIfChanged(surface.varyingGroup, true);
    }

    const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, maxRadius + extraRadius);
    surface.setBoundingSphere(sphere);

    return surface;
}

export function BlobSurfaceMeshVisual(materialId: number): UnitsVisual<BlobSurfaceMeshParams> {
    return UnitsMeshVisual<BlobSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(BlobSurfaceMeshParams),
        createGeometry: createBlobSurfaceMesh,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<BlobSurfaceMeshParams>, currentProps: PD.Values<BlobSurfaceMeshParams>) => {
            state.createGeometry = shouldUpdateGeometry(newProps, currentProps);
        }
    }, materialId);
}

//

async function createStructureBlobSurfaceMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<StructureBlobSurfaceMeshParams>, mesh?: Mesh): Promise<Mesh> {
    const { smoothness, radiusOffset } = props;
    const { transform, field, idField, radiusFactor, maxRadius } = await computeStructureBlobSurface(structure, theme.size, props).runInContext(ctx.runtime);

    const isoLevel = Math.exp(-smoothness) / radiusFactor;
    const surface = await computeMarchingCubesMesh({ isoLevel, scalarField: field, idField }, mesh).runAsChild(ctx.runtime);

    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) {
        Mesh.uniformTriangleGroup(surface);
        ValueCell.updateIfChanged(surface.varyingGroup, false);
    } else {
        ValueCell.updateIfChanged(surface.varyingGroup, true);
    }

    const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, maxRadius + extraRadius);
    surface.setBoundingSphere(sphere);

    return surface;
}

export function StructureBlobSurfaceMeshVisual(materialId: number): ComplexVisual<StructureBlobSurfaceMeshParams> {
    return ComplexMeshVisual<StructureBlobSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(StructureBlobSurfaceMeshParams),
        createGeometry: createStructureBlobSurfaceMesh,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<StructureBlobSurfaceMeshParams>, currentProps: PD.Values<StructureBlobSurfaceMeshParams>) => {
            state.createGeometry = shouldUpdateGeometry(newProps, currentProps);
        }
    }, materialId);
}
