/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from '../../../../mol-model/structure';
import { Task } from '../../../../mol-task';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { getUnitConformationAndRadius, getStructureConformationAndRadius, CommonSurfaceParams, ensureReasonableResolution } from './common';
import { computeBlobSurface, BlobSurfaceData } from '../../../../mol-math/geometry/blob-surface';
import { BaseGeometry } from '../../../../mol-geo/geometry/base';
import { SizeTheme } from '../../../../mol-theme/size';
import { PositionData } from '../../../../mol-math/geometry';
import { Boundary } from '../../../../mol-math/geometry/boundary';

export const BlobDensityParams = {
    blobSize: PD.Numeric(30, { min: 4, max: 200, step: 1 }, { description: 'Size of the spatial bins used to coarsen atoms into a small number of "blobs" before surfacing. Higher means fewer, bigger, blobbier blobs.' }),
    blobMethod: PD.MappedStatic('clustering', {
        grid: PD.Group({}),
        clustering: PD.Group({
            iterations: PD.Numeric(2, { min: 1, max: 20, step: 1 }, { description: 'Number of k-means refinement passes. More iterations better separate atoms into spatially coherent clusters (softens the grid-binning boundary-splitting artifact) at extra cost.' })
        })
    }, { description: 'How atoms are grouped into blobs before fitting a shape to each group. "grid" bins atoms on a fixed spatial grid (fast, single pass). "clustering" refines those grid bins further via k-means, which can better separate atoms that are close in space but land on opposite sides of a bin boundary.' }),
    resolution: PD.Numeric(1, { min: 0.1, max: 10, step: 0.1 }, { description: 'Grid resolution/cell spacing used to polygonize the union of blobs.', ...BaseGeometry.CustomQualityParamInfo }),
    resolutionFactor: PD.Numeric(16, { min: 0, max: 32, step: 1 }, { description: 'Scale the effective grid resolution proportionally to `blobSize`. Set to 0 to use `resolution` exactly as chosen, regardless of `blobSize` (can get expensive with large blobs).' }),
    blobShape: PD.MappedStatic('ellipsoid', {
        ellipsoid: PD.Group({}),
        sphericalHarmonics: PD.Group({
            degree: PD.Numeric(2, { min: 0, max: 12, step: 1 }, { description: 'Max degree of the fitted spherical-harmonics radial boundary. Higher degrees can hug non-ellipsoidal, even mildly concave, blob shapes more closely, at extra cost.' }),
            regularization: PD.Numeric(0.05, { min: 0, max: 1, step: 0.01 }, { description: 'Tikhonov regularization strength for the spherical-harmonics fit. Higher values relax degenerate/sparse groups (e.g. single-atom blobs, high degrees) toward a smooth mean-radius sphere instead of oscillating/overfitting; 0 disables it.' })
        })
    }, { description: 'Shape fitted to each group of atoms. "ellipsoid" fits a fixed quadratic boundary (fast). "sphericalHarmonics" fits an angularly-varying radial boundary that can better follow non-ellipsoidal blob shapes.' }),
    radiusOffset: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, { description: 'Extra/offset radius added to the atoms/coarse elements for blob calculation. Useful to create coarse, low resolution surfaces.' }),
    smoothness: PD.Numeric(1.5, { min: 1, max: 3, step: 0.1 }, { description: 'Smoothness of the blob surface, lower is smoother.' }),
    ...CommonSurfaceParams
};
export const DefaultBlobDensityProps = PD.getDefaultValues(BlobDensityParams);
export type BlobDensityProps = typeof DefaultBlobDensityProps

export type BlobDensityData = BlobSurfaceData

function getBlobDensityData(position: PositionData, boundary: Boundary, radius: (index: number) => number, props: BlobDensityProps): BlobDensityData {
    const { blobSize, blobMethod, blobShape, radiusOffset, smoothness, resolutionFactor } = props;
    const p = ensureReasonableResolution(boundary.box, props);
    // Scale the user's chosen resolution proportionally to `blobSize`
    const resolution = resolutionFactor > 0 && blobSize > resolutionFactor
        ? Math.min(p.resolution * (blobSize / resolutionFactor), 20)
        : p.resolution;
    return computeBlobSurface(position, boundary, radius, {
        blobSize,
        method: blobMethod.name,
        clusterIterations: blobMethod.name === 'clustering' ? blobMethod.params.iterations : 0,
        shape: blobShape.name === 'sphericalHarmonics' ? 'sh' : 'ellipsoid',
        shDegree: blobShape.name === 'sphericalHarmonics' ? blobShape.params.degree : 0,
        shRegularization: blobShape.name === 'sphericalHarmonics' ? blobShape.params.regularization : 0,
        resolution,
        radiusOffset,
        smoothness
    });
}

export function computeUnitBlobSurface(structure: Structure, unit: Unit, sizeTheme: SizeTheme<any>, props: BlobDensityProps) {
    const { position, boundary, radius } = getUnitConformationAndRadius(structure, unit, sizeTheme, props);
    return Task.create('Blob Surface', async () => {
        return getBlobDensityData(position, boundary, radius, props);
    });
}

export function computeStructureBlobSurface(structure: Structure, sizeTheme: SizeTheme<any>, props: BlobDensityProps) {
    const { position, boundary, radius } = getStructureConformationAndRadius(structure, sizeTheme, props);
    return Task.create('Blob Surface', async () => {
        return getBlobDensityData(position, boundary, radius, props);
    });
}
