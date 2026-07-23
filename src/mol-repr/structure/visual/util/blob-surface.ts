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
    resolution: PD.Numeric(1, { min: 0.1, max: 10, step: 0.1 }, { description: 'Base grid resolution/cell spacing used to polygonize the union of blobs. When `adjustResolution` is enabled this is the atomic-scale reference that is coarsened to match the (larger) blob feature size.', ...BaseGeometry.CustomQualityParamInfo }),
    adjustResolution: PD.Boolean(true, { description: 'Automatically coarsen `resolution` to match the blob feature size implied by `blobSize` (and, secondarily, `blobShape` and `blobMethod`). Blobs are much larger and smoother than individual atoms, so the atomic-scale resolution would otherwise produce a needlessly fine, expensive grid. Disable to use `resolution` exactly as chosen (can get expensive with large blobs).' }),
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

/**
 * Reference feature radius (~vdW/atomic radius) that the incoming quality `resolution` is
 * calibrated to resolve (see `getQualityProps`).
 */
const BlobReferenceRadius = 1.6;
/**
 * How much a blob's characteristic radius grows per unit `blobSize`. Atoms binned into cells of
 * side `blobSize` produce blobs whose principal semi-axis is ~`0.5 * blobSize` (a uniformly
 * filled cell has per-axis std-dev `blobSize / sqrt(12)`, and `PrincipalAxes` reports
 * `sqrt(3) * std-dev` as the semi-axis, i.e. `blobSize / 2`).
 */
const BlobRadiusPerSize = 0.5;
/**
 * Empirical fudge factor biasing the adapted resolution finer. The feature-size model keys off a
 * blob's broad radius, which over-estimates how coarse a grid can be and still look smooth
 * (visible surface detail lives at a finer scale than a blob's overall radius). Lower means
 * finer/more detail (and more expensive); `1` disables the bias (pure feature-size model).
 */
const BlobResolutionScale = 0.5;
/** Hard ceiling on the adjusted resolution (also re-guarded by box size in `ensureReasonableResolution`). */
const BlobMaxResolution = 20;

function getBlobDensityData(position: PositionData, boundary: Boundary, radius: (index: number) => number, props: BlobDensityProps): BlobDensityData {
    const { blobSize, blobMethod, blobShape, radiusOffset, smoothness, adjustResolution } = props;
    const p = ensureReasonableResolution(boundary.box, props);

    // The incoming `p.resolution` is calibrated to resolve atomic-scale features (radius
    // ~`BlobReferenceRadius`). Blobs are much larger and smoother, so `adjustResolution` coarsens
    // the grid to keep a roughly constant number of cells per blob feature instead of a needlessly
    // fine (and expensive) atomic-scale grid.
    let resolution = p.resolution;
    if (adjustResolution) {
        // Characteristic radius of a blob feature; `blobSize -> 0` recovers the atomic reference
        // radius (no coarsening), larger bins grow the feature linearly.
        const broadRadius = BlobRadiusPerSize * blobSize + BlobReferenceRadius;
        // Spherical-harmonics of degree `d` express angular detail down to ~`radius / (d/2)`
        // relative to the ellipsoid's fixed quadratic (degree-2) boundary, so higher degrees
        // resolve finer features and want a finer grid (less coarsening). Ellipsoid == degree 2.
        const degree = blobShape.name === 'sphericalHarmonics' ? Math.max(2, blobShape.params.degree) : 2;
        const featureRadius = broadRadius * (2 / degree);
        // Grid binning can split a coherent cluster across a bin boundary into smaller fragments
        // (slightly finer features); clustering merges them into rounder, slightly larger blobs.
        // Minor secondary effect.
        const methodFactor = blobMethod.name === 'clustering' ? 1 : 0.9;
        // Scale the atomic-calibrated resolution by the feature-size ratio (biased finer by
        // `BlobResolutionScale`), but never go finer than the base resolution and never coarser
        // than the hard cap.
        resolution = Math.min(Math.max(p.resolution * (featureRadius / BlobReferenceRadius) * methodFactor * BlobResolutionScale, p.resolution), BlobMaxResolution);
    }

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
