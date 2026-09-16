/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mesh } from '../../../../mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from '../../../../mol-geo/util/marching-cubes/algorithm';
import { Sphere3D } from '../../../../mol-math/geometry';
import { Structure } from '../../../../mol-model/structure';
import { Theme } from '../../../../mol-theme/theme';
import { ValueCell } from '../../../../mol-util/value-cell';
import { VisualContext } from '../../../visual';
import { DefaultBlobDensityProps, computeStructureBlobSurface, BlobDensityProps } from '../../../structure/visual/util/blob-surface';
import { getStructureTargetSizeTheme } from './structure';

/** The subset of representation props needed to build a blob surface target geometry. */
export interface BlobSurfaceTargetProps {
    blobSize: number
    blobMethod: BlobDensityProps['blobMethod']
    resolution: number
    adjustResolution: boolean
    blobShape: BlobDensityProps['blobShape']
    radiusOffset: number
    smoothness: number
}

/**
 * The geometry for a blob surface target is a coarse surface of the target's elements, built by
 * fitting a small number of blobs to them and polygonizing their union.
 */
export async function createBlobSurfaceTargetGeometry(ctx: VisualContext, target: Structure, theme: Theme, props: BlobSurfaceTargetProps, existing?: Mesh): Promise<Mesh> {
    const { blobSize, blobMethod, resolution, adjustResolution, blobShape, radiusOffset, smoothness } = props;

    const sizeTheme = getStructureTargetSizeTheme(theme);
    const { transform, field, idField, radiusFactor, maxRadius } = await computeStructureBlobSurface(target, sizeTheme, {
        ...DefaultBlobDensityProps,
        blobSize, blobMethod, resolution, adjustResolution, blobShape, radiusOffset, smoothness,
    }).runInContext(ctx.runtime);

    const isoLevel = Math.exp(-smoothness) / radiusFactor;
    const surface = await computeMarchingCubesMesh({ isoLevel, scalarField: field, idField }, existing).runAsChild(ctx.runtime);

    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) {
        Mesh.uniformTriangleGroup(surface);
        ValueCell.updateIfChanged(surface.varyingGroup, false);
    } else {
        ValueCell.updateIfChanged(surface.varyingGroup, true);
    }

    const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
    surface.setBoundingSphere(Sphere3D.expand(Sphere3D(), target.boundary.sphere, maxRadius + extraRadius));

    return surface;
}

/** Whether the blob surface geometry must be rebuilt because relevant props changed. */
export function blobSurfaceTargetGeometryPropsChanged(oldProps: BlobSurfaceTargetProps, newProps: BlobSurfaceTargetProps): boolean {
    return (
        newProps.blobSize !== oldProps.blobSize ||
        newProps.blobMethod.name !== oldProps.blobMethod.name ||
        (newProps.blobMethod.name === 'clustering' && oldProps.blobMethod.name === 'clustering' &&
            newProps.blobMethod.params.iterations !== oldProps.blobMethod.params.iterations) ||
        newProps.resolution !== oldProps.resolution ||
        newProps.adjustResolution !== oldProps.adjustResolution ||
        newProps.blobShape.name !== oldProps.blobShape.name ||
        (newProps.blobShape.name === 'sphericalHarmonics' && oldProps.blobShape.name === 'sphericalHarmonics' &&
            (newProps.blobShape.params.degree !== oldProps.blobShape.params.degree ||
                newProps.blobShape.params.regularization !== oldProps.blobShape.params.regularization)) ||
        newProps.radiusOffset !== oldProps.radiusOffset ||
        newProps.smoothness !== oldProps.smoothness
    );
}
