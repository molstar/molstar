/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Vec3, Mat4 } from '../../mol-math/linear-algebra';
import { Box3D } from '../../mol-math/geometry';
import { Grid, Volume } from '../../mol-model/volume';
import { RuntimeContext } from '../../mol-task';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { DirectVolume } from '../../mol-geo/geometry/direct-volume/direct-volume';
import { VisualContext } from '../visual';
import { Theme, ThemeRegistryContext } from '../../mol-theme/theme';
import { VolumeVisual, VolumeRepresentation, VolumeRepresentationProvider } from './representation';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { NullLocation } from '../../mol-model/location';
import { VisualUpdateState } from '../util';
import { RepresentationContext, RepresentationParamsGetter } from '../representation';
import { Interval } from '../../mol-data/int';
import { Loci, EmptyLoci } from '../../mol-model/loci';
import { PickingId } from '../../mol-geo/geometry/picking';
import { createVolumeTexture2d, createVolumeTexture3d, eachVolumeLoci, getVolumeTexture2dLayout } from './util';

function getBoundingBox(gridDimension: Vec3, transform: Mat4) {
    const bbox = Box3D();
    Box3D.add(bbox, gridDimension);
    Box3D.transform(bbox, bbox, transform);
    return bbox;
}

// 2d volume texture

export function createDirectVolume2d(ctx: RuntimeContext, webgl: WebGLContext, volume: Volume, directVolume?: DirectVolume) {
    const gridDimension = volume.grid.cells.space.dimensions as Vec3;
    const { width, height } = getVolumeTexture2dLayout(gridDimension);
    if (Math.max(width, height) > webgl.maxTextureSize / 2) {
        throw new Error('volume too large for direct-volume rendering');
    }

    const textureImage = createVolumeTexture2d(volume, 'normals');
    // debugTexture(createImageData(textureImage.array, textureImage.width, textureImage.height), 1/3)
    const transform = Grid.getGridToCartesianTransform(volume.grid);
    const bbox = getBoundingBox(gridDimension, transform);

    const texture = directVolume ? directVolume.gridTexture.ref.value : webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    texture.load(textureImage);

    const { unitToCartn, cellDim } = getUnitToCartn(volume.grid);
    return DirectVolume.create(bbox, gridDimension, transform, unitToCartn, cellDim, texture, volume.grid.stats, false, directVolume);
}

// 3d volume texture

function getUnitToCartn(grid: Grid) {
    if (grid.transform.kind === 'matrix') {
        return {
            unitToCartn: Mat4.mul(Mat4(),
                grid.transform.matrix,
                Mat4.fromScaling(Mat4(), grid.cells.space.dimensions as Vec3)
            ),
            cellDim: Mat4.getScaling(Vec3(), grid.transform.matrix)
        };
    }
    const box = grid.transform.fractionalBox;
    const size = Box3D.size(Vec3(), box);
    return {
        unitToCartn: Mat4.mul3(Mat4(),
            grid.transform.cell.fromFractional,
            Mat4.fromTranslation(Mat4(), box.min),
            Mat4.fromScaling(Mat4(), size)
        ),
        cellDim: Vec3.div(Vec3(), grid.transform.cell.size, grid.cells.space.dimensions as Vec3)
    };
}

export function createDirectVolume3d(ctx: RuntimeContext, webgl: WebGLContext, volume: Volume, directVolume?: DirectVolume) {
    const gridDimension = volume.grid.cells.space.dimensions as Vec3;
    if (Math.max(...gridDimension) > webgl.max3dTextureSize / 2) {
        throw new Error('volume too large for direct-volume rendering');
    }

    const textureVolume = createVolumeTexture3d(volume);
    const transform = Grid.getGridToCartesianTransform(volume.grid);
    const bbox = getBoundingBox(gridDimension, transform);

    const texture = directVolume ? directVolume.gridTexture.ref.value : webgl.resources.texture('volume-uint8', 'rgba', 'ubyte', 'linear');
    texture.load(textureVolume);

    const { unitToCartn, cellDim } = getUnitToCartn(volume.grid);
    return DirectVolume.create(bbox, gridDimension, transform, unitToCartn, cellDim, texture, volume.grid.stats, false, directVolume);
}

//

export async function createDirectVolume(ctx: VisualContext, volume: Volume, theme: Theme, props: PD.Values<DirectVolumeParams>, directVolume?: DirectVolume) {
    const { runtime, webgl } = ctx;
    if (webgl === undefined) throw new Error('DirectVolumeVisual requires `webgl` in props');

    return webgl.isWebGL2 ?
        createDirectVolume3d(runtime, webgl, volume, directVolume) :
        createDirectVolume2d(runtime, webgl, volume, directVolume);
}

function getLoci(volume: Volume, props: PD.Values<DirectVolumeParams>) {
    return props.renderMode.name === 'isosurface'
        ? Volume.Isosurface.Loci(volume, props.renderMode.params.isoValue)
        : Volume.Loci(volume);
}

export function getDirectVolumeLoci(pickingId: PickingId, volume: Volume, props: DirectVolumeProps, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        return Volume.Cell.Loci(volume, Interval.ofSingleton(groupId as Volume.CellIndex));
    }
    return EmptyLoci;
}

export function eachDirectVolume(loci: Loci, volume: Volume, props: DirectVolumeProps, apply: (interval: Interval) => boolean) {
    const isoValue = props.renderMode.name === 'isosurface'
        ? props.renderMode.params.isoValue : undefined;
    return eachVolumeLoci(loci, volume, isoValue, apply);
}

//

export const DirectVolumeParams = {
    ...DirectVolume.Params,
    quality: { ...DirectVolume.Params.quality, isEssential: false },
};
export type DirectVolumeParams = typeof DirectVolumeParams
export function getDirectVolumeParams(ctx: ThemeRegistryContext, volume: Volume) {
    const p = PD.clone(DirectVolumeParams);
    p.renderMode = DirectVolume.createRenderModeParam(volume.grid.stats);
    return p;
}
export type DirectVolumeProps = PD.Values<DirectVolumeParams>

export function DirectVolumeVisual(materialId: number): VolumeVisual<DirectVolumeParams> {
    return VolumeVisual<DirectVolume, DirectVolumeParams>({
        defaultProps: PD.getDefaultValues(DirectVolumeParams),
        createGeometry: createDirectVolume,
        createLocationIterator: (volume: Volume) => LocationIterator(volume.grid.cells.data.length, 1, 1, () => NullLocation),
        getLoci: getDirectVolumeLoci,
        eachLocation: eachDirectVolume,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<DirectVolumeParams>, currentProps: PD.Values<DirectVolumeParams>) => {
        },
        geometryUtils: DirectVolume.Utils,
        dispose: (geometry: DirectVolume) => {
            geometry.gridTexture.ref.value.destroy();
        }
    }, materialId);
}

export function DirectVolumeRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, DirectVolumeParams>): VolumeRepresentation<DirectVolumeParams> {
    return VolumeRepresentation('Direct Volume', ctx, getParams, DirectVolumeVisual, getLoci);
}

export const DirectVolumeRepresentationProvider = VolumeRepresentationProvider({
    name: 'direct-volume',
    label: 'Direct Volume',
    description: 'Direct rendering of volumetric data.',
    factory: DirectVolumeRepresentation,
    getParams: getDirectVolumeParams,
    defaultValues: PD.getDefaultValues(DirectVolumeParams),
    defaultColorTheme: { name: 'uniform' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (volume: Volume) => !Volume.isEmpty(volume)
});