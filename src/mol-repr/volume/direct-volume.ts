/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { VisualUpdateState } from '../util';
import { RepresentationContext, RepresentationParamsGetter } from '../representation';
import { Interval, OrderedSet } from '../../mol-data/int';
import { Loci, EmptyLoci } from '../../mol-model/loci';
import { PickingId } from '../../mol-geo/geometry/picking';
import { createVolumeCellLocationIterator, createVolumeTexture2d, createVolumeTexture3d, eachVolumeLoci, getVolumeTexture2dLayout } from './util';
import { Texture } from '../../mol-gl/webgl/texture';

function getBoundingBox(gridDimension: Vec3, transform: Mat4) {
    const bbox = Box3D();
    Box3D.add(bbox, gridDimension);
    Box3D.transform(bbox, bbox, transform);
    return bbox;
}

// 2d volume texture

export function createDirectVolume2d(ctx: RuntimeContext, webgl: WebGLContext, volume: Volume, props: PD.Values<DirectVolumeParams>, directVolume?: DirectVolume) {
    const gridDimension = volume.grid.cells.space.dimensions as Vec3;
    const { width, height } = getVolumeTexture2dLayout(gridDimension);
    if (Math.max(width, height) > webgl.maxTextureSize / 2) {
        throw new Error('volume too large for direct-volume rendering');
    }

    const dataType = props.dataType === 'halfFloat' && !webgl.extensions.textureHalfFloat ? 'float' : props.dataType;

    const textureImage = createVolumeTexture2d(volume, 'normals', 0, dataType);
    // debugTexture(createImageData(textureImage.array, textureImage.width, textureImage.height), 1/3)
    const transform = Grid.getGridToCartesianTransform(volume.grid);
    const bbox = getBoundingBox(gridDimension, transform);

    let texture: Texture;
    if (directVolume && directVolume.dataType.ref.value === dataType) {
        texture = directVolume.gridTexture.ref.value;
    } else {
        texture = dataType === 'byte'
            ? webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'linear')
            : dataType === 'halfFloat'
                ? webgl.resources.texture('image-float16', 'rgba', 'fp16', 'linear')
                : webgl.resources.texture('image-float32', 'rgba', 'float', 'linear');
    }
    texture.load(textureImage);

    const { unitToCartn, cellDim } = getUnitToCartn(volume.grid);
    const axisOrder = volume.grid.cells.space.axisOrderSlowToFast as Vec3;
    return DirectVolume.create(bbox, gridDimension, transform, unitToCartn, cellDim, texture, volume.grid.stats, false, axisOrder, dataType, directVolume);
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

export function createDirectVolume3d(ctx: RuntimeContext, webgl: WebGLContext, volume: Volume, props: PD.Values<DirectVolumeParams>, directVolume?: DirectVolume) {
    const gridDimension = volume.grid.cells.space.dimensions as Vec3;
    if (Math.max(...gridDimension) > webgl.max3dTextureSize / 2) {
        throw new Error('volume too large for direct-volume rendering');
    }

    const dataType = props.dataType === 'halfFloat' && !webgl.extensions.textureHalfFloat ? 'float' : props.dataType;

    const textureVolume = createVolumeTexture3d(volume, dataType);
    const transform = Grid.getGridToCartesianTransform(volume.grid);
    const bbox = getBoundingBox(gridDimension, transform);

    let texture: Texture;
    if (directVolume && directVolume.dataType.ref.value === dataType) {
        texture = directVolume.gridTexture.ref.value;
    } else {
        texture = dataType === 'byte'
            ? webgl.resources.texture('volume-uint8', 'rgba', 'ubyte', 'linear')
            : dataType === 'halfFloat'
                ? webgl.resources.texture('volume-float16', 'rgba', 'fp16', 'linear')
                : webgl.resources.texture('volume-float32', 'rgba', 'float', 'linear');
    }
    texture.load(textureVolume);

    const { unitToCartn, cellDim } = getUnitToCartn(volume.grid);
    const axisOrder = volume.grid.cells.space.axisOrderSlowToFast as Vec3;
    return DirectVolume.create(bbox, gridDimension, transform, unitToCartn, cellDim, texture, volume.grid.stats, false, axisOrder, dataType, directVolume);
}

//

export async function createDirectVolume(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: PD.Values<DirectVolumeParams>, directVolume?: DirectVolume) {
    const { runtime, webgl } = ctx;
    if (webgl === undefined) throw new Error('DirectVolumeVisual requires `webgl` in VisualContext');

    return webgl.isWebGL2 ?
        createDirectVolume3d(runtime, webgl, volume, props, directVolume) :
        createDirectVolume2d(runtime, webgl, volume, props, directVolume);
}

function getLoci(volume: Volume, props: PD.Values<DirectVolumeParams>) {
    const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
    return Volume.Loci(volume, instances);
}

export function getDirectVolumeLoci(pickingId: PickingId, volume: Volume, key: number, props: DirectVolumeProps, id: number) {
    const { objectId, groupId, instanceId } = pickingId;
    if (id === objectId) {
        const instances = OrderedSet.ofSingleton(instanceId as Volume.InstanceIndex);
        const indices = Interval.ofSingleton(groupId as Volume.CellIndex);
        return Volume.Cell.Loci(volume, [{ indices, instances }]);
    }
    return EmptyLoci;
}

export function eachDirectVolume(loci: Loci, volume: Volume, key: number, props: DirectVolumeProps, apply: (interval: Interval) => boolean) {
    return eachVolumeLoci(loci, volume, undefined, apply);
}

//

export const DirectVolumeParams = {
    ...DirectVolume.Params,
    quality: { ...DirectVolume.Params.quality, isEssential: false },
    dataType: PD.Select('byte', PD.arrayToOptions(['byte', 'float', 'halfFloat'] as const)),
};
export type DirectVolumeParams = typeof DirectVolumeParams
export function getDirectVolumeParams(ctx: ThemeRegistryContext, volume: Volume) {
    const params = PD.clone(DirectVolumeParams);
    params.controlPoints.getVolume = () => volume;
    return params;
}
export type DirectVolumeProps = PD.Values<DirectVolumeParams>

export function DirectVolumeVisual(materialId: number): VolumeVisual<DirectVolumeParams> {
    return VolumeVisual<DirectVolume, DirectVolumeParams>({
        defaultProps: PD.getDefaultValues(DirectVolumeParams),
        createGeometry: createDirectVolume,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getDirectVolumeLoci,
        eachLocation: eachDirectVolume,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<DirectVolumeParams>, currentProps: PD.Values<DirectVolumeParams>) => {
            state.createGeometry = newProps.dataType !== currentProps.dataType;
        },
        geometryUtils: DirectVolume.Utils,
        dispose: (geometry: DirectVolume) => {
            geometry.gridTexture.ref.value.destroy();
        },
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
    defaultColorTheme: { name: 'volume-value' },
    defaultSizeTheme: { name: 'uniform' },
    locationKinds: ['position-location', 'direct-location'],
    isApplicable: (volume: Volume) => !Volume.isEmpty(volume) && !Volume.Segmentation.get(volume)
});