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
import { eachVolumeLoci } from './util';

function getBoundingBox(gridDimension: Vec3, transform: Mat4) {
    const bbox = Box3D();
    Box3D.add(bbox, gridDimension);
    Box3D.transform(bbox, bbox, transform);
    return bbox;
}

// 2d volume texture

function getVolumeTexture2dLayout(dim: Vec3, maxTextureSize: number) {
    let width = 0;
    let height = dim[1];
    let rows = 1;
    let columns = dim[0];
    if (maxTextureSize < dim[0] * dim[2]) {
        columns =  Math.floor(maxTextureSize / dim[0]);
        rows = Math.ceil(dim[2] / columns);
        width = columns * dim[0];
        height *= rows;
    } else {
        width = dim[0] * dim[2];
    }
    width += columns; // horizontal padding
    height += rows; // vertical padding
    return { width, height, columns, rows };
}

function createVolumeTexture2d(volume: Volume, maxTextureSize: number) {
    const { cells: { space, data }, stats: { max, min } } = volume.grid;
    const dim = space.dimensions as Vec3;
    const { dataOffset: o } = space;
    const { width, height } = getVolumeTexture2dLayout(dim, maxTextureSize);

    const array = new Uint8Array(width * height * 4);
    const textureImage = { array, width, height };

    const diff = max - min;
    const [ xl, yl, zl ] = dim;
    const xlp = xl + 1; // horizontal padding
    const ylp = yl + 1; // vertical padding

    const n0 = Vec3();
    const n1 = Vec3();

    let i = 0;
    for (let z = 0; z < zl; ++z) {
        for (let y = 0; y < yl; ++y) {
            for (let x = 0; x < xl; ++x) {
                const column = Math.floor(((z * xlp) % width) / xlp);
                const row = Math.floor((z * xlp) / width);
                const px = column * xlp + x;
                const index = 4 * ((row * ylp * width) + (y * width) + px);
                const offset = o(x, y, z);

                Vec3.set(n0, data[o(x - 1, y, z)], data[o(x, y - 1, z)], data[o(x, y, z - 1)]);
                Vec3.set(n1, data[o(x + 1, y, z)], data[o(x, y + 1, z)], data[o(x, y, z + 1)]);
                Vec3.normalize(n0, Vec3.sub(n0, n0, n1));
                Vec3.addScalar(n0, Vec3.scale(n0, n0, 0.5), 0.5);
                Vec3.toArray(Vec3.scale(n0, n0, 255), array, i);

                array[index + 3] = ((data[offset] - min) / diff) * 255;
                i += 4;
            }
        }
    }

    return textureImage;
}

export function createDirectVolume2d(ctx: RuntimeContext, webgl: WebGLContext, volume: Volume, directVolume?: DirectVolume) {
    const gridDimension = volume.grid.cells.space.dimensions as Vec3;
    const textureImage = createVolumeTexture2d(volume, webgl.maxTextureSize);
    // debugTexture(createImageData(textureImage.array, textureImage.width, textureImage.height), 1/3)
    const transform = Grid.getGridToCartesianTransform(volume.grid);
    const bbox = getBoundingBox(gridDimension, transform);
    const dim = Vec3.create(gridDimension[0], gridDimension[1], gridDimension[2]);
    dim[0] += 1; // horizontal padding
    dim[1] += 1; // vertical padding

    const texture = directVolume ? directVolume.gridTexture.ref.value : webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    texture.load(textureImage);

    const { unitToCartn, cellDim } = getUnitToCartn(volume.grid);
    return DirectVolume.create(bbox, gridDimension, transform, unitToCartn, cellDim, texture, volume.grid.stats, false, directVolume);
}

// 3d volume texture

function createVolumeTexture3d(volume: Volume) {
    const { cells: { space, data }, stats: { max, min } } = volume.grid;
    const [ width, height, depth ] = space.dimensions as Vec3;
    const { dataOffset: o } = space;

    const array = new Uint8Array(width * height * depth * 4);
    const textureVolume = { array, width, height, depth };
    const diff = max - min;

    const n0 = Vec3();
    const n1 = Vec3();

    let i = 0;
    for (let z = 0; z < depth; ++z) {
        for (let y = 0; y < height; ++y) {
            for (let x = 0; x < width; ++x) {
                const offset = o(x, y, z);

                Vec3.set(n0, data[o(x - 1, y, z)], data[o(x, y - 1, z)], data[o(x, y, z - 1)]);
                Vec3.set(n1, data[o(x + 1, y, z)], data[o(x, y + 1, z)], data[o(x, y, z + 1)]);
                Vec3.normalize(n0, Vec3.sub(n0, n0, n1));
                Vec3.addScalar(n0, Vec3.scale(n0, n0, 0.5), 0.5);
                Vec3.toArray(Vec3.scale(n0, n0, 255), array, i);

                array[i + 3] = ((data[offset] - min) / diff) * 255;
                i += 4;
            }
        }
    }

    return textureVolume;
}

function getUnitToCartn(grid: Grid) {
    if (grid.transform.kind === 'matrix') {
        // TODO:
        return {
            unitToCartn: Mat4.mul(Mat4(),
                Grid.getGridToCartesianTransform(grid),
                Mat4.fromScaling(Mat4(), grid.cells.space.dimensions as Vec3)),
            cellDim: Vec3.create(1, 1, 1)
        };
    }
    const box = grid.transform.fractionalBox;
    const size = Box3D.size(Vec3(), box);
    return {
        unitToCartn: Mat4.mul3(Mat4(),
            grid.transform.cell.fromFractional,
            Mat4.fromTranslation(Mat4(), box.min),
            Mat4.fromScaling(Mat4(), size)),
        cellDim: Vec3.div(Vec3(), grid.transform.cell.size, grid.cells.space.dimensions as Vec3)
    };
}

export function createDirectVolume3d(ctx: RuntimeContext, webgl: WebGLContext, volume: Volume, directVolume?: DirectVolume) {
    const gridDimension = volume.grid.cells.space.dimensions as Vec3;
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
    return Volume.Loci(volume);
}

export function getDirectVolumeLoci(pickingId: PickingId, volume: Volume, props: DirectVolumeProps, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        return Volume.Cell.Loci(volume, Interval.ofSingleton(groupId as Volume.CellIndex));
    }
    return EmptyLoci;
}

export function eachDirectVolume(loci: Loci, volume: Volume, props: DirectVolumeProps, apply: (interval: Interval) => boolean) {
    return props.renderMode.name === 'isosurface'
        ? eachVolumeLoci(loci, volume, props.renderMode.params.isoValue, apply)
        : false;
}

//

export const DirectVolumeParams = {
    ...DirectVolume.Params,
    quality: { ...DirectVolume.Params.quality, isEssential: false },
};
export type DirectVolumeParams = typeof DirectVolumeParams
export function getDirectVolumeParams(ctx: ThemeRegistryContext, volume: Volume) {
    const p = PD.clone(DirectVolumeParams);
    p.renderMode = DirectVolume.createRenderModeParam(volume);
    return p;
}
export type DirectVolumeProps = PD.Values<DirectVolumeParams>

export function DirectVolumeVisual(materialId: number): VolumeVisual<DirectVolumeParams> {
    return VolumeVisual<DirectVolume, DirectVolumeParams>({
        defaultProps: PD.getDefaultValues(DirectVolumeParams),
        createGeometry: createDirectVolume,
        createLocationIterator: (volume: Volume) => LocationIterator(volume.grid.cells.data.length, 1, () => NullLocation),
        getLoci: getDirectVolumeLoci,
        eachLocation: eachDirectVolume,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<DirectVolumeParams>, currentProps: PD.Values<DirectVolumeParams>) => {
        },
        geometryUtils: DirectVolume.Utils
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