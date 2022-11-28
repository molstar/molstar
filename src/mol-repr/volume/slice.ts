/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Image } from '../../mol-geo/geometry/image/image';
import { ThemeRegistryContext, Theme } from '../../mol-theme/theme';
import { Grid, Volume } from '../../mol-model/volume';
import { VolumeVisual, VolumeRepresentation, VolumeRepresentationProvider } from './representation';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { VisualUpdateState } from '../util';
import { NullLocation } from '../../mol-model/location';
import { RepresentationContext, RepresentationParamsGetter } from '../representation';
import { VisualContext } from '../visual';
import { PickingId } from '../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../mol-model/loci';
import { Interval, SortedArray } from '../../mol-data/int';
import { transformPositionArray } from '../../mol-geo/util';
import { RenderableState } from '../../mol-gl/renderable';
import { Color } from '../../mol-util/color';
import { ColorTheme } from '../../mol-theme/color';
import { packIntToRGBArray } from '../../mol-util/number-packing';
import { eachVolumeLoci } from './util';

export async function createImage(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: PD.Values<SliceParams>, image?: Image) {
    const { dimension: { name: dim }, isoValue } = props;

    const { space, data } = volume.grid.cells;
    const { min, max } = volume.grid.stats;
    const isoVal = Volume.IsoValue.toAbsolute(isoValue, volume.grid.stats).absoluteValue;

    // TODO more color themes
    const color = 'color' in theme.color ? theme.color.color(NullLocation, false) : Color(0xffffff);
    const [r, g, b] = Color.toRgbNormalized(color);

    const {
        width, height,
        x, y, z,
        x0, y0, z0,
        nx, ny, nz
    } = getSliceInfo(volume.grid, props);

    const corners = new Float32Array(
        dim === 'x' ? [x, 0, 0, x, y, 0, x, 0, z, x, y, z] :
            dim === 'y' ? [0, y, 0, x, y, 0, 0, y, z, x, y, z] :
                [0, 0, z, 0, y, z, x, 0, z, x, y, z]
    );

    const imageArray = new Uint8Array(width * height * 4);
    const groupArray = getPackedGroupArray(volume.grid, props);

    let i = 0;
    for (let iy = y0; iy < ny; ++iy) {
        for (let ix = x0; ix < nx; ++ix) {
            for (let iz = z0; iz < nz; ++iz) {
                const val = space.get(data, ix, iy, iz);
                const normVal = (val - min) / (max - min);

                imageArray[i] = r * normVal * 2 * 255;
                imageArray[i + 1] = g * normVal * 2 * 255;
                imageArray[i + 2] = b * normVal * 2 * 255;
                imageArray[i + 3] = val >= isoVal ? 255 : 0;

                i += 4;
            }
        }
    }

    const imageTexture = { width, height, array: imageArray, flipY: true };
    const groupTexture = { width, height, array: groupArray, flipY: true };

    const transform = Grid.getGridToCartesianTransform(volume.grid);
    transformPositionArray(transform, corners, 0, 4);

    return Image.create(imageTexture, corners, groupTexture, image);
}

function getSliceInfo(grid: Grid, props: PD.Values<SliceParams>) {
    const { dimension: { name: dim, params: index } } = props;
    const { space } = grid.cells;

    let width, height;
    let x, y, z;
    let x0 = 0, y0 = 0, z0 = 0;
    let [nx, ny, nz] = space.dimensions;

    if (dim === 'x') {
        x = index, y = ny - 1, z = nz - 1;
        width = nz, height = ny;
        x0 = x, nx = x0 + 1;
    } else if (dim === 'y') {
        x = nx - 1, y = index, z = nz - 1;
        width = nz, height = nx;
        y0 = y, ny = y0 + 1;
    } else {
        x = nx - 1, y = ny - 1, z = index;
        width = nx, height = ny;
        z0 = z, nz = z0 + 1;
    }
    return {
        width, height,
        x, y, z,
        x0, y0, z0,
        nx, ny, nz
    };
}

function getPackedGroupArray(grid: Grid, props: PD.Values<SliceParams>) {
    const { space } = grid.cells;
    const { width, height, x0, y0, z0, nx, ny, nz } = getSliceInfo(grid, props);
    const groupArray = new Uint8Array(width * height * 4);

    let j = 0;
    for (let iy = y0; iy < ny; ++iy) {
        for (let ix = x0; ix < nx; ++ix) {
            for (let iz = z0; iz < nz; ++iz) {
                packIntToRGBArray(space.dataOffset(ix, iy, iz), groupArray, j);
                j += 4;
            }
        }
    }
    return groupArray;
}

function getGroupArray(grid: Grid, props: PD.Values<SliceParams>) {
    const { space } = grid.cells;
    const { width, height, x0, y0, z0, nx, ny, nz } = getSliceInfo(grid, props);
    const groupArray = new Uint32Array(width * height);

    let j = 0;
    for (let iy = y0; iy < ny; ++iy) {
        for (let ix = x0; ix < nx; ++ix) {
            for (let iz = z0; iz < nz; ++iz) {
                groupArray[j] = space.dataOffset(ix, iy, iz);
                j += 1;
            }
        }
    }
    return groupArray;
}

function getLoci(volume: Volume, props: PD.Values<SliceParams>) {
    // TODO cache somehow?
    const groupArray = getGroupArray(volume.grid, props);
    return Volume.Cell.Loci(volume, SortedArray.ofUnsortedArray(groupArray));
}

function getSliceLoci(pickingId: PickingId, volume: Volume, key: number, props: PD.Values<SliceParams>, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        const granularity = Volume.PickingGranularity.get(volume);
        if (granularity === 'volume') {
            return Volume.Loci(volume);
        } if (granularity === 'object') {
            return getLoci(volume, props);
        } else {
            return Volume.Cell.Loci(volume, Interval.ofSingleton(groupId as Volume.CellIndex));
        }
    }
    return EmptyLoci;
}

function eachSlice(loci: Loci, volume: Volume, key: number, props: PD.Values<SliceParams>, apply: (interval: Interval) => boolean) {
    return eachVolumeLoci(loci, volume, undefined, apply);
}

//

export const SliceParams = {
    ...Image.Params,
    quality: { ...Image.Params.quality, isEssential: false },
    dimension: PD.MappedStatic('x', {
        x: PD.Numeric(0, { min: 0, max: 0, step: 1 }),
        y: PD.Numeric(0, { min: 0, max: 0, step: 1 }),
        z: PD.Numeric(0, { min: 0, max: 0, step: 1 }),
    }, { isEssential: true }),
    isoValue: Volume.IsoValueParam,
};
export type SliceParams = typeof SliceParams
export function getSliceParams(ctx: ThemeRegistryContext, volume: Volume) {
    const p = PD.clone(SliceParams);
    const dim = volume.grid.cells.space.dimensions;
    p.dimension = PD.MappedStatic('x', {
        x: PD.Numeric(0, { min: 0, max: dim[0] - 1, step: 1 }),
        y: PD.Numeric(0, { min: 0, max: dim[1] - 1, step: 1 }),
        z: PD.Numeric(0, { min: 0, max: dim[2] - 1, step: 1 }),
    }, { isEssential: true });
    p.isoValue = Volume.createIsoValueParam(Volume.IsoValue.absolute(volume.grid.stats.min), volume.grid.stats);
    return p;
}

export function SliceVisual(materialId: number): VolumeVisual<SliceParams> {
    return VolumeVisual<Image, SliceParams>({
        defaultProps: PD.getDefaultValues(SliceParams),
        createGeometry: createImage,
        createLocationIterator: (volume: Volume) => LocationIterator(volume.grid.cells.data.length, 1, 1, () => NullLocation),
        getLoci: getSliceLoci,
        eachLocation: eachSlice,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<SliceParams>, currentProps: PD.Values<SliceParams>, newTheme: Theme, currentTheme: Theme) => {
            state.createGeometry = (
                newProps.dimension.name !== currentProps.dimension.name ||
                newProps.dimension.params !== currentProps.dimension.params ||
                !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats) ||
                !ColorTheme.areEqual(newTheme.color, currentTheme.color)
            );
        },
        geometryUtils: {
            ...Image.Utils,
            createRenderableState: (props: PD.Values<SliceParams>) => {
                const state = Image.Utils.createRenderableState(props);
                updateRenderableState(state, props);
                return state;
            },
            updateRenderableState
        }
    }, materialId);
}

function updateRenderableState(state: RenderableState, props: PD.Values<SliceParams>) {
    Image.Utils.updateRenderableState(state, props);
    state.opaque = false;
    state.writeDepth = true;
}

export function SliceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, SliceParams>): VolumeRepresentation<SliceParams> {
    return VolumeRepresentation('Slice', ctx, getParams, SliceVisual, getLoci);
}

export const SliceRepresentationProvider = VolumeRepresentationProvider({
    name: 'slice',
    label: 'Slice',
    description: 'Slice of volume rendered as image with interpolation.',
    factory: SliceRepresentation,
    getParams: getSliceParams,
    defaultValues: PD.getDefaultValues(SliceParams),
    defaultColorTheme: { name: 'uniform' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (volume: Volume) => !Volume.isEmpty(volume) && !Volume.Segmentation.get(volume)
});