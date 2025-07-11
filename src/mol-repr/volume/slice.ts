/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Image } from '../../mol-geo/geometry/image/image';
import { ThemeRegistryContext, Theme } from '../../mol-theme/theme';
import { Grid, Volume } from '../../mol-model/volume';
import { VolumeVisual, VolumeRepresentation, VolumeRepresentationProvider } from './representation';
import { PositionLocation } from '../../mol-geo/util/location-iterator';
import { VisualUpdateState } from '../util';
import { RepresentationContext, RepresentationParamsGetter } from '../representation';
import { VisualContext } from '../visual';
import { PickingId } from '../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../mol-model/loci';
import { Interval, OrderedSet, SortedArray } from '../../mol-data/int';
import { transformPositionArray } from '../../mol-geo/util';
import { Color } from '../../mol-util/color';
import { ColorTheme } from '../../mol-theme/color';
import { packIntToRGBArray } from '../../mol-util/number-packing';
import { createVolumeCellLocationIterator, eachVolumeLoci } from './util';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { Quat } from '../../mol-math/linear-algebra/3d/quat';
import { degToRad } from '../../mol-math/misc';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { clamp, normalize } from '../../mol-math/interpolate';
import { assertUnreachable } from '../../mol-util/type-helpers';

export const SliceParams = {
    ...Image.Params,
    quality: { ...Image.Params.quality, isEssential: false },
    dimension: PD.MappedStatic('x', {
        x: PD.Numeric(0, { min: 0, max: 0, step: 1 }, { immediateUpdate: true }),
        y: PD.Numeric(0, { min: 0, max: 0, step: 1 }, { immediateUpdate: true }),
        z: PD.Numeric(0, { min: 0, max: 0, step: 1 }, { immediateUpdate: true }),
        relativeX: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }, { immediateUpdate: true }),
        relativeY: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }, { immediateUpdate: true }),
        relativeZ: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }, { immediateUpdate: true }),
    }, { isEssential: true, hideIf: p => p.mode !== 'grid', description: 'Slice position in grid coordinates.' }),
    isoValue: Volume.IsoValueParam,
    mode: PD.Select('grid', PD.arrayToOptions(['grid', 'frame', 'plane'] as const), { description: 'Grid: slice through the volume along the grid axes in integer steps. Frame: slice through the volume along arbitrary axes in any step size. Plane: an arbitrary plane defined by point and normal.' }),
    offset: PD.Numeric(0, { min: -1, max: 1, step: 0.01 }, { isEssential: true, immediateUpdate: true, hideIf: p => p.mode !== 'frame', description: 'Relative offset from center.' }),
    axis: PD.Select('a', PD.arrayToOptions(['a', 'b', 'c'] as const), { isEssential: true, hideIf: p => p.mode !== 'frame', description: 'Axis of the frame.' }),
    rotation: PD.Group({
        axis: PD.Vec3(Vec3.create(1, 0, 0), {}, { description: 'Axis of rotation' }),
        angle: PD.Numeric(0, { min: -180, max: 180, step: 1 }, { immediateUpdate: true, description: 'Axis rotation angle in Degrees' }),
    }, { isExpanded: true, hideIf: p => p.mode !== 'frame' }),
    plane: PD.Group({
        point: PD.Vec3(Vec3.create(0, 0, 0), {}, { description: 'Plane point' }),
        normal: PD.Vec3(Vec3.create(1, 0, 0), {}, { description: 'Plane normal' }),
    }, { isExpanded: true, hideIf: p => p.mode !== 'plane' }),
};
export type SliceParams = typeof SliceParams
export type SliceProps = PD.Values<SliceParams>

export function getSliceParams(ctx: ThemeRegistryContext, volume: Volume) {
    const p = PD.clone(SliceParams);
    const dim = volume.grid.cells.space.dimensions;
    p.dimension = PD.MappedStatic('x', {
        x: PD.Numeric(0, { min: 0, max: dim[0] - 1, step: 1 }, { immediateUpdate: true }),
        y: PD.Numeric(0, { min: 0, max: dim[1] - 1, step: 1 }, { immediateUpdate: true }),
        z: PD.Numeric(0, { min: 0, max: dim[2] - 1, step: 1 }, { immediateUpdate: true }),
        relativeX: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }, { immediateUpdate: true }),
        relativeY: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }, { immediateUpdate: true }),
        relativeZ: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }, { immediateUpdate: true }),
    }, { isEssential: true, hideIf: p => p.mode !== 'grid' });
    p.isoValue = Volume.createIsoValueParam(Volume.IsoValue.absolute(volume.grid.stats.min), volume.grid.stats);
    return p;
}

export async function createImage(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: SliceProps, image?: Image): Promise<Image> {
    switch (props.mode) {
        case 'frame':
            return createFrameImage(ctx, volume, key, theme, props, image);
        case 'grid':
            return createGridImage(ctx, volume, key, theme, props, image);
        case 'plane':
            return createPlaneImage(ctx, volume, key, theme, props, image);
        default:
            assertUnreachable(props.mode);
    }
}

//

function getFrame(volume: Volume, props: SliceProps) {
    const { axis, rotation, mode } = props;

    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const cartnToGrid = Mat4.invert(Mat4(), gridToCartn);
    const [nx, ny, nz] = volume.grid.cells.space.dimensions;

    const a = nx - 1;
    const b = ny - 1;
    const c = nz - 1;

    const dirA = Vec3.create(a, 0, 0);
    const dirB = Vec3.create(0, b, 0);
    const dirC = Vec3.create(0, 0, c);

    const resolution = Math.max(a, b, c) / Math.max(nx, ny, nz);

    const min = Vec3.create(0, 0, 0);
    const max = Vec3.create(a, b, c);
    Vec3.transformMat4(min, min, gridToCartn);
    Vec3.transformMat4(max, max, gridToCartn);
    const center = Vec3.center(Vec3(), max, min);

    const size = Vec3();
    const major = Vec3();
    const minor = Vec3();
    const normal = Vec3();

    if (axis === 'c') {
        Vec3.set(size, a, b, c);
        Vec3.copy(major, dirA);
        Vec3.copy(minor, dirB);
        Vec3.copy(normal, dirC);
    } else if (axis === 'b') {
        Vec3.set(size, a, c, b);
        Vec3.copy(major, dirA);
        Vec3.copy(normal, dirB);
        Vec3.copy(minor, dirC);
    } else {
        Vec3.set(size, b, c, a);
        Vec3.copy(normal, dirA);
        Vec3.copy(major, dirB);
        Vec3.copy(minor, dirC);
    }

    if (rotation.angle !== 0) {
        const ra = Vec3();
        Vec3.scaleAndAdd(ra, ra, dirA, rotation.axis[0]);
        Vec3.scaleAndAdd(ra, ra, dirB, rotation.axis[1]);
        Vec3.scaleAndAdd(ra, ra, dirC, rotation.axis[2]);
        Vec3.normalize(ra, ra);

        const rm = Mat4.fromRotation(Mat4(), degToRad(rotation.angle), ra);
        Vec3.transformDirection(major, major, rm);
        Vec3.transformDirection(minor, minor, rm);
        Vec3.transformDirection(normal, normal, rm);
    }

    if (mode === 'frame') {
        const r = Vec3.distance(min, max);
        const s = Vec3.distance(min, max) * Math.SQRT2;
        Vec3.set(size, s, s, r);
    }

    Vec3.transformDirection(major, major, gridToCartn);
    Vec3.transformDirection(minor, minor, gridToCartn);
    Vec3.transformDirection(normal, normal, gridToCartn);

    const trim: Image.Trim = {
        type: 3,
        center: Vec3.create(a / 2, b / 2, c / 2),
        scale: Vec3.create(a, b, c),
        rotation: Quat.identity(),
        transform: cartnToGrid,
    };

    return { size, major, minor, normal, center, trim, resolution };
}

type SamplingInfo = {
    m: Mat4
    width: number
    height: number
};

function getSampledImage(volume: Volume, theme: Theme, info: SamplingInfo, isoValue: Volume.IsoValue, trim: Image.Trim, image?: Image): Image {
    const { m, width, height } = info;

    const { cells: { space }, stats } = volume.grid;
    const { min, max } = stats;

    const isUniform = theme.color.granularity === 'uniform';
    const color = 'color' in theme.color && theme.color.color
        ? theme.color.color
        : () => Color(0xffffff);

    const v = Vec3();

    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const cartnToGrid = Mat4.invert(Mat4(), gridToCartn);
    const [mx, my, mz] = space.dimensions;

    const imageArray = new Uint8Array(width * height * 4);
    const groupArray = new Uint8Array(width * height * 4);
    const valueArray = new Float32Array(width * height);

    const gridCoords = Vec3();

    const pl = PositionLocation(Vec3(), Vec3());
    const getTrilinearlyInterpolated = Grid.makeGetTrilinearlyInterpolated(volume.grid, 'none');

    let i = 0;
    for (let ih = 0; ih < height; ++ih) {
        for (let iw = 0; iw < width; ++iw) {
            const y = (clamp(iw + 0.5, 0, width - 1) / width) - 0.5;
            const x = (clamp(ih + 0.5, 0, height - 1) / height) - 0.5;
            Vec3.set(v, x, -y, 0);
            Vec3.transformMat4(v, v, m);

            Vec3.copy(gridCoords, v);
            Vec3.transformMat4(gridCoords, gridCoords, cartnToGrid);

            const ix = Math.trunc(gridCoords[0]);
            const iy = Math.trunc(gridCoords[1]);
            const iz = Math.trunc(gridCoords[2]);

            Vec3.copy(pl.position, v);
            const c = color(pl, false);
            Color.toArray(c, imageArray, i);

            const val = normalize(getTrilinearlyInterpolated(v), min, max);
            if (isUniform) {
                imageArray[i] *= val;
                imageArray[i + 1] *= val;
                imageArray[i + 2] *= val;
            }
            valueArray[i / 4] = val;

            if (ix >= 0 && ix < mx && iy >= 0 && iy < my && iz >= 0 && iz < mz) {
                packIntToRGBArray(space.dataOffset(ix, iy, iz), groupArray, i);
            }

            i += 4;
        }
    }

    const imageTexture = { width, height, array: imageArray, flipY: true };
    const groupTexture = { width, height, array: groupArray, flipY: true };
    const valueTexture = { width, height, array: valueArray, flipY: true };

    const corners = new Float32Array([
        -0.5, 0.5, 0,
        0.5, 0.5, 0,
        -0.5, -0.5, 0,
        0.5, -0.5, 0
    ]);
    transformPositionArray(m, corners, 0, 4);

    const isoLevel = clamp(normalize(Volume.IsoValue.toAbsolute(isoValue, stats).absoluteValue, min, max), 0, 1);

    const im = Image.create(imageTexture, corners, groupTexture, valueTexture, trim, isoLevel, image);
    im.setBoundingSphere(Volume.getBoundingSphere(volume));
    return im;
}

async function createFrameImage(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: SliceProps, image?: Image): Promise<Image> {
    const { offset, isoValue } = props;

    const { size, major, minor, normal, center, trim, resolution } = getFrame(volume, props);
    const scaleFactor = 1 / resolution;

    const scale = Vec3.create(size[0], size[1], 1);
    const offsetDir = Vec3.setMagnitude(Vec3(), normal, size[2] / 2);

    const width = Math.floor(size[1] * scaleFactor);
    const height = Math.floor(size[0] * scaleFactor);

    const m = Mat4.identity();
    const v = Vec3();
    const anchor = Vec3();

    Vec3.add(v, center, major);
    Mat4.targetTo(m, center, v, minor);
    Vec3.scaleAndAdd(anchor, center, offsetDir, offset);
    Mat4.setTranslation(m, anchor);
    Mat4.mul(m, m, Mat4.rotY90);
    Mat4.scale(m, m, scale);

    const info = { m, width, height };

    return getSampledImage(volume, theme, info, isoValue, trim, image);
}

async function createPlaneImage(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: SliceProps, image?: Image): Promise<Image> {
    const { plane: { point, normal }, isoValue } = props;

    const m = Mat4.fromPlane(Mat4(), normal, point);

    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const cartnToGrid = Mat4.invert(Mat4(), gridToCartn);
    const [mx, my, mz] = volume.grid.cells.space.dimensions;

    const a = mx - 1;
    const b = my - 1;
    const c = mz - 1;

    const resolution = Math.max(a, b, c) / Math.max(mx, my, mz);
    const scaleFactor = 1 / resolution;

    const s = Vec3.distance(Vec3.create(0, 0, 0), Vec3.create(a, b, c)) * Math.SQRT2;
    const size = Vec3.set(Vec3(), s, s, s);

    Mat4.mul(m, m, Mat4.rotY90);
    Mat4.scale(m, m, size);

    const width = Math.floor(size[1] * scaleFactor);
    const height = Math.floor(size[0] * scaleFactor);

    const trim: Image.Trim = {
        type: 3,
        center: Vec3.create(a / 2, b / 2, c / 2),
        scale: Vec3.create(a, b, c),
        rotation: Quat.identity(),
        transform: cartnToGrid,
    };

    const info: SamplingInfo = { m, width, height };

    return getSampledImage(volume, theme, info, isoValue, trim, image);
}

async function createGridImage(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: SliceProps, image?: Image): Promise<Image> {
    const { dimension: { name: dim }, isoValue } = props;

    const { cells: { space, data }, stats } = volume.grid;
    const { min, max } = stats;

    const isUniform = theme.color.granularity === 'uniform';
    const color = 'color' in theme.color && theme.color.color
        ? theme.color.color
        : () => Color(0xffffff);

    const {
        width, height,
        x, y, z,
        x0, y0, z0,
        nx, ny, nz
    } = getSliceInfo(volume.grid, props);

    const corners = new Float32Array(
        dim === 'x' || dim === 'relativeX' ? [x, 0, 0, x, y, 0, x, 0, z, x, y, z] :
            dim === 'y' || dim === 'relativeY' ? [0, y, 0, x, y, 0, 0, y, z, x, y, z] :
                [0, 0, z, 0, y, z, x, 0, z, x, y, z]
    );

    const imageArray = new Uint8Array(width * height * 4);
    const groupArray = getPackedGroupArray(volume.grid, props);
    const valueArray = new Float32Array(width * height);

    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const l = PositionLocation(Vec3(), Vec3());

    let i = 0;
    for (let iy = y0; iy < ny; ++iy) {
        for (let ix = x0; ix < nx; ++ix) {
            for (let iz = z0; iz < nz; ++iz) {
                Vec3.set(l.position, ix, iy, iz);
                Vec3.transformMat4(l.position, l.position, gridToCartn);
                Color.toArray(color(l, false), imageArray, i);

                const val = normalize(space.get(data, ix, iy, iz), min, max);
                if (isUniform) {
                    imageArray[i] *= val;
                    imageArray[i + 1] *= val;
                    imageArray[i + 2] *= val;
                }
                valueArray[i / 4] = val;

                i += 4;
            }
        }
    }

    const imageTexture = { width, height, array: imageArray, flipY: true };
    const groupTexture = { width, height, array: groupArray, flipY: true };
    const valueTexture = { width, height, array: valueArray, flipY: true };

    const transform = Grid.getGridToCartesianTransform(volume.grid);
    transformPositionArray(transform, corners, 0, 4);

    const trim = Image.createEmptyTrim();
    const isoLevel = clamp(normalize(Volume.IsoValue.toAbsolute(isoValue, stats).absoluteValue, min, max), 0, 1);

    const im = Image.create(imageTexture, corners, groupTexture, valueTexture, trim, isoLevel, image);
    im.setBoundingSphere(Volume.getBoundingSphere(volume));
    return im;
}

//

function getRelativeIndex(dim: number, index: number) {
    return clamp(Math.round((dim - 1) * index), 0, dim - 1);
}

function getSliceInfo(grid: Grid, props: SliceProps) {
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
    } else if (dim === 'z') {
        x = nx - 1, y = ny - 1, z = index;
        width = nx, height = ny;
        z0 = z, nz = z0 + 1;
    } else if (dim === 'relativeX') {
        x = getRelativeIndex(nx, index);
        y = ny - 1;
        z = nz - 1;
        width = nz, height = ny;
        x0 = x, nx = x0 + 1;
    } else if (dim === 'relativeY') {
        x = nx - 1;
        y = getRelativeIndex(ny, index);
        z = nz - 1;
        width = nz, height = nx;
        y0 = y, ny = y0 + 1;
    } else /* if (dim === 'relativeZ') */ {
        x = nx - 1;
        y = ny - 1;
        z = getRelativeIndex(nz, index);
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

function getPackedGroupArray(grid: Grid, props: SliceProps) {
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

function getGroupArray(grid: Grid, props: SliceProps) {
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

function getObjectLoci(volume: Volume, instances: OrderedSet<Volume.InstanceIndex>, props: SliceProps) {
    // TODO: cache somehow?
    if (props.mode === 'grid') {
        const groupArray = getGroupArray(volume.grid, props);
        const indices = SortedArray.ofUnsortedArray<Volume.CellIndex>(groupArray);
        return Volume.Cell.Loci(volume, [{ indices, instances }]);
    } else {
        // TODO: get exact groups
        return Volume.Loci(volume, instances);
    }
}

function getLoci(volume: Volume, props: SliceProps) {
    const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
    return getObjectLoci(volume, instances, props);
}

function getSliceLoci(pickingId: PickingId, volume: Volume, key: number, props: SliceProps, id: number) {
    const { objectId, groupId, instanceId } = pickingId;
    if (id === objectId) {
        const granularity = Volume.PickingGranularity.get(volume);
        const instances = OrderedSet.ofSingleton(instanceId as Volume.InstanceIndex);
        if (granularity === 'volume') {
            return Volume.Loci(volume, instances);
        } if (granularity === 'object') {
            return getObjectLoci(volume, instances, props);
        } else {
            const indices = Interval.ofSingleton(groupId as Volume.CellIndex);
            return Volume.Cell.Loci(volume, [{ indices, instances }]);
        }
    }
    return EmptyLoci;
}

function eachSlice(loci: Loci, volume: Volume, key: number, props: SliceProps, apply: (interval: Interval) => boolean) {
    return eachVolumeLoci(loci, volume, undefined, apply);
}

//

export function SliceVisual(materialId: number): VolumeVisual<SliceParams> {
    return VolumeVisual<Image, SliceParams>({
        defaultProps: PD.getDefaultValues(SliceParams),
        createGeometry: createImage,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getSliceLoci,
        eachLocation: eachSlice,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: SliceProps, currentProps: SliceProps, newTheme: Theme, currentTheme: Theme) => {
            state.createGeometry = (
                newProps.dimension.name !== currentProps.dimension.name ||
                newProps.dimension.params !== currentProps.dimension.params ||
                newProps.mode !== currentProps.mode ||
                !Vec3.equals(newProps.rotation.axis, currentProps.rotation.axis) ||
                newProps.rotation.angle !== currentProps.rotation.angle ||
                newProps.offset !== currentProps.offset ||
                newProps.axis !== currentProps.axis ||
                !Vec3.equals(newProps.plane.point, currentProps.plane.point) ||
                !Vec3.equals(newProps.plane.normal, currentProps.plane.normal) ||
                !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats) ||
                !ColorTheme.areEqual(newTheme.color, currentTheme.color)
            );
        },

        geometryUtils: Image.Utils
    }, materialId);
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