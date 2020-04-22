/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Image } from '../../mol-geo/geometry/image/image';
import { BaseGeometry } from '../../mol-geo/geometry/base';
import { ThemeRegistryContext, Theme } from '../../mol-theme/theme';
import { VolumeData } from '../../mol-model/volume';
import { VolumeVisual, VolumeRepresentation, VolumeRepresentationProvider } from './representation';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { VisualUpdateState } from '../util';
import { NullLocation } from '../../mol-model/location';
import { RepresentationContext, RepresentationParamsGetter } from '../representation';
import { VisualContext } from '../visual';
import { Volume } from '../../mol-model/volume/volume';
import { PickingId } from '../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../mol-model/loci';
import { Interval, OrderedSet } from '../../mol-data/int';
import { fillSerial } from '../../mol-util/array';

export async function createImage(ctx: VisualContext, volume: VolumeData, theme: Theme, props: PD.Values<SliceParams>, image?: Image) {
    const dim = parseInt(props.dimension.name.toString());
    // const index = props.dimension.params;

    const { space } = volume.data;

    let width: number, height: number;
    if (dim === 0) {
        width = space.dimensions[1];
        height = space.dimensions[2];
    } else if (dim === 1) {
        width = space.dimensions[0];
        height = space.dimensions[2];
    } else {
        width = space.dimensions[0];
        height = space.dimensions[1];
    }

    const n = width * height;

    // TODO fill with volume data values
    const imageTexture = { width, height, array: new Float32Array(n * 4) };

    for (let i = 0, il = n * 4; i < il; i += 4) {
        imageTexture.array[i] = 0;
        imageTexture.array[i + 1] = Math.random();
        imageTexture.array[i + 2] = Math.random();
        imageTexture.array[i + 3] = 1;
    }

    // TODO fill with linearized index into volume
    //      (to be used for picking which needs a volume location/loci)
    const groupTexture = { width, height, array: fillSerial(new Float32Array(n)) };

    // TODO four corners of a plane
    const corners = new Float32Array([
        0, 0, 0,
        0, 50, 0,
        0, 0, 50,
        0, 50, 50
    ]);

    return Image.create(imageTexture, corners, groupTexture, image);
}

function getLoci(volume: VolumeData, props: PD.Values<SliceParams>) {
    // TODO only slice indices
    return Volume.Loci(volume);
}

function getSliceLoci(pickingId: PickingId, volume: VolumeData, props: PD.Values<SliceParams>, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        return Volume.Cell.Loci(volume, Interval.ofSingleton(groupId as Volume.CellIndex));
    }
    return EmptyLoci;
}

function eachSlice(loci: Loci, volume: VolumeData, props: PD.Values<SliceParams>, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (Volume.isLoci(loci)) {
        if (!VolumeData.areEquivalent(loci.volume, volume)) return false;
        if (apply(Interval.ofLength(volume.data.data.length))) changed = true;
    } else if (Volume.Isosurface.isLoci(loci)) {
        if (!VolumeData.areEquivalent(loci.volume, volume)) return false;
        // TODO check isoValue
        if (apply(Interval.ofLength(volume.data.data.length))) changed = true;
    } else if (Volume.Cell.isLoci(loci)) {
        if (!VolumeData.areEquivalent(loci.volume, volume)) return false;
        if (Interval.is(loci.indices)) {
            if (apply(loci.indices)) changed = true;
        } else {
            OrderedSet.forEach(loci.indices, v => {
                if (apply(Interval.ofSingleton(v))) changed = true;
            });
        }
    }
    return changed;
}

//

export const SliceParams = {
    ...BaseGeometry.Params,
    ...Image.Params,
    dimension: PD.MappedStatic(0, {
        0: PD.Numeric(0, { min: 0, max: 0, step: 1 }),
        1: PD.Numeric(0, { min: 0, max: 0, step: 1 }),
        2: PD.Numeric(0, { min: 0, max: 0, step: 1 }),
    }, { isEssential: true }),
};
export type SliceParams = typeof SliceParams
export function getSliceParams(ctx: ThemeRegistryContext, volume: VolumeData) {
    const p = PD.clone(SliceParams);
    p.dimension = PD.MappedStatic(0, {
        0: PD.Numeric(0, { min: 0, max: volume.data.space.dimensions[0], step: 1 }),
        1: PD.Numeric(0, { min: 0, max: volume.data.space.dimensions[1], step: 1 }),
        2: PD.Numeric(0, { min: 0, max: volume.data.space.dimensions[2], step: 1 }),
    }, { isEssential: true });
    return p;
}

export function SliceVisual(materialId: number): VolumeVisual<SliceParams> {
    return VolumeVisual<Image, SliceParams>({
        defaultProps: PD.getDefaultValues(SliceParams),
        createGeometry: createImage,
        createLocationIterator: (volume: VolumeData) => LocationIterator(volume.data.data.length, 1, () => NullLocation),
        getLoci: getSliceLoci,
        eachLocation: eachSlice,
        setUpdateState: (state: VisualUpdateState, volume: VolumeData, newProps: PD.Values<SliceParams>, currentProps: PD.Values<SliceParams>) => {
            state.createGeometry = (
                newProps.dimension.name !== currentProps.dimension.name ||
                newProps.dimension.params !== currentProps.dimension.params
            );
        },
        geometryUtils: Image.Utils
    }, materialId);
}

export function SliceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<VolumeData, SliceParams>): VolumeRepresentation<SliceParams> {
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
    isApplicable: (volume: VolumeData) => volume.data.data.length > 0
});