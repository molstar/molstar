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
import { EmptyLoci } from '../../mol-model/loci';
import { NullLocation } from '../../mol-model/location';
import { RepresentationContext, RepresentationParamsGetter } from '../representation';
import { VisualContext } from '../visual';

export async function createImage(ctx: VisualContext, volume: VolumeData, theme: Theme, props: PD.Values<SliceParams>, image?: Image) {
    const { space } = volume.data;

    const width = space.dimensions[0], height = space.dimensions[1];
    const n = width * height;

    // TODO fill with volume data values
    const imageTexture = { width, height, array: new Float32Array(n * 4) };
    for (let i = 0, il = n * 4; i < il; ++i) {
        imageTexture.array[i] = Math.random();
    }

    // TODO fill with linearized index into volume
    //      (to be used for picking which needs a volume location/loci)
    const groupTexture = { width, height, array: new Float32Array(n * 1) };

    // TODO four corners of a plane
    const corners = new Float32Array([
        0, 0, 0,
        0, 50, 0,
        0, 0, 50,
        0, 50, 50
    ]);

    return Image.create(imageTexture, corners, groupTexture, image);
}

//

export const SliceParams = {
    ...BaseGeometry.Params,
    ...Image.Params
};
export type SliceParams = typeof SliceParams
export function getSliceParams(ctx: ThemeRegistryContext, volume: VolumeData) {
    return PD.clone(SliceParams);
}

export function SliceVisual(materialId: number): VolumeVisual<SliceParams> {
    return VolumeVisual<Image, SliceParams>({
        defaultProps: PD.getDefaultValues(SliceParams),
        createGeometry: createImage,
        createLocationIterator: (volume: VolumeData) => LocationIterator(1, 1, () => NullLocation),
        getLoci: () => EmptyLoci,
        eachLocation: () => false,
        setUpdateState: (state: VisualUpdateState, volume: VolumeData, newProps: PD.Values<SliceParams>, currentProps: PD.Values<SliceParams>) => {
        },
        geometryUtils: Image.Utils
    }, materialId);
}

export function SliceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<VolumeData, SliceParams>): VolumeRepresentation<SliceParams> {
    return VolumeRepresentation('Slice', ctx, getParams, SliceVisual);
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