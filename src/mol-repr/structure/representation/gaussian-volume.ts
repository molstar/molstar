/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { GaussianDensityVolumeParams, GaussianDensityVolumeVisual } from '../visual/gaussian-density-volume';
import { StructureRepresentation, StructureRepresentationProvider, ComplexRepresentation, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';

const GaussianVolumeVisuals = {
    'gaussian-volume': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, GaussianDensityVolumeParams>) => ComplexRepresentation('Gaussian volume', ctx, getParams, GaussianDensityVolumeVisual)
};

export const GaussianVolumeParams = {
    ...GaussianDensityVolumeParams,
};
export type GaussianVolumeParams = typeof GaussianVolumeParams
export function getGaussianVolumeParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(GaussianVolumeParams);
}

export type GaussianVolumeRepresentation = StructureRepresentation<GaussianVolumeParams>
export function GaussianVolumeRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, GaussianVolumeParams>): GaussianVolumeRepresentation {
    return Representation.createMulti('Gaussian Volume', ctx, getParams, StructureRepresentationStateBuilder, GaussianVolumeVisuals as unknown as Representation.Def<Structure, GaussianVolumeParams>);
}

export const GaussianVolumeRepresentationProvider = StructureRepresentationProvider({
    name: 'gaussian-volume',
    label: 'Gaussian Volume',
    description: 'Displays a gaussian molecular density using direct volume rendering.',
    factory: GaussianVolumeRepresentation,
    getParams: getGaussianVolumeParams,
    defaultValues: PD.getDefaultValues(GaussianVolumeParams),
    defaultColorTheme: { name: 'chain-id' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.elementCount > 0
});