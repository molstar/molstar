/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from 'mol-util/param-definition';
import { GaussianDensityVolumeParams, GaussianDensityVolumeVisual } from '../visual/gaussian-density-volume';
import { StructureRepresentation, StructureRepresentationProvider, ComplexRepresentation, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from 'mol-repr/representation';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { Structure } from 'mol-model/structure';

const MolecularVolumeVisuals = {
    'gaussian-volume': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, GaussianDensityVolumeParams>) => ComplexRepresentation('Gaussian volume', ctx, getParams, GaussianDensityVolumeVisual)
}

export const MolecularVolumeParams = {
    ...GaussianDensityVolumeParams,
}
export type MolecularVolumeParams = typeof MolecularVolumeParams
export function getMolecularVolumeParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(MolecularVolumeParams)
}

export type MolecularVolumeRepresentation = StructureRepresentation<MolecularVolumeParams>
export function MolecularVolumeRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MolecularVolumeParams>): MolecularVolumeRepresentation {
    return Representation.createMulti('Molecular Volume', ctx, getParams, StructureRepresentationStateBuilder, MolecularVolumeVisuals as unknown as Representation.Def<Structure, MolecularVolumeParams>)
}

export const MolecularVolumeRepresentationProvider: StructureRepresentationProvider<MolecularVolumeParams> = {
    label: 'Molecular Volume',
    description: 'Displays a gaussian molecular density using direct volume rendering.',
    factory: MolecularVolumeRepresentation,
    getParams: getMolecularVolumeParams,
    defaultValues: PD.getDefaultValues(MolecularVolumeParams),
    defaultColorTheme: 'polymer-id',
    defaultSizeTheme: 'uniform'
}