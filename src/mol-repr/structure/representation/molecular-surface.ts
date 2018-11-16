/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GaussianSurfaceVisual, GaussianSurfaceParams } from '../visual/gaussian-surface-mesh';
import { UnitsRepresentation } from '../units-representation';
import { GaussianWireframeVisual, GaussianWireframeParams } from '../visual/gaussian-surface-wireframe';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { GaussianDensityVolumeParams, GaussianDensityVolumeVisual } from '../visual/gaussian-density-volume';
import { StructureRepresentation, StructureRepresentationProvider } from '../representation';
import { Representation, RepresentationParamsGetter } from 'mol-repr/representation';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { Structure } from 'mol-model/structure';
import { BuiltInColorThemeOptions, getBuiltInColorThemeParams } from 'mol-theme/color';

const MolecularSurfaceVisuals = {
    'gaussian-surface': (getParams: RepresentationParamsGetter<Structure, GaussianSurfaceParams>) => UnitsRepresentation('Gaussian surface', getParams, GaussianSurfaceVisual),
    'gaussian-wireframe': (getParams: RepresentationParamsGetter<Structure, GaussianWireframeParams>) => UnitsRepresentation('Gaussian wireframe', getParams, GaussianWireframeVisual),
    'gaussian-volume': (getParams: RepresentationParamsGetter<Structure, GaussianDensityVolumeParams>) => UnitsRepresentation('Gaussian volume', getParams, GaussianDensityVolumeVisual)
}
type MolecularSurfaceVisualName = keyof typeof MolecularSurfaceVisuals
const MolecularSurfaceVisualOptions = Object.keys(MolecularSurfaceVisuals).map(name => [name, name] as [MolecularSurfaceVisualName, string])

export const MolecularSurfaceParams = {
    ...GaussianSurfaceParams,
    ...GaussianWireframeParams,
    ...GaussianDensityVolumeParams,
    colorTheme: PD.Mapped('polymer-index', BuiltInColorThemeOptions, getBuiltInColorThemeParams),
    visuals: PD.MultiSelect<MolecularSurfaceVisualName>(['gaussian-surface'], MolecularSurfaceVisualOptions),
}
PD.getDefaultValues(MolecularSurfaceParams).colorTheme.name
export type MolecularSurfaceParams = typeof MolecularSurfaceParams
export function getMolecularSurfaceParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(MolecularSurfaceParams)
}

export type MolecularSurfaceRepresentation = StructureRepresentation<MolecularSurfaceParams>
export function MolecularSurfaceRepresentation(getParams: RepresentationParamsGetter<Structure, MolecularSurfaceParams>): MolecularSurfaceRepresentation {
    return Representation.createMulti('Molecular Surface', getParams, MolecularSurfaceVisuals as unknown as Representation.Def<Structure, MolecularSurfaceParams>)
}

export const MolecularSurfaceRepresentationProvider: StructureRepresentationProvider<MolecularSurfaceParams> = {
    label: 'Molecular Surface',
    description: 'Displays a gaussian molecular surface.',
    factory: MolecularSurfaceRepresentation,
    getParams: getMolecularSurfaceParams,
    defaultValues: PD.getDefaultValues(MolecularSurfaceParams)
}



// export const MolecularSurfaceParams = {
//     ...GaussianSurfaceParams,
//     ...GaussianWireframeParams,
//     ...GaussianDensityVolumeParams,
// }
// export function getMolecularSurfaceParams(ctx: ThemeRegistryContext, structure: Structure) {
//     return MolecularSurfaceParams // TODO return copy
// }
// export type MolecularSurfaceProps = PD.DefaultValues<typeof MolecularSurfaceParams>

// export type MolecularSurfaceRepresentation = StructureRepresentation<MolecularSurfaceProps>

// export function MolecularSurfaceRepresentation(defaultProps: MolecularSurfaceProps): MolecularSurfaceRepresentation {
//     return Representation.createMulti('Molecular Surface', defaultProps, [
//         UnitsRepresentation('Gaussian surface', defaultProps, GaussianSurfaceVisual),
//         UnitsRepresentation('Gaussian wireframe', defaultProps, GaussianWireframeVisual),
//         UnitsRepresentation('Gaussian volume', defaultProps, GaussianDensityVolumeVisual)
//     ])
// }