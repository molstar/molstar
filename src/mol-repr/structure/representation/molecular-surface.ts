// /**
//  * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author Alexander Rose <alexander.rose@weirdbyte.de>
//  */

// import { GaussianSurfaceVisual, GaussianSurfaceParams } from '../visual/gaussian-surface-mesh';
// import { UnitsRepresentation } from '../units-representation';
// import { GaussianWireframeVisual, GaussianWireframeParams } from '../visual/gaussian-surface-wireframe';
// import { ParamDefinition as PD } from 'mol-util/param-definition';
// import { GaussianDensityVolumeParams, GaussianDensityVolumeVisual } from '../visual/gaussian-density-volume';
// import { StructureRepresentation } from '../representation';
// import { Representation } from 'mol-repr/representation';
// import { ThemeRegistryContext } from 'mol-theme/theme';
// import { Structure } from 'mol-model/structure';

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