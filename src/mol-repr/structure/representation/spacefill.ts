// /**
//  * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author Alexander Rose <alexander.rose@weirdbyte.de>
//  */

// import { ElementSphereVisual, ElementSphereParams } from '../visual/element-sphere';
// import { UnitsRepresentation } from '../units-representation';
// import { ParamDefinition as PD } from 'mol-util/param-definition';
// import { StructureRepresentation } from '../representation';
// import { Representation } from 'mol-repr/representation';
// import { ThemeRegistryContext } from 'mol-theme/theme';
// import { Structure } from 'mol-model/structure';

// export const SpacefillParams = {
//     ...ElementSphereParams,
// }
// export function getSpacefillParams(ctx: ThemeRegistryContext, structure: Structure) {
//     return SpacefillParams // TODO return copy
// }
// export type SpacefillProps = PD.DefaultValues<typeof SpacefillParams>

// export type SpacefillRepresentation = StructureRepresentation<SpacefillProps>

// export function SpacefillRepresentation(defaultProps: SpacefillProps): SpacefillRepresentation {
//     return Representation.createMulti('Spacefill', defaultProps, [
//         UnitsRepresentation('Sphere mesh', defaultProps, ElementSphereVisual)
//     ])
// }