// /**
//  * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author Alexander Rose <alexander.rose@weirdbyte.de>
//  */

// import { PolymerBackboneVisual, PolymerBackboneParams } from '../visual/polymer-backbone-cylinder';
// import { ParamDefinition as PD } from 'mol-util/param-definition';
// import { UnitsRepresentation } from '../units-representation';
// import { StructureRepresentation } from '../representation';
// import { Representation } from 'mol-repr/representation';
// import { ThemeRegistryContext } from 'mol-theme/theme';
// import { Structure } from 'mol-model/structure';

// export const BackboneParams = {
//     ...PolymerBackboneParams,
// }
// export function getBackboneParams(ctx: ThemeRegistryContext, structure: Structure) {
//     return BackboneParams // TODO return copy
// }
// export type BackboneProps = PD.DefaultValues<typeof BackboneParams>

// export type BackboneRepresentation = StructureRepresentation<BackboneProps>

// export function BackboneRepresentation(defaultProps: BackboneProps): BackboneRepresentation {
//     return Representation.createMulti('Backbone', defaultProps, [
//         UnitsRepresentation('Polymer backbone cylinder', defaultProps, PolymerBackboneVisual)
//     ])
// }