// /**
//  * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author Alexander Rose <alexander.rose@weirdbyte.de>
//  */

// import { CrossLinkRestraintVisual, CrossLinkRestraintParams } from '../visual/cross-link-restraint-cylinder';
// import { ParamDefinition as PD } from 'mol-util/param-definition';
// import { ComplexRepresentation } from '../complex-representation';
// import { StructureRepresentation } from '../representation';
// import { Representation } from 'mol-repr/representation';
// import { ThemeRegistryContext } from 'mol-theme/theme';
// import { Structure } from 'mol-model/structure';

// export const DistanceRestraintParams = {
//     ...CrossLinkRestraintParams,
// }
// export function getDistanceRestraintParams(ctx: ThemeRegistryContext, structure: Structure) {
//     return DistanceRestraintParams // TODO return copy
// }
// export type DistanceRestraintProps = PD.DefaultValues<typeof DistanceRestraintParams>

// export type DistanceRestraintRepresentation = StructureRepresentation<DistanceRestraintProps>

// export function DistanceRestraintRepresentation(defaultProps: DistanceRestraintProps): DistanceRestraintRepresentation {
//     return Representation.createMulti('Distance restraint', defaultProps, [
//         ComplexRepresentation('Cross-link restraint', defaultProps, CrossLinkRestraintVisual)
//     ])
// }