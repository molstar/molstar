// /**
//  * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author Alexander Rose <alexander.rose@weirdbyte.de>
//  */

// import { ElementSphereVisual, ElementSphereParams } from '../visual/element-sphere';
// import { IntraUnitLinkVisual, IntraUnitLinkParams } from '../visual/intra-unit-link-cylinder';
// import { InterUnitLinkVisual, InterUnitLinkParams } from '../visual/inter-unit-link-cylinder';
// import { ParamDefinition as PD } from 'mol-util/param-definition';
// import { UnitsRepresentation } from '../units-representation';
// import { ComplexRepresentation } from '../complex-representation';
// import { StructureRepresentation } from '../representation';
// import { Representation } from 'mol-repr/representation';
// import { ThemeRegistryContext } from 'mol-theme/theme';
// import { Structure } from 'mol-model/structure';

// export const BallAndStickParams = {
//     ...ElementSphereParams,
//     ...IntraUnitLinkParams,
//     ...InterUnitLinkParams,
//     // TODO
//     // unitKinds: PD.MultiSelect<UnitKind>('Unit Kind', '', ['atomic'], UnitKindOptions),
// }
// export function getBallAndStickParams(ctx: ThemeRegistryContext, structure: Structure) {
//     return BallAndStickParams // TODO return copy
// }
// export type BallAndStickProps = PD.DefaultValues<typeof BallAndStickParams>

// export type BallAndStickRepresentation = StructureRepresentation<BallAndStickProps>

// export function BallAndStickRepresentation(defaultProps: BallAndStickProps): BallAndStickRepresentation {
//     return Representation.createMulti('Ball & Stick', defaultProps, [
//         UnitsRepresentation('Element sphere mesh', defaultProps, ElementSphereVisual),
//         UnitsRepresentation('Intra-unit link cylinder', defaultProps, IntraUnitLinkVisual),
//         ComplexRepresentation('Inter-unit link cylinder', defaultProps, InterUnitLinkVisual)
//     ])
// }