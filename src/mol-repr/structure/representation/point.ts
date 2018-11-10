// /**
//  * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author Alexander Rose <alexander.rose@weirdbyte.de>
//  */

// import { ElementPointVisual, ElementPointParams } from '../visual/element-point';
// import { UnitsRepresentation } from '../units-representation';
// import { ParamDefinition as PD } from 'mol-util/param-definition';
// import { StructureRepresentation } from '../representation';
// import { Representation } from 'mol-repr/representation';
// import { ThemeRegistryContext } from 'mol-theme/theme';
// import { Structure } from 'mol-model/structure';

// export const PointParams = {
//     ...ElementPointParams,
// }
// export function getPointParams(ctx: ThemeRegistryContext, structure: Structure) {
//     return PointParams // TODO return copy
// }
// export type PointProps = PD.DefaultValues<typeof PointParams>

// export type PointRepresentation = StructureRepresentation<PointProps>

// export function PointRepresentation(defaultProps: PointProps): PointRepresentation {
//     return Representation.createMulti('Point', defaultProps, [
//         UnitsRepresentation('Point', defaultProps, ElementPointVisual)
//     ])
// }