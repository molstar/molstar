// /**
//  * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author Alexander Rose <alexander.rose@weirdbyte.de>
//  */

// import { CarbohydrateSymbolVisual, CarbohydrateSymbolParams } from '../visual/carbohydrate-symbol-mesh';
// import { CarbohydrateLinkVisual, CarbohydrateLinkParams } from '../visual/carbohydrate-link-cylinder';
// import { ParamDefinition as PD } from 'mol-util/param-definition';
// import { ComplexRepresentation } from '../complex-representation';
// import { StructureRepresentation } from '../representation';
// import { Representation } from 'mol-repr/representation';
// import { ThemeRegistryContext } from 'mol-theme/theme';
// import { Structure } from 'mol-model/structure';

// export const CarbohydrateParams = {
//     ...CarbohydrateSymbolParams,
//     ...CarbohydrateLinkParams,
// }
// export function getCarbohydrateParams(ctx: ThemeRegistryContext, structure: Structure) {
//     return CarbohydrateParams // TODO return copy
// }
// export type CarbohydrateProps = PD.DefaultValues<typeof CarbohydrateParams>

// export type CarbohydrateRepresentation = StructureRepresentation<CarbohydrateProps>

// export function CarbohydrateRepresentation(defaultProps: CarbohydrateProps): CarbohydrateRepresentation {
//     return Representation.createMulti('Carbohydrate', defaultProps, [
//         ComplexRepresentation('Carbohydrate symbol mesh', defaultProps, CarbohydrateSymbolVisual),
//         ComplexRepresentation('Carbohydrate link cylinder', defaultProps, CarbohydrateLinkVisual)
//     ])
// }