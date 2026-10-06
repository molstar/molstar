/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { type Structure, Model } from '@molstar/model/model/structure';
import { Representation, type RepresentationContext, type RepresentationParamsGetter } from '../../representation.js';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ComplexRepresentation } from '../complex-representation.js';
import {
  type StructureRepresentation,
  StructureRepresentationProvider,
  StructureRepresentationStateBuilder,
} from '../representation.js';
import { CarbohydrateLinkParams, CarbohydrateLinkVisual } from '../visual/carbohydrate-link-cylinder.js';
import { CarbohydrateSymbolParams, CarbohydrateSymbolVisual } from '../visual/carbohydrate-symbol-mesh.js';
import {
  CarbohydrateTerminalLinkParams,
  CarbohydrateTerminalLinkVisual,
} from '../visual/carbohydrate-terminal-link-cylinder.js';

const CarbohydrateVisuals = {
  'carbohydrate-symbol': (
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<Structure, CarbohydrateSymbolParams>,
  ) => ComplexRepresentation('Carbohydrate symbol mesh', ctx, getParams, CarbohydrateSymbolVisual),
  'carbohydrate-link': (
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<Structure, CarbohydrateLinkParams>,
  ) => ComplexRepresentation('Carbohydrate link cylinder', ctx, getParams, CarbohydrateLinkVisual),
  'carbohydrate-terminal-link': (
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<Structure, CarbohydrateTerminalLinkParams>,
  ) => ComplexRepresentation('Carbohydrate terminal link cylinder', ctx, getParams, CarbohydrateTerminalLinkVisual),
};

export const CarbohydrateParams = {
  ...CarbohydrateSymbolParams,
  ...CarbohydrateLinkParams,
  ...CarbohydrateTerminalLinkParams,
  visuals: PD.MultiSelect(
    ['carbohydrate-symbol', 'carbohydrate-link', 'carbohydrate-terminal-link'],
    PD.objectToOptions(CarbohydrateVisuals),
  ),
  bumpFrequency: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
  density: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
};
export type CarbohydrateParams = typeof CarbohydrateParams;
export function getCarbohydrateParams(ctx: ThemeRegistryContext, structure: Structure) {
  return CarbohydrateParams;
}

export type CarbohydrateRepresentation = StructureRepresentation<CarbohydrateParams>;
export function CarbohydrateRepresentation(
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<Structure, CarbohydrateParams>,
): CarbohydrateRepresentation {
  return Representation.createMulti(
    'Carbohydrate',
    ctx,
    getParams,
    StructureRepresentationStateBuilder,
    CarbohydrateVisuals as unknown as Representation.Def<Structure, CarbohydrateParams>,
  );
}

export const CarbohydrateRepresentationProvider = StructureRepresentationProvider({
  name: 'carbohydrate',
  label: 'Carbohydrate',
  description: 'Displays carbohydrate symbols (3D SNFG).',
  factory: CarbohydrateRepresentation,
  getParams: getCarbohydrateParams,
  defaultValues: PD.getDefaultValues(CarbohydrateParams),
  defaultColorTheme: { name: 'carbohydrate-symbol' },
  defaultSizeTheme: { name: 'uniform' },
  isApplicable: (structure: Structure) => {
    return structure.models.some((m) => Model.hasCarbohydrate(m));
  },
});
