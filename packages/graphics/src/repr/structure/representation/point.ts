/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import {
  ElementPointVisual,
  ElementPointParams,
  StructureElementPointVisual,
  type StructureElementPointParams,
} from '../visual/element-point.js';
import { UnitsRepresentation } from '../units-representation.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  ComplexRepresentation,
  type StructureRepresentation,
  StructureRepresentationProvider,
  StructureRepresentationStateBuilder,
} from '../representation.js';
import { Representation, type RepresentationParamsGetter, type RepresentationContext } from '../../representation.js';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import type { Structure } from '@molstar/model/model/structure';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';

const PointVisuals = {
  'element-point': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, ElementPointParams>) =>
    UnitsRepresentation('Element points', ctx, getParams, ElementPointVisual),
  'structure-element-point': (
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<Structure, StructureElementPointParams>,
  ) => ComplexRepresentation('Structure element points', ctx, getParams, StructureElementPointVisual),
};

export const PointParams = {
  ...ElementPointParams,
  density: PD.Numeric(0.1, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
  visuals: PD.MultiSelect(['element-point'], PD.objectToOptions(PointVisuals)),
};
export type PointParams = typeof PointParams;
export function getPointParams(ctx: ThemeRegistryContext, structure: Structure) {
  let params = PointParams;
  if (structure.unitSymmetryGroups.length > 10000) {
    params = PD.clone(params);
    params.visuals.defaultValue = ['structure-element-point'];
  }
  return params;
}

export type PointRepresentation = StructureRepresentation<PointParams>;
export function PointRepresentation(
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<Structure, PointParams>,
): PointRepresentation {
  return Representation.createMulti(
    'Point',
    ctx,
    getParams,
    StructureRepresentationStateBuilder,
    PointVisuals as unknown as Representation.Def<Structure, PointParams>,
  );
}

export const PointRepresentationProvider = StructureRepresentationProvider({
  name: 'point',
  label: 'Point',
  description: 'Displays elements (atoms, coarse spheres) as points.',
  factory: PointRepresentation,
  getParams: getPointParams,
  defaultValues: PD.getDefaultValues(PointParams),
  defaultColorTheme: { name: 'element-symbol' },
  defaultSizeTheme: { name: 'uniform' },
  isApplicable: (structure: Structure) => structure.elementCount > 0,
});
