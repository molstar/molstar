/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Structure } from '@molstar/model/model/structure';
import {
  Representation,
  type RepresentationContext,
  type RepresentationParamsGetter,
} from '@molstar/graphics/repr/representation';
import {
  ComplexRepresentation,
  type StructureRepresentation,
  StructureRepresentationProvider,
  StructureRepresentationStateBuilder,
} from '@molstar/graphics/repr/structure/representation';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { isMVSStructure } from '../is-mvs-model-prop.js';
import { CustomLabelTextParams, CustomLabelTextVisual } from './visual.js';

/** Components of "Custom Label" representation */
const CustomLabelVisuals = {
  'label-text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CustomLabelTextParams>) =>
    ComplexRepresentation('Label text', ctx, getParams, CustomLabelTextVisual),
};

/** Parameter definition for representation type "Custom Label" */
export type CustomLabelParams = typeof CustomLabelParams;
export const CustomLabelParams = {
  ...CustomLabelTextParams,
  visuals: PD.MultiSelect(['label-text'], PD.objectToOptions(CustomLabelVisuals)),
};

/** Parameter values for representation type "Custom Label" */
export type CustomLabelProps = PD.ValuesFor<CustomLabelParams>;

/** Structure representation type "Custom Label", allowing user-defined labels at at user-defined positions */
export type CustomLabelRepresentation = StructureRepresentation<CustomLabelParams>;
export function CustomLabelRepresentation(
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<Structure, CustomLabelParams>,
): CustomLabelRepresentation {
  const repr = Representation.createMulti(
    'Label',
    ctx,
    getParams,
    StructureRepresentationStateBuilder,
    CustomLabelVisuals as unknown as Representation.Def<Structure, CustomLabelParams>,
  );
  repr.setState({ pickable: false, markerActions: MarkerAction.None });
  return repr;
}

/** A thingy that is needed to register representation type "Custom Label", allowing user-defined labels at at user-defined positions */
export const CustomLabelRepresentationProvider = StructureRepresentationProvider({
  name: 'mvs-custom-label',
  label: 'MVS Custom Label',
  description: 'Displays labels with custom text',
  factory: CustomLabelRepresentation,
  getParams: () => CustomLabelParams,
  defaultValues: PD.getDefaultValues(CustomLabelParams),
  defaultColorTheme: { name: 'uniform' },
  defaultSizeTheme: { name: 'physical' },
  isApplicable: (structure: Structure) => structure.elementCount > 0 && isMVSStructure(structure),
});
