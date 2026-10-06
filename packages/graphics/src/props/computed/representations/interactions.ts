/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  Representation,
  type RepresentationParamsGetter,
  type RepresentationContext,
} from '@molstar/graphics/repr/representation';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import type { Structure } from '@molstar/model/model/structure';
import {
  UnitsRepresentation,
  type StructureRepresentation,
  StructureRepresentationStateBuilder,
  StructureRepresentationProvider,
  ComplexRepresentation,
} from '@molstar/graphics/repr/structure/representation';
import { InteractionsIntraUnitParams, InteractionsIntraUnitVisual } from './interactions-intra-unit-cylinder.js';
import { InteractionsProvider } from '@molstar/model/props/computed/interactions';
import { InteractionsInterUnitParams, InteractionsInterUnitVisual } from './interactions-inter-unit-cylinder.js';
import { BridgeParams, BridgeVisual } from './interactions-bridge-cylinder.js';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import { getUnitKindsParam } from '@molstar/graphics/repr/structure/params';

const InteractionsVisuals = {
  'intra-unit': (
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<Structure, InteractionsIntraUnitParams>,
  ) => UnitsRepresentation('Intra-unit interactions cylinder', ctx, getParams, InteractionsIntraUnitVisual),
  'inter-unit': (
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<Structure, InteractionsInterUnitParams>,
  ) => ComplexRepresentation('Inter-unit interactions cylinder', ctx, getParams, InteractionsInterUnitVisual),
  bridge: (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, BridgeParams>) =>
    ComplexRepresentation('Bridge cylinder', ctx, getParams, BridgeVisual),
};

export const InteractionsParams = {
  ...InteractionsIntraUnitParams,
  ...InteractionsInterUnitParams,
  ...BridgeParams,
  unitKinds: getUnitKindsParam(['atomic']),
  sizeFactor: PD.Numeric(0.2, { min: 0.01, max: 1, step: 0.01 }),
  visuals: PD.MultiSelect(['intra-unit', 'inter-unit', 'bridge'], PD.objectToOptions(InteractionsVisuals)),
};
export type InteractionsParams = typeof InteractionsParams;
export function getInteractionParams(ctx: ThemeRegistryContext, structure: Structure) {
  return PD.clone(InteractionsParams);
}

export type InteractionRepresentation = StructureRepresentation<InteractionsParams>;
export function InteractionRepresentation(
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<Structure, InteractionsParams>,
): InteractionRepresentation {
  return Representation.createMulti(
    'Interactions',
    ctx,
    getParams,
    StructureRepresentationStateBuilder,
    InteractionsVisuals as unknown as Representation.Def<Structure, InteractionsParams>,
  );
}

export const InteractionsRepresentationProvider = StructureRepresentationProvider({
  name: 'interactions',
  label: 'Non-covalent Interactions',
  description: 'Displays non-covalent interactions as dashed cylinders.',
  factory: InteractionRepresentation,
  getParams: getInteractionParams,
  defaultValues: PD.getDefaultValues(InteractionsParams),
  defaultColorTheme: { name: 'interaction-type' },
  defaultSizeTheme: { name: 'uniform' },
  isApplicable: (structure: Structure) => structure.elementCount > 0 && InteractionsProvider.isApplicable(structure),
  ensureCustomProperties: {
    attach: (ctx: CustomProperty.Context, structure: Structure) =>
      InteractionsProvider.attach(ctx, structure, void 0, true),
    detach: (data) => InteractionsProvider.ref(data, false),
  },
  getData: (structure: Structure, props: PD.Values<InteractionsParams>) => {
    return props.includeParent ? structure.asParent() : structure;
  },
  mustRecreate: (oldProps: PD.Values<InteractionsParams>, newProps: PD.Values<InteractionsParams>) => {
    return oldProps.includeParent !== newProps.includeParent;
  },
});
