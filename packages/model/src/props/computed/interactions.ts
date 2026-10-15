/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Structure } from '@molstar/model/model/structure';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { computeInteractions, type Interactions } from './interactions/interactions.js';
import { InteractionsParams as _InteractionsParams } from './interactions/params.js';
import { CustomStructureProperty } from '../common/custom-structure-property.js';
import type { CustomProperty } from '../common/custom-property.js';
import { CustomPropertyDescriptor } from '@molstar/model/model/custom-property';

export const InteractionsParams = {
  ..._InteractionsParams,
};
export type InteractionsParams = typeof InteractionsParams;
export type InteractionsProps = PD.Values<InteractionsParams>;

export type InteractionsValue = Interactions;

export const InteractionsProvider: CustomStructureProperty.Provider<InteractionsParams, InteractionsValue> =
  CustomStructureProperty.createProvider({
    label: 'Interactions',
    descriptor: CustomPropertyDescriptor({
      name: 'molstar_computed_interactions',
      // TODO `cifExport` and `symbol`
    }),
    type: 'local',
    defaultParams: InteractionsParams,
    getParams: (data: Structure) => InteractionsParams,
    isApplicable: (data: Structure) => !data.isCoarseGrained,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<InteractionsProps>) => {
      const p = { ...PD.getDefaultValues(InteractionsParams), ...props };
      return { value: await computeInteractions(ctx, data, p) };
    },
  });
