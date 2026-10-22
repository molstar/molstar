/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { Dnatco, DnatcoParams, type DnatcoSteps } from '@molstar/dnatco-extension/property';
import { CustomPropertyDescriptor } from '@molstar/model/model/custom-property';
import { Model } from '@molstar/model/model/structure';
import { CustomProperty } from '@molstar/model/props/common/custom-property';
import { CustomModelProperty } from '@molstar/model/props/common/custom-model-property';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

export const ConfalPyramidsParams = { ...DnatcoParams };
export type ConfalPyramidsParams = typeof ConfalPyramidsParams;
export type ConfalPyramidsProps = PD.Values<ConfalPyramidsParams>;

export const ConfalPyramidsProvider: CustomModelProperty.Provider<ConfalPyramidsParams, DnatcoSteps> =
  CustomModelProperty.createProvider({
    label: 'Confal Pyramids',
    descriptor: CustomPropertyDescriptor({
      name: 'confal_pyramids',
    }),
    type: 'static',
    defaultParams: ConfalPyramidsParams,
    getParams: (data: Model) => ConfalPyramidsParams,
    isApplicable: (data: Model) => Dnatco.isApplicable(data),
    obtain: async (ctx: CustomProperty.Context, data: Model, props: Partial<ConfalPyramidsProps>) => {
      const p = { ...PD.getDefaultValues(ConfalPyramidsParams), ...props };
      return Dnatco.fromCif(ctx, data, p);
    },
  });
