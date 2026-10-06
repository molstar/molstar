/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { CustomProperty } from '@molstar/model/props/common/custom-property';
import { CustomStructureProperty } from '@molstar/model/props/common/custom-structure-property';
import { CustomPropertyDescriptor } from '@molstar/model/model/custom-property';
import { Loci } from '@molstar/model/model/loci';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import type { LociLabelProvider } from '@molstar/plugin/state/manager/loci-label';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { filterDefined } from '@molstar/mvs/helpers/utils';
import { ElementSet, type Selector, SelectorParams } from './selector.js';

/** Parameter definition for custom structure property "CustomTooltips" */
export type CustomTooltipsParams = typeof CustomTooltipsParams;
export const CustomTooltipsParams = {
  tooltips: PD.ObjectList(
    {
      text: PD.Text('', { description: 'Text of the tooltip' }),
      selector: SelectorParams,
    },
    (obj) => obj.text,
  ),
};

/** Parameter values of custom structure property "CustomTooltips" */
export type CustomTooltipsProps = PD.Values<CustomTooltipsParams>;

/** Values of custom structure property "CustomTooltips" (and for its params at the same type) */
export type CustomTooltipsData = { selector: Selector; text: string; elementSet?: ElementSet }[];

/** Provider for custom structure property "CustomTooltips" */
export const CustomTooltipsProvider: CustomStructureProperty.Provider<CustomTooltipsParams, CustomTooltipsData> =
  CustomStructureProperty.createProvider({
    label: 'MVS Custom Tooltips',
    descriptor: CustomPropertyDescriptor<any, any>({
      name: 'mvs-custom-tooltips',
    }),
    type: 'local',
    defaultParams: CustomTooltipsParams,
    getParams: (data: Structure) => CustomTooltipsParams,
    isApplicable: (data: Structure) => data.root === data,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<CustomTooltipsProps>) => {
      const fullProps = { ...PD.getDefaultValues(CustomTooltipsParams), ...props };
      const value = fullProps.tooltips.map(
        (t) =>
          ({
            selector: t.selector,
            text: t.text,
          }) satisfies CustomTooltipsData[number],
      );
      return { value: value } satisfies CustomProperty.Data<CustomTooltipsData>;
    },
    isHidden: true,
  });

/** Label provider based on custom structure property "CustomTooltips" */
export const CustomTooltipsLabelProvider = {
  label: (loci: Loci): string | undefined => {
    switch (loci.kind) {
      case 'element-loci':
        if (!loci.structure.customPropertyDescriptors.hasReference(CustomTooltipsProvider.descriptor)) return undefined;
        const location = StructureElement.Loci.getFirstLocation(loci);
        if (!location) return undefined;
        const tooltipData = CustomTooltipsProvider.get(location.structure).value;
        if (!tooltipData || tooltipData.length === 0) return undefined;
        const texts = [];
        for (const tooltip of tooltipData) {
          const elements = (tooltip.elementSet ??= ElementSet.fromSelector(location.structure, tooltip.selector));
          if (ElementSet.has(elements, location)) texts.push(tooltip.text);
        }
        return filterDefined(texts).join(' | ');
      default:
        return undefined;
    }
  },
} satisfies LociLabelProvider;
