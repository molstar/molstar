/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { CustomStructureProperty } from '../../../mol-model-props/common/custom-structure-property';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';
import { Loci } from '../../../mol-model/loci';
import { Structure, StructureElement } from '../../../mol-model/structure';
import { LociLabelProvider } from '../../../mol-plugin-state/manager/loci-label';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { filterDefined } from '../helpers/utils';
import { ElementSet, Selector, SelectorParams } from './selector';


/** Parameter definition for custom structure property "CustomTooltips" */
export type CustomTooltipsParams = typeof CustomTooltipsParams
export const CustomTooltipsParams = {
    tooltips: PD.ObjectList(
        {
            text: PD.Text('', { description: 'Text of the tooltip' }),
            selector: SelectorParams,
        },
        obj => obj.text
    ),
};

/** Parameter values of custom structure property "CustomTooltips" */
export type CustomTooltipsProps = PD.Values<CustomTooltipsParams>

/** Values of custom structure property "CustomTooltips" (and for its params at the same type) */
export type CustomTooltipsData = { selector: Selector, text: string, elementSet?: ElementSet }[]


/** Provider for custom structure property "CustomTooltips" */
export const CustomTooltipsProvider: CustomStructureProperty.Provider<CustomTooltipsParams, CustomTooltipsData> = CustomStructureProperty.createProvider({
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
        const value = fullProps.tooltips.map(t => ({
            selector: t.selector,
            text: t.text,
        } satisfies CustomTooltipsData[number]));
        return { value: value } satisfies CustomProperty.Data<CustomTooltipsData>;
    },
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
                    const elements = tooltip.elementSet ??= ElementSet.fromSelector(location.structure, tooltip.selector);
                    if (ElementSet.has(elements, location)) texts.push(tooltip.text);
                }
                return filterDefined(texts).join(' | ');
            default:
                return undefined;
        }
    }
} satisfies LociLabelProvider;
