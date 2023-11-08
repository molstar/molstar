/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { AnnotationsProvider } from './annotation-prop';


/** Parameter definition for custom structure property "AnnotationTooltips" */
export const AnnotationTooltipsParams = {
    tooltips: PD.ObjectList(
        {
            annotationId: PD.Text('', { description: 'Reference to "Annotation" custom model property' }),
            fieldName: PD.Text('tooltip', { description: 'Annotation field (column) from which to take color values' }),
        },
        obj => `${obj.annotationId}:${obj.fieldName}`
    ),
};
export type AnnotationTooltipsParams = typeof AnnotationTooltipsParams

/** Values of custom structure property "AnnotationTooltips" (and for its params at the same type) */
export type AnnotationTooltipsProps = PD.Values<AnnotationTooltipsParams>


/** Provider for custom structure property "AnnotationTooltips" */
export const AnnotationTooltipsProvider: CustomStructureProperty.Provider<AnnotationTooltipsParams, AnnotationTooltipsProps> = CustomStructureProperty.createProvider({
    label: 'Annotation Tooltips',
    descriptor: CustomPropertyDescriptor<any, any>({
        name: 'mvs-annotation-tooltips',
    }),
    type: 'local',
    defaultParams: AnnotationTooltipsParams,
    getParams: (data: Structure) => AnnotationTooltipsParams,
    isApplicable: (data: Structure) => data.root === data,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<AnnotationTooltipsProps>) => {
        const fullProps = { ...PD.getDefaultValues(AnnotationTooltipsParams), ...props };
        return { value: fullProps } satisfies CustomProperty.Data<AnnotationTooltipsProps>;
    },
});


/** Label provider based on data from "Annotation" custom model property */
export const AnnotationTooltipsLabelProvider = {
    label: (loci: Loci): string | undefined => {
        switch (loci.kind) {
            case 'element-loci':
                if (!loci.structure.customPropertyDescriptors.hasReference(AnnotationTooltipsProvider.descriptor)) return undefined;
                const location = StructureElement.Loci.getFirstLocation(loci);
                if (!location) return undefined;
                const tooltipProps = AnnotationTooltipsProvider.get(location.structure).value;
                if (!tooltipProps || tooltipProps.tooltips.length === 0) return undefined;
                const annotations = AnnotationsProvider.get(location.unit.model).value;
                const texts = tooltipProps.tooltips.map(p => annotations?.getAnnotation(p.annotationId)?.getValueForLocation(location, p.fieldName));
                return filterDefined(texts).join(' | ');
            default:
                return undefined;
        }
    }
} satisfies LociLabelProvider;
