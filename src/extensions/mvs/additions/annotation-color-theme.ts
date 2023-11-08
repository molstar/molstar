/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Location } from '../../../mol-model/location';
import { Bond, StructureElement } from '../../../mol-model/structure';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { ColorNames } from '../../../mol-util/color/names';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { decodeColor } from '../helpers/utils';
import { getAnnotationForStructure } from './annotation-prop';


/** Parameter definition for color theme "Annotation" */
export const AnnotationColorThemeParams = {
    annotationId: PD.Text('', { description: 'Reference to "Annotation" custom model property' }),
    fieldName: PD.Text('color', { description: 'Annotation field (column) from which to take color values' }),
    background: PD.Color(ColorNames.gainsboro, { description: 'Color for elements without annotation' }),
};
export type AnnotationColorThemeParams = typeof AnnotationColorThemeParams

/** Parameter values for color theme "Annotation" */
export type AnnotationColorThemeProps = PD.Values<AnnotationColorThemeParams>


/** Return color theme that assigns colors based on an annotation file.
 * The annotation file itself is handled by a custom model property (`AnnotationsProvider`),
 * the color theme then just uses this property. */
export function AnnotationColorTheme(ctx: ThemeDataContext, props: AnnotationColorThemeProps): ColorTheme<AnnotationColorThemeParams> {
    let color: LocationColor = () => props.background;

    if (ctx.structure && !ctx.structure.isEmpty) {
        const { annotation } = getAnnotationForStructure(ctx.structure, props.annotationId);
        if (annotation) {
            const colorForStructureElementLocation = (location: StructureElement.Location) => {
                // if (annot.getAnnotationForLocation(location)?.color !== annot.getAnnotationForLocation_Reference(location)?.color) throw new Error('AssertionError');
                return decodeColor(annotation?.getValueForLocation(location, props.fieldName)) ?? props.background;
            };
            const auxLocation = StructureElement.Location.create(ctx.structure);

            color = (location: Location) => {
                if (StructureElement.Location.is(location)) {
                    return colorForStructureElementLocation(location);
                } else if (Bond.isLocation(location)) {
                    // this will be applied for each bond twice, to get color of each half (a* refers to the adjacent atom, b* to the opposite atom)
                    auxLocation.unit = location.aUnit;
                    auxLocation.element = location.aUnit.elements[location.aIndex];
                    return colorForStructureElementLocation(auxLocation);
                }
                return props.background;
            };
        } else {
            console.error(`Annotation source "${props.annotationId}" not present`);
        }
    }

    return {
        factory: AnnotationColorTheme,
        granularity: 'group',
        preferSmoothing: true,
        color: color,
        props: props,
        description: 'Assigns colors based on custom annotation data.',
    };
}


/** A thingy that is needed to register color theme "Annotation" */
export const AnnotationColorThemeProvider: ColorTheme.Provider<AnnotationColorThemeParams, 'mvs-annotation'> = {
    name: 'mvs-annotation',
    label: 'Annotation',
    category: ColorTheme.Category.Misc,
    factory: AnnotationColorTheme,
    getParams: ctx => AnnotationColorThemeParams,
    defaultValues: PD.getDefaultValues(AnnotationColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => true,
};
