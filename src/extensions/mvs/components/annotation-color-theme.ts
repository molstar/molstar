/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Location } from '../../../mol-model/location';
import { Bond, StructureElement } from '../../../mol-model/structure';
import type { ColorTheme, LocationColor } from '../../../mol-theme/color';
import type { ThemeDataContext } from '../../../mol-theme/theme';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { decodeColor } from '../helpers/utils';
import { getMVSAnnotationForStructure, MVSAnnotation } from './annotation-prop';
import { isMVSStructure } from './is-mvs-model-prop';


export const MVSCategoricalPaletteParams = {
    colors: PD.MappedStatic('list', {
        list: PD.ColorList('category-10', { description: 'TODO', presetKind: 'set' }),
        mapping: PD.ObjectList({ // TODO consider having PD.ColorMapping for this, with UI and presets, like PD.ColorList
            value: PD.Text(),
            color: PD.Color(ColorNames.white),
        }, e => `${e.value}: ${Color.toHexStyle(e.color)}`, { description: 'TODO' }),
    }),
    repeatColorList: PD.Boolean(false, { hideIf: g => g.colors.name !== 'list', description: 'Repeat color list once all colors are depleted (only applies if `color` is a list or a named palette).' }),
    sort: PD.Select('none', [['none', 'None'], ['lexical', 'Lexical'], ['numeric', 'Numeric']] as const, { hideIf: g => g.colors.name !== 'list', description: 'Sort real values before assigning colors from a list or named palette.' }),
    /** Sort direction (only applies if `sort` is provided). */
    sortDirection: PD.Select('ascending', [['ascending', 'Ascending'], ['descending', 'Descending']] as const, { hideIf: g => g.colors.name !== 'list' || g.sort === 'none', description: 'Sort direction (only applies if `sort` is provided)' }),
    setMissingColor: PD.Boolean(false, { description: 'Allow setting a color for missing values.' }),
    missingColor: PD.Color(ColorNames.white, { hideIf: g => !g.setMissingColor, description: 'Color to use when (a) `color` is a dictionary and given key is not present, or (b) `color` is a list or a named palette and there are more real values than listed values and `repeat_color_list` is not true.' }),
    // ... TODO adjust descriptions
};
export type MVSCategoricalPaletteParams = typeof MVSCategoricalPaletteParams
export type MVSCategoricalPaletteProps = PD.Values<MVSCategoricalPaletteParams>


/** Parameter definition for color theme "MVS Annotation" */
export const MVSAnnotationColorThemeParams = {
    annotationId: PD.Text('', { description: 'Reference to "Annotation" custom model property' }),
    fieldName: PD.Text('color', { description: 'Annotation field (column) from which to take color values' }),
    background: PD.Color(ColorNames.gainsboro, { description: 'Color for elements without annotation' }),
    palette: PD.MappedStatic('direct', {
        'direct': PD.EmptyGroup(),
        'categorical': PD.Group(MVSCategoricalPaletteParams),
    }),
};
export type MVSAnnotationColorThemeParams = typeof MVSAnnotationColorThemeParams

/** Parameter values for color theme "MVS Annotation" */
export type MVSAnnotationColorThemeProps = PD.Values<MVSAnnotationColorThemeParams>

type x = Extract<MVSCategoricalPaletteProps['colors'], { name: 'list' }>['params']

/** Return color theme that assigns colors based on an annotation file.
 * The annotation file itself is handled by a custom model property (`MVSAnnotationsProvider`),
 * the color theme then just uses this property. */
export function MVSAnnotationColorTheme(ctx: ThemeDataContext, props: MVSAnnotationColorThemeProps): ColorTheme<MVSAnnotationColorThemeParams> {
    let color: LocationColor = () => props.background;

    if (ctx.structure && !ctx.structure.isEmpty) {
        const { annotation } = getMVSAnnotationForStructure(ctx.structure, props.annotationId);
        if (annotation) {
            const paletteFunction = makePaletteFunction(props.palette, annotation, props.fieldName);

            const colorForStructureElementLocation = (location: StructureElement.Location) => {
                // if (annot.getAnnotationForLocation(location)?.color !== annot.getAnnotationForLocation_Reference(location)?.color) throw new Error('AssertionError');
                const annotValue = annotation?.getValueForLocation(location, props.fieldName);
                const color = annotValue !== undefined ? paletteFunction(annotValue) : undefined;
                return color ?? props.background;
                // return decodeColor(annotation?.getValueForLocation(location, props.fieldName)) ?? props.background;
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
        factory: MVSAnnotationColorTheme,
        granularity: 'group',
        preferSmoothing: true,
        color: color,
        props: props,
        description: 'Assigns colors based on custom MolViewSpec annotation data.',
    };
}


/** A thingy that is needed to register color theme "MVS Annotation" */
export const MVSAnnotationColorThemeProvider: ColorTheme.Provider<MVSAnnotationColorThemeParams, 'mvs-annotation'> = {
    name: 'mvs-annotation',
    label: 'MVS Annotation',
    category: 'Miscellaneous', // ColorTheme.Category.Misc can cause webpack build error due to import ordering
    factory: MVSAnnotationColorTheme,
    getParams: ctx => MVSAnnotationColorThemeParams,
    defaultValues: PD.getDefaultValues(MVSAnnotationColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && isMVSStructure(ctx.structure),
};


function makePaletteFunction(props: MVSAnnotationColorThemeProps['palette'], annotation: MVSAnnotation, fieldName: string): (v: string) => Color | undefined {
    if (props.name === 'direct') return decodeColor;

    if (props.name === 'categorical') {
        const colorMap: { [value: string]: Color } = {};
        if (props.params.colors.name === 'mapping') {
            for (const { value, color } of props.params.colors.params) {
                colorMap[value] = color;
            }
        } else if (props.params.colors.name === 'list') {
            const values = annotation.getDistinctValuesInField(fieldName);
            if (props.params.sort === 'lexical') values.sort();
            else if (props.params.sort === 'numeric') values.sort((a, b) => Number.parseFloat(a) - Number.parseFloat(b));
            if (props.params.sortDirection === 'descending') values.reverse();

            const colorList = props.params.colors.params.colors.map(entry => typeof entry === 'number' ? entry : entry[0]);
            let next = 0;
            for (const value of values) {
                colorMap[value] = colorList[next++];
                if (next >= colorList.length && props.params.repeatColorList) next = 0; // else will get index-out-of-range and assign undefined
            }
        }
        const missingColor = props.params.setMissingColor ? props.params.missingColor : undefined;
        console.log('makePaletteFunction categorical', props, missingColor, colorMap);
        return (value: string) => colorMap[value] ?? missingColor;
    }

    throw new Error(`NotImplementedError: makePaletteFunction for ${(props as any).name}`);
}
