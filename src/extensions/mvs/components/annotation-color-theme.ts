/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { SortedArray } from '../../../mol-data/int';
import { Location } from '../../../mol-model/location';
import { Bond, StructureElement } from '../../../mol-model/structure';
import type { ColorTheme, LocationColor } from '../../../mol-theme/color';
import type { ThemeDataContext } from '../../../mol-theme/theme';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { MaybeFloatParamDefinition } from '../helpers/param-definition';
import { SplitColor } from '../helpers/utils';
import { getMVSAnnotationForStructure, MVSAnnotation } from './annotation-prop';
import { isMVSStructure } from './is-mvs-model-prop';
import { PDOptionalSplitColor, PDSplitColor, SplitColorProp } from './split-uniform-color-theme';


function fmtFloat(x: number): string {
    if (x === Infinity) return '\u221e';
    if (x === -Infinity) return '-\u221e';
    return x.toString();
}

export const MVSCategoricalPaletteParams = {
    colors: PD.MappedStatic('list', {
        list: PD.ObjectList({
            color: PDSplitColor(),
        }, e => SplitColorProp.toString(e.color), { description: 'List of colors.' }),
        dictionary: PD.ObjectList({
            value: PD.Text(),
            color: PDSplitColor(),
        }, e => `${e.value}: ${SplitColorProp.toString(e.color)}`, { description: 'Mapping of annotation values to colors.' }),
    }),
    repeatColorList: PD.Boolean(false, { hideIf: g => g.colors.name !== 'list', description: 'Repeat color list once all colors are depleted (only applies if `colors` is a list).' }),
    sort: PD.Select('none', [['none', 'None'], ['lexical', 'Lexical'], ['numeric', 'Numeric']] as const, { hideIf: g => g.colors.name !== 'list', description: 'Sort actual annotation values before assigning colors from a list (none = take values in order of their first occurrence).' }),
    sortDirection: PD.Select('ascending', [['ascending', 'Ascending'], ['descending', 'Descending']] as const, { hideIf: g => g.colors.name !== 'list', description: 'Sort direction.' }),
    caseInsensitive: PD.Boolean(false, { description: 'Treat annotation values as case-insensitive strings.' }),
    missingColor: PDOptionalSplitColor({ description: 'Color to use when (a) `colors` is a dictionary and given key is not present, or (b) `colors` is a list and there are more actual annotation values than listed colors and `repeat_color_list` is not true.' }),
};
export type MVSCategoricalPaletteParams = typeof MVSCategoricalPaletteParams
export type MVSCategoricalPaletteProps = PD.Values<MVSCategoricalPaletteParams>

export const MVSDiscretePaletteParams = {
    colors: PD.ObjectList({
        color: PDOptionalSplitColor(),
        fromValue: PD.Numeric(-Infinity),
        toValue: PD.Numeric(Infinity),
    }, e => `${SplitColorProp.toString(e.color)} [${fmtFloat(e.fromValue)}, ${fmtFloat(e.toValue)}]`, { description: 'Mapping of annotation value ranges to colors.' }),
    mode: PD.Select('normalized', [['normalized', 'Normalized'], ['absolute', 'Absolute']] as const, { description: 'Defines whether the annotation values should be normalized before assigning color based on checkpoints in `colors` (`x_normalized = (x - x_min) / (x_max - x_min)`, where `[x_min, x_max]` are either `value_domain` if provided, or the lowest and the highest value encountered in the annotation).' }),
    xMin: MaybeFloatParamDefinition({ hideIf: g => g.mode !== 'normalized', placeholder: 'auto', description: 'Defines `x_min` for normalization of annotation values. If not provided, minimum of the actual values will be used. Only used when `mode` is `"normalized"' }),
    xMax: MaybeFloatParamDefinition({ hideIf: g => g.mode !== 'normalized', placeholder: 'auto', description: 'Defines `x_max` for normalization of annotation values. If not provided, maximum of the actual values will be used. Only used when `mode` is `"normalized"' }),
};
export type MVSDiscretePaletteParams = typeof MVSDiscretePaletteParams
export type MVSDiscretePaletteProps = PD.Values<MVSDiscretePaletteParams>

export const MVSContinuousPaletteParams = {
    colors: PD.ObjectList({
        color: PDSplitColor(),
        checkpoint: PD.Numeric(0),
    }, e => `${SplitColorProp.toString(e.color)} [${fmtFloat(e.checkpoint)}]`, { description: 'List of colors with checkpoints.' }),
    mode: PD.Select('normalized', [['normalized', 'Normalized'], ['absolute', 'Absolute']] as const, { description: 'Defines whether the annotation values should be normalized before assigning color based on checkpoints in `colors` (`x_normalized = (x - x_min) / (x_max - x_min)`, where `[x_min, x_max]` are either `value_domain` if provided, or the lowest and the highest value encountered in the annotation).' }),
    xMin: MaybeFloatParamDefinition({ hideIf: g => g.mode !== 'normalized', placeholder: 'auto', description: 'Defines `x_min` for normalization of annotation values. If not provided, minimum of the actual values will be used. Only used when `mode` is `"normalized"' }),
    xMax: MaybeFloatParamDefinition({ hideIf: g => g.mode !== 'normalized', placeholder: 'auto', description: 'Defines `x_max` for normalization of annotation values. If not provided, maximum of the actual values will be used. Only used when `mode` is `"normalized"' }),
    underflowColor: PDOptionalSplitColor({ description: 'Color for values below the lowest checkpoint.' }),
    overflowColor: PDOptionalSplitColor({ description: 'Color for values above the highest checkpoint.' }),
};
export type MVSContinuousPaletteParams = typeof MVSContinuousPaletteParams
export type MVSContinuousPaletteProps = PD.Values<MVSContinuousPaletteParams>


/** Parameter definition for color theme "MVS Annotation" */
export const MVSAnnotationColorThemeParams = {
    annotationId: PD.Text('', { description: 'Reference to "Annotation" custom model property' }),
    fieldName: PD.Text('color', { description: 'Annotation field (column) from which to take color values' }),
    background: PD.Color(ColorNames.gainsboro, { description: 'Color for elements without annotation' }),
    palette: PD.MappedStatic('direct', {
        'direct': PD.EmptyGroup(),
        'categorical': PD.Group(MVSCategoricalPaletteParams),
        'discrete': PD.Group(MVSDiscretePaletteParams),
        'continuous': PD.Group(MVSContinuousPaletteParams),
    }),
};
export type MVSAnnotationColorThemeParams = typeof MVSAnnotationColorThemeParams

/** Parameter values for color theme "MVS Annotation" */
export type MVSAnnotationColorThemeProps = PD.Values<MVSAnnotationColorThemeParams>

/** Return color theme that assigns colors based on an annotation file.
 * The annotation file itself is handled by a custom model property (`MVSAnnotationsProvider`),
 * the color theme then just uses this property. */
export function MVSAnnotationColorTheme(ctx: ThemeDataContext, props: MVSAnnotationColorThemeProps): ColorTheme<MVSAnnotationColorThemeParams> {
    let color: LocationColor = () => props.background;

    if (ctx.structure && !ctx.structure.isEmpty) {
        const { annotation } = getMVSAnnotationForStructure(ctx.structure, props.annotationId);
        if (annotation) {
            const paletteFunction = makePaletteFunction(props.palette, annotation, props.fieldName);

            const colorForStructureElementLocation = (location: StructureElement.Location, isSecondary: boolean) => {
                const annotValue = annotation?.getValueForLocation(location, props.fieldName);
                const color = annotValue !== undefined ? paletteFunction(annotValue, isSecondary) : undefined;
                return color ?? props.background;
            };
            const auxLocation = StructureElement.Location.create(ctx.structure);

            color = (location: Location, isSecondary) => {
                if (StructureElement.Location.is(location)) {
                    return colorForStructureElementLocation(location, isSecondary);
                } else if (Bond.isLocation(location)) {
                    // this will be applied for each bond twice, to get color of each half (a* refers to the adjacent atom, b* to the opposite atom)
                    auxLocation.unit = location.aUnit;
                    auxLocation.element = location.aUnit.elements[location.aIndex];
                    return colorForStructureElementLocation(auxLocation, isSecondary);
                }
                return props.background;
            };
        } else {
            console.error(`Annotation source "${props.annotationId}" not present`);
        }
    }

    return {
        factory: MVSAnnotationColorTheme,
        granularity: 'groupInstance',
        preferSmoothing: false,
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

type PaletteFunction = (value: string, isSecondary: boolean) => Color | undefined;

function makePaletteFunction(props: MVSAnnotationColorThemeProps['palette'], annotation: MVSAnnotation, fieldName: string): PaletteFunction {
    if (props.name === 'direct') return paletteFunctionDirect;
    if (props.name === 'categorical') return makePaletteFunctionCategorical(props.params, annotation, fieldName);
    if (props.name === 'discrete') return makePaletteFunctionDiscrete(props.params as MVSDiscretePaletteProps, annotation, fieldName);
    if (props.name === 'continuous') return makePaletteFunctionContinuous(props.params as MVSContinuousPaletteProps, annotation, fieldName);
    throw new Error(`NotImplementedError: makePaletteFunction for ${(props as any).name}`);
}

const _colors: [Color, Color] = [ColorNames.black, ColorNames.black];

const paletteFunctionDirect: PaletteFunction = (value, isSecondary) => {
    SplitColor.decode(value, _colors);
    return _colors[isSecondary ? 1 : 0];
    // if (isSecondary) {
    //     return colors.color2 ?? colors.color1;
    // } else {
    //     return colors.color1;
    // }
};

function makePaletteFunctionCategorical(props: MVSCategoricalPaletteProps, annotation: MVSAnnotation, fieldName: string): PaletteFunction {
    const colorMap: { [value: string]: [Color, Color] | undefined } = {};
    if (props.colors.name === 'dictionary') {
        for (const { value, color } of props.colors.params) {
            const key = props.caseInsensitive ? value.toUpperCase() : value;
            colorMap[key] = SplitColorProp.toTuple(color);
        }
    } else if (props.colors.name === 'list') {
        const values = annotation.getDistinctValuesInField(fieldName, props.caseInsensitive);
        if (props.sort === 'lexical') values.sort();
        else if (props.sort === 'numeric') values.sort((a, b) => Number.parseFloat(a) - Number.parseFloat(b));
        if (props.sortDirection === 'descending') values.reverse();

        const colorList = props.colors.params.map(item => SplitColorProp.toTuple(item.color));
        let next = 0;
        for (const value of values) {
            colorMap[value] = colorList[next++];
            if (next >= colorList.length) {
                if (props.repeatColorList) next = 0;
                else break;
            }
        }
    }
    const missingColor = SplitColorProp.toTuple(props.missingColor);
    if (props.caseInsensitive) {
        return (value: string, isSecondary) => (colorMap[value.toUpperCase()] ?? missingColor)[isSecondary ? 1 : 0];
    } else {
        return (value: string, isSecondary) => (colorMap[value] ?? missingColor)[isSecondary ? 1 : 0];
    }
}

function makePaletteFunctionDiscrete(props: MVSDiscretePaletteProps, annotation: MVSAnnotation, fieldName: string): PaletteFunction {
    if (props.colors.length === 0) return () => undefined;

    const bins = props.colors.map(item => ({ ...item, color: SplitColorProp.toTuple(item.color) }));
    const scale = makeNumericPaletteScale(props, annotation, fieldName);

    return (value: string, isSecondary: boolean) => {
        const xAbs = parseFloat(value);
        if (isNaN(xAbs)) return undefined;
        const x = scale(xAbs);

        for (let i = bins.length - 1; i >= 0; i--) {
            const { color, fromValue, toValue } = bins[i];
            if (fromValue <= x && x <= toValue) return color[isSecondary ? 1 : 0];
        }
    };
}

function makePaletteFunctionContinuous(props: MVSContinuousPaletteProps, annotation: MVSAnnotation, fieldName: string): PaletteFunction {
    const { colors, checkpoints } = makeContinuousPaletteCheckpoints(props);
    if (colors.length === 0) return () => undefined;

    const scale = makeNumericPaletteScale(props, annotation, fieldName);
    const underflowColor = SplitColorProp.toTuple(props.underflowColor);
    const overflowColor = SplitColorProp.toTuple(props.overflowColor);

    return (value: string, isSecondary: boolean) => {
        const secFlag = isSecondary ? 1 : 0;
        const xAbs = parseFloat(value);
        if (isNaN(xAbs)) return undefined;
        const x = scale(xAbs);
        const gteIdx = SortedArray.findPredecessorIndex(checkpoints, x); // Index of the first greater or equal checkpoint
        if (gteIdx === 0) {
            if (x === checkpoints[0]) return colors[0][secFlag];
            else return underflowColor[secFlag];
        }
        if (gteIdx === checkpoints.length) {
            return overflowColor[secFlag];
        }
        const q = (x - checkpoints[gteIdx - 1]) / (checkpoints[gteIdx] - checkpoints[gteIdx - 1]);
        return Color.interpolate(colors[gteIdx - 1][secFlag], colors[gteIdx][secFlag], q);
    };
}

function makeNumericPaletteScale(props: MVSContinuousPaletteProps | MVSDiscretePaletteProps, annotation: MVSAnnotation, fieldName: string): (x: number) => number {
    if (props.mode === 'normalized') {
        // Mode normalized
        let xMin = props.xMin;
        let xMax = props.xMax;
        if (xMin === null || xMax === null) {
            const values = annotation.getDistinctValuesInField(fieldName, false).map(parseFloat).filter(x => !isNaN(x));
            if (values.length > 0) {
                xMin ??= values.reduce((a, b) => a < b ? a : b); // xMin ??= min(values)
                xMax ??= values.reduce((a, b) => a > b ? a : b); // xMax ??= max(values)
            } else {
                xMin ??= 0;
                xMax ??= 1;
            }
        }
        if (xMin === xMax) {
            return x => (x < xMin ? -0.5 : x === xMin ? 0.5 : 1.5);
        } else {
            return x => (x - xMin) / (xMax - xMin);
        }
    } else {
        // Mode absolute
        return x => x;
    }
}

export function makeContinuousPaletteCheckpoints(props: MVSContinuousPaletteProps) {
    const sorted = props.colors.sort((a, b) => a.checkpoint - b.checkpoint);
    const colors = sorted.map(t => SplitColorProp.toTuple(t.color));
    const checkpoints = SortedArray.ofSortedArray(sorted.map(t => t.checkpoint));
    return { colors, checkpoints };
}
