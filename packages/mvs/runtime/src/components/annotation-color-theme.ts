/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Location } from '@molstar/model/model/location';
import { Bond, StructureElement } from '@molstar/model/model/structure';
import type { ColorTheme, LocationColor } from '@molstar/graphics/theme/color';
import type { ThemeDataContext } from '@molstar/graphics/theme/theme';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { makeValuePaletteFunction } from '@molstar/core/util/color/value-palette';
import type {
  CategoricalValuePaletteProps,
  ContinuousValuePaletteProps,
  DiscreteValuePaletteProps,
} from '@molstar/core/util/color/value-palette';
import {
  CategoricalValuePaletteParams,
  ContinuousValuePaletteParams,
  DiscreteValuePaletteParams,
  ValuePaletteParam,
} from '@molstar/core/util/color/value-palette-params';
import { getMVSAnnotationForStructure } from './annotation-prop.js';
import { isMVSStructure } from './is-mvs-model-prop.js';

export const MVSCategoricalPaletteParams = CategoricalValuePaletteParams;
export type MVSCategoricalPaletteParams = CategoricalValuePaletteParams;
export type MVSCategoricalPaletteProps = CategoricalValuePaletteProps;

export const MVSDiscretePaletteParams = DiscreteValuePaletteParams;
export type MVSDiscretePaletteParams = DiscreteValuePaletteParams;
export type MVSDiscretePaletteProps = DiscreteValuePaletteProps;

export const MVSContinuousPaletteParams = ContinuousValuePaletteParams;
export type MVSContinuousPaletteParams = ContinuousValuePaletteParams;
export type MVSContinuousPaletteProps = ContinuousValuePaletteProps;

/** Parameter definition for color theme "MVS Annotation" */
export const MVSAnnotationColorThemeParams = {
  annotationId: PD.Text('', { description: 'Reference to "Annotation" custom model property' }),
  fieldName: PD.Text('color', { description: 'Annotation field (column) from which to take color values' }),
  background: PD.Color(ColorNames.gainsboro, { description: 'Color for elements without annotation' }),
  palette: ValuePaletteParam(),
};
export type MVSAnnotationColorThemeParams = typeof MVSAnnotationColorThemeParams;

/** Parameter values for color theme "MVS Annotation" */
export type MVSAnnotationColorThemeProps = PD.Values<MVSAnnotationColorThemeParams>;

/** Return color theme that assigns colors based on an annotation file.
 * The annotation file itself is handled by a custom model property (`MVSAnnotationsProvider`),
 * the color theme then just uses this property. */
export function MVSAnnotationColorTheme(
  ctx: ThemeDataContext,
  props: MVSAnnotationColorThemeProps,
): ColorTheme<MVSAnnotationColorThemeParams> {
  let color: LocationColor = () => props.background;

  if (ctx.structure && !ctx.structure.isEmpty) {
    const { annotation } = getMVSAnnotationForStructure(ctx.structure, props.annotationId);
    if (annotation) {
      // The field's values only when the palette needs them (a list, a normalized domain without both ends).
      const paletteFunction = makeValuePaletteFunction(props.palette, () =>
        annotation.getDistinctValuesInField(props.fieldName, false),
      );

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
  getParams: (ctx) => MVSAnnotationColorThemeParams,
  defaultValues: PD.getDefaultValues(MVSAnnotationColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && isMVSStructure(ctx.structure),
};
