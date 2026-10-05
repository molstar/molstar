/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorScale } from '@molstar/core/util/color';
import { StructureElement, type Unit, Bond, type ElementIndex } from '@molstar/model/model/structure';
import type { Location } from '@molstar/model/model/location';
import type { ColorTheme } from '../color.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import { AtomPartialCharge } from '@molstar/model/formats/structure/property/partial-charge';
import { ColorThemeCategory } from './categories.js';

const DefaultPartialChargeColor = Color(0xffff99);
const Description = `Assigns a color based on the partial charge of an atom.`;

export const PartialChargeColorThemeParams = {
  domain: PD.Interval([-1, 1]),
  list: PD.ColorList('red-white-blue', { presetKind: 'scale' }),
};
export type PartialChargeColorThemeParams = typeof PartialChargeColorThemeParams;
export function getPartialChargeColorThemeParams(ctx: ThemeDataContext) {
  return PartialChargeColorThemeParams; // TODO return copy
}

function getPartialCharge(unit: Unit, element: ElementIndex) {
  return AtomPartialCharge.Provider.get(unit.model)?.data.value(element);
}

export function PartialChargeColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<PartialChargeColorThemeParams>,
): ColorTheme<PartialChargeColorThemeParams> {
  const scale = ColorScale.create({
    domain: props.domain,
    listOrName: props.list.colors,
  });

  function color(location: Location): Color {
    if (StructureElement.Location.is(location)) {
      const q = getPartialCharge(location.unit, location.element);
      return q !== undefined ? scale.color(q) : DefaultPartialChargeColor;
    } else if (Bond.isLocation(location)) {
      const q = getPartialCharge(location.aUnit, location.aUnit.elements[location.aIndex]);
      return q !== undefined ? scale.color(q) : DefaultPartialChargeColor;
    }
    return DefaultPartialChargeColor;
  }

  return {
    factory: PartialChargeColorTheme,
    granularity: 'group',
    preferSmoothing: true,
    color,
    props,
    description: Description,
    legend: scale ? scale.legend : undefined,
  };
}

export const PartialChargeColorThemeProvider: ColorTheme.Provider<PartialChargeColorThemeParams, 'partial-charge'> = {
  name: 'partial-charge',
  label: 'Partial Charge',
  category: ColorThemeCategory.Atom,
  factory: PartialChargeColorTheme,
  getParams: getPartialChargeColorThemeParams,
  defaultValues: PD.getDefaultValues(PartialChargeColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) =>
    !!ctx.structure && ctx.structure.models.some((m) => AtomPartialCharge.Provider.get(m) !== undefined),
};
