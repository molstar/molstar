/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorScale } from '@molstar/core/util/color';
import { StructureElement, Unit, Bond, type ElementIndex } from '@molstar/model/model/structure';
import type { Location } from '@molstar/model/model/location';
import type { ColorTheme } from '../color.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import { ColorThemeCategory } from './categories.js';

const DefaultOccupancyColor = Color(0xcccccc);
const Description = `Assigns a color based on the occupancy of an atom.`;

export const OccupancyColorThemeParams = {
  domain: PD.Interval([0, 1]),
  list: PD.ColorList('purples', { presetKind: 'scale' }),
};
export type OccupancyColorThemeParams = typeof OccupancyColorThemeParams;
export function getOccupancyColorThemeParams(ctx: ThemeDataContext) {
  return OccupancyColorThemeParams; // TODO return copy
}

export function getOccupancy(unit: Unit, element: ElementIndex): number {
  if (Unit.isAtomic(unit)) {
    return unit.model.atomicConformation.occupancy.value(element);
  } else {
    return 0;
  }
}

export function OccupancyColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<OccupancyColorThemeParams>,
): ColorTheme<OccupancyColorThemeParams> {
  const scale = ColorScale.create({
    reverse: false,
    domain: props.domain,
    listOrName: props.list.colors,
  });

  function color(location: Location): Color {
    if (StructureElement.Location.is(location)) {
      return scale.color(getOccupancy(location.unit, location.element));
    } else if (Bond.isLocation(location)) {
      return scale.color(getOccupancy(location.aUnit, location.aUnit.elements[location.aIndex]));
    }
    return DefaultOccupancyColor;
  }

  return {
    factory: OccupancyColorTheme,
    granularity: 'group',
    preferSmoothing: true,
    color,
    props,
    description: Description,
    legend: scale ? scale.legend : undefined,
  };
}

export const OccupancyColorThemeProvider: ColorTheme.Provider<OccupancyColorThemeParams, 'occupancy'> = {
  name: 'occupancy',
  label: 'Occupancy',
  category: ColorThemeCategory.Atom,
  factory: OccupancyColorTheme,
  getParams: getOccupancyColorThemeParams,
  defaultValues: PD.getDefaultValues(OccupancyColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) =>
    !!ctx.structure && ctx.structure.models.some((m) => m.atomicConformation.occupancy.isDefined),
};
