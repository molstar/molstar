/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '@molstar/core/util/color';
import type { Location } from '@molstar/model/model/location';
import { StructureElement, Bond } from '@molstar/model/model/structure';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import type { ColorTheme, LocationColor } from '../color.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import { getPaletteParams, getPalette } from '@molstar/core/util/color/palette';
import type { TableLegend, ScaleLegend } from '@molstar/core/util/legend';
import { ColorThemeCategory } from './categories.js';

const DefaultColor = Color(0xcccccc);
const Description =
  'Gives every element (atom or coarse sphere/gaussian) a unique color based on the position (index) of the element in the list of elements in the structure.';

export const ElementIndexColorThemeParams = {
  ...getPaletteParams({ type: 'colors', colorList: 'red-yellow-blue' }),
};
export type ElementIndexColorThemeParams = typeof ElementIndexColorThemeParams;
export function getElementIndexColorThemeParams(ctx: ThemeDataContext) {
  return ElementIndexColorThemeParams; // TODO return copy
}

export function ElementIndexColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<ElementIndexColorThemeParams>,
): ColorTheme<ElementIndexColorThemeParams> {
  let color: LocationColor;
  let legend: ScaleLegend | TableLegend | undefined;

  if (ctx.structure) {
    const { units } = ctx.structure.root;
    const unitCount = units.length;
    const cummulativeElementCount = new Int32Array(unitCount);
    const unitIdIndex = new Map<number, number>();
    // per unit its first element where its elements are consecutive (a chain mostly), else -1: the element's index in
    // the unit by subtraction, not by a search per location (each atom of each copy of an assembly)
    const consecutiveStart = new Int32Array(unitCount);

    let elementCount = 0;
    for (let i = 0; i < unitCount; ++i) {
      const { elements } = units[i];
      cummulativeElementCount[i] = elementCount;
      elementCount += elements.length;
      unitIdIndex.set(units[i].id, i);
      const n = elements.length;
      const first = n > 0 ? elements[0] : 0;
      consecutiveStart[i] = n > 0 && elements[n - 1] - first + 1 === n ? first : -1;
    }

    const palette = getPalette(elementCount, props);
    legend = palette.legend;

    // the last unit's index: consecutive locations are mostly of one unit
    let lastUnitId = -1;
    let lastUnitIndex = -1;
    const getUnitIndex = (unitId: number) => {
      if (unitId !== lastUnitId) {
        lastUnitId = unitId;
        lastUnitIndex = unitIdIndex.get(unitId)!;
      }
      return lastUnitIndex;
    };
    const getUnitElementIndex = (unitIndex: number, element: ElementIndex) => {
      const start = consecutiveStart[unitIndex];
      return start >= 0 ? element - start : OrderedSet.findPredecessorIndex(units[unitIndex].elements, element);
    };

    color = (location: Location): Color => {
      if (StructureElement.Location.is(location)) {
        const unitIndex = getUnitIndex(location.unit.id);
        return palette.color(cummulativeElementCount[unitIndex] + getUnitElementIndex(unitIndex, location.element));
      } else if (Bond.isLocation(location)) {
        const unitIndex = getUnitIndex(location.aUnit.id);
        const element = location.aUnit.elements[location.aIndex];
        return palette.color(cummulativeElementCount[unitIndex] + getUnitElementIndex(unitIndex, element));
      }
      return DefaultColor;
    };
  } else {
    color = () => DefaultColor;
  }

  return {
    factory: ElementIndexColorTheme,
    granularity: 'groupInstance',
    preferSmoothing: true,
    color,
    props,
    description: Description,
    legend,
  };
}

export const ElementIndexColorThemeProvider: ColorTheme.Provider<ElementIndexColorThemeParams, 'element-index'> = {
  name: 'element-index',
  label: 'Element Index',
  category: ColorThemeCategory.Atom,
  factory: ElementIndexColorTheme,
  getParams: getElementIndexColorThemeParams,
  defaultValues: PD.getDefaultValues(ElementIndexColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
};
