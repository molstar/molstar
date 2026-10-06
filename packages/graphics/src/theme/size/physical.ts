/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Unit, Bond, type ElementIndex } from '@molstar/model/model/structure';
import type { Location } from '@molstar/model/model/location';
import type { SizeTheme } from '../size.js';
import { VdwRadius } from '@molstar/model/model/structure/model/properties/atomic';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';

const DefaultSize = 1;
const Description = 'Assigns a physical size, i.e. vdW radius for atoms or given radius for coarse spheres.';

export const PhysicalSizeThemeParams = {
  scale: PD.Numeric(1, { min: 0.1, max: 5, step: 0.1 }),
};
export type PhysicalSizeThemeParams = typeof PhysicalSizeThemeParams;
export function getPhysicalSizeThemeParams(ctx: ThemeDataContext) {
  return PhysicalSizeThemeParams; // TODO return copy
}

export function getPhysicalRadius(unit: Unit, element: ElementIndex): number {
  if (Unit.isAtomic(unit)) {
    return VdwRadius(unit.model.atomicHierarchy.atoms.type_symbol.value(element));
  } else if (Unit.isSpheres(unit)) {
    return unit.model.coarseConformation.spheres.radius[element];
  } else {
    return 0;
  }
}

/**
 * Create attribute data with the physical size of an element,
 * i.e. vdw for atoms and radius for coarse spheres
 */
export function PhysicalSizeTheme(
  ctx: ThemeDataContext,
  props: PD.Values<PhysicalSizeThemeParams>,
): SizeTheme<PhysicalSizeThemeParams> {
  const scale = props.scale === void 0 ? 1 : props.scale;

  function size(location: Location): number {
    let size: number;
    if (StructureElement.Location.is(location)) {
      size = scale * getPhysicalRadius(location.unit, location.element);
    } else if (Bond.isLocation(location)) {
      size =
        scale *
        Math.min(
          getPhysicalRadius(location.aUnit, location.aUnit.elements[location.aIndex]),
          getPhysicalRadius(location.bUnit, location.bUnit.elements[location.bIndex]),
        );
    } else {
      size = scale * DefaultSize;
    }
    return size;
  }

  return {
    factory: PhysicalSizeTheme,
    granularity: 'group',
    size,
    props,
    description: Description,
  };
}

export const PhysicalSizeThemeProvider: SizeTheme.Provider<PhysicalSizeThemeParams, 'physical'> = {
  name: 'physical',
  label: 'Physical',
  category: '',
  factory: PhysicalSizeTheme,
  getParams: getPhysicalSizeThemeParams,
  defaultValues: PD.getDefaultValues(PhysicalSizeThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
};
