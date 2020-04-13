/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Unit, Bond, ElementIndex } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { SizeTheme } from '../size';
import { VdwRadius } from '../../mol-model/structure/model/properties/atomic';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';

const DefaultSize = 1;
const Description = 'Assigns a physical size, i.e. vdW radius for atoms or given radius for coarse spheres.';

export const PhysicalSizeThemeParams = {};
export type PhysicalSizeThemeParams = typeof PhysicalSizeThemeParams
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
export function PhysicalSizeTheme(ctx: ThemeDataContext, props: PD.Values<PhysicalSizeThemeParams>): SizeTheme<PhysicalSizeThemeParams> {
    function size(location: Location): number {
        let size: number;
        if (StructureElement.Location.is(location)) {
            size = getPhysicalRadius(location.unit, location.element);
        } else if (Bond.isLocation(location)) {
            size = getPhysicalRadius(location.aUnit, location.aUnit.elements[location.aIndex]);
        } else {
            size = DefaultSize;
        }
        return size;
    }

    return {
        factory: PhysicalSizeTheme,
        granularity: 'group',
        size,
        props,
        description: Description
    };
}

export const PhysicalSizeThemeProvider: SizeTheme.Provider<PhysicalSizeThemeParams, 'physical'> = {
    name: 'physical',
    label: 'Physical',
    category: '',
    factory: PhysicalSizeTheme,
    getParams: getPhysicalSizeThemeParams,
    defaultValues: PD.getDefaultValues(PhysicalSizeThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};