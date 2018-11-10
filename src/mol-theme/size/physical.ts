/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Unit, Link, ElementIndex } from 'mol-model/structure';
import { Location } from 'mol-model/location';
import { SizeTheme } from '../size';
import { VdwRadius } from 'mol-model/structure/model/properties/atomic';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from 'mol-theme/theme';

const DefaultSize = 1
const Description = 'Assigns a physical size.'

export const PhysicalSizeThemeParams = {}
export function getPhysicalSizeThemeParams(ctx: ThemeDataContext) {
    return PhysicalSizeThemeParams // TODO return copy
}
export type PhysicalSizeThemeProps = PD.DefaultValues<typeof PhysicalSizeThemeParams>

export function getPhysicalRadius(unit: Unit, element: ElementIndex): number {
    if (Unit.isAtomic(unit)) {
        return VdwRadius(unit.model.atomicHierarchy.atoms.type_symbol.value(element))
    } else if (Unit.isSpheres(unit)) {
        return unit.model.coarseConformation.spheres.radius[element]
    } else {
        return 0
    }
}

/**
 * Create attribute data with the physical size of an element,
 * i.e. vdw for atoms and radius for coarse spheres
 */
export function PhysicalSizeTheme(ctx: ThemeDataContext, props: PhysicalSizeThemeProps): SizeTheme<PhysicalSizeThemeProps> {
    function size(location: Location): number {
        let size: number
        if (StructureElement.isLocation(location)) {
            size = getPhysicalRadius(location.unit, location.element)
        } else if (Link.isLocation(location)) {
            size = getPhysicalRadius(location.aUnit, location.aUnit.elements[location.aIndex])
        } else {
            size = DefaultSize
        }
        return size
    }

    return {
        granularity: 'group',
        size,
        props,
        description: Description
    }
}

export const PhysicalSizeThemeProvider: SizeTheme.Provider<typeof PhysicalSizeThemeParams> = {
    factory: PhysicalSizeTheme, params: getPhysicalSizeThemeParams
}