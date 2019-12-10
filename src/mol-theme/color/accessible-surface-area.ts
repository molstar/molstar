/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ThemeDataContext } from '../theme';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { ColorScale, Color } from '../../mol-util/color';
import { Unit, StructureElement } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { ComputedAccessibleSurfaceArea } from '../../mol-model-props/computed/accessible-surface-area';
import { ColorListName, ColorListOptionsScale } from '../../mol-util/color/lists';
import { VdWLookup } from '../../mol-model-props/computed/accessible-surface-area/shrake-rupley/common';

const DefaultColor = Color(0xFFFFFF)
const Description = 'Assigns a color based on the relative accessible surface area of a residue.'

export const AccessibleSurfaceAreaColorThemeParams = {
    list: PD.ColorList<ColorListName>('rainbow', ColorListOptionsScale)
}
export type AccessibleSurfaceAreaColorThemeParams = typeof AccessibleSurfaceAreaColorThemeParams
export function getAccessibleSurfaceAreaColorThemeParams(ctx: ThemeDataContext) {
    return AccessibleSurfaceAreaColorThemeParams // TODO return copy
}
export function AccessibleSurfaceAreaColorTheme(ctx: ThemeDataContext, props: PD.Values<AccessibleSurfaceAreaColorThemeParams>): ColorTheme<AccessibleSurfaceAreaColorThemeParams> {
    let color: LocationColor

    const scale = ColorScale.create({
        listOrName: props.list,
        minLabel: 'buried',
        maxLabel: 'exposed',
        domain: [0.0, 1.0]
    })

    const accessibleSurfaceArea = ctx.structure ? ComputedAccessibleSurfaceArea.get(ctx.structure)!.asa : undefined

    if (accessibleSurfaceArea) {
        color = (location: Location): Color => {
            if (StructureElement.Location.is(location)) {
                if (Unit.isAtomic(location.unit)) {
                    const value = accessibleSurfaceArea.relativeAccessibleSurfaceArea[location.unit.residueIndex[location.element]];
                    return value !== VdWLookup[0] /* signals missing value */ ? scale.color(value) : DefaultColor;
                }
            }

            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        factory: AccessibleSurfaceAreaColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: scale ? scale.legend : undefined
    }
}

export const AccessibleSurfaceAreaColorThemeProvider: ColorTheme.Provider<AccessibleSurfaceAreaColorThemeParams> = {
    label: 'Accessible Surface Area',
    factory: AccessibleSurfaceAreaColorTheme,
    getParams: getAccessibleSurfaceAreaColorThemeParams,
    defaultValues: PD.getDefaultValues(AccessibleSurfaceAreaColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
}