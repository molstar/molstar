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
import { Unit, StructureElement, StructureProperties } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { AccessibleSurfaceAreaProvider } from '../../mol-model-props/computed/accessible-surface-area';
import { ColorListName, ColorListOptionsScale } from '../../mol-util/color/lists';
import { AccessibleSurfaceArea } from '../../mol-model-props/computed/accessible-surface-area/shrake-rupley';
import { CustomProperty } from '../../mol-model-props/common/custom-property';

const DefaultColor = Color(0xFAFAFA)
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

    const { label_comp_id } = StructureProperties.residue
    const accessibleSurfaceArea = ctx.structure && AccessibleSurfaceAreaProvider.get(ctx.structure)
    const contextHash = accessibleSurfaceArea?.version

    if (accessibleSurfaceArea?.value && ctx.structure) {
        const { getSerialIndex } = ctx.structure.root.serialMapping
        const { area, serialResidueIndex } = accessibleSurfaceArea.value

        color = (location: Location): Color => {
            if (StructureElement.Location.is(location) && Unit.isAtomic(location.unit)) {
                const rSI = serialResidueIndex[getSerialIndex(location.unit, location.element)]
                return rSI === -1 ? DefaultColor : scale.color(AccessibleSurfaceArea.normalize(label_comp_id(location), area[rSI]))
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
        contextHash,
        description: Description,
        legend: scale ? scale.legend : undefined
    }
}

export const AccessibleSurfaceAreaColorThemeProvider: ColorTheme.Provider<AccessibleSurfaceAreaColorThemeParams> = {
    label: 'Accessible Surface Area',
    factory: AccessibleSurfaceAreaColorTheme,
    getParams: getAccessibleSurfaceAreaColorThemeParams,
    defaultValues: PD.getDefaultValues(AccessibleSurfaceAreaColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
    ensureCustomProperties: (ctx: CustomProperty.Context, data: ThemeDataContext) => {
        return data.structure ? AccessibleSurfaceAreaProvider.attach(ctx, data.structure) : Promise.resolve()
    }
}