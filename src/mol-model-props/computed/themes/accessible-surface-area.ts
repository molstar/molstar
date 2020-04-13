/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Color, ColorScale } from '../../../mol-util/color';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { StructureElement, Unit } from '../../../mol-model/structure';
import { AccessibleSurfaceAreaProvider } from '../accessible-surface-area';
import { AccessibleSurfaceArea } from '../accessible-surface-area/shrake-rupley';
import { CustomProperty } from '../../common/custom-property';
import { Location } from '../../../mol-model/location';

const DefaultColor = Color(0xFAFAFA);
const Description = 'Assigns a color based on the relative accessible surface area of a residue.';

export const AccessibleSurfaceAreaColorThemeParams = {
    list: PD.ColorList('rainbow', { presetKind: 'scale' })
};
export type AccessibleSurfaceAreaColorThemeParams = typeof AccessibleSurfaceAreaColorThemeParams
export function getAccessibleSurfaceAreaColorThemeParams(ctx: ThemeDataContext) {
    return AccessibleSurfaceAreaColorThemeParams; // TODO return copy
}
export function AccessibleSurfaceAreaColorTheme(ctx: ThemeDataContext, props: PD.Values<AccessibleSurfaceAreaColorThemeParams>): ColorTheme<AccessibleSurfaceAreaColorThemeParams> {
    let color: LocationColor;

    const scale = ColorScale.create({
        listOrName: props.list.colors,
        minLabel: 'buried',
        maxLabel: 'exposed',
        domain: [0.0, 1.0]
    });

    const accessibleSurfaceArea = ctx.structure && AccessibleSurfaceAreaProvider.get(ctx.structure);
    const contextHash = accessibleSurfaceArea?.version;

    if (accessibleSurfaceArea?.value && ctx.structure) {
        const asa = accessibleSurfaceArea.value;

        color = (location: Location): Color => {
            if (StructureElement.Location.is(location) && Unit.isAtomic(location.unit)) {
                const value = AccessibleSurfaceArea.getNormalizedValue(location, asa);
                return value === -1 ? DefaultColor : scale.color(value);
            }
            return DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: AccessibleSurfaceAreaColorTheme,
        granularity: 'group',
        color,
        props,
        contextHash,
        description: Description,
        legend: scale ? scale.legend : undefined
    };
}

export const AccessibleSurfaceAreaColorThemeProvider: ColorTheme.Provider<AccessibleSurfaceAreaColorThemeParams, 'accessible-surface-area'> = {
    name: 'accessible-surface-area',
    label: 'Accessible Surface Area',
    category: ColorTheme.Category.Residue,
    factory: AccessibleSurfaceAreaColorTheme,
    getParams: getAccessibleSurfaceAreaColorThemeParams,
    defaultValues: PD.getDefaultValues(AccessibleSurfaceAreaColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? AccessibleSurfaceAreaProvider.attach(ctx, data.structure, void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && data.structure.customPropertyDescriptors.reference(AccessibleSurfaceAreaProvider.descriptor, false)
    }
};