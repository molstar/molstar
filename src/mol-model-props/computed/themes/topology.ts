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
import { CustomProperty } from '../../common/custom-property';
import { Location } from '../../../mol-model/location';

const DefaultColor = Color(0xFAFAFA);
const Description = 'Assigns a color based on the membrane topology of a residue.';

export const TopologyColorThemeParams = {
    list: PD.ColorList('rainbow', { presetKind: 'scale' }) // change to binary
};
export type TopologyColorThemeParams = typeof TopologyColorThemeParams
export function getTopologyColorThemeParams(ctx: ThemeDataContext) {
    return TopologyColorThemeParams; // TODO return copy
}
export function TopologyColorTheme(ctx: ThemeDataContext, props: PD.Values<TopologyColorThemeParams>): ColorTheme<TopologyColorThemeParams> {
    let color: LocationColor;

    const scale = ColorScale.create({
        listOrName: props.list.colors,
        minLabel: 'membrane',
        maxLabel: 'non-membrane',
        domain: [0.0, 1.0]
    }); // prolly not needed

    // const accessibleSurfaceArea = ctx.structure && AccessibleSurfaceAreaProvider.get(ctx.structure);
    // const contextHash = accessibleSurfaceArea?.version;

    if (/*accessibleSurfaceArea?.value &&*/ ctx.structure) {
        // const asa = accessibleSurfaceArea.value;

        color = (location: Location): Color => {
            // if (StructureElement.Location.is(location) && Unit.isAtomic(location.unit)) {
                // const value = AccessibleSurfaceArea.getNormalizedValue(location, asa);
                // return value === -1 ? DefaultColor : scale.color(value);
            // }
            return DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: TopologyColorTheme,
        granularity: 'group',
        color,
        props,
        /*contextHash*/'',
        description: Description,
        legend: scale ? scale.legend : undefined
    };
}

export const TopologyColorThemeProvider: ColorTheme.Provider<TopologyColorThemeParams, 'topology'> = {
    name: 'topology',
    label: 'Membrane Topology',
    category: ColorTheme.Category.Residue,
    factory: TopologyColorTheme,
    getParams: getTopologyColorThemeParams,
    defaultValues: PD.getDefaultValues(TopologyColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? TopologyProvider.attach(ctx, data.structure, void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && data.structure.customPropertyDescriptors.reference(TopologyProvider.descriptor, false)
    }
};