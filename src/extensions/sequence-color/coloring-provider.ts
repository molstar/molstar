/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { ColorTypeLocation } from '../../mol-geo/geometry/color-data';
import { Structure, StructureElement } from '../../mol-model/structure';
import { ColorTheme, LocationColor } from '../../mol-theme/color';
import { ThemeDataContext } from '../../mol-theme/theme';
import { Color } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ElementSet } from '../mvs/components/selector';
import { SequenceColorProperty } from './prop';


// export type FooProvider = ColorTheme.Provider<any, string, ColorTypeLocation>;


/** Get color assigned to a loci in custom structure property "SequenceColor" */
function sequenceColorForLoci(loci: StructureElement.Loci): Color | undefined {
    const colorData = SequenceColorProperty.Provider.get(loci.structure).value;
    if (!colorData || colorData.items.length === 0) return undefined;
    const location = StructureElement.Loci.getFirstLocation(loci);
    if (location === undefined) return undefined;
    return sequenceColorForLocation(colorData, location);
}

function sequenceColorForLocation(colorData: SequenceColorProperty.Data, location: StructureElement.Location): Color | undefined {
    const unitCache = colorData.colorCache[location.unit.id] ??= {};
    if (!(location.element in unitCache)) { // not using ??= pattern here, because cache may contain undefineds
        unitCache[location.element] = findSequenceColorForLocation(colorData, location);
    }
    return unitCache[location.element];
}

function findSequenceColorForLocation(colorData: SequenceColorProperty.Data, location: StructureElement.Location): Color | undefined {
    for (let i = colorData.items.length - 1; i >= 0; i--) { // last color matters
        const item = colorData.items[i];
        const elements = item.elementSet ??= ElementSet.fromSelector(location.structure, item.selector);
        if (ElementSet.has(elements, location)) {
            return item.color;
        }
    }
    return undefined;
}


export const CustomSequenceColorThemeParams = {};
export type CustomSequenceColorThemeParams = typeof CustomSequenceColorThemeParams
export type CustomSequenceColorThemeProps = PD.Values<CustomSequenceColorThemeParams>

const NoColor = Color(-1);

export function CustomSequenceColorTheme(ctx: ThemeDataContext, props: CustomSequenceColorThemeProps): ColorTheme<CustomSequenceColorThemeProps> { // TODO custom prop value as param?
    let color: LocationColor = () => NoColor;

    if (ctx.structure && !ctx.structure.isEmpty) {
        const colorData = SequenceColorProperty.Provider.get(ctx.structure).value;
        if (colorData?.items.length) {
            color = (location) => {
                if (StructureElement.Location.is(location)) return sequenceColorForLocation(colorData, location) ?? NoColor;
                return NoColor;
            }
        }
    }

    return {
        factory: CustomSequenceColorTheme,
        granularity: 'groupInstance',
        color: color,
        props: props,
        description: 'Assigns colors based on custom structure property `SequenceColorProperty`.',
    };
}

/** A thingy that is needed to register color theme "MVS Annotation" */
export const CustomSequenceColorThemeProvider: ColorTheme.Provider<CustomSequenceColorThemeParams, 'custom-sequence-color'> = {
    name: 'custom-sequence-color',
    label: 'Custom Sequence Color',
    category: 'Miscellaneous',
    factory: CustomSequenceColorTheme,
    getParams: ctx => CustomSequenceColorThemeParams,
    defaultValues: PD.getDefaultValues(CustomSequenceColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
};
