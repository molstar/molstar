/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { StructureElement, Unit, ElementIndex, Bond } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { ColorTheme, LocationColor } from '../../mol-theme/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { ColorNames } from '../../mol-util/color/names';
import { MmcifFormat } from '../../mol-model-formats/structure/mmcif';
import { defaults } from '../../mol-util';

const DefaultPetworldColor = Color(0xEEEE00);
const Description = 'Petworld coloring.';

export const PetworldColorThemeParams = {
    value: PD.Color(Color(0xCCCCCC)),
    saturation: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
};
export type PetworldColorThemeParams = typeof PetworldColorThemeParams
export function getPetworldColorThemeParams(ctx: ThemeDataContext) {
    return PetworldColorThemeParams; // TODO return copy
}

function getAtomicCompId(unit: Unit.Atomic, element: ElementIndex) {
    return unit.model.atomicHierarchy.atoms.label_comp_id.value(element);
}

const lipids = ['PEA', 'PSE', 'CLR', 'PCH'];

export function PetworldColorTheme(ctx: ThemeDataContext, props: PD.Values<PetworldColorThemeParams>): ColorTheme<PetworldColorThemeParams> {
    let value = defaults(props.value, DefaultPetworldColor);
    value = Color.saturate(value, props.saturation);
    value = Color.lighten(value, props.lightness);

    let lipid = ColorNames.lightgrey;
    lipid = Color.saturate(lipid, props.saturation);
    lipid = Color.lighten(lipid, props.lightness);

    let color: LocationColor;

    const source = ctx.structure?.model.sourceData;
    if (MmcifFormat.is(source)) {
        color = (location: Location): Color => {
            if (StructureElement.Location.is(location)) {
                if (Unit.isAtomic(location.unit)) {
                    const compId = getAtomicCompId(location.unit, location.element);
                    if (lipids.includes(compId)) {
                        return lipid;
                    }
                }
                return value;
            } else if (Bond.isLocation(location)) {
                if (Unit.isAtomic(location.aUnit)) {
                    const compId = getAtomicCompId(location.aUnit, location.aUnit.elements[location.aIndex]);
                    if (lipids.includes(compId)) {
                        return lipid;
                    }
                }
                return value;
            }
            return DefaultPetworldColor;
        };
    } else {
        color = () => DefaultPetworldColor;
    }

    return {
        factory: PetworldColorTheme,
        granularity: 'group',
        preferSmoothing: false,
        color,
        props,
        description: Description,
    };
}

export const PetworldColorThemeProvider: ColorTheme.Provider<PetworldColorThemeParams, 'petworld'> = {
    name: 'petworld',
    label: 'Petworld',
    category: ColorTheme.Category.Misc,
    factory: PetworldColorTheme,
    getParams: getPetworldColorThemeParams,
    defaultValues: PD.getDefaultValues(PetworldColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};