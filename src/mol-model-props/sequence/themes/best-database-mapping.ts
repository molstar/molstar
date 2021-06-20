/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Location } from '../../../mol-model/location';
import { Bond, StructureElement, Unit } from '../../../mol-model/structure';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { Color } from '../../../mol-util/color';
import { getPalette, getPaletteParams } from '../../../mol-util/color/palette';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { CustomProperty } from '../../common/custom-property';
import { BestDatabaseSequenceMapping } from '../best-database-mapping';

const DefaultColor = Color(0xFAFAFA);
const Description = 'Assigns a color based on best dababase sequence mapping.';

// same colors for same accessions
const globalAccessionMap = new Map<string, number>();

export const BestDatabaseSequenceMappingColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: 'set-1' }),
};
export type BestDatabaseSequenceMappingColorThemeParams = typeof BestDatabaseSequenceMappingColorThemeParams
export function getBestDatabaseSequenceMappingColorThemeParams(ctx: ThemeDataContext) {
    return BestDatabaseSequenceMappingColorThemeParams; // TODO return copy
}
export function BestDatabaseSequenceMappingColorTheme(ctx: ThemeDataContext, props: PD.Values<BestDatabaseSequenceMappingColorThemeParams>): ColorTheme<BestDatabaseSequenceMappingColorThemeParams> {
    let color: LocationColor;

    if (ctx.structure) {
        for (const m of ctx.structure.models) {
            const mapping = BestDatabaseSequenceMapping.Provider.get(m).value;
            if (!mapping) continue;
            for (const acc of mapping.accession) {
                if (!acc || globalAccessionMap.has(acc)) continue;
                globalAccessionMap.set(acc, globalAccessionMap.size);
            }
        }

        const l = StructureElement.Location.create(ctx.structure);
        const palette = getPalette(globalAccessionMap.size + 1, props, { valueLabel: i => `${i}` });
        const colorMap = new Map<string, Color>();

        const getColor = (location: StructureElement.Location) => {
            const key = BestDatabaseSequenceMapping.getKey(location);
            if (!key) return DefaultColor;

            if (colorMap.has(key)) return colorMap.get(key)!;

            const color = palette.color(globalAccessionMap.get(key)!);
            colorMap.set(key, color);
            return color;
        };

        color = (location: Location): Color => {
            if (StructureElement.Location.is(location) && Unit.isAtomic(location.unit)) {
                return getColor(location);
            } else if (Bond.isLocation(location)) {
                l.unit = location.aUnit;
                l.element = location.aUnit.elements[location.aIndex];
                return getColor(l);
            }
            return DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: BestDatabaseSequenceMappingColorTheme,
        granularity: 'group',
        preferSmoothing: true,
        color,
        props,
        description: Description,
    };
}

export const BestDatabaseSequenceMappingColorThemeProvider: ColorTheme.Provider<BestDatabaseSequenceMappingColorThemeParams, 'best-sequence-database-mapping'> = {
    name: 'best-sequence-database-mapping',
    label: 'Best Database Sequence Mapping',
    category: ColorTheme.Category.Residue,
    factory: BestDatabaseSequenceMappingColorTheme,
    getParams: getBestDatabaseSequenceMappingColorThemeParams,
    defaultValues: PD.getDefaultValues(BestDatabaseSequenceMappingColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure?.models.some(m => BestDatabaseSequenceMapping.Provider.isApplicable(m)),
    ensureCustomProperties: {
        attach: async (ctx: CustomProperty.Context, data: ThemeDataContext) => {
            if (!data.structure) return;

            for (const m of data.structure.models) {
                await BestDatabaseSequenceMapping.Provider.attach(ctx, m, void 0, true);
            }
        },
        detach: (data) => {
            if (!data.structure) return;
            for (const m of data.structure.models) {
                BestDatabaseSequenceMapping.Provider.ref(m, false);
            }
        }
    }
};