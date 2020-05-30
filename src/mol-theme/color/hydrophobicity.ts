/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorScale } from '../../mol-util/color';
import { StructureElement, Unit, Bond, ElementIndex } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { ColorTheme } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { ResidueHydrophobicity } from '../../mol-model/structure/model/types';

const Description = 'Assigns a color to every amino acid according to the "Experimentally determined hydrophobicity scale for proteins at membrane interfaces" by Wimely and White (doi:10.1038/nsb1096-842).';

export const HydrophobicityColorThemeParams = {
    list: PD.ColorList('red-yellow-green', { presetKind: 'scale' }),
    scale: PD.Select('DGwif', [['DGwif', 'DG water-membrane'], ['DGwoct', 'DG water-octanol'], ['Oct-IF', 'DG difference']] as const)
};
export type HydrophobicityColorThemeParams = typeof HydrophobicityColorThemeParams
export function getHydrophobicityColorThemeParams(ctx: ThemeDataContext) {
    return HydrophobicityColorThemeParams; // TODO return copy
}

const scaleIndexMap = { 'DGwif': 0, 'DGwoct': 1, 'Oct-IF': 2 };

export function hydrophobicity(compId: string, scaleIndex: number): number {
    const c = (ResidueHydrophobicity as { [k: string]: number[] })[compId];
    return c === undefined ? 0 : c[scaleIndex];
}

function getAtomicCompId(unit: Unit.Atomic, element: ElementIndex) {
    return unit.model.atomicHierarchy.atoms.label_comp_id.value(unit.residueIndex[element]);
}

function getCoarseCompId(unit: Unit.Spheres | Unit.Gaussians, element: ElementIndex) {
    const seqIdBegin = unit.coarseElements.seq_id_begin.value(element);
    const seqIdEnd = unit.coarseElements.seq_id_end.value(element);
    if (seqIdBegin === seqIdEnd) {
        const entityKey = unit.coarseElements.entityKey[element];
        const seq = unit.model.sequence.byEntityKey[entityKey].sequence;
        return seq.compId.value(seqIdBegin - 1); // 1-indexed
    }
}

export function HydrophobicityColorTheme(ctx: ThemeDataContext, props: PD.Values<HydrophobicityColorThemeParams>): ColorTheme<HydrophobicityColorThemeParams> {
    const scaleIndex = scaleIndexMap[props.scale];

    // get domain
    let min = Infinity;
    let max = -Infinity;
    for (const name in ResidueHydrophobicity) {
        const val = (ResidueHydrophobicity as { [k: string]: number[] })[name][scaleIndex];
        min = Math.min(min, val);
        max = Math.max(max, val);
    }

    const scale = ColorScale.create({
        listOrName: props.list.colors,
        domain: [ max, min ],
        minLabel: 'Hydrophobic',
        maxLabel: 'Hydrophilic'
    });

    function color(location: Location): Color {
        let compId: string | undefined;
        if (StructureElement.Location.is(location)) {
            if (Unit.isAtomic(location.unit)) {
                compId = getAtomicCompId(location.unit, location.element);
            } else {
                compId = getCoarseCompId(location.unit, location.element);
            }
        } else if (Bond.isLocation(location)) {
            if (Unit.isAtomic(location.aUnit)) {
                compId = getAtomicCompId(location.aUnit, location.aUnit.elements[location.aIndex]);
            } else {
                compId = getCoarseCompId(location.aUnit, location.aUnit.elements[location.aIndex]);
            }
        }
        return scale.color(compId ? hydrophobicity(compId, scaleIndex) : 0);
    }

    return {
        factory: HydrophobicityColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: scale ? scale.legend : undefined
    };
}

export const HydrophobicityColorThemeProvider: ColorTheme.Provider<HydrophobicityColorThemeParams, 'hydrophobicity'> = {
    name: 'hydrophobicity',
    label: 'Hydrophobicity',
    category: ColorTheme.Category.Residue,
    factory: HydrophobicityColorTheme,
    getParams: getHydrophobicityColorThemeParams,
    defaultValues: PD.getDefaultValues(HydrophobicityColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};