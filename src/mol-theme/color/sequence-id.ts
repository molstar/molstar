/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, Bond, ElementIndex } from '../../mol-model/structure';

import { ColorScale, Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { ColorTheme } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every polymer residue a color based on its `seq_id` value.';

export const SequenceIdColorThemeParams = {
    list: PD.ColorList('turbo', { presetKind: 'scale' }),
};
export type SequenceIdColorThemeParams = typeof SequenceIdColorThemeParams
export function getSequenceIdColorThemeParams(ctx: ThemeDataContext) {
    return SequenceIdColorThemeParams; // TODO return copy
}

function getSeqId(unit: Unit, element: ElementIndex): number {
    const { model } = unit;
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            const residueIndex = model.atomicHierarchy.residueAtomSegments.index[element];
            return model.atomicHierarchy.residues.label_seq_id.value(residueIndex);
        case Unit.Kind.Spheres:
            return Math.round(
                (model.coarseHierarchy.spheres.seq_id_begin.value(element) +
                    model.coarseHierarchy.spheres.seq_id_end.value(element)) / 2
            );
        case Unit.Kind.Gaussians:
            return Math.round(
                (model.coarseHierarchy.gaussians.seq_id_begin.value(element) +
                    model.coarseHierarchy.gaussians.seq_id_end.value(element)) / 2
            );
    }
}

function getSequenceLength(unit: Unit, element: ElementIndex) {
    const { model } = unit;
    let entityId = '';
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            const chainIndex = model.atomicHierarchy.chainAtomSegments.index[element];
            entityId = model.atomicHierarchy.chains.label_entity_id.value(chainIndex);
            break;
        case Unit.Kind.Spheres:
            entityId = model.coarseHierarchy.spheres.entity_id.value(element);
            break;
        case Unit.Kind.Gaussians:
            entityId = model.coarseHierarchy.gaussians.entity_id.value(element);
            break;
    }
    if (entityId === '') return 0;
    const entityIndex = model.entities.getEntityIndex(entityId);
    if (entityIndex === -1) return 0;
    return model.sequence.byEntityKey[entityIndex].sequence.length;
}

export function SequenceIdColorTheme(ctx: ThemeDataContext, props: PD.Values<SequenceIdColorThemeParams>): ColorTheme<SequenceIdColorThemeParams> {
    const scale = ColorScale.create({
        listOrName: props.list.colors,
        minLabel: 'Start',
        maxLabel: 'End',
    });
    const color = (location: Location): Color => {
        if (StructureElement.Location.is(location)) {
            const { unit, element } = location;
            const seq_id = getSeqId(unit, element);
            if (seq_id > 0) {
                scale.setDomain(0, getSequenceLength(unit, element) - 1);
                return scale.color(seq_id);
            }
        } else if (Bond.isLocation(location)) {
            const { aUnit, aIndex } = location;
            const seq_id = getSeqId(aUnit, aUnit.elements[aIndex]);
            if (seq_id > 0) {
                scale.setDomain(0, getSequenceLength(aUnit, aUnit.elements[aIndex]) - 1);
                return scale.color(seq_id);
            }
        }
        return DefaultColor;
    };

    return {
        factory: SequenceIdColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: scale ? scale.legend : undefined
    };
}

export const SequenceIdColorThemeProvider: ColorTheme.Provider<SequenceIdColorThemeParams, 'sequence-id'> = {
    name: 'sequence-id',
    label: 'Sequence Id',
    category: ColorTheme.Category.Residue,
    factory: SequenceIdColorTheme,
    getParams: getSequenceIdColorThemeParams,
    defaultValues: PD.getDefaultValues(SequenceIdColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};