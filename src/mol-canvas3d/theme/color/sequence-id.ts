/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, Link, ElementIndex } from 'mol-model/structure';

import { ColorScale, Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { ColorThemeProps, ColorTheme } from '../color';
import { ColorOther } from 'mol-util/color/tables';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives every polymer residue a color based on its `seq_id` value.'

function getSeqId(unit: Unit, element: ElementIndex): number {
    const { model } = unit
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            const residueIndex = model.atomicHierarchy.residueAtomSegments.index[element]
            return model.atomicHierarchy.residues.label_seq_id.value(residueIndex)
        case Unit.Kind.Spheres:
            return Math.round(
                (model.coarseHierarchy.spheres.seq_id_begin.value(element) +
                    model.coarseHierarchy.spheres.seq_id_end.value(element)) / 2
            )
        case Unit.Kind.Gaussians:
            return Math.round(
                (model.coarseHierarchy.gaussians.seq_id_begin.value(element) +
                    model.coarseHierarchy.gaussians.seq_id_end.value(element)) / 2
            )
    }
}

function getSequenceLength(unit: Unit, element: ElementIndex) {
    const { model } = unit
    let entityId = ''
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            const chainIndex = model.atomicHierarchy.chainAtomSegments.index[element]
            entityId = model.atomicHierarchy.chains.label_entity_id.value(chainIndex)
            break
        case Unit.Kind.Spheres:
            entityId = model.coarseHierarchy.spheres.entity_id.value(element)
            break
        case Unit.Kind.Gaussians:
            entityId = model.coarseHierarchy.gaussians.entity_id.value(element)
            break
    }
    if (entityId === '') return 0
    const entityIndex = model.entities.getEntityIndex(entityId)
    if (entityIndex === -1) return 0
    return model.sequence.byEntityKey[entityIndex].sequence.sequence.length
}

export function SequenceIdColorTheme(props: ColorThemeProps): ColorTheme {
    const p = {
        ...props,
        colors: ColorOther.rainbow,
        minLabel: 'Start',
        maxLabel: 'End',
    }

    const scale = ColorScale.create(p)
    const color = (location: Location): Color => {
        if (StructureElement.isLocation(location)) {
            const { unit, element } = location
            const seq_id = getSeqId(unit, element)
            if (seq_id > 0) {
                scale.setDomain(0, getSequenceLength(unit, element) - 1)
                return scale.color(seq_id)
            }
        } else if (Link.isLocation(location)) {
            const { aUnit, aIndex } = location
            const seq_id = getSeqId(aUnit, aUnit.elements[aIndex])
            if (seq_id > 0) {
                scale.setDomain(0, getSequenceLength(aUnit, aUnit.elements[aIndex]) - 1)
                return scale.color(seq_id)
            }
        }
        return DefaultColor
    }

    return {
        granularity: 'group',
        color,
        description: Description,
        // legend: scale ? TableLegend(table) : undefined
        legend: scale ? scale.legend : undefined
    }
}