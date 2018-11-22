/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit, StructureElement, StructureProperties as Props, Link } from 'mol-model/structure';
import { Loci } from 'mol-model/loci';
import { OrderedSet } from 'mol-data/int';

// for `labelFirst`, don't create right away to avoid problems with circular dependencies/imports
let elementLocA: StructureElement
let elementLocB: StructureElement

function setElementLocation(loc: StructureElement, unit: Unit, index: StructureElement.UnitIndex) {
    loc.unit = unit
    loc.element = unit.elements[index]
}

export function labelFirst(loci: Loci): string {
    switch (loci.kind) {
        case 'structure-loci':
            return loci.structure.models.map(m => m.label).join(', ')
        case 'element-loci':
            const e = loci.elements[0]
            if (e) {
                const el = e.unit.elements[OrderedSet.getAt(e.indices, 0)];
                return elementLabel(StructureElement.create(e.unit, el))
            } else {
                return 'Unknown'
            }
        case 'link-loci':
            const link = loci.links[0]
            return link ? linkLabel(link) : 'Unknown'
        case 'group-loci':
            const g = loci.groups[0]
            if (g) {
                return loci.shape.labels.ref.value[OrderedSet.getAt(g.ids, 0)]
            } else {
                return 'Unknown'
            }
        case 'every-loci':
            return 'Everything'
        case 'empty-loci':
            return 'Nothing'
    }
}

export function linkLabel(link: Link.Location) {
    if (!elementLocA) elementLocA = StructureElement.create()
    if (!elementLocB) elementLocB = StructureElement.create()
    setElementLocation(elementLocA, link.aUnit, link.aIndex)
    setElementLocation(elementLocB, link.bUnit, link.bIndex)
    return `${elementLabel(elementLocA)} - ${elementLabel(elementLocB)}`
}

export function elementLabel(location: StructureElement) {
    const model = location.unit.model.label
    const instance = location.unit.conformation.operator.name
    let label = ''

    if (Unit.isAtomic(location.unit)) {
        const asym_id = Props.chain.auth_asym_id(location)
        const seq_id = Props.residue.auth_seq_id(location)
        const comp_id = Props.residue.auth_comp_id(location)
        const atom_id = Props.atom.auth_atom_id(location)
        const alt_id = Props.atom.label_alt_id(location)
        label = `[${comp_id}]${seq_id}:${asym_id}.${atom_id}${alt_id ? `%${alt_id}` : ''}`
    } else if (Unit.isCoarse(location.unit)) {
        const asym_id = Props.coarse.asym_id(location)
        const seq_id_begin = Props.coarse.seq_id_begin(location)
        const seq_id_end = Props.coarse.seq_id_end(location)
        if (seq_id_begin === seq_id_end) {
            const entityIndex = Props.coarse.entityKey(location)
            const seq = location.unit.model.sequence.byEntityKey[entityIndex]
            const comp_id = seq.compId.value(seq_id_begin - 1) // 1-indexed
            label = `[${comp_id}]${seq_id_begin}:${asym_id}`
        } else {
            label = `${seq_id_begin}-${seq_id_end}:${asym_id}`
        }
    } else {
        label = 'unknown'
    }

    return `${model} ${instance} ${label}`
}