/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit, Element, Queries } from 'mol-model/structure';
import { Loci } from 'mol-model/loci';

const elementLocA = Element.Location()
const elementLocB = Element.Location()

function setElementLocation(loc: Element.Location, unit: Unit, index: number) {
    loc.unit = unit
    loc.element = unit.elements[index]
}

export function labelFirst(loci: Loci): string {
    switch (loci.kind) {
        case 'element-loci':
            const e = loci.elements[0]
            if (e && e.indices[0] !== undefined) {
                return elementLabel(Element.Location(e.unit, e.indices[0]))
            } else {
                return 'Unknown'
            }
        case 'link-loci':
            const bond = loci.links[0]
            if (bond) {
                setElementLocation(elementLocA, bond.aUnit, bond.aIndex)
                setElementLocation(elementLocB, bond.bUnit, bond.bIndex)
                return `${elementLabel(elementLocA)} - ${elementLabel(elementLocB)}`
            } else {
                return 'Unknown'
            }
        case 'every-loci':
            return 'Everything'
        case 'empty-loci':
            return 'Nothing'
    }
}

export function elementLabel(loc: Element.Location) {
    const model = loc.unit.model.label
    const instance = loc.unit.conformation.operator.name
    let element = ''

    if (Unit.isAtomic(loc.unit)) {
        const asym_id = Queries.props.chain.auth_asym_id(loc)
        const seq_id = Queries.props.residue.auth_seq_id(loc)
        const comp_id = Queries.props.residue.auth_comp_id(loc)
        const atom_id = Queries.props.atom.auth_atom_id(loc)
        element = `[${comp_id}]${seq_id}:${asym_id}.${atom_id}`
    } else if (Unit.isCoarse(loc.unit)) {
        const asym_id = Queries.props.coarse.asym_id(loc)
        const seq_id_begin = Queries.props.coarse.seq_id_begin(loc)
        const seq_id_end = Queries.props.coarse.seq_id_end(loc)
        if (seq_id_begin === seq_id_end) {
            const entityKey = Queries.props.coarse.entityKey(loc)
            const seq = loc.unit.model.sequence.byEntityKey[entityKey]
            const comp_id = seq.compId.value(seq_id_begin)
            element = `[${comp_id}]${seq_id_begin}:${asym_id}`
        } else {
            element = `${seq_id_begin}-${seq_id_end}:${asym_id}`
        }
    } else {
        element = 'unknown'
    }

    return `${model} ${instance} ${element}`
}