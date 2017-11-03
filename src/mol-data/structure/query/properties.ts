/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Atom } from '../structure'

const constant = {
    true: Atom.property(l => true),
    false: Atom.property(l => false),
    zero: Atom.property(l => 0)
}

const atom = {
    x: Atom.property(l => l.unit.x(l.atom)),
    y: Atom.property(l => l.unit.y(l.atom)),
    z: Atom.property(l => l.unit.z(l.atom)),

    type_symbol: Atom.property(l => l.unit.hierarchy.atoms.type_symbol.value(l.atom))
}

const residue = {
    key: Atom.property(l => l.unit.hierarchy.residueKey.value(l.unit.residueIndex[l.atom])),

    auth_seq_id: Atom.property(l => l.unit.hierarchy.residues.auth_seq_id.value(l.unit.residueIndex[l.atom])),
    auth_comp_id: Atom.property(l => l.unit.hierarchy.residues.auth_comp_id.value(l.unit.residueIndex[l.atom]))
}

const chain = {
    auth_asym_id: Atom.property(l => l.unit.hierarchy.chains.auth_asym_id.value(l.unit.chainIndex[l.atom]))
}

const Properties = {
    constant,
    atom,
    residue,
    chain
}

type Properties = typeof Properties
export default Properties