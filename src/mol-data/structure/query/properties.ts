/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Atom } from '../structure'

export const constant = {
    true: Atom.property(l => true),
    false: Atom.property(l => false),
    zero: Atom.property(l => 0)
}

export const atom = {
    type_symbol: Atom.property(l => l.unit.hierarchy.atoms.type_symbol.value(l.atom))
}

export const residue = {
    auth_seq_id: Atom.property(l => l.unit.hierarchy.residues.auth_seq_id.value(l.unit.residueIndex[l.atom])),
    auth_comp_id: Atom.property(l => l.unit.hierarchy.residues.auth_comp_id.value(l.unit.residueIndex[l.atom]))
}

export const chain = {
    auth_asym_id: Atom.property(l => l.unit.hierarchy.chains.auth_asym_id.value(l.unit.chainIndex[l.atom]))
}