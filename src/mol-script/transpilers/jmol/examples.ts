/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * Adapted from MolQL project
 */

export const examples = [{
    name: 'Residue 50 or 135',
    value: '50 or 135'
}, {
    name: 'Atoms with no covalent bonds',
    value: 'bondcount = 0'
}, {
    name: 'All 3-10 helices',
    value: 'substructure = "helix310"'
}, {
    name: 'Metal atoms',
    value: 'metal'
}, {
    name: 'Atoms invloved in aromatic bonds',
    value: 'isAromatic'
}, {
    name: 'Pyrimidine residues',
    value: 'pyrimidine'
}];
