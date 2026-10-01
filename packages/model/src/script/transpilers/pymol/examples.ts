/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 *
 * Adapted from MolQL project
 */

export const examples = [{
    name: 'ALA residues',
    value: 'resn ALA'
}, {
    name: 'Atoms named "C", "O", "N", or "CA"',
    value: 'name c+o+n+ca'
}, {
    name: 'Residues with helix or sheet secondary structure',
    value: 'ss h+s'
}, {
    name: 'C-alpha atoms of residues 100 to 180 in chain A',
    value: 'A/100-180/CA'
}, {
    name: 'Residues 100 to 180',
    value: 'resi 100-180'
}, {
    name: 'Atoms that are 1 ang + vdw radius away from polymer',
    value: 'polymer gap 1'
}, {
    name: 'Residues within 4 ang of HEM',
    value: 'byres resn HEM around 4'
}, {
    name: 'HEM and residues within 4 ang',
    value: 'byres resn HEM expand 4'
}, {
    name: 'Solvent close (2.5 ang) to polymer',
    value: 'solvent NEAR_TO 2.5 OF polymer'
}, {
    name: 'Cystein residues within 3 ang of HEM',
    value: 'byres resn CYS WITHIN 3 OF resn HEM'
}, {
    name: 'Solvent atoms 4 ang away from oxygen',
    value: 'solvent beyond 4 of (name O and not solvent)'
}, {
    name: 'All rings in PHE',
    value: 'byring resn PHE'
}, {
    name: 'CYS and all bound residues',
    value: 'byres BOUND_TO resn CYS'
}, {
    name: 'HEM and atoms up to 7 bonds away',
    value: 'resn HEM extend 7'
}, {
    name: 'Atoms with alternate location A or none',
    value: 'alt A+""'
}];
