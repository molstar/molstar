/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 *
 * Adapted from MolQL project
 */

export const examples = [{
    name: 'All water residues',
    value: 'water'
}, {
    name: 'All C-alpha atoms',
    value: 'name CA'
}, {
    name: 'Residue 35',
    value: 'resid 35'
}, {
    name: 'C-alpha atoms of ALA',
    value: 'name CA and resname ALA'
}, {
    name: 'Backbone atoms',
    value: 'backbone'
}, {
    name: 'Non-protein atoms',
    value: 'not protein'
}, {
    name: 'Protein backbone or hydrogen atoms',
    value: 'protein (backbone or name H)'
}, {
    name: 'Atoms heavier than 20',
    value: 'mass > 20'
}, {
    name: 'Atoms with two bonds',
    value: 'numbonds = 2'
}, {
    name: 'Atoms with an absolute charge greater 1',
    value: 'abs(charge) > 1'
}, {
    name: 'Atoms with an x coordinate between -25 and -20',
    value: 'x < -20 and x > -25'
}, {
    name: 'Helices',
    value: 'structure H'
}, {
    name: 'Atoms with name "A 1"',
    value: "name 'A 1'"
}, {
    name: 'Atoms with name "A *"',
    value: "name 'A *'"
}, {
    name: 'Atoms with names starting with C',
    value: 'name "C.*"'
}, {
    name: 'Atoms within 10 ang of [25, 15, 10]',
    value: 'sqr(x+25)+sqr(y+15)+sqr(z+10) <= sqr(10)'
}, {
    name: 'Atoms within 5 ang of iron atoms',
    value: 'within 5 of name FE'
}, {
    name: 'Atoms around 10 ang of HEM residue',
    value: 'exwithin 10 of resname HEM'
}, {
    name: 'ALA residues within 15 ang of HEM',
    value: 'resname ALA within 15 of resname HEM'
}, {
    name: 'All groups that include an iron atom',
    value: 'same resid as name FE'
}, {
    name: 'Atoms with mass between 12 and 17.5',
    value: 'mass 12 to 17.5'
}, {
    name: 'Residues 60, 80, 90 and 142',
    value: 'resid 60 80 90 142'
}/* , {
    name: 'Residues ala, arg, asn, asp, cys, and tyr',
    value: 'resname ALA to CYS TYR'
}*/];
