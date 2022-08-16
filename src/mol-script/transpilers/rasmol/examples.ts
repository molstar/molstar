/*
 * Copyright (c) 2017 MolQL contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 */

export const rasmolSelectionsExamples = [{
    name: 'Residue 50 or 135',
    value: '(50) or (135)'
}, {
    name: 'Residue 50 or 135',
    value: '(50,135)'
}, {
    name: 'Residue 20-22 or 8 or 2-16',
    value: '(20-22,8,2-16)'
}, {
    name: 'Residue with the index more than 10',
    value: 'resno > 10'
}, {
    name: 'Atoms assinged as backbone of protein or nucleic acid',
    value: 'backbone'
}, {
    name: 'ALA, LEU, or aromatic residue within 10 angstrom from HEM ',
    value: '( [ALA] or [LEU] or aromatic ) and within( 10.0 ,[HEM])'
}, {
    name: 'Residue within 10 angstrom from HEM but not HEM',
    value: 'not [HEM] and within(5, [HEM])'
}, {
    name: 'Residue within 10 angstrom from HEM but not HEM',
    value: '( within(5, [HEM]) ) and not [HEM]'
}, {
    name: 'C-beta atom ALA10 in chain A or all C-alpha atoms or all CG atoms from LYS',
    value: '[ALA]10:A.CB or :.CA or [LYS].CG'
}, {
    name: 'Pyrimidine residues',
    value: 'pyrimidine'
}];
