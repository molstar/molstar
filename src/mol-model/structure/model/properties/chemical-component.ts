/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MoleculeType, ComponentType } from '../types'

export interface ChemicalComponent {
    id: string
    type: ComponentType
    moleculeType: MoleculeType
    name: string
    synonyms: string[]
    formula: string
    formulaWeight: number
}

export type ChemicalComponentMap = ReadonlyMap<string, ChemicalComponent>

const CommonChemicalComponents: ChemicalComponent[] = [
    {
        id: 'ALA',
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'ALANINE',
        synonyms: [],
        formula: 'C3 H7 N O2',
        formulaWeight: 89.093
    },
    { 
        id: 'ARG',
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'ARGININE',
        synonyms: [],
        formula: 'C6 H15 N4 O2 1',
        formulaWeight: 175.209
    },
    { 
        id: 'ASN', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'ASPARAGINE',
        synonyms: [],
        formula: 'C4 H8 N2 O3',
        formulaWeight: 132.118
    },
    { 
        id: 'ASP', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'ASPARTIC ACID',
        synonyms: [],
        formula: 'C4 H7 N O4',
        formulaWeight: 133.103
    },
    { 
        id: 'CYS', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'CYSTEINE',
        synonyms: [],
        formula: 'C3 H7 N O2 S',
        formulaWeight: 121.158
    },
    { 
        id: 'GLN', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'GLUTAMINE',
        synonyms: [],
        formula: 'C5 H10 N2 O3',
        formulaWeight: 146.144
    },
    { 
        id: 'GLU', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'GLUTAMIC ACID',
        synonyms: [],
        formula: 'C5 H9 N O4',
        formulaWeight: 147.129
    },
    { 
        id: 'GLY', 
        type: ComponentType['peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'GLYCINE',
        synonyms: [],
        formula: 'C2 H5 N O2',
        formulaWeight: 75.067
    },
    { 
        id: 'HIS', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'HISTIDINE',
        synonyms: [],
        formula: 'C6 H10 N3 O2 1',
        formulaWeight: 156.162
    },
    { 
        id: 'ILE', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'ISOLEUCINE',
        synonyms: [],
        formula: 'C6 H13 N O2',
        formulaWeight: 131.173
    },
    { 
        id: 'LEU', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'LEUCINE',
        synonyms: [],
        formula: 'C6 H13 N O2',
        formulaWeight: 131.173
    },
    { 
        id: 'LYS', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'LYSINE',
        synonyms: [],
        formula: 'C6 H15 N2 O2 1',
        formulaWeight: 147.195
    },
    { 
        id: 'MET', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'METHIONINE',
        synonyms: [],
        formula: 'C5 H11 N O2 S',
        formulaWeight: 149.211 
    },
    { 
        id: 'PHE', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'PHENYLALANINE',
        synonyms: [],
        formula: 'C9 H11 N O2',
        formulaWeight: 165.19 
    },
    { 
        id: 'PRO', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'PROLINE',
        synonyms: [],
        formula: 'C5 H9 N O2',
        formulaWeight: 115.13
    },
    { // 'O' as per IUPAC definition
        id: 'PYL', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'PYRROLYSINE',
        synonyms: [],
        formula: 'C12 H21 N3 O3',
        formulaWeight: 255.31
    },
    { // 'U' as per IUPAC definition
        id: 'SEC', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'SELENOCYSTEINE',
        synonyms: [],
        formula: 'C3 H7 N O2 Se',
        formulaWeight: 168.05
    },
    {
        id: 'SER', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'SERINE',
        synonyms: [],
        formula: 'C3 H7 N O3',
        formulaWeight: 105.09
    },
    {
        id: 'THR', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'THREONINE',
        synonyms: [],
        formula: 'C4 H9 N O3',
        formulaWeight: 119.12
    },
    {
        id: 'TRP', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'TRYPTOPHAN',
        synonyms: [],
        formula: 'C11 H12 N2 O2',
        formulaWeight: 204.22
    },
    {
        id: 'TYR', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'TYROSINE',
        synonyms: [],
        formula: 'C9 H11 N O3',
        formulaWeight: 181.19
    },
    {
        id: 'VAL', 
        type: ComponentType['L-peptide linking'],
        moleculeType: MoleculeType.protein,
        name: 'VALINE',
        synonyms: [],
        formula: 'C5 H11 N O2',
        formulaWeight: 117.15
    }
]
export const CommonChemicalComponentMap = new Map()
for (let i = 0, il = CommonChemicalComponents.length; i < il; ++i) {
    CommonChemicalComponentMap.set(CommonChemicalComponents[i].id, CommonChemicalComponents[i])
}