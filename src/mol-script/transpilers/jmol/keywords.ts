/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Adapted from MolQL project
 */


import { MolScriptBuilder } from '../../../mol-script/language/builder';
const B = MolScriptBuilder;
import * as h from '../helper';
import { KeywordDict } from '../types';

const ResDict = {
    acidic: ['ASP', 'GLU'],
    aliphatic: ['ALA', 'GLY', 'ILE', 'LEU', 'VAL'],
    amino: ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ASX', 'GLX', 'UNK'],
    aromatic: ['HIS', 'PHE', 'TRP', 'TYR'],
    basic: ['ARG', 'HIS', 'LYS'],
    buried: ['ALA', 'CYS', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'VAL'],
    cg: ['CYT', 'C', 'GUA', 'G'],
    cyclic: ['HIS', 'PHE', 'PRO', 'TRP', 'TYR'],
    hydrophobic: ['ALA', 'GLY', 'ILE', 'LEU', 'MET', 'PHE', 'PRO', 'TRP', 'TYR', 'VAL'],
    large: ['ARG', 'GLU', 'GLN', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'TRP', 'TYR'],
    medium: ['ASN', 'ASP', 'CYS', 'PRO', 'THR', 'VAL'],
    small: ['ALA', 'GLY', 'SER'],

    nucleic: ['G', 'C', 'A', 'T', 'U', 'I', 'DG', 'DC', 'DA', 'DT', 'DU', 'DI', '+G', '+C', '+A', '+T', '+U', '+I']
};

const Backbone = {
    nucleic: ['P', "O3'", "O5'", "C5'", "C4'", "C3'", 'OP1', 'OP2', 'O3*', 'O5*', 'C5*', 'C4*', 'C3*',
        "C2'", "C1'", "O4'", "O2'"],
    protein: ['C', 'N', 'CA']
};

function nucleicExpr() {
    return B.struct.combinator.merge([
        B.struct.generator.atomGroups({
            'residue-test': B.core.set.has([
                B.set(...ResDict.nucleic),
                B.ammp('label_comp_id')
            ])
        }),
        B.struct.filter.pick({
            0: B.struct.generator.atomGroups({
                'group-by': B.ammp('residueKey')
            }),
            test: B.core.logic.and([
                B.core.rel.eq([B.struct.atomSet.atomCount(), 1]),
                B.core.rel.eq([B.ammp('label_atom_id'), B.atomName('P')]),
            ])
        }),
        B.struct.filter.pick({
            0: B.struct.generator.atomGroups({
                'group-by': B.ammp('residueKey')
            }),
            test: B.core.logic.or([
                B.core.set.isSubset([
                    h.atomNameSet(["C1'", "C2'", "O3'", "C3'", "C4'", "C5'", "O5'"]),
                    B.ammpSet('label_atom_id')
                ]),
                B.core.set.isSubset([
                    h.atomNameSet(['C1*', 'C2*', 'O3*', 'C3*', 'C4*', 'C5*', 'O5*']),
                    B.ammpSet('label_atom_id')
                ])
            ])
        })
    ]);
}

// TODO: improve, see keywords.protein['@desc'] below
function proteinExpr() {
    return B.struct.generator.atomGroups({
        'residue-test': B.core.set.has([
            B.set(...ResDict.amino),
            B.ammp('label_comp_id')
        ])
    });
}

// TODO: improve, see keywords.backbone['@desc'] below
function backboneExpr() {
    return B.struct.combinator.merge([
        B.struct.modifier.intersectBy({
            0: B.struct.generator.atomGroups({
                'residue-test': B.core.set.has([
                    B.core.type.set(ResDict.amino),
                    B.ammp('label_comp_id')
                ])
            }),
            by: B.struct.generator.atomGroups({
                'atom-test': B.core.set.has([
                    B.core.type.set(Backbone.protein),
                    B.ammp('label_atom_id')
                ])
            })
        }),
        B.struct.modifier.intersectBy({
            0: B.struct.generator.atomGroups({
                'residue-test': B.core.set.has([
                    B.core.type.set(ResDict.nucleic),
                    B.ammp('label_comp_id')
                ])
            }),
            by: B.struct.generator.atomGroups({
                'atom-test': B.core.set.has([
                    B.core.type.set(Backbone.nucleic),
                    B.ammp('label_atom_id')
                ])
            })
        }),
    ]);
}

export const keywords: KeywordDict = {
    // general terms
    all: {
        '@desc': 'all atoms; same as *',
        abbr: ['*'],
        map: () => B.struct.generator.all()
    },
    bonded: {
        '@desc': 'covalently bonded',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.rel.gr([B.struct.atomProperty.core.bondCount({
                flags: B.struct.type.bondFlags(['covalent', 'metallic', 'sulfide'])
            }), 0])
        })
    },
    clickable: {
        '@desc': 'actually visible -- having some visible aspect such as wireframe, spacefill, or a label showing, or the alpha-carbon or phosphorus atom in a biomolecule that is rendered with only cartoon, rocket, or other biomolecule-specific shape.'
    },
    connected: {
        '@desc': 'bonded in any way, including hydrogen bonds',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.rel.gr([B.struct.atomProperty.core.bondCount({
                flags: B.struct.type.bondFlags()
            }), 0])
        })
    },
    displayed: {
        '@desc': 'displayed using the display or hide command; not necessarily visible'
    },
    hidden: {
        '@desc': 'hidden using the display or hide command'
    },
    none: {
        '@desc': 'no atoms',
        map: () => B.struct.generator.empty()
    },
    selected: {
        '@desc': 'atoms that have been selected; defaults to all when a file is first loaded'
    },
    thisModel: {
        '@desc': 'atoms in the current frame set, as defined by frame, model, or animation commands. If more than one model is in this set, "thisModel" refers to all of them, regardless of atom displayed/hidden status.'
    },
    visible: {
        '@desc': 'visible in any way, including PDB residue atoms for which a cartoon or other such rendering makes their group visible, even if they themselves are not visible.'
    },
    subset: {
        '@desc': 'the currently defined subset. Note that if a subset is currently defined, then select/display all is the same as select/display subset, restrict none is the same as restrict not subset. In addition, select not subset selects nothing.'
    },
    specialPosition: {
        '@desc': 'atoms in crystal structures that are at special positions - that is, for which there is more than one operator that leads to them.'
    },
    unitcell: {
        '@desc': 'atoms within the current unitcell, which may be offset. This includes atoms on the faces and at the vertices of the unitcell.'
    },
    polyhedra: {
        '@desc': 'all central atoms for which polyhedra have been created. See also polyhera(n), below. (Jmol 14.4)'
    },
    nonmetal: {
        '@desc': '_H,_He,_B,_C,_N,_O,_F,_Ne,_Si,_P,_S,_Cl,_Ar,_As,_Se,_Br,_Kr,_Te,_I,_Xe,_At,_Rn',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.set.has([
                B.set(...['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Si', 'P', 'S', 'Cl', 'Ar', 'As', 'Se', 'Br', 'Kr', 'Te', 'I', 'Xe', 'At', 'Rn'].map(B.es)),
                B.acp('elementSymbol')
            ])
        })
    },
    metal: {
        '@desc': '!nonmetal',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.logic.not([
                B.core.set.has([
                    B.set(...['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Si', 'P', 'S', 'Cl', 'Ar', 'As', 'Se', 'Br', 'Kr', 'Te', 'I', 'Xe', 'At', 'Rn'].map(B.es)),
                    B.acp('elementSymbol')
                ])
            ])
        })
    },
    alkaliMetal: {
        '@desc': '_Li,_Na,_K,_Rb,_Cs,_Fr',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.set.has([
                B.set(...['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr'].map(B.es)),
                B.acp('elementSymbol')
            ])
        })
    },
    alkalineEarth: {
        '@desc': '_Be,_Mg,_Ca,_Sr,_Ba,_Ra',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.set.has([
                B.set(...['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'].map(B.es)),
                B.acp('elementSymbol')
            ])
        })
    },
    nobleGas: {
        '@desc': '_He,_Ne,_Ar,_Kr,_Xe,_Rn',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.set.has([
                B.set(...['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn'].map(B.es)),
                B.acp('elementSymbol')
            ])
        })
    },
    metalloid: {
        '@desc': '_B,_Si,_Ge,_As,_Sb,_Te',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.set.has([
                B.set(...['B', 'Si', 'Ge', 'As', 'Sb', 'Te'].map(B.es)),
                B.acp('elementSymbol')
            ])
        })
    },
    transitionMetal: {
        '@desc': '(includes La and Ac) elemno>=21 and elemno<=30, elemno=57, elemno=89, elemno>=39 and elemno<=48, elemno>=72 and elemno<=80, elemno>=104 and elemno<=112',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.logic.or([
                B.core.rel.inRange([B.acp('atomicNumber'), 21, 30]),
                B.core.rel.inRange([B.acp('atomicNumber'), 39, 48]),
                B.core.rel.inRange([B.acp('atomicNumber'), 72, 80]),
                B.core.rel.inRange([B.acp('atomicNumber'), 104, 112]),
                B.core.set.has([B.set(57, 89), B.acp('atomicNumber')])
            ])
        })
    },
    lanthanide: {
        '@desc': '(does not include La) elemno>57 and elemno<=71',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.rel.inRange([B.acp('atomicNumber'), 57, 71])
        })
    },
    actinide: {
        '@desc': '(does not include Ac) elemno>89 and elemno<=103',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.rel.inRange([B.acp('atomicNumber'), 89, 103])
        })
    },
    isaromatic: {
        '@desc': 'atoms connected with the AROMATIC, AROMATICSINGLE, or AROMATICDOUBLE bond types',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.rel.gr([
                B.struct.atomProperty.core.bondCount({
                    flags: B.struct.type.bondFlags(['aromatic'])
                }),
                0
            ])
        })
    },

    carbohydrate: {
        '@desc': ''
    },
    ions: {
        '@desc': '(specifically the PDB designations "PO4" and "SO4")'
    },
    ligand: {
        '@desc': '(originally "hetero and not solvent"; changed to "!(protein,nucleic,water,UREA)" for Jmol 12.2)'
    },
    nucleic: {
        '@desc': 'any group that (a) has one of the following group names: G, C, A, T, U, I, DG, DC, DA, DT, DU, DI, +G, +C, +A, +T, +U, +I; or (b) can be identified as a group that is only one atom, with name "P"; or (c) has all of the following atoms (prime, \', can replace * here): C1*, C2*, C3*, O3*, C4*, C5*, and O5*.',
        map: () => nucleicExpr()
    },
    purine: {
        '@desc': 'any nucleic group that (a) has one of the following group names: A, G, I, DA, DG, DI, +A, +G, or +I; or (b) also has atoms N7, C8, and N9.',
        map: () => B.struct.modifier.intersectBy({
            0: nucleicExpr(),
            by: B.struct.combinator.merge([
                B.struct.generator.atomGroups({
                    'residue-test': B.core.set.has([
                        B.set(...['A', 'G', 'I', 'DA', 'DG', 'DI', '+A', '+G', '+I']),
                        B.ammp('label_comp_id')
                    ])
                }),
                B.struct.filter.pick({
                    0: B.struct.generator.atomGroups({
                        'group-by': B.ammp('residueKey')
                    }),
                    test: B.core.set.isSubset([
                        h.atomNameSet(['N7', 'C8', 'N9']),
                        B.ammpSet('label_atom_id')
                    ])
                })
            ])
        })
    },
    pyrimidine: {
        '@desc': 'any nucleic group that (a) has one of the following group names: C, T, U, DC, DT, DU, +C, +T, +U; or (b) also has atom O2.',
        map: () => B.struct.modifier.intersectBy({
            0: nucleicExpr(),
            by: B.struct.combinator.merge([
                B.struct.generator.atomGroups({
                    'residue-test': B.core.set.has([
                        B.set(...['C', 'T', 'U', 'DC', 'DT', 'DU', '+C', '+T', '+U']),
                        B.ammp('label_comp_id')
                    ])
                }),
                B.struct.filter.pick({
                    0: B.struct.generator.atomGroups({
                        'group-by': B.ammp('residueKey')
                    }),
                    test: B.core.logic.or([
                        B.core.set.has([
                            B.ammpSet('label_atom_id'),
                            B.atomName('O2*')
                        ]),
                        B.core.set.has([
                            B.ammpSet('label_atom_id'),
                            B.atomName("O2'")
                        ])
                    ])
                })
            ])
        })
    },
    dna: {
        '@desc': 'any nucleic group that (a) has one of the following group names: DG, DC, DA, DT, DU, DI, T, +G, +C, +A, +T; or (b) has neither atom O2* or O2\'.',
        map: () => B.struct.modifier.intersectBy({
            0: nucleicExpr(),
            by: B.struct.combinator.merge([
                B.struct.generator.atomGroups({
                    'residue-test': B.core.set.has([
                        B.set(...['DG', 'DC', 'DA', 'DT', 'DU', 'DI', 'T', '+G', '+C', '+A', '+T']),
                        B.ammp('label_comp_id')
                    ])
                }),
                B.struct.filter.pick({
                    0: B.struct.generator.atomGroups({
                        'group-by': B.ammp('residueKey')
                    }),
                    test: B.core.logic.not([
                        B.core.logic.or([
                            B.core.set.has([
                                B.ammpSet('label_atom_id'),
                                B.atomName('O2*')
                            ]),
                            B.core.set.has([
                                B.ammpSet('label_atom_id'),
                                B.atomName("O2'")
                            ])
                        ])
                    ])
                })
            ])
        })
    },
    rna: {
        '@desc': 'any nucleic group that (a) has one of the following group names: G, C, A, U, I, +U, +I; or (b) has atom O2* or O2\'.',
        map: () => B.struct.modifier.intersectBy({
            0: nucleicExpr(),
            by: B.struct.combinator.merge([
                B.struct.generator.atomGroups({
                    'residue-test': B.core.set.has([
                        B.set(...['G', 'C', 'A', 'U', 'I', '+U', '+I']),
                        B.ammp('label_comp_id')
                    ])
                }),
                B.struct.filter.pick({
                    0: B.struct.generator.atomGroups({
                        'group-by': B.ammp('residueKey')
                    }),
                    test: B.core.logic.or([
                        B.core.set.has([
                            B.ammpSet('label_atom_id'),
                            B.atomName('O2*')
                        ]),
                        B.core.set.has([
                            B.ammpSet('label_atom_id'),
                            B.atomName("O2'")
                        ])
                    ])
                })
            ])
        })
    },
    protein: {
        '@desc': 'defined as a group that (a) has one of the following group names: ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, ASX, GLX, or UNK; or (b) contains PDB atom designations [C, O, CA, and N] bonded correctly; or (c) does not contain "O" but contains [C, CA, and N] bonded correctly; or (d) has only one atom, which has name CA and does not have the group name CA (indicating a calcium atom).',
        map: () => proteinExpr()
    },
    acidic: {
        '@desc': 'ASP GLU',
        map: () => h.resnameExpr(ResDict.acidic)
    },
    acyclic: {
        '@desc': 'amino and not cyclic',
        map: () => B.struct.modifier.intersectBy({
            0: h.resnameExpr(ResDict.amino),
            by: h.invertExpr(h.resnameExpr(ResDict.cyclic))
        })
    },
    aliphatic: {
        '@desc': 'ALA GLY ILE LEU VAL',
        map: () => h.resnameExpr(ResDict.aliphatic)
    },
    amino: {
        '@desc': 'all twenty standard amino acids, plus ASX, GLX, UNK',
        map: () => h.resnameExpr(ResDict.amino)
    },
    aromatic: {
        '@desc': 'HIS PHE TRP TYR (see also "isaromatic" for aromatic bonds)',
        map: () => h.resnameExpr(ResDict.aromatic)
    },
    basic: {
        '@desc': 'ARG HIS LYS',
        map: () => h.resnameExpr(ResDict.basic)
    },
    buried: {
        '@desc': 'ALA CYS ILE LEU MET PHE TRP VAL',
        map: () => h.resnameExpr(ResDict.buried)
    },
    charged: {
        '@desc': 'same as acidic or basic -- ASP GLU, ARG HIS LYS',
        map: () => h.resnameExpr(ResDict.acidic.concat(ResDict.basic))
    },
    cyclic: {
        '@desc': 'HIS PHE PRO TRP TYR',
        map: () => h.resnameExpr(ResDict.cyclic)
    },
    helix: {
        '@desc': 'secondary structure-related.',
        map: () => B.struct.generator.atomGroups({
            'residue-test': B.core.flags.hasAny([
                B.struct.type.secondaryStructureFlags(['helix']),
                B.ammp('secondaryStructureFlags')
            ])
        })
    },
    helixalpha: {
        '@desc': 'secondary structure-related.',
        map: () => B.struct.generator.atomGroups({
            'residue-test': B.core.flags.hasAny([
                B.struct.type.secondaryStructureFlags(['alpha']),
                B.ammp('secondaryStructureFlags')
            ])
        })
    },
    helix310: {
        '@desc': 'secondary structure-related.',
        map: () => B.struct.generator.atomGroups({
            'residue-test': B.core.flags.hasAny([
                B.struct.type.secondaryStructureFlags(['3-10']),
                B.ammp('secondaryStructureFlags')
            ])
        })
    },
    helixpi: {
        '@desc': 'secondary structure-related.',
        map: () => B.struct.generator.atomGroups({
            'residue-test': B.core.flags.hasAny([
                B.struct.type.secondaryStructureFlags(['pi']),
                B.ammp('secondaryStructureFlags')
            ])
        })
    },
    hetero: {
        '@desc': 'PDB atoms designated as HETATM',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.ammp('isHet')
        })
    },
    hydrophobic: {
        '@desc': 'ALA GLY ILE LEU MET PHE PRO TRP TYR VAL',
        map: () => h.resnameExpr(ResDict.hydrophobic)
    },
    large: {
        '@desc': 'ARG GLU GLN HIS ILE LEU LYS MET PHE TRP TYR',
        map: () => h.resnameExpr(ResDict.large)
    },
    medium: {
        '@desc': 'ASN ASP CYS PRO THR VAL',
        map: () => h.resnameExpr(ResDict.medium)
    },
    negative: {
        '@desc': 'same as acidic -- ASP GLU',
        map: () => h.resnameExpr(ResDict.acidic)
    },
    neutral: {
        '@desc': 'amino and not (acidic or basic)',
        map: () => B.struct.modifier.intersectBy({
            0: h.resnameExpr(ResDict.amino),
            by: h.invertExpr(h.resnameExpr(ResDict.acidic.concat(ResDict.basic)))
        })
    },
    polar: {
        '@desc': 'amino and not hydrophobic',
        map: () => B.struct.modifier.intersectBy({
            0: h.resnameExpr(ResDict.amino),
            by: h.invertExpr(h.resnameExpr(ResDict.hydrophobic))
        })
    },
    positive: {
        '@desc': 'same as basic -- ARG HIS LYS',
        map: () => h.resnameExpr(ResDict.basic)
    },
    sheet: {
        '@desc': 'secondary structure-related',
        map: () => B.struct.generator.atomGroups({
            'residue-test': B.core.flags.hasAny([
                B.struct.type.secondaryStructureFlags(['sheet']),
                B.ammp('secondaryStructureFlags')
            ])
        })
    },
    small: {
        '@desc': 'ALA GLY SER',
        map: () => h.resnameExpr(ResDict.small)
    },
    surface: {
        '@desc': 'amino and not buried',
        map: () => B.struct.modifier.intersectBy({
            0: h.resnameExpr(ResDict.amino),
            by: h.invertExpr(h.resnameExpr(ResDict.buried))
        })
    },
    turn: {
        '@desc': 'secondary structure-related',
        map: () => B.struct.generator.atomGroups({
            'residue-test': B.core.flags.hasAny([
                B.struct.type.secondaryStructureFlags(['turn']),
                B.ammp('secondaryStructureFlags')
            ])
        })
    },
    alpha: {
        '@desc': '(*.CA)',
        map: () => B.struct.generator.atomGroups({
            'atom-test': B.core.rel.eq([
                B.atomName('CA'),
                B.ammp('label_atom_id')
            ])
        })
    },
    base: {
        '@desc': '(nucleic bases)'
    },
    backbone: {
        '@desc': '(*.C, *.CA, *.N, and all nucleic other than the bases themselves)',
        abbr: ['mainchain'],
        map: () => backboneExpr()
    },
    sidechain: {
        '@desc': '((protein or nucleic) and not backbone)'
    },
    spine: {
        '@desc': '(*.CA, *.N, *.C for proteins; *.P, *.O3\', *.O5\', *.C3\', *.C4\', *.C5 for nucleic acids)'
    },
    leadatom: {
        '@desc': '(*.CA, *.P, and terminal *.O5\')'
    },
    solvent: {
        '@desc': 'PDB "HOH", water, also the connected set of H-O-H in any model'
    },
};


