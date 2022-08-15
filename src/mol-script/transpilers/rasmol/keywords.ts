/**
 * Copyright (c) 2017-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
**/


import { MolScriptBuilder } from '../../../mol-script/language/builder';
const B = MolScriptBuilder;
import * as h from '../helper';
import { KeywordDict } from '../types';

const Backbone = {
    nucleic: ['P', "O3'", "O5'", "C5'", "C4'", "C3'", 'OP1', 'OP2', 'O3*', 'O5*', 'C5*', 'C4*', 'C3*'],
    protein: ['C', 'N', 'CA', 'O']
};


function nucleicExpr() {
    return B.struct.combinator.merge([
        B.struct.generator.atomGroups({
	  'residue-test': B.core.set.has([
                B.core.type.set(['G', 'C', 'A', 'T', 'U', 'I', 'DG', 'DC', 'DA', 'DT', 'DU', 'DI', '+G', '+C', '+A', '+T', '+U', '+I']),
                B.ammp('label_comp_id')
	  ])
        }),
        B.struct.filter.pick({
	    0: B.struct.generator.atomGroups({
                'group-by': B.ammp('residueKey')
	    }),
	    test: B.core.logic.and([
                B.core.set.isSubset([
                // B.core.type.set([ 'P', 'O1P', 'O2P' ]),
                    h.atomNameSet(['P']),
                    B.ammpSet('label_atom_id')
                ]),
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



const ResDict = {
    aliphatic: ['ALA', 'GLY', 'ILE', 'LEU', 'VAL'],
    amino: ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ASX', 'GLX', 'UNK'],
    acidic: ['ASP', 'GLU'],
    aromatic: ['HIS', 'PHE', 'TRP', 'TYR'],
    basic: ['ARG', 'HIS', 'LYS'],
    buried: ['ALA', 'CYS', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'VAL'],
    cg: ['CYT', 'C', 'GUA', 'G'],
    cyclic: ['HIS', 'PHE', 'PRO', 'TRP', 'TYR'],
    hydrophobic: ['ALA', 'GLY', 'ILE', 'LEU', 'MET', 'PHE', 'PRO', 'TRP', 'TYR', 'VAL'],
    large: ['ARG', 'GLU', 'GLN', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'TRP', 'TYR'],
    medium: ['ASN', 'ASP', 'CYS', 'PRO', 'THR', 'VAL'],
    small: ['ALA', 'GLY', 'SER'],
    nucleic: ['A', 'C', 'T', 'G', 'U', 'DA', 'DC', 'DT', 'DG', 'DU'],
    protein: ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'CYX', 'GLN', 'GLU', 'GLY', 'HIS', 'HID', 'HIE', 'HIP', 'ILE', 'LEU', 'LYS', 'MET', 'MSE', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'],
    solvent: ['HOH', 'WAT', 'H20', 'TIP', 'SOL']
};


export const keywords: KeywordDict = {
    // general terms
    all: {
        '@desc': 'all atoms; same as *',
        abbr: ['*'],
        map: () => B.struct.generator.all()
    },
    none: {
        '@desc': 'no atoms',
        map: () => B.struct.generator.empty()
    },
    selected: {
        '@desc': 'atoms that have been selected; defaults to all when a file is first loaded'
    },
    unitcell: {
        '@desc': 'atoms within the current unitcell, which may be offset. This includes atoms on the faces and at the vertices of the unitcell.'
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
        '@desc': 'defined as a group that (a) has one of the following group names: ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU}, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, ASX, GLX, or UNK; or (b) contains PDB atom designations [C, O, CA, and N] bonded correctly; or (c) does not contain "O" but contains [C, CA, and N] bonded correctly; or (d) has only one atom, which has name CA and does not have the group name CA (indicating a calcium atom).',
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
    solvent: {
        '@desc': 'PDB "HOH", water, also the connected set of H-O-H in any model',
	    map: () => h.resnameExpr(ResDict.solvent)
    },
};

function backboneExpr() {
    return B.struct.combinator.merge([
        B.struct.generator.queryInSelection({
            0: proteinExpr(),
            query: B.struct.generator.atomGroups({
                'atom-test': B.core.set.has([
                    h.atomNameSet(Backbone.protein),
                    B.ammp('label_atom_id')
                ])
            })
        }),
        B.struct.generator.queryInSelection({
            0: nucleicExpr(),
            query: B.struct.generator.atomGroups({
                'atom-test': B.core.set.has([
                    h.atomNameSet(Backbone.nucleic),
                    B.ammp('label_atom_id')
                ])
            })
        })
    ]);
}

function proteinExpr() {
    return B.struct.filter.pick({
        0: B.struct.generator.atomGroups({
            'group-by': B.ammp('residueKey')
        }),
        test: B.core.set.isSubset([
            h.atomNameSet(['C', 'N', 'CA', 'O']),
            B.ammpSet('label_atom_id')
        ])
    });
}

