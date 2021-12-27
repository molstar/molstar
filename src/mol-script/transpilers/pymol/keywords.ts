/*
 * Copyright (c) 2017 MolQL contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as h from '../helper'
import { KeywordDict } from '../types'
import B from '../../molql/builder'

const ResDict = {
  nucleic: ['A', 'C', 'T', 'G', 'U', 'DA', 'DC', 'DT', 'DG', 'DU'],
  protein: ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'CYX', 'GLN', 'GLU', 'GLY', 'HIS', 'HID', 'HIE', 'HIP', 'ILE', 'LEU', 'LYS', 'MET', 'MSE', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'],
  solvent: ['HOH', 'WAT', 'H20', 'TIP', 'SOL']
}

const keywords: KeywordDict = {
  all: {
    '@desc': 'All atoms currently loaded into PyMOL',
    abbr: ['*'],
    map: () => B.struct.generator.atomGroups()
  },
  none: {
    '@desc': 'No atoms (empty selection)',
    map: () => B.struct.generator.empty()
  },
  hydrogens: {
    '@desc': 'All hydrogen atoms currently loaded into PyMOL',
    abbr: ['hydro', 'h.'],
    map: () => B.struct.generator.atomGroups({
      'atom-test': B.core.rel.eq([
        B.acp('elementSymbol'),
        B.es('H')
      ])
    })
  },
  hetatm: {
    '@desc': 'All atoms loaded from Protein Data Bank HETATM records',
    abbr: ['het'],
    map: () => B.struct.generator.atomGroups({
      'atom-test': B.core.rel.eq([B.ammp('isHet'), true])
    })
  },
  visible: {
    '@desc': 'All atoms in enabled objects with at least one visible representation',
    abbr: ['v.']
  },
  polymer: {
    '@desc': 'All atoms on the polymer (not het). Finds atoms with residue identifiers matching a known polymer, such a peptide and DNA.',
    abbr: ['pol.'],
    map: () => B.struct.generator.atomGroups({
      'residue-test': B.core.set.has([
        B.core.type.set(ResDict.nucleic.concat(ResDict.protein)),
        B.ammp('label_comp_id')
      ])
    })
  },
  backbone: {
    '@desc': 'Polymer backbone atoms (new in PyMOL 1.6.1)',
    abbr: ['bb.']
  },
  sidechain: {
    '@desc': 'Polymer non-backbone atoms (new in PyMOL 1.6.1)',
    abbr: ['sc.']
  },
  present: {
    '@desc': 'All atoms with defined coordinates in the current state (used in creating movies)',
    abbr: ['pr.']
  },
  center: {
    '@desc': 'Pseudo-atom at the center of the scene'
  },
  origin: {
    '@desc': 'Pseudo-atom at the origin of rotation',
  },
  enabled: {
    '@desc': 'All enabled objects or selections from the object list.',
  },
  masked: {
    '@desc': 'All masked atoms.',
    abbr: ['msk.']
  },
  protected: {
    '@desc': 'All protected atoms.',
    abbr: ['pr.']
  },
  bonded: {
    '@desc': 'All bonded atoms',
    map: () => B.struct.generator.atomGroups({
      'atom-test': B.core.rel.gr([B.struct.atomProperty.core.bondCount({
        flags: B.struct.type.bondFlags(['covalent', 'metallic', 'sulfide'])
      }), 0])
    })
  },
  donors: {
    '@desc': 'All hydrogen bond donor atoms.',
    abbr: ['don.']
  },
  acceptors: {
    '@desc': 'All hydrogen bond acceptor atoms.',
    abbr: ['acc.']
  },
  fixed: {
    '@desc': 'All fixed atoms.',
    abbr: ['fxd.']
  },
  restrained: {
    '@desc': 'All restrained atoms.',
    abbr: ['rst.']
  },
  organic: {
    '@desc': 'All atoms in non-polymer organic compounds (e.g. ligands, buffers). Finds carbon-containing molecules that do not match known polymers.',
    abbr: ['org.'],
    map: () => h.asAtoms(B.struct.modifier.expandProperty({
      '0': B.struct.modifier.union([
        B.struct.generator.queryInSelection({
          '0': B.struct.generator.atomGroups({
            'residue-test': B.core.logic.not([
              B.core.set.has([
                B.core.type.set(ResDict.nucleic.concat(ResDict.protein)),
                B.ammp('label_comp_id')
              ])
            ])
          }),
          query: B.struct.generator.atomGroups({
            'atom-test': B.core.rel.eq([
              B.es('C'),
              B.acp('elementSymbol')
            ])
          })
        })
      ]),
      property: B.ammp('residueKey')
    }))
  },
  inorganic: {
    '@desc': 'All non-polymer inorganic atoms/ions. Finds atoms in molecules that do not contain carbon and do not match any known solvent residues.',
    abbr: ['ino.'],
    map: () => h.asAtoms(B.struct.modifier.expandProperty({
      '0': B.struct.modifier.union([
        B.struct.filter.pick({
          '0': B.struct.generator.atomGroups({
            'residue-test': B.core.logic.not([
              B.core.set.has([
                B.core.type.set(ResDict.nucleic.concat(ResDict.protein).concat(ResDict.solvent)),
                B.ammp('label_comp_id')
              ])
            ]),
            'group-by': B.ammp('residueKey')
          }),
          test: B.core.logic.not([
            B.core.set.has([
              B.struct.atomSet.propertySet([ B.acp('elementSymbol') ]),
              B.es('C')
            ])
          ])
        })
      ]),
      property: B.ammp('residueKey')
    }))
  },
  solvent: {
    '@desc': 'All water molecules. The hardcoded solvent residue identifiers are currently: HOH, WAT, H20, TIP, SOL.',
    abbr: ['sol.'],
    map: () => B.struct.generator.atomGroups({
      'residue-test': B.core.set.has([
        B.core.type.set(ResDict.solvent),
        B.ammp('label_comp_id')
      ])
    })
  },
  guide: {
    '@desc': 'All protein CA and nucleic acid C4*/C4',
    map: () => B.struct.combinator.merge([
      B.struct.generator.atomGroups({
        'atom-test': B.core.rel.eq([
          B.atomName('CA'),
          B.ammp('label_atom_id')
        ]),
        'residue-test': B.core.set.has([
          B.core.type.set(ResDict.protein),
          B.ammp('label_comp_id')
        ])
      }),
      B.struct.generator.atomGroups({
        'atom-test': B.core.set.has([
          h.atomNameSet(['C4*', 'C4']),
          B.ammp('label_atom_id')
        ]),
        'residue-test': B.core.set.has([
          B.core.type.set(ResDict.nucleic),
          B.ammp('label_comp_id')
        ])
      })
    ]),
  },
  metals: {
    '@desc': 'All metal atoms (new in PyMOL 1.6.1)'
  }
}

export default keywords
