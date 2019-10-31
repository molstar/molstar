/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { StateSelection } from '../../mol-state';
import { PluginStateObject } from '../state/objects';
import { QueryContext, StructureSelection, StructureQuery, StructureElement } from '../../mol-model/structure';
import { compile } from '../../mol-script/runtime/query/compiler';
import { Loci } from '../../mol-model/loci';
import { PluginContext } from '../context';
import Expression from '../../mol-script/language/expression';
import { LinkType } from '../../mol-model/structure/model/types';

export interface StructureSelectionQuery {
    label: string
    query: StructureQuery
    expression: Expression
    description: string
}

export function StructureSelectionQuery(label: string, expression: Expression, description = ''): StructureSelectionQuery {
    return { label, expression, query: compile<StructureSelection>(expression), description }
}

const all = StructureSelectionQuery('All', MS.struct.generator.all())

const polymer = StructureSelectionQuery('Polymer', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer'])
    })
]))

const trace = StructureSelectionQuery('Trace', MS.struct.modifier.union([
    MS.struct.combinator.merge([
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
                'chain-test': MS.core.set.has([
                    MS.set('sphere', 'gaussian'), MS.ammp('objectPrimitive')
                ])
            })
        ]),
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
                'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
                'atom-test': MS.core.set.has([MS.set('CA', 'P'), MS.ammp('label_atom_id')])
            })
        ])
    ])
]))

const protein = StructureSelectionQuery('Protein', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.logic.and([
            MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
            MS.core.str.match([
                MS.re('(polypeptide|cyclic-pseudo-peptide)', 'i'),
                MS.ammp('entitySubtype')
            ])
        ])
    })
]))

const nucleic = StructureSelectionQuery('Nucleic', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.logic.and([
            MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
            MS.core.str.match([
                MS.re('(nucleotide|peptide nucleic acid)', 'i'),
                MS.ammp('entitySubtype')
            ])
        ])
    })
]))

const proteinAndNucleic = StructureSelectionQuery('Protein | Nucleic', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.logic.and([
            MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
            MS.core.str.match([
                MS.re('(polypeptide|cyclic-pseudo-peptide|nucleotide|peptide nucleic acid)', 'i'),
                MS.ammp('entitySubtype')
            ])
        ])
    })
]))

const water = StructureSelectionQuery('Water', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'water'])
    })
]))

const branched = StructureSelectionQuery('Carbohydrate', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.logic.or([
            MS.core.rel.eq([MS.ammp('entityType'), 'branched']),
            MS.core.logic.and([
                MS.core.rel.eq([MS.ammp('entityType'), 'non-polymer']),
                MS.core.str.match([
                    MS.re('oligosaccharide', 'i'),
                    MS.ammp('entitySubtype')
                ])
            ])
        ])
    })
]))

const branchedPlusConnected = StructureSelectionQuery('Carbohydrate with Connected', MS.struct.modifier.union([
    MS.struct.modifier.includeConnected({
        0: branched.expression, 'layer-count': 1, 'as-whole-residues': true
    })
]))

const branchedConnectedOnly = StructureSelectionQuery('Connected to Carbohydrate', MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: branchedPlusConnected.expression,
        by: branched.expression
    })
]))

const ligand = StructureSelectionQuery('Ligand', MS.struct.modifier.union([
    MS.struct.combinator.merge([
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'entity-test': MS.core.logic.and([
                    MS.core.rel.eq([MS.ammp('entityType'), 'non-polymer']),
                    MS.core.logic.not([MS.core.str.match([
                        MS.re('oligosaccharide', 'i'),
                        MS.ammp('entitySubtype')
                    ])])
                ]),
                'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
                'residue-test': MS.core.logic.not([
                    MS.core.str.match([MS.re('saccharide', 'i'), MS.ammp('chemCompType')])
                ])
            })
        ]),
        // this is to get non-polymer components in polymer entities,
        // e.g. PXZ in 4HIV or generally ACE
        //
        // one option to optimize this is to expose `_entity_poly.nstd_monomer` in molql
        // and only check those entities residue by residue
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
                'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
                'residue-test': MS.core.str.match([MS.re('non-polymer', 'i'), MS.ammp('chemCompType')])
            })
        ])
    ]),
]))

// don't include branched entities as they have their own link representation
const ligandPlusConnected = StructureSelectionQuery('Ligand with Connected', MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: MS.struct.modifier.union([
            MS.struct.modifier.includeConnected({
                0: ligand.expression,
                'layer-count': 1,
                'as-whole-residues': true,
                'link-test': MS.core.flags.hasAny([
                    MS.struct.linkProperty.flags(),
                    MS.core.type.bitflags([
                        LinkType.Flag.Covalent | LinkType.Flag.MetallicCoordination
                    ])
                ])
            })
        ]),
        by: branched.expression
    })
]))

const ligandConnectedOnly = StructureSelectionQuery('Connected to Ligand', MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: ligandPlusConnected.expression,
        by: ligand.expression
    })
]))

// residues connected to ligands or branched entities
const connectedOnly = StructureSelectionQuery('Connected to Ligand | Carbohydrate', MS.struct.modifier.union([
    MS.struct.combinator.merge([
        branchedConnectedOnly.expression,
        ligandConnectedOnly.expression
    ]),
]))

const modified = StructureSelectionQuery('Modified Residues', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
        'residue-test': MS.ammp('isModified')
    })
]))

const nonStandardPolymer = StructureSelectionQuery('Non-standard Residues in Polymers', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
        'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
        'residue-test': MS.ammp('isNonStandard')
    })
]))

const coarse = StructureSelectionQuery('Coarse Elements', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'chain-test': MS.core.set.has([
            MS.set('sphere', 'gaussian'), MS.ammp('objectPrimitive')
        ])
    })
]))

const surroundings = StructureSelectionQuery('Surrounding Residues (5 \u212B)', MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: MS.struct.modifier.includeSurroundings({
            0: MS.internal.generator.current(),
            radius: 5,
            'as-whole-residues': true
        }),
        by: MS.internal.generator.current()
    })
]), 'Select residues within 5 \u212B of the current selection.')

const complement = StructureSelectionQuery('Inverse / Complement', MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: MS.struct.generator.all(),
        by: MS.internal.generator.current()
    })
]), 'Select everything not in the current selection.')

export const StructureSelectionQueries = {
    all,
    polymer,
    trace,
    protein,
    nucleic,
    proteinAndNucleic,
    water,
    branched,
    branchedPlusConnected,
    branchedConnectedOnly,
    ligand,
    ligandPlusConnected,
    ligandConnectedOnly,
    connectedOnly,
    modified,
    nonStandardPolymer,
    coarse,
    surroundings,
    complement,
}

//

export type SelectionModifier = 'add' | 'remove' | 'only'

export class StructureSelectionHelper {
    private get structures() {
        return this.plugin.state.dataState.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure)).map(s => s.obj!.data)
    }

    private _set(modifier: SelectionModifier, loci: Loci) {
        switch (modifier) {
            case 'add':
                this.plugin.interactivity.lociSelects.select({ loci })
                break
            case 'remove':
                this.plugin.interactivity.lociSelects.deselect({ loci })
                break
            case 'only':
                this.plugin.interactivity.lociSelects.selectOnly({ loci })
                break
        }
    }

    set(modifier: SelectionModifier, query: StructureQuery) {
        for (const s of this.structures) {
            const current = this.plugin.helpers.structureSelectionManager.get(s)
            const currentSelection = Loci.isEmpty(current)
                ? StructureSelection.Empty(s)
                : StructureSelection.Singletons(s, StructureElement.Loci.toStructure(current))

            const result = query(new QueryContext(s, { currentSelection }))
            const loci = StructureSelection.toLociWithSourceUnits(result)
            this._set(modifier, loci)
        }
    }

    constructor(private plugin: PluginContext) {

    }
}