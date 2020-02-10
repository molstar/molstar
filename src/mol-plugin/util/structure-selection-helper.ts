/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { StateSelection, StateBuilder } from '../../mol-state';
import { PluginStateObject } from '../state/objects';
import { QueryContext, StructureSelection, StructureQuery, StructureElement } from '../../mol-model/structure';
import { compile } from '../../mol-script/runtime/query/compiler';
import { Loci } from '../../mol-model/loci';
import { PluginContext } from '../context';
import Expression from '../../mol-script/language/expression';
import { BondType, ProteinBackboneAtoms, NucleicBackboneAtoms, SecondaryStructureType } from '../../mol-model/structure/model/types';
import { StateTransforms } from '../state/transforms';

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

// TODO maybe pre-calculate atom properties like backbone/sidechain
const backbone = StructureSelectionQuery('Backbone', MS.struct.modifier.union([
    MS.struct.combinator.merge([
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'entity-test': MS.core.logic.and([
                    MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
                    MS.core.str.match([
                        MS.re('(polypeptide|cyclic-pseudo-peptide)', 'i'),
                        MS.ammp('entitySubtype')
                    ])
                ]),
                'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
                'atom-test': MS.core.set.has([MS.set(...Array.from(ProteinBackboneAtoms.values())), MS.ammp('label_atom_id')])
            })
        ]),
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'entity-test': MS.core.logic.and([
                    MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
                    MS.core.str.match([
                        MS.re('(nucleotide|peptide nucleic acid)', 'i'),
                        MS.ammp('entitySubtype')
                    ])
                ]),
                'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
                'atom-test': MS.core.set.has([MS.set(...Array.from(NucleicBackboneAtoms.values())), MS.ammp('label_atom_id')])
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

const proteinOrNucleic = StructureSelectionQuery('Protein or Nucleic', MS.struct.modifier.union([
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

const helix = StructureSelectionQuery('Helix', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.logic.and([
            MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
            MS.core.str.match([
                MS.re('(polypeptide|cyclic-pseudo-peptide)', 'i'),
                MS.ammp('entitySubtype')
            ])
        ]),
        'residue-test': MS.core.flags.hasAny([
            MS.ammp('secondaryStructureFlags'),
            MS.core.type.bitflags([SecondaryStructureType.Flag.Helix])
        ])
    })
]))

const beta = StructureSelectionQuery('Beta Strand/Sheet', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.logic.and([
            MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
            MS.core.str.match([
                MS.re('(polypeptide|cyclic-pseudo-peptide)', 'i'),
                MS.ammp('entitySubtype')
            ])
        ]),
        'residue-test': MS.core.flags.hasAny([
            MS.ammp('secondaryStructureFlags'),
            MS.core.type.bitflags([SecondaryStructureType.Flag.Beta])
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
        // this is to get non-polymer and peptide terminus components in polymer entities,
        // - non-polymer, e.g. PXZ in 4HIV or generally ACE
        // - carboxy terminus, e.g. FC0 in 4BP9, or ETA in 6DDE
        // - amino terminus, e.g. ARF in 3K4V, or 4MM in 3EGV
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
                'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
                'residue-test': MS.core.str.match([
                    MS.re('non-polymer|(amino|carboxy) terminus', 'i'),
                    MS.ammp('chemCompType')
                ])
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
                'bond-test': MS.core.flags.hasAny([
                    MS.struct.bondProperty.flags(),
                    MS.core.type.bitflags([
                        BondType.Flag.Covalent | BondType.Flag.MetallicCoordination
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
const connectedOnly = StructureSelectionQuery('Connected to Ligand or Carbohydrate', MS.struct.modifier.union([
    MS.struct.combinator.merge([
        branchedConnectedOnly.expression,
        ligandConnectedOnly.expression
    ]),
]))

const disulfideBridges = StructureSelectionQuery('Disulfide Bridges', MS.struct.modifier.union([
    MS.struct.modifier.wholeResidues([
        MS.struct.modifier.union([
            MS.struct.generator.bondedAtomicPairs({
                0: MS.core.flags.hasAny([
                    MS.struct.bondProperty.flags(),
                    MS.core.type.bitflags([BondType.Flag.Disulfide])
                ])
            })
        ])
    ])
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

const ring = StructureSelectionQuery('Rings in Residues', MS.struct.modifier.union([
    MS.struct.generator.rings()
]))

const aromaticRing = StructureSelectionQuery('Aromatic Rings in Residues', MS.struct.modifier.union([
    MS.struct.generator.rings({ 'only-aromatic': true })
]))

const surroundings = StructureSelectionQuery('Surrounding Residues (5 \u212B) of Selection', MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: MS.struct.modifier.includeSurroundings({
            0: MS.internal.generator.current(),
            radius: 5,
            'as-whole-residues': true
        }),
        by: MS.internal.generator.current()
    })
]), 'Select residues within 5 \u212B of the current selection.')

const complement = StructureSelectionQuery('Inverse / Complement of Selection', MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: MS.struct.generator.all(),
        by: MS.internal.generator.current()
    })
]), 'Select everything not in the current selection.')

const bonded = StructureSelectionQuery('Residues Bonded to Selection', MS.struct.modifier.union([
    MS.struct.modifier.includeConnected({
        0: MS.internal.generator.current(), 'layer-count': 1, 'as-whole-residues': true
    })
]), 'Select residues covalently bonded to current selection.')

export const StructureSelectionQueries = {
    all,
    polymer,
    trace,
    backbone,
    protein,
    nucleic,
    proteinOrNucleic,
    helix,
    beta,
    water,
    branched,
    branchedPlusConnected,
    branchedConnectedOnly,
    ligand,
    ligandPlusConnected,
    ligandConnectedOnly,
    connectedOnly,
    disulfideBridges,
    modified,
    nonStandardPolymer,
    coarse,
    ring,
    aromaticRing,
    surroundings,
    complement,
    bonded,
}

export function applyBuiltInSelection(to: StateBuilder.To<PluginStateObject.Molecule.Structure>, query: keyof typeof StructureSelectionQueries, customTag?: string) {
    return to.apply(StateTransforms.Model.StructureSelectionFromExpression,
        { expression: StructureSelectionQueries[query].expression, label: StructureSelectionQueries[query].label },
        { tags: customTag ? [query, customTag] : [query] });
}

//

export type SelectionModifier = 'add' | 'remove' | 'only'

export class StructureSelectionHelper {
    private get structures() {
        return this.plugin.state.dataState.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure)).map(s => s.obj!.data)
    }

    private _set(modifier: SelectionModifier, loci: Loci, applyGranularity = true) {
        switch (modifier) {
            case 'add':
                this.plugin.interactivity.lociSelects.select({ loci }, applyGranularity)
                break
            case 'remove':
                this.plugin.interactivity.lociSelects.deselect({ loci }, applyGranularity)
                break
            case 'only':
                this.plugin.interactivity.lociSelects.selectOnly({ loci }, applyGranularity)
                break
        }
    }

    set(modifier: SelectionModifier, query: StructureQuery, applyGranularity = true) {
        for (const s of this.structures) {
            const current = this.plugin.helpers.structureSelectionManager.get(s)
            const currentSelection = Loci.isEmpty(current)
                ? StructureSelection.Empty(s)
                : StructureSelection.Singletons(s, StructureElement.Loci.toStructure(current))

            const result = query(new QueryContext(s, { currentSelection }))
            const loci = StructureSelection.toLociWithSourceUnits(result)
            this._set(modifier, loci, applyGranularity)
        }
    }

    constructor(private plugin: PluginContext) {

    }
}