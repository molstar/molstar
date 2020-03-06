/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { StateSelection, StateBuilder } from '../../mol-state';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { QueryContext, StructureSelection, StructureQuery, StructureElement, Structure } from '../../mol-model/structure';
import { compile } from '../../mol-script/runtime/query/compiler';
import { Loci } from '../../mol-model/loci';
import { PluginContext } from '../context';
import Expression from '../../mol-script/language/expression';
import { BondType, ProteinBackboneAtoms, NucleicBackboneAtoms, SecondaryStructureType } from '../../mol-model/structure/model/types';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { SetUtils } from '../../mol-util/set';
import { ValidationReport, ValidationReportProvider } from '../../mol-model-props/rcsb/validation-report';
import { CustomProperty } from '../../mol-model-props/common/custom-property';
import { Task } from '../../mol-task';
import { AccessibleSurfaceAreaSymbols, AccessibleSurfaceAreaProvider } from '../../mol-model-props/computed/accessible-surface-area';
import { stringToWords } from '../../mol-util/string';

export enum StructureSelectionCategory {
    Type = 'Type',
    Structure = 'Structure Property',
    Atom = 'Atom Property',
    Bond = 'Bond Property',
    Residue = 'Residue Property',
    AminoAcid = 'Amino Acid',
    NucleicBase = 'Nucleic Base',
    Manipulate = 'Manipulate Selection',
    Validation = 'Validation',
    Misc = 'Miscellaneous',
    Internal = 'Internal',
}

export { StructureSelectionQuery }

interface StructureSelectionQuery {
    readonly label: string
    readonly expression: Expression
    readonly description: string
    readonly category: string
    readonly isHidden: boolean
    readonly query: StructureQuery
    readonly ensureCustomProperties?: (ctx: CustomProperty.Context, structure: Structure) => Promise<void>
}

interface StructureSelectionQueryProps {
    description?: string,
    category?: string
    isHidden?: boolean
    ensureCustomProperties?: (ctx: CustomProperty.Context, structure: Structure) => Promise<void>
}

function StructureSelectionQuery(label: string, expression: Expression, props: StructureSelectionQueryProps = {}): StructureSelectionQuery {
    let _query: StructureQuery
    return {
        label,
        expression,
        description: props.description || '',
        category: props.category ?? StructureSelectionCategory.Misc,
        isHidden: !!props.isHidden,
        get query() {
            if (!_query) _query = compile<StructureSelection>(expression)
            return _query
        },
        ensureCustomProperties: props.ensureCustomProperties
    }
}

const all = StructureSelectionQuery('All', MS.struct.generator.all(), { category: '' })

const polymer = StructureSelectionQuery('Polymer', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.logic.and([
            MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
            MS.core.str.match([
                MS.re('(polypeptide|cyclic-pseudo-peptide|nucleotide|peptide nucleic acid)', 'i'),
                MS.ammp('entitySubtype')
            ])
        ])
    })
]), { category: StructureSelectionCategory.Type })

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
]), { category: StructureSelectionCategory.Structure })

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
                'atom-test': MS.core.set.has([MS.set(...SetUtils.toArray(ProteinBackboneAtoms)), MS.ammp('label_atom_id')])
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
                'atom-test': MS.core.set.has([MS.set(...SetUtils.toArray(NucleicBackboneAtoms)), MS.ammp('label_atom_id')])
            })
        ])
    ])
]), { category: StructureSelectionCategory.Structure })

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
]), { category: StructureSelectionCategory.Type })

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
]), { category: StructureSelectionCategory.Type })

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
]), { category: StructureSelectionCategory.Residue })

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
]), { category: StructureSelectionCategory.Residue })

const water = StructureSelectionQuery('Water', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'water'])
    })
]), { category: StructureSelectionCategory.Type })

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
]), { category: StructureSelectionCategory.Type })

const branchedPlusConnected = StructureSelectionQuery('Carbohydrate with Connected', MS.struct.modifier.union([
    MS.struct.modifier.includeConnected({
        0: branched.expression, 'layer-count': 1, 'as-whole-residues': true
    })
]), { category: StructureSelectionCategory.Internal })

const branchedConnectedOnly = StructureSelectionQuery('Connected to Carbohydrate', MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: branchedPlusConnected.expression,
        by: branched.expression
    })
]), { category: StructureSelectionCategory.Internal })

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
]), { category: StructureSelectionCategory.Type })

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
]), { category: StructureSelectionCategory.Internal })

const ligandConnectedOnly = StructureSelectionQuery('Connected to Ligand', MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: ligandPlusConnected.expression,
        by: ligand.expression
    })
]), { category: StructureSelectionCategory.Internal })

// residues connected to ligands or branched entities
const connectedOnly = StructureSelectionQuery('Connected to Ligand or Carbohydrate', MS.struct.modifier.union([
    MS.struct.combinator.merge([
        branchedConnectedOnly.expression,
        ligandConnectedOnly.expression
    ]),
]), { category: StructureSelectionCategory.Internal })

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
]), { category: StructureSelectionCategory.Bond })

const nonStandardPolymer = StructureSelectionQuery('Non-standard Residues in Polymers', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
        'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
        'residue-test': MS.ammp('isNonStandard')
    })
]), { category: StructureSelectionCategory.Residue })

const coarse = StructureSelectionQuery('Coarse Elements', MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'chain-test': MS.core.set.has([
            MS.set('sphere', 'gaussian'), MS.ammp('objectPrimitive')
        ])
    })
]), { category: StructureSelectionCategory.Type })

const ring = StructureSelectionQuery('Rings in Residues', MS.struct.modifier.union([
    MS.struct.generator.rings()
]), { category: StructureSelectionCategory.Residue })

const aromaticRing = StructureSelectionQuery('Aromatic Rings in Residues', MS.struct.modifier.union([
    MS.struct.generator.rings({ 'only-aromatic': true })
]), { category: StructureSelectionCategory.Residue })

const surroundings = StructureSelectionQuery('Surrounding Residues (5 \u212B) of Selection', MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: MS.struct.modifier.includeSurroundings({
            0: MS.internal.generator.current(),
            radius: 5,
            'as-whole-residues': true
        }),
        by: MS.internal.generator.current()
    })
]), {
    description: 'Select residues within 5 \u212B of the current selection.',
    category: StructureSelectionCategory.Manipulate
})

const complement = StructureSelectionQuery('Inverse / Complement of Selection', MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: MS.struct.generator.all(),
        by: MS.internal.generator.current()
    })
]), {
    description: 'Select everything not in the current selection.',
    category: StructureSelectionCategory.Manipulate
})

const bonded = StructureSelectionQuery('Residues Bonded to Selection', MS.struct.modifier.union([
    MS.struct.modifier.includeConnected({
        0: MS.internal.generator.current(), 'layer-count': 1, 'as-whole-residues': true
    })
]), {
    description: 'Select residues covalently bonded to current selection.',
    category: StructureSelectionCategory.Manipulate
})

const hasClash = StructureSelectionQuery('Residues with Clashes', MS.struct.modifier.union([
    MS.struct.modifier.wholeResidues([
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
                'atom-test': ValidationReport.symbols.hasClash.symbol(),
            })
        ])
    ])
]), {
    description: 'Select residues with clashes in the wwPDB validation report.',
    category: StructureSelectionCategory.Residue,
    ensureCustomProperties: (ctx, structure) => {
        return ValidationReportProvider.attach(ctx, structure.models[0])
    }
})

const isBuried = StructureSelectionQuery('Buried Protein Residues', MS.struct.modifier.union([
    MS.struct.modifier.wholeResidues([
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
                'residue-test': AccessibleSurfaceAreaSymbols.isBuried.symbol(),
            })
        ])
    ])
]), {
    description: 'Select buried protein residues.',
    category: StructureSelectionCategory.Residue,
    ensureCustomProperties: (ctx, structure) => {
        return AccessibleSurfaceAreaProvider.attach(ctx, structure)
    }
})

const isAccessible = StructureSelectionQuery('Accessible Protein Residues', MS.struct.modifier.union([
    MS.struct.modifier.wholeResidues([
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
                'residue-test': AccessibleSurfaceAreaSymbols.isAccessible.symbol(),
            })
        ])
    ])
]), {
    description: 'Select accessible protein residues.',
    category: StructureSelectionCategory.Residue,
    ensureCustomProperties: (ctx, structure) => {
        return AccessibleSurfaceAreaProvider.attach(ctx, structure)
    }
})

const StandardAminoAcids = [
    [['HIS'], 'HISTIDINE'],
    [['ARG'], 'ARGININE'],
    [['LYS'], 'LYSINE'],
    [['ILE'], 'ISOLEUCINE'],
    [['PHE'], 'PHENYLALANINE'],
    [['LEU'], 'LEUCINE'],
    [['TRP'], 'TRYPTOPHAN'],
    [['ALA'], 'ALANINE'],
    [['MET'], 'METHIONINE'],
    [['CYS'], 'CYSTEINE'],
    [['ASN'], 'ASPARAGINE'],
    [['VAL'], 'VALINE'],
    [['GLY'], 'GLYCINE'],
    [['SER'], 'SERINE'],
    [['GLN'], 'GLUTAMINE'],
    [['TYR'], 'TYROSINE'],
    [['ASP'], 'ASPARTIC ACID'],
    [['GLU'], 'GLUTAMIC ACID'],
    [['THR'], 'THREONINE'],
    [['SEC'], 'SELENOCYSTEINE'],
    [['PYL'], 'PYRROLYSINE'],
].sort((a, b) => a[1] < b[1] ? -1 : a[1] > b[1] ? 1 : 0) as [string[], string][]

const StandardNucleicBases = [
    [['A', 'DA'], 'ADENOSINE'],
    [['C', 'DC'], 'CYTIDINE'],
    [['T', 'DT'], 'THYMIDINE'],
    [['G', 'DG'], 'GUANOSINE'],
    [['I', 'DI'], 'INOSINE'],
    [['U', 'DU'], 'URIDINE'],
].sort((a, b) => a[1] < b[1] ? -1 : a[1] > b[1] ? 1 : 0) as [string[], string][]

function ResidueQuery([names, label]: [string[], string], category: string) {
    return StructureSelectionQuery(`${stringToWords(label)} (${names.join(', ')})`, MS.struct.modifier.union([
        MS.struct.generator.atomGroups({
            'residue-test': MS.core.set.has([MS.set(...names), MS.ammp('auth_comp_id')])
        })
    ]), { category })
}

export const StructureSelectionQueries = {
    all,
    polymer,
    trace,
    backbone,
    protein,
    nucleic,
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
    nonStandardPolymer,
    coarse,
    ring,
    aromaticRing,
    surroundings,
    complement,
    bonded,

    hasClash,
    isBuried,
    isAccessible
}

export const StructureSelectionQueryList = [
    ...Object.values(StructureSelectionQueries),
    ...StandardAminoAcids.map(v => ResidueQuery(v, StructureSelectionCategory.AminoAcid)),
    ...StandardNucleicBases.map(v => ResidueQuery(v, StructureSelectionCategory.NucleicBase)),
]

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

    async set(modifier: SelectionModifier, selectionQuery: StructureSelectionQuery, applyGranularity = true) {
        this.plugin.runTask(Task.create('Structure Selection', async runtime => {
            const ctx = { fetch: this.plugin.fetch, runtime }
            for (const s of this.structures) {
                const current = this.plugin.helpers.structureSelectionManager.get(s)
                const currentSelection = Loci.isEmpty(current)
                    ? StructureSelection.Empty(s)
                    : StructureSelection.Singletons(s, StructureElement.Loci.toStructure(current))

                if (selectionQuery.ensureCustomProperties) {
                    await selectionQuery.ensureCustomProperties(ctx, s)
                }

                const result = selectionQuery.query(new QueryContext(s, { currentSelection }))
                const loci = StructureSelection.toLociWithSourceUnits(result)
                this._set(modifier, loci, applyGranularity)
            }
        }))
    }

    constructor(private plugin: PluginContext) {

    }
}