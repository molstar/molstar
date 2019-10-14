/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { StateSelection } from '../../mol-state';
import { PluginStateObject } from '../state/objects';
import { QueryContext, StructureSelection, StructureQuery } from '../../mol-model/structure';
import { compile } from '../../mol-script/runtime/query/compiler';
import { Loci } from '../../mol-model/loci';
import { PluginContext } from '../context';
import Expression from '../../mol-script/language/expression';

const all = MS.struct.generator.all()

const polymer = MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer'])
    })
])

const trace = MS.struct.modifier.union([
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
])

const protein = MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.logic.and([
            MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
            MS.core.str.match([
                MS.re('(polypeptide|cyclic-pseudo-peptide)', 'i'),
                MS.ammp('entitySubtype')
            ])
        ])
    })
])

const nucleic = MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.logic.and([
            MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
            MS.core.str.match([
                MS.re('(nucleotide|peptide nucleic acid)', 'i'),
                MS.ammp('entitySubtype')
            ])
        ])
    })
])

const proteinAndNucleic = MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.logic.and([
            MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
            MS.core.str.match([
                MS.re('(polypeptide|cyclic-pseudo-peptide|nucleotide|peptide nucleic acid)', 'i'),
                MS.ammp('entitySubtype')
            ])
        ])
    })
])

const water = MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'water'])
    })
])

const branched = MS.struct.modifier.union([
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
])

const branchedPlusConnected = MS.struct.modifier.union([
    MS.struct.modifier.includeConnected({
        0: branched, 'layer-count': 1, 'as-whole-residues': true
    })
])

const branchedConnectedOnly = MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: branchedPlusConnected,
        by: branched
    })
])

const ligand = MS.struct.modifier.union([
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
        // this is to get non-polymer components in polymer entities, e.g. PXZ in 4HIV
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
])

// don't include branched entities as they have their own link representation
const ligandPlusConnected = MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: MS.struct.modifier.union([
            MS.struct.modifier.includeConnected({
                0: ligand,
                'layer-count': 1,
                'as-whole-residues': true
            })
        ]),
        by: branched
    })
])

const ligandConnectedOnly = MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: ligandPlusConnected,
        by: ligand
    })
])

// residues connected to ligands or branched entities
const connectedOnly = MS.struct.modifier.union([
    MS.struct.combinator.merge([
        branchedConnectedOnly,
        ligandConnectedOnly
    ]),
])

const modified = MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
        'residue-test': MS.ammp('isModified')
    })
])

const coarse = MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'chain-test': MS.core.set.has([
            MS.set('sphere', 'gaussian'), MS.ammp('objectPrimitive')
        ])
    })
])

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
    coarse,
}

export const CompiledStructureSelectionQueries = (function () {
    const ret: { [K in keyof typeof StructureSelectionQueries]: StructureQuery } = Object.create(null);
    for (const k of Object.keys(StructureSelectionQueries)) {
        (ret as any)[k] = compile<StructureSelection>((StructureSelectionQueries as any)[k]);
    }
    return ret;
})();

//

export type SelectionModifier = 'add' | 'remove' | 'only'

export class StructureSelectionHelper {
    private get structures() {
        const state = this.plugin.state.dataState
        return state.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure))
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

    set(modifier: SelectionModifier, query: Expression) {
        const compiled = compile<StructureSelection>(query)

        for (const so of this.structures) {
            const s = so.obj!.data
            const result = compiled(new QueryContext(s))
            const loci = StructureSelection.toLoci2(result)
            this._set(modifier, loci)
        }
    }

    constructor(private plugin: PluginContext) {

    }
}