/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { StateSelection } from '../../mol-state';
import { PluginStateObject } from '../state/objects';
import { QueryContext, StructureSelection } from '../../mol-model/structure';
import { compile } from '../../mol-script/runtime/query/compiler';
import { Loci } from '../../mol-model/loci';
import { PluginContext } from '../context';
import Expression from '../../mol-script/language/expression';

const all = MS.struct.generator.all()

const polymers = MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer'])
    })
])

const backboneTrace = MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'atom-test': MS.core.logic.or([
            MS.core.rel.eq([MS.ammp('label_atom_id'), 'CA']),
            MS.core.rel.eq([MS.ammp('label_atom_id'), 'P'])
        ])
    })
])

const water = MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'water'])
    })
])

const branched = MS.struct.modifier.union([
    MS.struct.combinator.merge([
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'branched'])
            })
        ]),
        MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'non-polymer']),
                'residue-test': MS.core.str.match([MS.re('saccharide', 'i'), MS.ammp('chemCompType')])
            })
        ])
    ])
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

const ligands = MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
        'entity-test': MS.core.logic.and([
            MS.core.rel.neq([MS.ammp('entityType'), 'branched']),
            MS.core.rel.eq([MS.ammp('entityType'), 'non-polymer'])
        ]),
        'residue-test': MS.core.logic.not([
            MS.core.str.match([MS.re('saccharide', 'i'), MS.ammp('chemCompType')])
        ])
    })
])

const ligandsPlusConnected = MS.struct.modifier.union([
    MS.struct.modifier.includeConnected({
        0: ligands, 'layer-count': 1, 'as-whole-residues': true
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
    polymers,
    backboneTrace,
    water,
    branched,
    branchedPlusConnected,
    branchedConnectedOnly,
    ligands,
    ligandsPlusConnected,
    coarse,
}

//

type SelectionModifier = 'add' | 'remove' | 'only'

export class StructureSelectionHelper {
    private get structures() {
        const state = this.plugin.state.dataState
        return state.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure))
    }

    private _set(modifier: SelectionModifier, loci: Loci) {
        switch (modifier) {
            case 'add':
                this.plugin.interactivity.lociSelections.add({ loci })
                break
            case 'remove':
                this.plugin.interactivity.lociSelections.remove({ loci })
                break
            case 'only':
                this.plugin.interactivity.lociSelections.only({ loci })
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