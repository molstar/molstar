/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject } from '../../mol-plugin/state/objects';
import { StateTransforms } from '../../mol-plugin/state/transforms';
import { StateTransformer, StateSelection, StateObjectCell, StateTransform } from '../../mol-state';
import { StructureElement, Structure, StructureSelection, QueryContext } from '../../mol-model/structure';
import { PluginContext } from '../context';
import { StructureRepresentation3DHelpers } from '../state/transforms/representation';
import Expression from '../../mol-script/language/expression';
import { compile } from '../../mol-script/runtime/query/compiler';

type StructureTransform = StateObjectCell<PluginStateObject.Molecule.Structure, StateTransform<StateTransformer<any, PluginStateObject.Molecule.Structure, any>>>
const RepresentationManagerTag = 'representation-controls'

function getRepresentationManagerTag(type: string) {
    return `${RepresentationManagerTag}-${type}`
}

function getCombinedLoci(mode: SelectionModifier, loci: StructureElement.Loci, currentLoci: StructureElement.Loci): StructureElement.Loci {
    switch (mode) {
        case 'add': return StructureElement.Loci.union(loci, currentLoci)
        case 'remove': return StructureElement.Loci.subtract(currentLoci, loci)
        case 'only': return loci
    }
}

type SelectionModifier = 'add' | 'remove' | 'only'

export class StructureRepresentationHelper {
    private async _set(modifier: SelectionModifier, type: string, loci: StructureElement.Loci, structure: StructureTransform) {
        const state = this.plugin.state.dataState
        const update = state.build();
        const s = structure.obj!.data

        const selections = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure, structure.transform.ref).withTag(getRepresentationManagerTag(type)));

        if (selections.length > 0) {
            const currentLoci = StructureElement.Query.toLoci(selections[0].params!.values.query, s)
            const combinedLoci = getCombinedLoci(modifier, loci, currentLoci)

            update.to(selections[0]).update({
                ...selections[0].params!.values,
                query: StructureElement.Query.fromLoci(combinedLoci)
            })
        } else {
            const combinedLoci = getCombinedLoci(modifier, loci, StructureElement.Loci(s, []))

            update.to(structure.transform.ref)
                .apply(
                    StateTransforms.Model.LociStructureSelection,
                    {
                        query: StructureElement.Query.fromLoci(combinedLoci),
                        label: type
                    },
                    { tags: [ RepresentationManagerTag, getRepresentationManagerTag(type) ] }
                )
                .apply(
                    StateTransforms.Representation.StructureRepresentation3D,
                    StructureRepresentation3DHelpers.getDefaultParams(this.plugin, type as any, s)
                )
        }

        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }));
    }

    async set(modifier: SelectionModifier, type: string, lociGetter: (structure: Structure) => StructureElement.Loci) {
        const state = this.plugin.state.dataState;
        const structures = state.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure));

        for (const structure of structures) {
            const s = structure.obj!.data
            const loci = lociGetter(s)
            await this._set(modifier, type, loci, structure)
        }
    }

    async setFromExpression(modifier: SelectionModifier, type: string, expression: Expression) {
        return this.set(modifier, type, (structure) => {
            const compiled = compile<StructureSelection>(expression)
            const result = compiled(new QueryContext(structure))
            return StructureSelection.toLoci2(result)
        })
    }

    constructor(private plugin: PluginContext) {

    }
}