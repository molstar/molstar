/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject } from '../../mol-plugin/state/objects';
import { StateTransforms } from '../../mol-plugin/state/transforms';
import { StateTransformer, StateSelection, StateObjectCell, StateTransform } from '../../mol-state';
import { StructureElement } from '../../mol-model/structure';
import { isEmptyLoci } from '../../mol-model/loci';
import { PluginContext } from '../context';
import { parseMolScript } from '../../mol-script/language/parser';
import { StructureRepresentation3DHelpers } from '../state/transforms/representation';
import Expression from '../../mol-script/language/expression';
import { formatMolScript } from '../../mol-script/language/expression-formatter';
import { MolScriptBuilder as MS } from '../../mol-script/language/builder';

type StructureTransform = StateObjectCell<PluginStateObject.Molecule.Structure, StateTransform<StateTransformer<any, PluginStateObject.Molecule.Structure, any>>>
const RepresentationManagerTag = 'representation-controls'

function getRepresentationManagerTag(type: string) {
    return `${RepresentationManagerTag}-${type}`
}

function getCombinedExpression(modifier: SelectionModifier, expression: Expression, currentExpression: Expression): Expression {
    switch (modifier) {
        case 'add': return MS.struct.combinator.merge([ currentExpression, expression ])
        case 'remove': return MS.struct.modifier.exceptBy({ 0: currentExpression, by: expression })
        case 'only': return expression
        case 'all': return MS.struct.generator.all()
    }
}

type SelectionModifier = 'add' | 'remove' | 'only' | 'all'

export class StructureRepresentationHelper {
    async set(modifier: SelectionModifier, type: string, expression: Expression, structure: StructureTransform) {
        const state = this.plugin.state.dataState
        const update = state.build();
        const s = structure.obj!.data

        const selections = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure, structure.transform.ref).withTag(getRepresentationManagerTag(type)));

        if (selections.length > 0) {
            const parsedExpressions = parseMolScript(selections[0].params!.values.query.expression)
            if (parsedExpressions.length === 0) return
            const currentExpression = parsedExpressions[0]
            const combinedExpression = getCombinedExpression(modifier, expression, currentExpression)

            update.to(selections[0]).update({
                ...selections[0].params!.values,
                query: { language: 'mol-script', expression: formatMolScript(combinedExpression) }
            })
        } else {
            const combinedExpression = getCombinedExpression(modifier, expression, MS.struct.generator.empty())

            update.to(structure.transform.ref)
                .apply(
                    StateTransforms.Model.UserStructureSelection,
                    {
                        query: { language: 'mol-script', expression: formatMolScript(combinedExpression) },
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

    async setSelected(modifier: SelectionModifier, type: string) {
        const state = this.plugin.state.dataState;
        const structures = state.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure));

        for (const structure of structures) {
            const s = structure.obj!.data
            const _loci = this.plugin.helpers.structureSelectionManager.get(s)
            const loci = isEmptyLoci(_loci) ? StructureElement.Loci(s, []) : _loci
            const expression = StructureElement.Loci.toScriptExpression(loci)

            await this.set(modifier, type, expression, structure)
        }
    }

    async hideAll(type: string) {
        const state = this.plugin.state.dataState;
        const update = state.build();

        state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure).withTag(getRepresentationManagerTag(type))).forEach(structure => update.delete(structure.transform.ref));

        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }));
    }

    constructor(private plugin: PluginContext) {

    }
}