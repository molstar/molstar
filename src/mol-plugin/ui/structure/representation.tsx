/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { PluginStateObject } from '../../../mol-plugin/state/objects';
import { StateTransforms } from '../../../mol-plugin/state/transforms';
import { StateTransformer, StateSelection, StateObjectCell, StateTransform, StateBuilder } from '../../../mol-state';
import { ParamDefinition as PD} from '../../../mol-util/param-definition';
import { ParameterControls } from '../controls/parameters';
import { StructureElement, QueryContext, StructureSelection } from '../../../mol-model/structure';
import { isEmptyLoci } from '../../../mol-model/loci';
import { PluginContext } from '../../context';
import { getExpression } from './util';
import { parseMolScript } from '../../../mol-script/language/parser';
import { transpileMolScript } from '../../../mol-script/script/mol-script/symbols';
import { compile } from '../../../mol-script/runtime/query/compiler';
import { StructureRepresentation3DHelpers } from '../../state/transforms/representation';

type RepresentationEachStructureCallback = (update: StateBuilder.Root, structure: StateObjectCell<PluginStateObject.Molecule.Structure, StateTransform<StateTransformer<any, PluginStateObject.Molecule.Structure, any>>>) => void
const RepresentationManagerTag = 'representation-controls'

function getRepresentationManagerTag(type: string) {
    return `${RepresentationManagerTag}-${type}`
}

function getCombinedLoci(mode: 'add' | 'remove' | 'only' | 'all', loci: StructureElement.Loci, currentLoci: StructureElement.Loci): StructureElement.Loci {
    switch (mode) {
        case 'add': return StructureElement.Loci.union(loci, currentLoci)
        case 'remove': return StructureElement.Loci.subtract(currentLoci, loci)
        case 'only': return loci
        case 'all': return StructureElement.Loci.all(loci.structure)
    }
}

export class RepresentationControls extends PluginUIComponent<{}, { params: PD.Values<ReturnType<typeof RepresentationControls.getParams>> }> {
    state = { params: PD.getDefaultValues(RepresentationControls.getParams(this.plugin)) }

    static getParams = (plugin: PluginContext) => {
        const { types } = plugin.structureRepresentation.registry
        return {
            type: PD.Select(types[0][0], types)
        }
    }

    componentDidMount() {

    }

    private async eachStructure(callback: RepresentationEachStructureCallback) {
        const state = this.plugin.state.dataState;
        const structures = state.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure));

        const update = state.build();
        for (const s of structures) {
            callback(update, s)
        }

        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }));
    }

    set = async (mode: 'add' | 'remove' | 'only' | 'all') => {
        const state = this.plugin.state.dataState
        const { type } = this.state.params

        await this.eachStructure((update, structure) => {
            const s = structure.obj!.data
            const _loci = this.plugin.helpers.structureSelection.get(s)
            const loci = isEmptyLoci(_loci) ? StructureElement.Loci(s, []) : _loci

            const selections = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure, structure.transform.ref).withTag(getRepresentationManagerTag(type)));

            if (selections.length > 0) {
                const parsed = parseMolScript(selections[0].params!.values.query.expression)
                if (parsed.length === 0) return

                const query = transpileMolScript(parsed[0])
                const compiled = compile(query)
                const result = compiled(new QueryContext(structure.obj!.data))
                const currentLoci = StructureSelection.toLoci2(result)

                const combinedLoci = getCombinedLoci(mode, loci, currentLoci)

                update.to(selections[0]).update({
                    ...selections[0].params!.values,
                    query: { language: 'mol-script', expression: getExpression(combinedLoci) }
                })
            } else {
                const combinedLoci = getCombinedLoci(mode, loci, StructureElement.Loci(loci.structure, []))

                update.to(structure.transform.ref)
                    .apply(
                        StateTransforms.Model.UserStructureSelection,
                        {
                            query: { language: 'mol-script', expression: getExpression(combinedLoci) },
                            label: type
                        },
                        { tags: [ RepresentationManagerTag, getRepresentationManagerTag(type) ] }
                    )
                    .apply(
                        StateTransforms.Representation.StructureRepresentation3D,
                        StructureRepresentation3DHelpers.getDefaultParams(this.plugin, type as any, s)
                    )
            }
        })
    }

    show = async () => { this.set('add') }
    hide = async () => { this.set('remove') }
    only = async () => { this.set('only') }
    showAll = async () => { this.set('all') }

    hideAll = async () => {
        const { type } = this.state.params
        const state = this.plugin.state.dataState;
        const update = state.build();

        state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure).withTag(getRepresentationManagerTag(type))).forEach(structure => update.delete(structure.transform.ref));

        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }));
    }

    render() {
        return <div className='msp-transform-wrapper'>
            <div className='msp-transform-header'>
                <button className='msp-btn msp-btn-block'>Current Selection Representation</button>
            </div>
            <div>
                <ParameterControls params={RepresentationControls.getParams(this.plugin)} values={this.state.params} onChange={p => {
                    const params = { ...this.state.params, [p.name]: p.value };
                    this.setState({ params });
                }}/>

                <div className='msp-btn-row-group'>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={this.show}>Show</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={this.hide}>Hide</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={this.only}>Only</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={this.showAll}>Show All</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={this.hideAll}>Hide All</button>
                </div>
            </div>
        </div>
    }
}