/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { PluginStateObject } from '../../../mol-plugin/state/objects';
import { StateTransforms } from '../../../mol-plugin/state/transforms';
import { StateSelection, StateObjectCell, StateTransform, StateBuilder } from '../../../mol-state';
import { ParamDefinition as PD} from '../../../mol-util/param-definition';
import { ColorNames } from '../../../mol-util/color/tables';
import { ParameterControls } from '../controls/parameters';
import { Structure } from '../../../mol-model/structure';
import { isEmptyLoci } from '../../../mol-model/loci';
import { PluginContext } from '../../context';
import { getExpression } from './util';


type OverpaintEachReprCallback = (update: StateBuilder.Root, repr: StateObjectCell<PluginStateObject.Molecule.Structure.Representation3D, StateTransform<typeof StateTransforms.Representation.StructureRepresentation3D>>, rootStructure: Structure, overpaint?: StateObjectCell<any, StateTransform<typeof StateTransforms.Representation.OverpaintStructureRepresentation3D>>) => void
const OverpaintManagerTag = 'overpaint-controls'

export class StructureOverpaintControls extends PluginUIComponent<{}, { params: PD.Values<ReturnType<typeof StructureOverpaintControls.getParams>> }> {
    state = { params: PD.getDefaultValues(StructureOverpaintControls.getParams(this.plugin)) }

    static getParams = (plugin: PluginContext) => {
        const { types } = plugin.structureRepresentation.registry
        return {
            color: PD.Color(ColorNames.cyan),
            type: PD.MultiSelect(types.map(t => t[0]), types)
        }
    }

    componentDidMount() {

    }

    private async eachRepr(callback: OverpaintEachReprCallback) {
        const state = this.plugin.state.dataState;
        const reprs = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure.Representation3D));

        const update = state.build();
        for (const r of reprs) {
            const overpaint = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.OverpaintStructureRepresentation3D, r.transform.ref).withTag(OverpaintManagerTag));

            const structure = r.obj!.data.source.data
            const rootStructure = structure.parent || structure

            callback(update, r, rootStructure, overpaint[0])
        }

        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }));
    }

    set = async (clear: boolean) => {
        await this.eachRepr((update, repr, rootStructure, overpaint) => {
            if (!this.state.params.type.includes(repr.params!.values.type.name)) return

            const loci = this.plugin.helpers.structureSelection.get(rootStructure)
            if (isEmptyLoci(loci) || loci.elements.length === 0) return
            const expression = getExpression(loci)

            const layer = {
                script: { language: 'mol-script', expression },
                color: this.state.params.color,
                clear
            }

            if (overpaint) {
                update.to(overpaint).update({ layers: [ ...overpaint.params!.values.layers, layer ], alpha: 1 })
            } else {
                update.to(repr.transform.ref)
                    .apply(StateTransforms.Representation.OverpaintStructureRepresentation3D, { layers: [ layer ], alpha: 1 }, { tags: OverpaintManagerTag });
            }
        })
    }

    add = async () => {
        this.set(false)
    }

    clear = async () => {
        this.set(true)
    }

    clearAll = async () => {
        await this.eachRepr((update, repr, rootStructure, overpaint) => {
            if (overpaint) update.delete(overpaint.transform.ref)
        })
    }

    render() {
        return <div className='msp-transform-wrapper'>
            <div className='msp-transform-header'>
                <button className='msp-btn msp-btn-block'>Current Selection Overpaint</button>
            </div>
            <div>
                <ParameterControls params={StructureOverpaintControls.getParams(this.plugin)} values={this.state.params} onChange={p => {
                    const params = { ...this.state.params, [p.name]: p.value };
                    this.setState({ params });
                }}/>

                <div className='msp-btn-row-group'>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={this.add}>Add</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={this.clear}>Clear</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={this.clearAll}>Clear All</button>
                </div>
            </div>
        </div>
    }
}