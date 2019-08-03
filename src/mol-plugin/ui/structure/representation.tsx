/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { ParamDefinition as PD} from '../../../mol-util/param-definition';
import { ParameterControls } from '../controls/parameters';
import { PluginContext } from '../../context';

export class StructureRepresentationControls extends PluginUIComponent<{}, { params: PD.Values<ReturnType<typeof StructureRepresentationControls.getParams>> }> {
    state = { params: PD.getDefaultValues(StructureRepresentationControls.getParams(this.plugin)) }

    static getParams = (plugin: PluginContext) => {
        const { types } = plugin.structureRepresentation.registry
        return {
            type: PD.Select(types[0][0], types)
        }
    }

    componentDidMount() {

    }

    set = (mode: 'add' | 'remove' | 'only' | 'all') => {
        this.plugin.helpers.structureRepresentation.setSelected(mode, this.state.params.type)
    }

    show = () => { this.set('add') }
    hide = () => { this.set('remove') }
    only = () => { this.set('only') }
    showAll = () => { this.set('all') }

    hideAll = () => {
        this.plugin.helpers.structureRepresentation.hideAll(this.state.params.type)
    }

    render() {
        return <div className='msp-transform-wrapper'>
            <div className='msp-transform-header'>
                <button className='msp-btn msp-btn-block'>Current Selection Representation</button>
            </div>
            <div>
                <ParameterControls params={StructureRepresentationControls.getParams(this.plugin)} values={this.state.params} onChange={p => {
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