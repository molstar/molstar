/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { ParamDefinition as PD} from '../../../mol-util/param-definition';
import { ColorNames } from '../../../mol-util/color/tables';
import { ParameterControls } from '../controls/parameters';
import { PluginContext } from '../../context';
import { Color } from '../../../mol-util/color';

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

    set = (color: Color | -1) => {
        this.plugin.helpers.structureOverpaint.set(color, this.state.params.type)
    }

    add = () => {
        this.set(this.state.params.color)
    }

    clear = () => {
        this.set(-1)
    }

    clearAll = () => {
        this.plugin.helpers.structureOverpaint.clearAll()
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