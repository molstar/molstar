/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { formatStructureSelectionStats } from '../../util/structure-element-selection';
import { StructureSelectionQueries } from '../../util/structure-selection-helper';
import { ButtonSelect, Options } from '../controls/common';
import { PluginCommands } from '../../command';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Interactivity } from '../../util/interactivity';
import { ParameterControls } from '../controls/parameters';
import { camelCaseToWords } from '../../../mol-util/string';

type SelectionModifier = 'add' | 'remove' | 'only'

const StructureSelectionParams = {
    granularity: Interactivity.Params.granularity,
}

export class StructureSelectionControls extends PluginUIComponent<{}, {}> {
    state = {}

    componentDidMount() {
        this.subscribe(this.plugin.events.interactivity.selectionUpdated, () => {
            this.forceUpdate()
        });

        this.subscribe(this.plugin.events.interactivity.propsUpdated, () => {
            this.forceUpdate()
        });
    }

    get stats() {
        return formatStructureSelectionStats(this.plugin.helpers.structureSelectionManager.stats)
    }

    setProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        if (p.name === 'granularity') {
            PluginCommands.Interactivity.SetProps.dispatch(this.plugin, { props: { granularity: p.value } });
        }
    }

    get values () {
        return {
            granularity: this.plugin.interactivity.props.granularity
        }
    }

    set = (modifier: SelectionModifier, value: string) => {
        const query = StructureSelectionQueries[value as keyof typeof StructureSelectionQueries]
        this.plugin.helpers.structureSelection.set(modifier, query)
    }

    add = (value: string) => this.set('add', value)
    remove = (value: string) => this.set('remove', value)
    only = (value: string) => this.set('only', value)

    render() {
        const queries = Object.keys(StructureSelectionQueries).map(name => {
            return [name, camelCaseToWords(name)] as [string, string]
        })

        return <div className='msp-transform-wrapper'>
            <div className='msp-transform-header'>
                <button className='msp-btn msp-btn-block'>Selection</button>
            </div>
            <div>
                <div className='msp-control-row msp-row-text'>
                    <div>{this.stats}</div>
                </div>
                <ParameterControls params={StructureSelectionParams} values={this.values} onChange={this.setProps} />
                <div className='msp-control-row'>
                    <div className='msp-select-row'>
                        <ButtonSelect label='Add' onChange={this.add}>
                            <optgroup label='Add'>
                                {Options(queries)}
                            </optgroup>
                        </ButtonSelect>
                        <ButtonSelect label='Remove' onChange={this.remove}>
                            <optgroup label='Remove'>
                                {Options(queries)}
                            </optgroup>
                        </ButtonSelect>
                        <ButtonSelect label='Only' onChange={this.only}>
                            <optgroup label='Only'>
                                {Options(queries)}
                            </optgroup>
                        </ButtonSelect>
                    </div>
                </div>
            </div>
        </div>
    }
}