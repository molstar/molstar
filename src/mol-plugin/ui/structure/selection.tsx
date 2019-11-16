/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { CollapsableControls, CollapsableState } from '../base';
import { StructureSelectionQueries, SelectionModifier } from '../../util/structure-selection-helper';
import { ButtonSelect, Options } from '../controls/common';
import { PluginCommands } from '../../command';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Interactivity } from '../../util/interactivity';
import { ParameterControls } from '../controls/parameters';
import { stripTags } from '../../../mol-util/string';

const SSQ = StructureSelectionQueries
const DefaultQueries: (keyof typeof SSQ)[] = [
    'all', 'polymer', 'trace', 'protein', 'nucleic', 'water', 'branched', 'ligand', 'nonStandardPolymer',
    'surroundings', 'complement', 'bonded'
]

const StructureSelectionParams = {
    granularity: Interactivity.Params.granularity,
}

interface StructureSelectionControlsState extends CollapsableState {
    minRadius: number,
    extraRadius: number,
    durationMs: number,

    isDisabled: boolean
}

export class StructureSelectionControls<P, S extends StructureSelectionControlsState> extends CollapsableControls<P, S> {
    componentDidMount() {
        this.subscribe(this.plugin.events.interactivity.selectionUpdated, () => {
            this.forceUpdate()
        });

        this.subscribe(this.plugin.events.interactivity.propsUpdated, () => {
            this.forceUpdate()
        });

        this.subscribe(this.plugin.state.dataState.events.isUpdating, v => this.setState({ isDisabled: v }))
    }

    get stats() {
        const stats = this.plugin.helpers.structureSelectionManager.stats
        if (stats.structureCount === 0 || stats.elementCount === 0) {
            return 'Selected nothing'
        } else {
            return `Selected ${stripTags(stats.label)}`
        }
    }

    focus = () => {
        const { extraRadius, minRadius, durationMs } = this.state
        if (this.plugin.helpers.structureSelectionManager.stats.elementCount === 0) return
        const { sphere } = this.plugin.helpers.structureSelectionManager.getBoundary();
        const radius = Math.max(sphere.radius + extraRadius, minRadius);
        this.plugin.canvas3d.camera.focus(sphere.center, radius, durationMs);
    }

    setProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        if (p.name === 'granularity') {
            PluginCommands.Interactivity.SetProps.dispatch(this.plugin, { props: { granularity: p.value } });
        }
    }

    get values () {
        return {
            granularity: this.plugin.interactivity.props.granularity,
        }
    }

    set = (modifier: SelectionModifier, value: string) => {
        const query = SSQ[value as keyof typeof SSQ]
        this.plugin.helpers.structureSelection.set(modifier, query.query, false)
    }

    add = (value: string) => this.set('add', value)
    remove = (value: string) => this.set('remove', value)
    only = (value: string) => this.set('only', value)

    defaultState() {
        return {
            isCollapsed: false,
            header: 'Selection',

            minRadius: 8,
            extraRadius: 4,
            durationMs: 250,

            isDisabled: false
        } as S
    }

    renderControls() {
        const queries = Object.keys(StructureSelectionQueries)
            .map(name => [name, SSQ[name as keyof typeof SSQ].label] as [string, string])
            .filter(pair => DefaultQueries.includes(pair[0] as keyof typeof SSQ))

        return <div>
            <div className='msp-control-row msp-row-text'>
                <button className='msp-btn msp-btn-block' onClick={this.focus}>
                    <span className={`msp-icon msp-icon-focus-on-visual`} style={{ position: 'absolute', left: '5px' }} />
                    {this.stats}
                </button>
            </div>
            <ParameterControls params={StructureSelectionParams} values={this.values} onChange={this.setProps} isDisabled={this.state.isDisabled} />
            <div className='msp-control-row'>
                <div className='msp-select-row'>
                    <ButtonSelect label='Select' onChange={this.add} disabled={this.state.isDisabled}>
                        <optgroup label='Select'>
                            {Options(queries)}
                        </optgroup>
                    </ButtonSelect>
                    <ButtonSelect label='Deselect' onChange={this.remove} disabled={this.state.isDisabled}>
                        <optgroup label='Deselect'>
                            {Options(queries)}
                        </optgroup>
                    </ButtonSelect>
                    <ButtonSelect label='Only' onChange={this.only} disabled={this.state.isDisabled}>
                        <optgroup label='Only'>
                            {Options(queries)}
                        </optgroup>
                    </ButtonSelect>
                </div>
            </div>
        </div>
    }
}