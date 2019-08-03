/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { StructureSelection, QueryFn, Queries as _Queries } from '../../../mol-model/structure';
import { formatStructureSelectionStats } from '../../util/structure-element-selection';
import { StructureSelectionQueries } from '../../util/structure-selection-helper';

export class StructureSelectionControls extends PluginUIComponent<{}, {}> {
    state = {}

    componentDidMount() {
        this.subscribe(this.plugin.events.interactivity.selectionUpdated, () => {
            this.forceUpdate()
        });
    }

    get stats() {
        return formatStructureSelectionStats(this.plugin.helpers.structureSelectionManager.stats)
    }

    select = (query: QueryFn<StructureSelection>) => {
        this.plugin.helpers.structureSelection.select(query)
    }

    clear = () => {
        this.plugin.helpers.structureSelection.clearSelection()
    }

    render() {
        return <div className='msp-transform-wrapper'>
            <div className='msp-transform-header'>
                <button className='msp-btn msp-btn-block'>Current Selection</button>
            </div>
            <div>
                <div className='msp-control-row msp-row-text'>
                    <div>{this.stats}</div>
                </div>
                <div className='msp-btn-row-group'>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.select(StructureSelectionQueries.all())}>All</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.clear()}>None</button>
                </div>
                <div className='msp-btn-row-group'>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.select(StructureSelectionQueries.polymers())}>Polymers</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.select(StructureSelectionQueries.ligands())}>Ligands</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.select(StructureSelectionQueries.water())}>Water</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.select(StructureSelectionQueries.coarse())}>Coarse</button>
                </div>
            </div>
        </div>
    }
}