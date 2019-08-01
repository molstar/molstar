/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { MolScriptBuilder as MS } from '../../../mol-script/language/builder';
import { StateSelection } from '../../../mol-state';
import { PluginStateObject } from '../../state/objects';
import { QueryContext, StructureSelection, QueryFn, Queries as _Queries } from '../../../mol-model/structure';
import { compile } from '../../../mol-script/runtime/query/compiler';
import { ButtonsType } from '../../../mol-util/input/input-observer';
import { EmptyLoci } from '../../../mol-model/loci';
import { formatStructureSelectionStats } from '../../util/structure-element-selection';

const Queries = {
    all: () => compile<StructureSelection>(MS.struct.generator.all()),
    polymers: () => _Queries.internal.atomicSequence(),
    water: () => _Queries.internal.water(),
    ligands: () => _Queries.internal.atomicHet(),
    coarse: () => _Queries.internal.spheres(),
}

export class StructureSelectionControls extends PluginUIComponent<{}, {}> {
    state = {}

    componentDidMount() {
        this.subscribe(this.plugin.events.interactivity.selectionUpdated, () => {
            this.forceUpdate()
        });
    }

    get stats() {
        return formatStructureSelectionStats(this.plugin.helpers.structureSelection.stats)
    }

    select = (query: QueryFn<StructureSelection>) => {
        const state = this.plugin.state.dataState
        const structures = state.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure))
        const { structureSelection } = this.plugin.helpers

        structureSelection.clear()
        for (const so of structures) {
            const s = so.obj!.data
            const result = query(new QueryContext(s))
            const loci = StructureSelection.toLoci2(result)

            // TODO use better API when available
            this.plugin.interactivity.lociSelections.apply({
                current: { loci },
                buttons: ButtonsType.Flag.Secondary,
                modifiers: { shift: false, alt: false, control: true, meta: false }
            })
        }
    }

    clear = () => {
        // TODO use better API when available
        this.plugin.interactivity.lociSelections.apply({
            current: { loci: EmptyLoci },
            buttons: ButtonsType.Flag.Secondary,
            modifiers: { shift: false, alt: false, control: true, meta: false }
        })
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
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.select(Queries.all())}>All</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.clear()}>None</button>
                </div>
                <div className='msp-btn-row-group'>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.select(Queries.polymers())}>Polymers</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.select(Queries.ligands())}>Ligands</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.select(Queries.water())}>Water</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.select(Queries.coarse())}>Coarse</button>
                </div>
            </div>
        </div>
    }
}