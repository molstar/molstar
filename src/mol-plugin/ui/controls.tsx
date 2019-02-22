/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginCommands } from 'mol-plugin/command';
import { UpdateTrajectory } from 'mol-plugin/state/actions/structure';
import { PluginUIComponent } from './base';
import { LociLabelEntry } from 'mol-plugin/util/loci-label-manager';

export class Controls extends PluginUIComponent<{ }, { }> {
    render() {
        return <>

        </>;
    }
}

export class TrajectoryControls extends PluginUIComponent {
    render() {
        return <div>
            <button className='msp-btn msp-btn-link' onClick={() => PluginCommands.State.ApplyAction.dispatch(this.plugin, {
                state: this.plugin.state.dataState,
                action: UpdateTrajectory.create({ action: 'advance', by: -1 })
            })} title='Previou Model'>◀</button>
            <button className='msp-btn msp-btn-link' onClick={() => PluginCommands.State.ApplyAction.dispatch(this.plugin, {
                state: this.plugin.state.dataState,
                action: UpdateTrajectory.create({ action: 'reset' })
            })} title='First Model'>↻</button>
            <button className='msp-btn msp-btn-link' onClick={() => PluginCommands.State.ApplyAction.dispatch(this.plugin, {
                state: this.plugin.state.dataState,
                action: UpdateTrajectory.create({ action: 'advance', by: +1 })
            })} title='Next Model'>►</button><br />
        </div>
    }
}

export class LociLabelControl extends PluginUIComponent<{}, { entries: ReadonlyArray<LociLabelEntry> }> {
    state = { entries: [] }

    componentDidMount() {
        this.subscribe(this.plugin.events.labels.highlight, e => this.setState({ entries: e.entries }));
    }

    render() {
        return <div style={{ textAlign: 'right' }}>
            {this.state.entries.map((e, i) => <div key={'' + i}>{e}</div>)}
        </div>
    }
}