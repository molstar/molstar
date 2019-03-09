/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from 'mol-plugin/command';
import * as React from 'react';
import { PluginUIComponent } from './base';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ParameterControls } from './controls/parameters';
import { Icon } from './controls/common';

export class CameraSnapshots extends PluginUIComponent<{ }, { }> {
    render() {
        return <div>
            <div className='msp-section-header'><Icon name='code' /> Camera Snapshots</div>
            <CameraSnapshotControls />
            <CameraSnapshotList />
        </div>;
    }
}

class CameraSnapshotControls extends PluginUIComponent<{ }, { name: string, description: string }> {
    static Params = {
        name: PD.Text(),
        description: PD.Text()
    }
    state = PD.getDefaultValues(CameraSnapshotControls.Params);

    add = () => {
        PluginCommands.Camera.Snapshots.Add.dispatch(this.plugin, this.state);
        this.setState({ name: '', description: '' })
    }

    clear = () => {
        PluginCommands.Camera.Snapshots.Clear.dispatch(this.plugin, {});
    }

    render() {
        return <div>
            <ParameterControls params={CameraSnapshotControls.Params} values={this.state} onEnter={this.add} onChange={p => this.setState({ [p.name]: p.value } as any)}  />

            <div className='msp-btn-row-group'>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.add}>Add</button>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.clear}>Clear</button>
            </div>
        </div>;
    }
}

class CameraSnapshotList extends PluginUIComponent<{ }, { }> {
    componentDidMount() {
        this.subscribe(this.plugin.events.state.cameraSnapshots.changed, () => this.forceUpdate());
    }

    apply(id: string) {
        return () => PluginCommands.Camera.Snapshots.Apply.dispatch(this.plugin, { id });
    }

    remove(id: string) {
        return () => {
            PluginCommands.Camera.Snapshots.Remove.dispatch(this.plugin, { id });
        }
    }

    render() {
        return <ul style={{ listStyle: 'none' }} className='msp-state-list'>
            {this.plugin.state.cameraSnapshots.state.entries.valueSeq().map(e =><li key={e!.id}>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.apply(e!.id)}>{e!.name || e!.timestamp} <small>{e!.description}</small></button>
                <button onClick={this.remove(e!.id)} className='msp-btn msp-btn-link msp-state-list-remove-button'>
                    <span className='msp-icon msp-icon-remove' />
                </button>
            </li>)}
        </ul>;
    }
}