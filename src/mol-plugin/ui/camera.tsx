/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from 'mol-plugin/command';
import * as React from 'react';
import { PluginComponent } from './base';

export class CameraSnapshots extends PluginComponent<{ }, { }> {
    render() {
        return <div>
            <h3>Camera Snapshots</h3>
            <CameraSnapshotControls />
            <CameraSnapshotList />
        </div>;
    }
}

class CameraSnapshotControls extends PluginComponent<{ }, { name: string, description: string }> {
    state = { name: '', description: '' };

    add = () => {
        PluginCommands.Camera.Snapshots.Add.dispatch(this.plugin, this.state);
        this.setState({ name: '', description: '' })
    }

    clear = () => {
        PluginCommands.Camera.Snapshots.Clear.dispatch(this.plugin, {});
    }

    render() {
        return <div>
            <input type='text' value={this.state.name} placeholder='Name...' style={{ width: '33%', display: 'block', float: 'left' }} onChange={e => this.setState({ name: e.target.value })} />
            <input type='text' value={this.state.description} placeholder='Description...' style={{ width: '67%', display: 'block' }} onChange={e => this.setState({ description: e.target.value })} />
            <button style={{ float: 'right' }} onClick={this.clear}>Clear</button>
            <button onClick={this.add}>Add</button>
        </div>;
    }
}

class CameraSnapshotList extends PluginComponent<{ }, { }> {
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
        return <ul style={{ listStyle: 'none' }}>
            {this.plugin.state.cameraSnapshots.entries.valueSeq().map(e =><li key={e!.id}>
                <button onClick={this.apply(e!.id)}>Set</button>
                &nbsp;{e!.name} <small>{e!.description}</small>
                <button onClick={this.remove(e!.id)} style={{ float: 'right' }}>X</button>
            </li>)}
        </ul>;
    }
}