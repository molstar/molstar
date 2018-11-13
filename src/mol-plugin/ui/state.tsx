/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from 'mol-plugin/command';
import * as React from 'react';
import { PluginComponent } from './base';
import { shallowEqual } from 'mol-util';
import { List } from 'immutable';
import { LogEntry } from 'mol-util/log-entry';

export class StateSnapshots extends PluginComponent<{ }, { serverUrl: string }> {
    state = { serverUrl: 'http://webchem.ncbr.muni.cz/molstar-state' }

    updateServerUrl = (serverUrl: string) => { this.setState({ serverUrl }) };

    render() {
        return <div>
            <h3>State Snapshots</h3>
            <StateSnapshotControls serverUrl={this.state.serverUrl} serverChanged={this.updateServerUrl} />
            <b>Local</b>
            <LocalStateSnapshotList />
            <RemoteStateSnapshotList serverUrl={this.state.serverUrl} />
        </div>;
    }
}

class StateSnapshotControls extends PluginComponent<{ serverUrl: string, serverChanged: (url: string) => void }, { name: string, description: string, serverUrl: string, isUploading: boolean }> {
    state = { name: '', description: '', serverUrl: this.props.serverUrl, isUploading: false };

    add = () => {
        PluginCommands.State.Snapshots.Add.dispatch(this.plugin, { name: this.state.name, description: this.state.description });
        this.setState({ name: '', description: '' })
    }

    clear = () => {
        PluginCommands.State.Snapshots.Clear.dispatch(this.plugin, {});
    }

    shouldComponentUpdate(nextProps: { serverUrl: string, serverChanged: (url: string) => void }, nextState: { name: string, description: string, serverUrl: string, isUploading: boolean }) {
        return !shallowEqual(this.props, nextProps) || !shallowEqual(this.state, nextState);
    }

    upload = async () => {
        this.setState({ isUploading: true });
        await PluginCommands.State.Snapshots.Upload.dispatch(this.plugin, { name: this.state.name, description: this.state.description, serverUrl: this.state.serverUrl });
        this.setState({ isUploading: false });
    }

    render() {
        return <div>
            <input type='text' value={this.state.name} placeholder='Name...' style={{ width: '33%', display: 'block', float: 'left' }} onChange={e => this.setState({ name: e.target.value })} />
            <input type='text' value={this.state.description} placeholder='Description...' style={{ width: '67%', display: 'block' }} onChange={e => this.setState({ description: e.target.value })} />
            <input type='text' value={this.state.serverUrl} placeholder='Server URL...' style={{ width: '100%', display: 'block' }} onChange={e => {
                this.setState({ serverUrl: e.target.value });
                this.props.serverChanged(e.target.value);
            }} />
            <button style={{ float: 'right' }} onClick={this.clear}>Clear</button>
            <button onClick={this.add}>Add Local</button>
            <button onClick={this.upload} disabled={this.state.isUploading}>Upload</button>
        </div>;
    }
}

class LocalStateSnapshotList extends PluginComponent<{ }, { }> {
    componentDidMount() {
        this.subscribe(this.plugin.events.state.snapshots.changed, () => this.forceUpdate());
    }

    apply(id: string) {
        return () => PluginCommands.State.Snapshots.Apply.dispatch(this.plugin, { id });
    }

    remove(id: string) {
        return () => {
            PluginCommands.State.Snapshots.Remove.dispatch(this.plugin, { id });
        }
    }

    render() {
        return <ul style={{ listStyle: 'none' }}>
            {this.plugin.state.snapshots.entries.valueSeq().map(e =><li key={e!.id}>
                <button onClick={this.apply(e!.id)}>Set</button>
                &nbsp;{e!.name} <small>{e!.description}</small>
                <button onClick={this.remove(e!.id)} style={{ float: 'right' }}>X</button>
            </li>)}
        </ul>;
    }
}

type RemoteEntry = { url: string, timestamp: number, id: string, name: string, description: string }
class RemoteStateSnapshotList extends PluginComponent<{ serverUrl: string }, { entries: List<RemoteEntry>, isFetching: boolean }> {
    state = { entries: List<RemoteEntry>(), isFetching: false };

    componentDidMount() {
        this.subscribe(this.plugin.events.state.snapshots.changed, () => this.forceUpdate());
        this.refresh();
    }

    refresh = async () => {
        try {
            this.setState({ isFetching: true });
            const req = await fetch(`${this.props.serverUrl}/list`);
            const json: RemoteEntry[] = await req.json();
            this.setState({ entries: List<RemoteEntry>(json.map((e: RemoteEntry) => ({ ...e, url: `${this.props.serverUrl}/get/${e.id}` }))), isFetching: false })
        } catch (e) {
            this.plugin.log(LogEntry.error('Fetching Remote Snapshots: ' + e));
            this.setState({ entries: List<RemoteEntry>(), isFetching: false })
        }
    }

    fetch(url: string) {
        return () => PluginCommands.State.Snapshots.Fetch.dispatch(this.plugin, { url });
    }

    render() {
        return <div>
            <b>Remote</b> <button onClick={this.refresh} disabled={this.state.isFetching}>Refresh</button>
            <ul style={{ listStyle: 'none' }}>
                {this.state.entries.valueSeq().map(e =><li key={e!.id}>
                    <button onClick={this.fetch(e!.url)} disabled={this.state.isFetching}>Fetch</button>
                    &nbsp;{e!.name} <small>{e!.description}</small>
                </li>)}
            </ul>
        </div>;
    }
}