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
import { ParameterControls } from './controls/parameters';
import { ParamDefinition as PD} from 'mol-util/param-definition';
import { Subject } from 'rxjs';

export class StateSnapshots extends PluginComponent<{ }, { serverUrl: string }> {
    state = { serverUrl: 'https://webchem.ncbr.muni.cz/molstar-state' }

    updateServerUrl = (serverUrl: string) => { this.setState({ serverUrl }) };

    render() {
        return <div>
            <div className='msp-section-header'>State Snapshots</div>
            <StateSnapshotControls serverUrl={this.state.serverUrl} serverChanged={this.updateServerUrl} />
            <LocalStateSnapshotList />
            <RemoteStateSnapshotList serverUrl={this.state.serverUrl} />
        </div>;
    }
}

// TODO: this is not nice: device some custom event system.
const UploadedEvent = new Subject();

class StateSnapshotControls extends PluginComponent<{ serverUrl: string, serverChanged: (url: string) => void }, { name: string, description: string, serverUrl: string, isUploading: boolean }> {
    state = { name: '', description: '', serverUrl: this.props.serverUrl, isUploading: false };

    static Params = {
        name: PD.Text(),
        description: PD.Text(),
        serverUrl: PD.Text()
    }

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
        this.plugin.log.message('Snapshot uploaded.');
        UploadedEvent.next();
    }

    render() {
        return <div>
            <ParameterControls params={StateSnapshotControls.Params} values={this.state} onEnter={this.add} onChange={p => {
                this.setState({ [p.name]: p.value } as any);
                if (p.name === 'serverUrl') this.props.serverChanged(p.value);
            }}/>

            <div className='msp-btn-row-group'>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.add}>Add Local</button>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.upload} disabled={this.state.isUploading}>Upload</button>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.clear}>Clear</button>
            </div>
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
        return <ul style={{ listStyle: 'none' }} className='msp-state-list'>
            {this.plugin.state.snapshots.entries.valueSeq().map(e =><li key={e!.id}>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.apply(e!.id)}>{e!.name || e!.timestamp} <small>{e!.description}</small></button>
                <button onClick={this.remove(e!.id)} className='msp-btn msp-btn-link msp-state-list-remove-button'>
                    <span className='msp-icon msp-icon-remove' />
                </button>
            </li>)}
        </ul>;
    }
}

type RemoteEntry = { url: string, removeUrl: string, timestamp: number, id: string, name: string, description: string }
class RemoteStateSnapshotList extends PluginComponent<{ serverUrl: string }, { entries: List<RemoteEntry>, isFetching: boolean }> {
    state = { entries: List<RemoteEntry>(), isFetching: false };

    componentDidMount() {
        this.subscribe(this.plugin.events.state.snapshots.changed, () => this.forceUpdate());
        this.refresh();
        this.subscribe(UploadedEvent, this.refresh);
    }

    refresh = async () => {
        try {
            this.setState({ isFetching: true });
            const req = await fetch(`${this.props.serverUrl}/list`);
            const json: RemoteEntry[] = await req.json();
            this.setState({
                entries: List<RemoteEntry>(json.map((e: RemoteEntry) => ({
                    ...e,
                    url: `${this.props.serverUrl}/get/${e.id}`,
                    removeUrl: `${this.props.serverUrl}/remove/${e.id}`
                }))),
                isFetching: false })
        } catch (e) {
            this.plugin.log.error('Fetching Remote Snapshots: ' + e);
            this.setState({ entries: List<RemoteEntry>(), isFetching: false })
        }
    }

    fetch(url: string) {
        return () => PluginCommands.State.Snapshots.Fetch.dispatch(this.plugin, { url });
    }

    remove(url: string) {
        return async () => {
            this.setState({ entries: List() });
            try {
                await fetch(url);
            } catch { }
            this.refresh();
        }
    }

    render() {
        return <div>
            <button title='Click to Refresh' style={{fontWeight: 'bold'}} className='msp-btn msp-btn-block msp-form-control msp-section-header' onClick={this.refresh} disabled={this.state.isFetching}>↻ Remote Snapshots</button>

            <ul style={{ listStyle: 'none' }} className='msp-state-list'>
                {this.state.entries.valueSeq().map(e =><li key={e!.id}>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={this.fetch(e!.url)}>{e!.name || e!.timestamp} <small>{e!.description}</small></button>
                    <button onClick={this.remove(e!.removeUrl)} className='msp-btn msp-btn-link msp-state-list-remove-button'>
                        <span className='msp-icon msp-icon-remove' />
                    </button>
                </li>)}
            </ul>
        </div>;
    }
}