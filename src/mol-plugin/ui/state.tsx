/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from 'mol-plugin/command';
import * as React from 'react';
import { PluginUIComponent, PurePluginUIComponent } from './base';
import { shallowEqual } from 'mol-util';
import { OrderedMap } from 'immutable';
import { ParameterControls } from './controls/parameters';
import { ParamDefinition as PD} from 'mol-util/param-definition';
import { PluginState } from 'mol-plugin/state';
import { urlCombine } from 'mol-util/url';

export class StateSnapshots extends PluginUIComponent<{ }> {

    render() {
        return <div>
            <div className='msp-section-header'>State</div>
            <StateSnapshotControls />
            <LocalStateSnapshotList />
            <RemoteStateSnapshots />
        </div>;
    }
}

class StateSnapshotControls extends PluginUIComponent<
    { },
    { params: PD.Values<typeof StateSnapshotControls.Params> }> {

    state = { params: PD.getDefaultValues(StateSnapshotControls.Params) };

    static Params = {
        name: PD.Text(),
        options: PD.Group({
            description: PD.Text(),
            ...PluginState.GetSnapshotParams
        })
    };

    add = () => {
        PluginCommands.State.Snapshots.Add.dispatch(this.plugin, { name: this.state.params.name, description: this.state.params.options.description });
        this.setState({
            params: {
                name: '',
                options: {
                    ...this.state.params.options,
                    description: ''
                }
            }
        });
    }

    clear = () => {
        PluginCommands.State.Snapshots.Clear.dispatch(this.plugin, {});
    }

    shouldComponentUpdate(nextProps: any, nextState: any) {
        return !shallowEqual(this.props, nextProps) || !shallowEqual(this.state, nextState);
    }

    downloadToFile = () => {
        PluginCommands.State.Snapshots.DownloadToFile.dispatch(this.plugin, { name: this.state.params.name });
    }

    open = (e: React.ChangeEvent<HTMLInputElement>) => {
        if (!e.target.files || !e.target.files![0]) return;

        PluginCommands.State.Snapshots.OpenFile.dispatch(this.plugin, { file: e.target.files![0] });
    }

    render() {
        // TODO: proper styling
        return <div>
            <div className='msp-btn-row-group' style={{ marginBottom: '10px' }}>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.downloadToFile}>Download JSON</button>
                <div className='msp-btn msp-btn-block msp-btn-action msp-loader-msp-btn-file'>
                    {'Open JSON'} <input onChange={this.open} type='file' multiple={false} accept='.json' />
                </div>
            </div>

            <ParameterControls params={StateSnapshotControls.Params} values={this.state.params} onEnter={this.add} onChange={p => {
                this.setState({ params: { ...this.state.params, [p.name]: p.value } } as any);
            }}/>

            <div className='msp-btn-row-group'>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.add}>Save</button>
                {/* <button className='msp-btn msp-btn-block msp-form-control' onClick={this.upload} disabled={this.state.isUploading}>Upload</button> */}
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.clear}>Clear</button>
            </div>
        </div>;
    }
}

class LocalStateSnapshotList extends PluginUIComponent<{ }, { }> {
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
        const current = this.plugin.state.snapshots.state.current;
        return <ul style={{ listStyle: 'none' }} className='msp-state-list'>
            {this.plugin.state.snapshots.state.entries.valueSeq().map(e =><li key={e!.snapshot.id}>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.apply(e!.snapshot.id)}>
                    <span style={{ fontWeight: e!.snapshot.id === current ? 'bold' : void 0}}>{e!.name || new Date(e!.timestamp).toLocaleString()}</span> <small>{e!.description}</small>
                </button>
                <button onClick={this.remove(e!.snapshot.id)} className='msp-btn msp-btn-link msp-state-list-remove-button'>
                    <span className='msp-icon msp-icon-remove' />
                </button>
            </li>)}
        </ul>;
    }
}

type RemoteEntry = { url: string, removeUrl: string, timestamp: number, id: string, name: string, description: string }
class RemoteStateSnapshots extends PluginUIComponent<
    { },
    { params: PD.Values<typeof RemoteStateSnapshots.Params>, entries: OrderedMap<string, RemoteEntry>, isBusy: boolean }> {

    state = { params: PD.getDefaultValues(RemoteStateSnapshots.Params), entries: OrderedMap<string, RemoteEntry>(), isBusy: false };

    static Params = {
        name: PD.Text(),
        options: PD.Group({
            description: PD.Text(),
            serverUrl: PD.Text('https://webchem.ncbr.muni.cz/molstar-state')
        })
    };

    componentDidMount() {
        this.refresh();
        // this.subscribe(UploadedEvent, this.refresh);
    }

    serverUrl(q?: string) {
        if (!q) return this.state.params.options.serverUrl;
        return urlCombine(this.state.params.options.serverUrl, q);
    }

    refresh = async () => {
        try {
            this.setState({ isBusy: true });
            const json = await this.plugin.runTask<RemoteEntry[]>(this.plugin.fetch({ url: this.serverUrl('list'), type: 'json'  }));
            const entries = OrderedMap<string, RemoteEntry>().asMutable();
            for (const e of json) {
                entries.set(e.id, {
                    ...e,
                    url: this.serverUrl(`get/${e.id}`),
                    removeUrl: this.serverUrl(`remove/${e.id}`)
                });
            }

            this.setState({ entries: entries.asImmutable(), isBusy: false })
        } catch (e) {
            this.plugin.log.error('Fetching Remote Snapshots: ' + e);
            this.setState({ entries: OrderedMap(), isBusy: false })
        }
    }

    upload = async () => {
        this.setState({ isBusy: true });
        if (this.plugin.state.snapshots.state.entries.size === 0) {
            await PluginCommands.State.Snapshots.Add.dispatch(this.plugin, { name: this.state.params.name, description: this.state.params.options.description });
        }

        await PluginCommands.State.Snapshots.Upload.dispatch(this.plugin, {
            name: this.state.params.name,
            description: this.state.params.options.description,
            serverUrl: this.state.params.options.serverUrl
        });
        this.setState({ isBusy: false });
        this.plugin.log.message('Snapshot uploaded.');
        this.refresh();
    }

    fetch = async (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        const entry = this.state.entries.get(id);
        if (!entry) return;

        this.setState({ isBusy: true });
        try {
            await PluginCommands.State.Snapshots.Fetch.dispatch(this.plugin, { url: entry.url });
        } finally {
            this.setState({ isBusy: false });
        }
    }

    remove = async (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        const entry = this.state.entries.get(id);
        if (!entry) return;
        this.setState({ entries: this.state.entries.remove(id) });

        try {
            await fetch(entry.removeUrl);
        } catch { }
    }

    render() {
        return <div>
            <div className='msp-section-header'>Remote State</div>

            <ParameterControls params={RemoteStateSnapshots.Params} values={this.state.params} onEnter={this.upload} onChange={p => {
                this.setState({ params: { ...this.state.params, [p.name]: p.value } } as any);
            }} isDisabled={this.state.isBusy}/>

            <div className='msp-btn-row-group'>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.upload} disabled={this.state.isBusy}>Upload</button>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.refresh} disabled={this.state.isBusy}>Refresh</button>
            </div>

            <RemoteStateSnapshotList entries={this.state.entries} isBusy={this.state.isBusy} serverUrl={this.state.params.options.serverUrl}
                fetch={this.fetch} remove={this.remove} />
        </div>;
    }
}

class RemoteStateSnapshotList extends PurePluginUIComponent<
    { entries: OrderedMap<string, RemoteEntry>, serverUrl: string, isBusy: boolean, fetch: (e: React.MouseEvent<HTMLElement>) => void, remove: (e: React.MouseEvent<HTMLElement>) => void },
    { }> {

    open = async (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        const entry = this.props.entries.get(id);
        if (!entry) return;

        e.preventDefault();
        let url = `${window.location}`, qi = url.indexOf('?');
        if (qi > 0) url = url.substr(0, qi);

        window.open(`${url}?snapshot-url=${encodeURIComponent(entry.url)}`, '_blank');
    }

    render() {
        return <ul style={{ listStyle: 'none' }} className='msp-state-list'>
            {this.props.entries.valueSeq().map(e =><li key={e!.id}>
                <button data-id={e!.id} className='msp-btn msp-btn-block msp-form-control' onClick={this.props.fetch}
                    disabled={this.props.isBusy} onContextMenu={this.open} title='Click to download, right-click to open in a new tab.'>
                    {e!.name || new Date(e!.timestamp).toLocaleString()} <small>{e!.description}</small>
                </button>
                <button data-id={e!.id} onClick={this.props.remove} className='msp-btn msp-btn-link msp-state-list-remove-button' disabled={this.props.isBusy}>
                    <span className='msp-icon msp-icon-remove' />
                </button>
            </li>)}
        </ul>;
    }
}
