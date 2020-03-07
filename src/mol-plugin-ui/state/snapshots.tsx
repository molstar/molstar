/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from '../../mol-plugin/commands';
import * as React from 'react';
import { PluginUIComponent, PurePluginUIComponent } from '../base';
import { shallowEqualObjects } from '../../mol-util';
import { OrderedMap } from 'immutable';
import { ParameterControls } from '../controls/parameters';
import { ParamDefinition as PD} from '../../mol-util/param-definition';
import { PluginState } from '../../mol-plugin/state';
import { urlCombine } from '../../mol-util/url';
import { IconButton, SectionHeader } from '../controls/common';
import { formatTimespan } from '../../mol-util/now';
import { PluginConfig } from '../../mol-plugin/config';
import { Icon } from '../controls/icons';

export class StateSnapshots extends PluginUIComponent<{ }> {
    downloadToFile = () => {
        PluginCommands.State.Snapshots.DownloadToFile(this.plugin, { });
    }

    open = (e: React.ChangeEvent<HTMLInputElement>) => {
        if (!e.target.files || !e.target.files![0]) return;

        PluginCommands.State.Snapshots.OpenFile(this.plugin, { file: e.target.files![0] });
    }

    render() {
        return <div>
            <SectionHeader icon='floppy' title='Plugin State' />
            <LocalStateSnapshots />
            <LocalStateSnapshotList />
            {this.plugin.spec.components?.remoteState !== 'none' && <RemoteStateSnapshots />}

            <div className='msp-btn-row-group' style={{ marginTop: '10px' }}>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.downloadToFile}>Download JSON</button>
                <div className='msp-btn msp-btn-block msp-btn-action msp-loader-msp-btn-file'>
                    {'Open JSON'} <input onChange={this.open} type='file' multiple={false} accept='.json' />
                </div>
            </div>
        </div>;
    }
}

class LocalStateSnapshots extends PluginUIComponent<
    { },
    { params: PD.Values<typeof LocalStateSnapshots.Params> }> {

    state = { params: PD.getDefaultValues(LocalStateSnapshots.Params) };

    static Params = {
        name: PD.Text(),
        options: PD.Group({
            description: PD.Text(),
            ...PluginState.GetSnapshotParams
        })
    };

    add = () => {
        PluginCommands.State.Snapshots.Add(this.plugin, {
            name: this.state.params.name,
            description: this.state.params.options.description,
            params: this.state.params.options
        });
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
        PluginCommands.State.Snapshots.Clear(this.plugin, {});
    }

    shouldComponentUpdate(nextProps: any, nextState: any) {
        return !shallowEqualObjects(this.props, nextProps) || !shallowEqualObjects(this.state, nextState);
    }

    render() {
        // TODO: proper styling
        return <div>
            <ParameterControls params={LocalStateSnapshots.Params} values={this.state.params} onEnter={this.add} onChange={p => {
                const params = { ...this.state.params, [p.name]: p.value };
                this.setState({ params } as any);
                this.plugin.state.snapshots.currentGetSnapshotParams = params.options;
            }}/>

            <div className='msp-btn-row-group'>
                {/* <button className='msp-btn msp-btn-block msp-form-control' onClick={this.add}><Icon name='floppy' /> Save</button> */}
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

    apply = (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        PluginCommands.State.Snapshots.Apply(this.plugin, { id });
    }

    remove = (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        PluginCommands.State.Snapshots.Remove(this.plugin, { id });
    }

    moveUp = (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        PluginCommands.State.Snapshots.Move(this.plugin, { id, dir: -1 });
    }

    moveDown = (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        PluginCommands.State.Snapshots.Move(this.plugin, { id, dir: 1 });
    }

    replace = (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        PluginCommands.State.Snapshots.Replace(this.plugin, { id, params: this.plugin.state.snapshots.currentGetSnapshotParams });
    }

    render() {
        const current = this.plugin.state.snapshots.state.current;
        return <ul style={{ listStyle: 'none' }} className='msp-state-list'>
            {this.plugin.state.snapshots.state.entries.map(e => <li key={e!.snapshot.id}>
                <button data-id={e!.snapshot.id} className='msp-btn msp-btn-block msp-form-control' onClick={this.apply}>
                    <span style={{ fontWeight: e!.snapshot.id === current ? 'bold' : void 0}}>
                        {e!.name || new Date(e!.timestamp).toLocaleString()}</span> <small>
                        {`${e!.snapshot.durationInMs ? formatTimespan(e!.snapshot.durationInMs, false) + `${e!.description ? ', ' : ''}` : ''}${e!.description ? e!.description : ''}`}
                    </small>
                </button>
                <div>
                    <IconButton data-id={e!.snapshot.id} icon='up-thin' title='Move Up' onClick={this.moveUp} isSmall={true} />
                    <IconButton data-id={e!.snapshot.id} icon='down-thin' title='Move Down' onClick={this.moveDown} isSmall={true} />
                    <IconButton data-id={e!.snapshot.id} icon='switch' title='Replace' onClick={this.replace} isSmall={true} />
                    <IconButton data-id={e!.snapshot.id} icon='remove' title='Remove' onClick={this.remove} isSmall={true} />
                </div>
            </li>)}
        </ul>;
    }
}

export type RemoteEntry = { url: string, removeUrl: string, timestamp: number, id: string, name: string, description: string, isSticky?: boolean }
export class RemoteStateSnapshots extends PluginUIComponent<
    { listOnly?: boolean },
    { params: PD.Values<RemoteStateSnapshots['Params']>, entries: OrderedMap<string, RemoteEntry>, isBusy: boolean }> {

    Params = {
        name: PD.Text(),
        options: PD.Group({
            description: PD.Text(),
            playOnLoad: PD.Boolean(false),
            serverUrl: PD.Text(this.plugin.config.get(PluginConfig.State.CurrentServer))
        })
    };

    state = { params: PD.getDefaultValues(this.Params), entries: OrderedMap<string, RemoteEntry>(), isBusy: false };

    ListOnlyParams = {
        options: PD.Group({
            serverUrl: PD.Text(this.plugin.config.get(PluginConfig.State.CurrentServer))
        }, { isFlat: true })
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
            this.plugin.config.set(PluginConfig.State.CurrentServer, this.state.params.options.serverUrl);

            const json = (await this.plugin.runTask<RemoteEntry[]>(this.plugin.fetch({ url: this.serverUrl('list'), type: 'json'  }))) || [];

            json.sort((a, b) => {
                if (a.isSticky === b.isSticky) return a.timestamp - b.timestamp;
                return a.isSticky ? -1 : 1;
            });

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
        this.plugin.config.set(PluginConfig.State.CurrentServer, this.state.params.options.serverUrl);

        if (this.plugin.state.snapshots.state.entries.size === 0) {
            await PluginCommands.State.Snapshots.Add(this.plugin, {
                name: this.state.params.name,
                description: this.state.params.options.description,
                params: this.plugin.state.snapshots.currentGetSnapshotParams
            });
        }

        await PluginCommands.State.Snapshots.Upload(this.plugin, {
            name: this.state.params.name,
            description: this.state.params.options.description,
            playOnLoad: this.state.params.options.playOnLoad,
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
            await PluginCommands.State.Snapshots.Fetch(this.plugin, { url: entry.url });
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
        return <>
            <SectionHeader title='Remote States' />

            {!this.props.listOnly && <>
                <ParameterControls params={this.Params} values={this.state.params} onEnter={this.upload} onChange={p => {
                    this.setState({ params: { ...this.state.params, [p.name]: p.value } } as any);
                }} isDisabled={this.state.isBusy}/>
                <div className='msp-btn-row-group'>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={this.upload} disabled={this.state.isBusy}><Icon name='upload' /> Upload</button>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={this.refresh} disabled={this.state.isBusy}>Refresh</button>
                </div>
            </>}

            <RemoteStateSnapshotList entries={this.state.entries} isBusy={this.state.isBusy} serverUrl={this.state.params.options.serverUrl}
                fetch={this.fetch} remove={this.props.listOnly ? void 0 : this.remove} />

            {this.props.listOnly && <>
                <ParameterControls params={this.ListOnlyParams} values={this.state.params} onEnter={this.upload} onChange={p => {
                    this.setState({ params: { ...this.state.params, [p.name]: p.value } } as any);
                }} isDisabled={this.state.isBusy}/>
                <div className='msp-btn-row-group'>
                    <button className='msp-btn msp-btn-block msp-form-control' onClick={this.refresh} disabled={this.state.isBusy}>Refresh</button>
                </div>
            </>}
        </>;
    }
}

class RemoteStateSnapshotList extends PurePluginUIComponent<
    { entries: OrderedMap<string, RemoteEntry>, serverUrl: string, isBusy: boolean, fetch: (e: React.MouseEvent<HTMLElement>) => void, remove?: (e: React.MouseEvent<HTMLElement>) => void },
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
                {!e!.isSticky && this.props.remove && <div>
                    <IconButton data-id={e!.id} icon='remove' title='Remove' onClick={this.props.remove} disabled={this.props.isBusy} />
                </div>}
            </li>)}
        </ul>;
    }
}
