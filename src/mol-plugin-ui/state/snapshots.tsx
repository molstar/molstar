/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedMap } from 'immutable';
import * as React from 'react';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginState } from '../../mol-plugin/state';
import { shallowEqualObjects } from '../../mol-util';
import { formatTimespan } from '../../mol-util/now';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { urlCombine } from '../../mol-util/url';
import { PluginUIComponent, PurePluginUIComponent } from '../base';
import { Button, ExpandGroup, IconButton, SectionHeader } from '../controls/common';
import { Icon, SaveOutlinedSvg, GetAppSvg, OpenInBrowserSvg, WarningSvg, DeleteOutlinedSvg, AddSvg, ArrowUpwardSvg, SwapHorizSvg, ArrowDownwardSvg, RefreshSvg, CloudUploadSvg } from '../controls/icons';
import { ParameterControls } from '../controls/parameters';

export class StateSnapshots extends PluginUIComponent<{}> {
    render() {
        return <div>
            <SectionHeader icon={SaveOutlinedSvg} title='Plugin State' />

            <div style={{ marginBottom: '10px' }}>
                <ExpandGroup header='Save Options' initiallyExpanded={false}>
                    <LocalStateSnapshotParams />
                </ExpandGroup>
            </div>

            <LocalStateSnapshots />
            <LocalStateSnapshotList />

            <SectionHeader title='Save as File' accent='blue' />
            <StateExportImportControls />

            {this.plugin.spec.components?.remoteState !== 'none' && <RemoteStateSnapshots />}

        </div>;
    }
}

export class StateExportImportControls extends PluginUIComponent<{ onAction?: () => void }> {
    downloadToFileJson = () => {
        this.props.onAction?.();
        PluginCommands.State.Snapshots.DownloadToFile(this.plugin, { type: 'json' });
    }

    downloadToFileZip = () => {
        this.props.onAction?.();
        PluginCommands.State.Snapshots.DownloadToFile(this.plugin, { type: 'zip' });
    }

    open = (e: React.ChangeEvent<HTMLInputElement>) => {
        if (!e.target.files || !e.target.files[0]) {
            this.plugin.log.error('No state file selected');
            return;
        }

        this.props.onAction?.();
        PluginCommands.State.Snapshots.OpenFile(this.plugin, { file: e.target.files[0] });
    }

    render() {
        return <>
            <div className='msp-flex-row'>
                <Button icon={GetAppSvg} onClick={this.downloadToFileJson} title='Save the state description. Input data are loaded using the provided sources. Does not work if local files are used as input.'>
                    State
                </Button>
                <Button icon={GetAppSvg} onClick={this.downloadToFileZip} title='Save the state including the input data.'>
                    Session
                </Button>
                <div className='msp-btn msp-btn-block msp-btn-action msp-loader-msp-btn-file'>
                    <Icon svg={OpenInBrowserSvg} inline /> Open <input onChange={this.open} type='file' multiple={false} accept='.molx,.molj' />
                </div>
            </div>
            <div className='msp-help-text' style={{ padding: '10px'}}>
                <Icon svg={WarningSvg} /> This is an experimental feature and stored states/sessions might not be openable in a future version.
            </div>
        </>;
    }
}

export class LocalStateSnapshotParams extends PluginUIComponent {
    componentDidMount() {
        this.subscribe(this.plugin.state.snapshotParams, () => this.forceUpdate());
    }

    render() {
        return <ParameterControls params={PluginState.SnapshotParams} values={this.plugin.state.snapshotParams.value} onChangeValues={this.plugin.state.setSnapshotParams} />;
    }
}

export class LocalStateSnapshots extends PluginUIComponent<
{},
{ params: PD.Values<typeof LocalStateSnapshots.Params> }> {
    state = { params: PD.getDefaultValues(LocalStateSnapshots.Params) };

    static Params = {
        name: PD.Text(),
        description: PD.Text()
    };

    add = () => {
        PluginCommands.State.Snapshots.Add(this.plugin, {
            name: this.state.params.name,
            description: this.state.params.description
        });
    }

    updateParams = (params: PD.Values<typeof LocalStateSnapshots.Params>) => this.setState({ params });

    clear = () => {
        PluginCommands.State.Snapshots.Clear(this.plugin, {});
    }

    shouldComponentUpdate(nextProps: any, nextState: any) {
        return !shallowEqualObjects(this.props, nextProps) || !shallowEqualObjects(this.state, nextState);
    }

    render() {
        return <div>
            <ParameterControls params={LocalStateSnapshots.Params} values={this.state.params} onEnter={this.add} onChangeValues={this.updateParams} />
            <div className='msp-flex-row'>
                <IconButton onClick={this.clear} svg={DeleteOutlinedSvg} title='Remove All' />
                <Button onClick={this.add} icon={AddSvg} style={{ textAlign: 'right' }} commit>Add</Button>
            </div>
        </div>;
    }
}

export class LocalStateSnapshotList extends PluginUIComponent<{}, {}> {
    componentDidMount() {
        this.subscribe(this.plugin.managers.snapshot.events.changed, () => this.forceUpdate());
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
        PluginCommands.State.Snapshots.Replace(this.plugin, { id });
    }

    render() {
        const current = this.plugin.managers.snapshot.state.current;
        return <ul style={{ listStyle: 'none', marginTop: '10px' }} className='msp-state-list'>
            {this.plugin.managers.snapshot.state.entries.map(e => <li key={e!.snapshot.id} className='msp-flex-row'>
                <Button data-id={e!.snapshot.id} onClick={this.apply} className='msp-no-overflow'>
                    <span style={{ fontWeight: e!.snapshot.id === current ? 'bold' : void 0 }}>
                        {e!.name || new Date(e!.timestamp).toLocaleString()}</span> <small>
                        {`${e!.snapshot.durationInMs ? formatTimespan(e!.snapshot.durationInMs, false) + `${e!.description ? ', ' : ''}` : ''}${e!.description ? e!.description : ''}`}
                    </small>
                </Button>
                <IconButton svg={ArrowUpwardSvg} data-id={e!.snapshot.id} title='Move Up' onClick={this.moveUp} flex='20px' />
                <IconButton svg={ArrowDownwardSvg} data-id={e!.snapshot.id} title='Move Down' onClick={this.moveDown} flex='20px' />
                <IconButton svg={SwapHorizSvg} data-id={e!.snapshot.id} title='Replace' onClick={this.replace} flex='20px' />
                <IconButton svg={DeleteOutlinedSvg} data-id={e!.snapshot.id} title='Remove' onClick={this.remove} flex='20px' />
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

    private _mounted = false;
    componentDidMount() {
        this.refresh();
        // TODO: solve this by using "PluginComponent" with behaviors intead
        this._mounted = true;
        // this.subscribe(UploadedEvent, this.refresh);
    }

    componentWillUnmount() {
        this._mounted = false;
    }

    serverUrl(q?: string) {
        if (!q) return this.state.params.options.serverUrl;
        return urlCombine(this.state.params.options.serverUrl, q);
    }

    refresh = async () => {
        try {
            this.setState({ isBusy: true });
            this.plugin.config.set(PluginConfig.State.CurrentServer, this.state.params.options.serverUrl);

            const json = (await this.plugin.runTask<RemoteEntry[]>(this.plugin.fetch({ url: this.serverUrl('list'), type: 'json' }))) || [];

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

            if (this._mounted) this.setState({ entries: entries.asImmutable(), isBusy: false });
        } catch (e) {
            this.plugin.log.error('Fetching Remote Snapshots: ' + e);
            if (this._mounted) this.setState({ entries: OrderedMap(), isBusy: false });
        }
    }

    upload = async () => {
        this.setState({ isBusy: true });
        this.plugin.config.set(PluginConfig.State.CurrentServer, this.state.params.options.serverUrl);

        await PluginCommands.State.Snapshots.Upload(this.plugin, {
            name: this.state.params.name,
            description: this.state.params.options.description,
            playOnLoad: this.state.params.options.playOnLoad,
            serverUrl: this.state.params.options.serverUrl
        });

        this.plugin.log.message('Snapshot uploaded.');

        if (this._mounted) {
            this.setState({ isBusy: false });
            this.refresh();
        }
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
            if (this._mounted) this.setState({ isBusy: false });
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
            <SectionHeader title='Remote States' accent='blue' />

            {!this.props.listOnly && <>
                <ParameterControls params={this.Params} values={this.state.params} onEnter={this.upload} onChange={p => {
                    this.setState({ params: { ...this.state.params, [p.name]: p.value } } as any);
                }} isDisabled={this.state.isBusy} />
                <div className='msp-flex-row'>
                    <IconButton onClick={this.refresh} disabled={this.state.isBusy} svg={RefreshSvg} />
                    <Button icon={CloudUploadSvg} onClick={this.upload} disabled={this.state.isBusy} commit>Upload</Button>
                </div>
            </>}

            <RemoteStateSnapshotList entries={this.state.entries} isBusy={this.state.isBusy} serverUrl={this.state.params.options.serverUrl}
                fetch={this.fetch} remove={this.props.listOnly ? void 0 : this.remove} />

            {this.props.listOnly && <div style={{ marginTop: '10px' }}>
                <ParameterControls params={this.ListOnlyParams} values={this.state.params} onEnter={this.upload} onChange={p => {
                    this.setState({ params: { ...this.state.params, [p.name]: p.value } } as any);
                }} isDisabled={this.state.isBusy} />
                <div className='msp-flex-row'>
                    <Button onClick={this.refresh} disabled={this.state.isBusy} icon={RefreshSvg}>Refresh</Button>
                </div>
            </div>}
        </>;
    }
}

class RemoteStateSnapshotList extends PurePluginUIComponent<
{ entries: OrderedMap<string, RemoteEntry>, serverUrl: string, isBusy: boolean, fetch: (e: React.MouseEvent<HTMLElement>) => void, remove?: (e: React.MouseEvent<HTMLElement>) => void },
{}> {

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
        return <ul style={{ listStyle: 'none', marginTop: '10px' }} className='msp-state-list'>
            {this.props.entries.valueSeq().map(e => <li key={e!.id} className='msp-flex-row'>
                <Button data-id={e!.id} onClick={this.props.fetch}
                    disabled={this.props.isBusy} onContextMenu={this.open} title='Click to download, right-click to open in a new tab.'>
                    {e!.name || new Date(e!.timestamp).toLocaleString()} <small>{e!.description}</small>
                </Button>
                {!e!.isSticky && this.props.remove && <IconButton svg={DeleteOutlinedSvg} data-id={e!.id} title='Remove' onClick={this.props.remove} disabled={this.props.isBusy} small />}
            </li>)}
        </ul>;
    }
}
