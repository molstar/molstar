/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
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
import { Button, ControlRow, ExpandGroup, IconButton, SectionHeader } from '../controls/common';
import { Icon, SaveOutlinedSvg, GetAppSvg, OpenInBrowserSvg, WarningSvg, DeleteOutlinedSvg, AddSvg, ArrowUpwardSvg, SwapHorizSvg, ArrowDownwardSvg, RefreshSvg, CloudUploadSvg, CheckSvg, TuneSvg } from '../controls/icons';
import { ParamHelp, ParameterControls, ToggleParamHelpButton } from '../controls/parameters';
import { PluginStateSnapshotManager } from '../../mol-plugin-state/manager/snapshots';
import { PluginContext } from '../../mol-plugin/context';

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
    };

    downloadToFileZip = () => {
        this.props.onAction?.();
        PluginCommands.State.Snapshots.DownloadToFile(this.plugin, { type: 'zip' });
    };

    open = (e: React.ChangeEvent<HTMLInputElement>) => {
        if (!e.target.files || !e.target.files[0]) {
            this.plugin.log.error('No state file selected');
            return;
        }

        this.props.onAction?.();
        PluginCommands.State.Snapshots.OpenFile(this.plugin, { file: e.target.files[0] });
    };

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
            <div className='msp-help-text' style={{ padding: '10px' }}>
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
    };

    updateParams = (params: PD.Values<typeof LocalStateSnapshots.Params>) => this.setState({ params });

    clear = () => {
        PluginCommands.State.Snapshots.Clear(this.plugin, {});
    };

    shouldComponentUpdate(nextProps: any, nextState: any) {
        return !shallowEqualObjects(this.props, nextProps) || !shallowEqualObjects(this.state, nextState);
    }

    render() {
        return <div>
            <AddSnapshot parent={this} />
        </div>;
    }
}

function invalidateSnapshotKey(plugin: PluginContext, key: string | undefined, currentId?: string) {
    if (!key) return false;
    return plugin.managers.snapshot.state.entries.some(e => (!currentId || e.snapshot.id !== currentId) && e.key === key);
}

function AddSnapshot({ parent }: { parent: LocalStateSnapshots }) {
    const [state, setState] = React.useState<PluginStateSnapshotManager.EntryParams>({ key: '', name: '', description: '' });

    const add = () => {
        PluginCommands.State.Snapshots.Add(parent.plugin, {
            key: state.key,
            name: state.name,
            description: state.description,
        });
        setState({ key: '', name: '', description: '' });
    };

    const keyExists = invalidateSnapshotKey(parent.plugin, state.key);

    return <>
        <EditSnapshotParams state={state} setState={setState} apply={add} />
        <div className='msp-flex-row'>
            <IconButton onClick={parent.clear} svg={DeleteOutlinedSvg} title='Remove All' />
            <Button onClick={add} icon={keyExists ? undefined : AddSvg} style={{ textAlign: 'right' }} commit={keyExists ? 'off' : 'on'} disabled={keyExists}>
                {keyExists
                    ? 'Key must be unique'
                    : 'Add'}
            </Button>
        </div>
    </>;
}


export class LocalStateSnapshotList extends PluginUIComponent<{}, { editingId?: string }> {
    state = { editingId: undefined as string | undefined };

    componentDidMount() {
        this.subscribe(this.plugin.managers.snapshot.events.changed, () => this.forceUpdate());
    }

    edit = (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        const current = this.state.editingId;
        this.setState({ editingId: id === current ? undefined : id });
    };

    doneEdit = () => this.setState({ editingId: undefined });

    apply = (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        PluginCommands.State.Snapshots.Apply(this.plugin, { id });
    };

    remove = (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        PluginCommands.State.Snapshots.Remove(this.plugin, { id });
    };

    moveUp = (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        PluginCommands.State.Snapshots.Move(this.plugin, { id, dir: -1 });
    };

    moveDown = (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        PluginCommands.State.Snapshots.Move(this.plugin, { id, dir: 1 });
    };

    replace = (e: React.MouseEvent<HTMLElement>) => {
        // TODO: add option change name/description
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        PluginCommands.State.Snapshots.Replace(this.plugin, { id });
    };

    render() {
        const current = this.plugin.managers.snapshot.state.current;
        const items: JSX.Element[] = [];
        this.plugin.managers.snapshot.state.entries.forEach(e => {
            items.push(<li key={e!.snapshot.id} className='msp-flex-row'>
                <Button data-id={e!.snapshot.id} onClick={this.apply} className='msp-no-overflow'>
                    <span style={{ fontWeight: e!.snapshot.id === current ? 'bold' : void 0 }}>
                        {!!e!.key && `[${e.key}] `}
                        {e!.name || new Date(e!.timestamp).toLocaleString()}</span> <small>
                        {`${e!.snapshot.durationInMs ? formatTimespan(e!.snapshot.durationInMs, false) : ''}`}
                    </small>
                </Button>
                <IconButton svg={TuneSvg} data-id={e!.snapshot.id} title='Edit' onClick={this.edit} flex='28px' />
                <IconButton svg={ArrowUpwardSvg} data-id={e!.snapshot.id} title='Move Up' onClick={this.moveUp} flex='20px' />
                <IconButton svg={ArrowDownwardSvg} data-id={e!.snapshot.id} title='Move Down' onClick={this.moveDown} flex='20px' />
                <IconButton svg={SwapHorizSvg} data-id={e!.snapshot.id} title='Replace' onClick={this.replace} flex='20px' />
                <IconButton svg={DeleteOutlinedSvg} data-id={e!.snapshot.id} title='Remove' onClick={this.remove} flex='20px' />
            </li>);
            if (this.state.editingId === e!.snapshot.id) {
                items.push(<EditSnapshot key={`${e!.snapshot.id}-edit`} entry={e} plugin={this.plugin} done={this.doneEdit} />);
            }
            const image = e.image && this.plugin.managers.asset.get(e.image)?.file;
            if (image) {
                items.push(<li key={`${e!.snapshot.id}-image`} className='msp-state-image-row'>
                    <Button data-id={e!.snapshot.id} onClick={this.apply}>
                        <img draggable={false} src={URL.createObjectURL(image)}/>
                    </Button>
                </li>);
            }
        });
        return <>
            <ul style={{ listStyle: 'none', marginTop: '10px' }} className='msp-state-list'>
                {items}
            </ul>
        </>;
    }
}

function EditSnapshotParams({ state, setState, apply }: { state: PluginStateSnapshotManager.EntryParams, setState: (s: PluginStateSnapshotManager.EntryParams) => any, apply: (s: PluginStateSnapshotManager.EntryParams) => any }) {
    const keyRef = React.useRef<HTMLElement>();
    const descRef = React.useRef<HTMLElement>();
    const [showKeyHelp, setShowKeyHelp] = React.useState(false);

    return <>
        <ControlRow
            label='Name'
            control={ <input type='text'
                value={state.name}
                placeholder='Name'
                onChange={e => setState({ ...state, name: e.target.value })}
                onKeyUp={e => {
                    if (e.key === 'Enter') keyRef.current?.focus();
                }}
            />}
        />
        <ControlRow
            label={<>
                Key
                <ToggleParamHelpButton show={showKeyHelp} toggle={() => setShowKeyHelp(prev => !prev)} />
            </>}
            control={ <input type='text'
                ref={keyRef as any}
                value={state.key}
                placeholder='Key (optional)'
                onChange={e => setState({ ...state, key: e.target.value })}
                onKeyUp={e => {
                    if (e.key === 'Enter') descRef.current?.focus();
                }}
            />}
        />
        {showKeyHelp && <div className='msp-control-offset'>
            <ParamHelp description='Optional snapshot key used to activate snapshots from descriptions, labels, etc.' />
        </div>}
        <div className='msp-flex-row msp-text-area-wrapper' style={{ marginBottom: 1 }}>
            <textarea
                ref={descRef as any}
                // NOTE: curly brackets are required to support \n in the placeholder, do not remove
                placeholder={'Markdown Description\n\n- Use [title](#key) to link to a snapshot'}
                className='msp-form-control'
                value={state.description}
                onChange={e => setState({ ...state, description: e.target.value })}
                onKeyUp={e => {
                    if (e.key === 'Enter' && e.ctrlKey) apply(state);
                }}
            />
        </div>
    </>;
}

function EditSnapshot({ entry, plugin, done }: { entry: PluginStateSnapshotManager.Entry, plugin: PluginContext, done: () => any }) {
    const [state, setState] = React.useState<PluginStateSnapshotManager.EntryParams>({ key: entry.key ?? '', name: entry.name ?? '', description: entry.description ?? '', descriptionFormat: entry.descriptionFormat });

    const apply = () => {
        plugin.managers.snapshot.update(entry, state);
        done();
    };

    const keyExists = invalidateSnapshotKey(plugin, state.key, entry.snapshot.id);

    return <>
        <EditSnapshotParams state={state} setState={setState} apply={apply} />
        <div className='msp-flex-row' style={{ marginBottom: 1 }}>
            <Button onClick={apply} icon={keyExists ? undefined : CheckSvg} style={{ textAlign: 'right' }} commit={keyExists ? 'off' : 'on'} disabled={keyExists}>
                {keyExists
                    ? 'Key must be unique'
                    : 'Apply'}
            </Button>
        </div>
    </>;
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
        super.componentWillUnmount();
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
            console.error(e);
            this.plugin.log.error('Error fetching remote snapshots');
            if (this._mounted) this.setState({ entries: OrderedMap(), isBusy: false });
        }
    };

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
    };


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
    };

    remove = async (e: React.MouseEvent<HTMLElement>) => {
        const id = e.currentTarget.getAttribute('data-id');
        if (!id) return;
        const entry = this.state.entries.get(id);
        if (!entry) return;
        this.setState({ entries: this.state.entries.remove(id) });

        try {
            await fetch(entry.removeUrl);
        } catch (e) {
            console.error(e);
        }
    };

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
        let url = `${window.location}`;
        const qi = url.indexOf('?');
        if (qi > 0) url = url.substr(0, qi);

        window.open(`${url}?snapshot-url=${encodeURIComponent(entry.url)}`, '_blank');
    };

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
