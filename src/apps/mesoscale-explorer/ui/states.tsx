/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { LoadCellPackModel, LoadCellPackModelParams } from '../../../extensions/cellpack/model';
import { LoadPetworldModel } from '../../../extensions/petworld/model';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { Button, ExpandGroup } from '../../../mol-plugin-ui/controls/common';
import { GetAppSvg, Icon, OpenInBrowserSvg } from '../../../mol-plugin-ui/controls/icons';
import { ApplyActionControl } from '../../../mol-plugin-ui/state/apply-action';
import { LocalStateSnapshotList, LocalStateSnapshotParams, LocalStateSnapshots } from '../../../mol-plugin-ui/state/snapshots';
import { PluginCommands } from '../../../mol-plugin/commands';
import { PluginContext } from '../../../mol-plugin/context';
import { StateAction, StateTransform } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

export const LoadModel = StateAction.build({
    display: { name: 'Load', description: 'Open or download a model' },
    params: {
        type: PD.Select('cellpack', PD.arrayToOptions(['cellpack', 'petworld'] as const), { label: 'Type' }),
        files: PD.FileList({ accept: '.cif,.bcif', multiple: true, description: 'Cellpack- or Petworld-style cif file.', label: 'File(s)' }),
    },
    from: PluginStateObject.Root
})(({ params }, ctx: PluginContext) => Task.create('Model Loader', async taskCtx => {
    if (params.files === null) {
        ctx.log.error('No file(s) selected');
        return;
    }

    PluginCommands.State.RemoveObject(ctx, { state: ctx.state.data, ref: StateTransform.RootRef });

    if (params.type === 'cellpack') {
        for (const file of params.files) {
            try {
                await PluginCommands.State.ApplyAction(ctx, {
                    state: ctx.state.data,
                    action: LoadCellPackModel.create({
                        ...PD.getDefaultValues(LoadCellPackModelParams),
                        source: {
                            name: 'cif',
                            params: file
                        },
                    }),
                    ref: ctx.state.data.tree.root.ref
                });
            } catch (e) {
                console.error(e);
                ctx.log.error(`Error opening file '${file.name}'`);
            }
        }
    } else if (params.type === 'petworld') {
        for (const file of params.files) {
            try {
                await PluginCommands.State.ApplyAction(ctx, {
                    state: ctx.state.data,
                    action: LoadPetworldModel.create({
                        cif: file
                    }),
                    ref: ctx.state.data.tree.root.ref
                });
            } catch (e) {
                console.error(e);
                ctx.log.error(`Error opening file '${file.name}'`);
            }
        }
    } else {
        ctx.log.error('Unsupported type');
    }
}));

export class LoaderControls extends PluginUIComponent {
    componentDidMount() {

    }

    render() {
        return <div style={{ margin: '5px' }}>
            <ApplyActionControl state={this.plugin.state.data} action={LoadModel} nodeRef={this.plugin.state.data.tree.root.ref} applyLabel={'Load'} hideHeader />
        </div>;
    }
}

export class SessionControls extends PluginUIComponent {
    downloadToFileZip = () => {
        PluginCommands.State.Snapshots.DownloadToFile(this.plugin, { type: 'zip' });
    };

    open = (e: React.ChangeEvent<HTMLInputElement>) => {
        if (!e.target.files || !e.target.files[0]) {
            this.plugin.log.error('No state file selected');
            return;
        }
        PluginCommands.State.Snapshots.OpenFile(this.plugin, { file: e.target.files[0] });
    };

    render() {
        return <div style={{ margin: '5px' }}>
            <div className='msp-flex-row'>
                <Button icon={GetAppSvg} onClick={this.downloadToFileZip} title='Save the state.'>
                    Save
                </Button>
                <div className='msp-btn msp-btn-block msp-btn-action msp-loader-msp-btn-file'>
                    <Icon svg={OpenInBrowserSvg} inline /> Open <input onChange={this.open} type='file' multiple={false} accept='.molx,.molj' />
                </div>
            </div>
        </div>;
    }
}

export class SnapshotControls extends PluginUIComponent<{}> {
    render() {
        return <div style={{ margin: '5px' }}>
            <div style={{ marginBottom: '10px' }}>
                <LocalStateSnapshotList />
            </div>
            <div style={{ marginBottom: '10px' }}>
                <LocalStateSnapshots />
            </div>

            <div style={{ marginBottom: '10px' }}>
                <ExpandGroup header='Snapshot Options' initiallyExpanded={false}>
                    <LocalStateSnapshotParams />
                </ExpandGroup>
            </div>
        </div>;
    }
}
