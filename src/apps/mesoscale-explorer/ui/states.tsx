/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';
import { MmcifProvider } from '../../../mol-plugin-state/formats/trajectory';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { Button, ExpandGroup } from '../../../mol-plugin-ui/controls/common';
import { GetAppSvg, Icon, OpenInBrowserSvg } from '../../../mol-plugin-ui/controls/icons';
import { ApplyActionControl } from '../../../mol-plugin-ui/state/apply-action';
import { LocalStateSnapshotList, LocalStateSnapshotParams, LocalStateSnapshots } from '../../../mol-plugin-ui/state/snapshots';
import { PluginCommands } from '../../../mol-plugin/commands';
import { PluginContext } from '../../../mol-plugin/context';
import { StateAction, StateObjectRef, StateTransform } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { Color } from '../../../mol-util/color/color';
import { getFileNameInfo } from '../../../mol-util/file-info';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { MesoscaleExplorerState } from '../app';
import { createCellpackHierarchy } from '../data/cellpack/preset';
import { createPetworldHierarchy } from '../data/petworld/preset';

function adjustPluginProps(ctx: PluginContext) {
    ctx.managers.interactivity.setProps({ granularity: 'chain' });
    ctx.canvas3d?.setProps({
        multiSample: { mode: 'off' },
        cameraClipping: { far: false, minNear: 50 },
        renderer: {
            colorMarker: true,
            highlightColor: Color(0xffffff),
            highlightStrength: 0,
            selectStrength: 0,
            dimColor: Color(0xffffff),
            dimStrength: 0,
            markerPriority: 2,
            interiorColorFlag: false,
            interiorDarkening: 0.15,
            exposure: 1.1,
        },
        marking: {
            enabled: false,
            highlightEdgeColor: Color(0x394e65),
            selectEdgeStrength: 0,
            ghostEdgeStrength: 1,
            innerEdgeFactor: 2.5,
        },
        postprocessing: {
            occlusion: {
                name: 'on',
                params: {
                    samples: 32,
                    multiScale: {
                        name: 'on',
                        params: {
                            levels: [
                                { radius: 2, bias: 1.0 },
                                { radius: 5, bias: 1.0 },
                                { radius: 8, bias: 1.0 },
                                { radius: 11, bias: 1.0 },
                            ],
                            nearThreshold: 10,
                            farThreshold: 1500,
                        }
                    },
                    radius: 5,
                    bias: 1,
                    blurKernelSize: 11,
                    resolutionScale: 1,
                    color: Color(0x000000),
                }
            },
            shadow: {
                name: 'on',
                params: {
                    bias: 0.6,
                    maxDistance: 80,
                    steps: 3,
                    tolerance: 1.0,
                }
            },
            outline: {
                name: 'on',
                params: {
                    scale: 1,
                    threshold: 0.15,
                    color: Color(0x000000),
                    includeTransparent: true,
                }
            }
        }
    });
}

async function createHierarchy(ctx: PluginContext, ref: string) {
    const parsed = await MmcifProvider.parse(ctx, ref);

    const tr = StateObjectRef.resolveAndCheck(ctx.state.data, parsed.trajectory)?.obj?.data;
    if (!tr) throw new Error('no trajectory');

    if (!MmcifFormat.is(tr.representative.sourceData)) {
        throw new Error('not mmcif');
    }

    const { frame } = tr.representative.sourceData.data;
    if (frame.categories.pdbx_model) {
        await createPetworldHierarchy(ctx, parsed.trajectory);
    } else {
        await createCellpackHierarchy(ctx, parsed.trajectory);
    }
}

export const LoadExample = StateAction.build({
    display: { name: 'Load', description: 'Load an example' },
    params: (a, ctx: PluginContext) => {
        const entries = (ctx.customState as MesoscaleExplorerState).examples || [];
        return {
            entry: PD.Select(0, entries.map((s, i) => [i, s.label])),
        };
    },
    from: PluginStateObject.Root
})(({ params }, ctx: PluginContext) => Task.create('Loading example...', async taskCtx => {
    const entries = (ctx.customState as MesoscaleExplorerState).examples || [];
    const { url, type } = entries[params.entry];
    console.time('LoadExample');
    if (type === 'molx' || type === 'molj') {
        await PluginCommands.State.Snapshots.OpenUrl(ctx, { url, type });
    } else {
        PluginCommands.State.Snapshots.Clear(ctx);
        PluginCommands.State.RemoveObject(ctx, { state: ctx.state.data, ref: StateTransform.RootRef });
        adjustPluginProps(ctx);

        const isBinary = type === 'bcif';
        const data = await ctx.builders.data.download({ url, isBinary });
        await createHierarchy(ctx, data.ref);
    }
    console.timeEnd('LoadExample');
}));

export const LoadModel = StateAction.build({
    display: { name: 'Load', description: 'Load a model' },
    params: {
        files: PD.FileList({ accept: '.cif,.bcif,.cif.gz,.bcif.gz,.cif.zip,.bcif.zip', multiple: true, description: 'Cellpack- or Petworld-style cif file.', label: 'File(s)' }),
    },
    from: PluginStateObject.Root
})(({ params }, ctx: PluginContext) => Task.create('Loading model...', async taskCtx => {
    if (params.files === null) {
        ctx.log.error('No file(s) selected');
        return;
    }

    console.time('LoadModel');
    PluginCommands.State.Snapshots.Clear(ctx);
    PluginCommands.State.RemoveObject(ctx, { state: ctx.state.data, ref: StateTransform.RootRef });
    adjustPluginProps(ctx);

    for (const file of params.files) {
        try {
            const info = getFileNameInfo(file.file!.name);
            const isBinary = ctx.dataFormats.binaryExtensions.has(info.ext);
            const { data } = await ctx.builders.data.readFile({ file, isBinary });
            await createHierarchy(ctx, data.ref);
        } catch (e) {
            console.error(e);
            ctx.log.error(`Error opening file '${file.name}'`);
        }
    }
    console.timeEnd('LoadModel');
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

export class ExampleControls extends PluginUIComponent {
    componentDidMount() {

    }

    render() {
        return <div style={{ margin: '5px' }}>
            <ApplyActionControl state={this.plugin.state.data} action={LoadExample} nodeRef={this.plugin.state.data.tree.root.ref} applyLabel={'Load'} hideHeader />
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
                <Button icon={GetAppSvg} onClick={this.downloadToFileZip} title='Download the state.'>
                    Download
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