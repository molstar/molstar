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
import { ColorNames } from '../../../mol-util/color/names';
import { getFileInfo } from '../../../mol-util/file-info';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { createCellpackHierarchy } from '../data/cellpack/preset';
import { createPetworldHierarchy } from '../data/petworld/preset';

function adjustPluginProps(ctx: PluginContext) {
    ctx.managers.interactivity.setProps({ granularity: 'chain' });
    ctx.canvas3d?.setProps({
        multiSample: { mode: 'off' },
        cameraClipping: { far: false },
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
                    levels: [
                        { radius: 1, bias: 0.4 },
                        { radius: 5, bias: 0.6 },
                        { radius: 8, bias: 1.0 },
                    ],
                    distanceFactor: 10,
                    blurKernelSize: 11,
                    resolutionScale: 1,
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
                    threshold: 0.33,
                    color: ColorNames.black,
                    includeTransparent: true,
                }
            }
        }
    });
}

export const LoadModel = StateAction.build({
    display: { name: 'Load', description: 'Open or download a model' },
    params: {
        files: PD.FileList({ accept: '.cif,.bcif', multiple: true, description: 'Cellpack- or Petworld-style cif file.', label: 'File(s)' }),
    },
    from: PluginStateObject.Root
})(({ params }, ctx: PluginContext) => Task.create('Model Loader', async taskCtx => {
    if (params.files === null) {
        ctx.log.error('No file(s) selected');
        return;
    }

    PluginCommands.State.RemoveObject(ctx, { state: ctx.state.data, ref: StateTransform.RootRef });

    adjustPluginProps(ctx);

    for (const file of params.files) {
        try {
            const info = getFileInfo(file.file!);
            const isBinary = ctx.dataFormats.binaryExtensions.has(info.ext);
            const { data } = await ctx.builders.data.readFile({ file, isBinary });
            const parsed = await MmcifProvider.parse(ctx, data);

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
        } catch (e) {
            console.error(e);
            ctx.log.error(`Error opening file '${file.name}'`);
        }
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
