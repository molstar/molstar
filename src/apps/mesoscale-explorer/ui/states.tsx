/**
 * Copyright (c) 2022-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Markdown from 'react-markdown';
import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';
import { MmcifProvider } from '../../../mol-plugin-state/formats/trajectory';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { Button, ExpandGroup, IconButton } from '../../../mol-plugin-ui/controls/common';
import { GetAppSvg, HelpOutlineSvg, TourSvg, Icon, OpenInBrowserSvg } from '../../../mol-plugin-ui/controls/icons';
import { ApplyActionControl } from '../../../mol-plugin-ui/state/apply-action';
import { LocalStateSnapshotList, LocalStateSnapshotParams, LocalStateSnapshots } from '../../../mol-plugin-ui/state/snapshots';
import { PluginCommands } from '../../../mol-plugin/commands';
import { PluginContext } from '../../../mol-plugin/context';
import { StateAction, StateObjectRef, StateTransform } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { Color } from '../../../mol-util/color/color';
import { getFileNameInfo } from '../../../mol-util/file-info';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ExampleEntry, MesoscaleExplorerState } from '../app';
import { createCellpackHierarchy } from '../data/cellpack/preset';
import { createGenericHierarchy } from '../data/generic/preset';
import { createMmcifHierarchy } from '../data/mmcif/preset';
import { createPetworldHierarchy } from '../data/petworld/preset';
import { MesoscaleState, MesoscaleStateObject, setGraphicsCanvas3DProps } from '../data/state';
import { driver } from 'driver.js';
import 'driver.js/dist/driver.css';

const driverObj = driver();

function adjustPluginProps(ctx: PluginContext) {
    ctx.managers.interactivity.setProps({ granularity: 'chain' });
    ctx.canvas3d?.setProps({
        multiSample: { mode: 'off' },
        cameraClipping: { far: false, minNear: 50 },
        sceneRadiusFactor: 2,
        renderer: {
            colorMarker: true,
            highlightColor: Color(0xffffff),
            highlightStrength: 0,
            selectColor: Color(0xffffff),
            selectStrength: 0,
            dimColor: Color(0xffffff),
            dimStrength: 1,
            markerPriority: 2,
            interiorColorFlag: false,
            interiorDarkening: 0.15,
            exposure: 1.1,
            xrayEdgeFalloff: 3,
        },
        marking: {
            enabled: true,
            highlightEdgeColor: Color(0x999999),
            selectEdgeColor: Color(0xffff00),
            highlightEdgeStrength: 1,
            selectEdgeStrength: 1,
            ghostEdgeStrength: 1,
            innerEdgeFactor: 2.5,
            edgeScale: 2,
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
                    blurDepthBias: 0.5,
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
                    includeTransparent: false,
                }
            }
        }
    });

    const { graphics } = MesoscaleState.get(ctx);
    setGraphicsCanvas3DProps(ctx, graphics);
}

async function createHierarchy(ctx: PluginContext, ref: string) {
    const parsed = await MmcifProvider.parse(ctx, ref);

    const tr = StateObjectRef.resolveAndCheck(ctx.state.data, parsed.trajectory)?.obj?.data;
    if (!tr) throw new Error('no trajectory');

    if (!MmcifFormat.is(tr.representative.sourceData)) {
        throw new Error('not mmcif');
    }

    const { frame, db } = tr.representative.sourceData.data;

    let hasCellpackAssemblyMethodDetails = false;
    const { method_details } = db.pdbx_struct_assembly;
    for (let i = 0, il = method_details.rowCount; i < il; ++i) {
        if (method_details.value(i).toUpperCase() === 'CELLPACK') {
            hasCellpackAssemblyMethodDetails = true;
            break;
        }
    }

    if (frame.categories.pdbx_model) {
        await createPetworldHierarchy(ctx, parsed.trajectory);
    } else if (
        frame.header.toUpperCase().includes('CELLPACK') ||
        hasCellpackAssemblyMethodDetails
    ) {
        await createCellpackHierarchy(ctx, parsed.trajectory);
    } else {
        await createMmcifHierarchy(ctx, parsed.trajectory);
    }
}

async function reset(ctx: PluginContext) {
    const customState = ctx.customState as MesoscaleExplorerState;
    delete customState.stateRef;
    customState.stateCache = {};
    ctx.managers.asset.clear();

    await PluginCommands.State.Snapshots.Clear(ctx);
    await PluginCommands.State.RemoveObject(ctx, { state: ctx.state.data, ref: StateTransform.RootRef });

    await MesoscaleState.init(ctx);
    adjustPluginProps(ctx);
}

export async function loadExampleEntry(ctx: PluginContext, entry: ExampleEntry) {
    const { url, type } = entry;
    await loadUrl(ctx, url, type);
    MesoscaleState.set(ctx, {
        description: entry.description || entry.label,
        link: entry.link,
    });
}

export async function loadUrl(ctx: PluginContext, url: string, type: 'molx' | 'molj' | 'cif' | 'bcif') {
    if (type === 'molx' || type === 'molj') {
        await PluginCommands.State.Snapshots.OpenUrl(ctx, { url, type });
    } else {
        await reset(ctx);
        const isBinary = type === 'bcif';
        const data = await ctx.builders.data.download({ url, isBinary });
        await createHierarchy(ctx, data.ref);
    }
}

export async function loadPdb(ctx: PluginContext, id: string) {
    await reset(ctx);
    const url = `https://models.rcsb.org/${id.toUpperCase()}.bcif`;
    const data = await ctx.builders.data.download({ url, isBinary: true });
    await createHierarchy(ctx, data.ref);
}

export async function loadPdbDev(ctx: PluginContext, id: string) {
    await reset(ctx);
    const nId = id.toUpperCase().startsWith('PDBDEV_') ? id : `PDBDEV_${id.padStart(8, '0')}`;
    const url = `https://pdb-dev.wwpdb.org/bcif/${nId.toUpperCase()}.bcif`;
    const data = await ctx.builders.data.download({ url, isBinary: true });
    await createHierarchy(ctx, data.ref);
}

//

export const LoadDatabase = StateAction.build({
    display: { name: 'Database', description: 'Load from Database' },
    params: (a, ctx: PluginContext) => {
        return {
            source: PD.Select('pdb', PD.objectToOptions({ pdb: 'PDB', pdbDev: 'PDB-Dev' })),
            entry: PD.Text(''),
        };
    },
    from: PluginStateObject.Root
})(({ params }, ctx: PluginContext) => Task.create('Loading from database...', async taskCtx => {
    if (params.source === 'pdb') {
        await loadPdb(ctx, params.entry);
    } else if (params.source === 'pdbDev') {
        await loadPdbDev(ctx, params.entry);
    }
}));

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
    await loadExampleEntry(ctx, entries[params.entry]);
}));

export const LoadModel = StateAction.build({
    display: { name: 'Load', description: 'Load a model' },
    params: {
        files: PD.FileList({ accept: '.cif,.bcif,.cif.gz,.bcif.gz,.zip', multiple: true, description: 'mmCIF or Cellpack- or Petworld-style cif file.', label: 'File(s)' }),
    },
    from: PluginStateObject.Root
})(({ params }, ctx: PluginContext) => Task.create('Loading model...', async taskCtx => {
    if (params.files === null || params.files.length === 0) {
        ctx.log.error('No file(s) selected');
        return;
    }

    await reset(ctx);

    const firstFile = params.files[0];
    const firstInfo = getFileNameInfo(firstFile.file!.name);

    if (firstInfo.name.endsWith('zip')) {
        try {
            await createGenericHierarchy(ctx, firstFile);
        } catch (e) {
            console.error(e);
            ctx.log.error(`Error opening file '${firstFile.name}'`);
        }
    } else {
        for (const file of params.files) {
            try {
                const info = getFileNameInfo(file.file!.name);
                if (!['cif', 'bcif'].includes(info.ext)) continue;

                const isBinary = ctx.dataFormats.binaryExtensions.has(info.ext);
                const { data } = await ctx.builders.data.readFile({ file, isBinary });
                await createHierarchy(ctx, data.ref);
            } catch (e) {
                console.error(e);
                ctx.log.error(`Error opening file '${file.name}'`);
            }
        }
    }
}));

//

export class DatabaseControls extends PluginUIComponent {
    componentDidMount() {

    }

    render() {
        return <div style={{ margin: '5px' }}>
            <ApplyActionControl state={this.plugin.state.data} action={LoadDatabase} nodeRef={this.plugin.state.data.tree.root.ref} applyLabel={'Load'} hideHeader />
        </div>;
    }
}

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

export async function openState(ctx: PluginContext, file: File) {
    const customState = ctx.customState as MesoscaleExplorerState;
    delete customState.stateRef;
    customState.stateCache = {};
    ctx.managers.asset.clear();

    await PluginCommands.State.Snapshots.Clear(ctx);
    await PluginCommands.State.Snapshots.OpenFile(ctx, { file });

    const cell = ctx.state.data.selectQ(q => q.ofType(MesoscaleStateObject))[0];
    if (!cell) throw new Error('Missing MesoscaleState');

    customState.stateRef = cell.transform.ref;
    customState.graphicsMode = cell.obj?.data.graphics || customState.graphicsMode;
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

        openState(this.plugin, e.target.files[0]);
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

export class ExplorerInfo extends PluginUIComponent<{}, { isDisabled: boolean, showHelp: boolean }> {
    state = {
        isDisabled: false,
        showHelp: false
    };

    componentDidMount() {
        this.subscribe(this.plugin.state.data.behaviors.isUpdating, v => {
            this.setState({ isDisabled: v });
        });

        this.subscribe(this.plugin.state.events.cell.stateUpdated, e => {
            if (!this.state.isDisabled && MesoscaleState.has(this.plugin) && MesoscaleState.ref(this.plugin) === e.ref) {
                this.forceUpdate();
            }
        });
    }

    setupDriver = () => {
        // // setup the tour of the interface
        driverObj.setSteps([
            // left panel
            { element: '#explorerinfo', popover: { title: 'Explorer Header Info', description: 'Here is the header of explorer showing the version and access to the documentation and this tour. You can right and left keyboard key to navigate the tour.', side: 'left', align: 'start' } },
            { element: '#database', popover: { title: 'Import from PDB', description: 'Here you can load structure directly from PDB and PDB-DEV', side: 'bottom', align: 'start' } },
            { element: '#loader', popover: { title: 'Importing from file', description: 'Here is the local file loader (.molx, .molj, .zip, .cif and .bcif)', side: 'bottom', align: 'start' } },
            { element: '#example', popover: { title: 'Example', description: 'You can select one of the example model or tour we provide', side: 'left', align: 'start' } },
            { element: '#session', popover: { title: 'Session', description: 'Download the current session in .molx format', side: 'top', align: 'start' } },
            { element: '#snaplist', popover: { title: 'Snapshot List', description: 'Show the current list of snapshot, you can reorder them, and edit the snapshot title, key and description. You cannot edit the state of the snapshot.', side: 'right', align: 'start' } },
            { element: '#snap', popover: { title: 'Snapshot Add', description: 'Save the current state (e.g. everything you see, camera, color, visiblity, etc..) in a snapshot with optional title, key and description', side: 'right', align: 'start' } },
            { element: '#snapoption', popover: { title: 'Snapshot options', description: 'Theses options olso get saved in the snapshot, set them prior to add a snapshot to see their effect when playing the animation', side: 'right', align: 'start' } },
            { element: '#exportanimation', popover: { title: 'Export Animation', description: 'Theses options olso get saved in the snapshot, set them prior to add a snapshot to see their effect when playing the animation', side: 'right', align: 'start' } },
            // { element: '#viewportsettings', popover: { title: 'Viewport settings', description: 'Theses options olso get saved in the snapshot, set them prior to add a snapshot to see their effect when playing the animation', side: 'right', align: 'start' } },
            // { element: '#behaviorsettings', popover: { title: 'Behavior settings', description: 'Theses options olso get saved in the snapshot, set them prior to add a snapshot to see their effect when playing the animation', side: 'right', align: 'start' } },
            // viewport
            { element: '#snapinfo', popover: { title: 'Snapshot Description', description: 'Save the current state (e.g. everything you see, camera, color, visiblity, etc..) in a snapshot with optional title, key and description', side: 'right', align: 'start' } },
            // { element: '#snapinfoctrl', popover: { title: 'Snapshot Description Control', description: 'Theses options olso get saved in the snapshot, set them prior to add a snapshot to see their effect when playing the animation', side: 'right', align: 'start' } },
            // { element: '#canvainfo', popover: { title: 'Selection Description', description: 'Save the current state (e.g. everything you see, camera, color, visiblity, etc..) in a snapshot with optional title, key and description', side: 'right', align: 'start' } },
            // right panel
            { element: '#modelinfo', popover: { title: 'Model informations', description: 'Theses options olso get saved in the snapshot, set them prior to add a snapshot to see their effect when playing the animation', side: 'right', align: 'start' } },
            { element: '#selestyle', popover: { title: 'Seletion style', description: 'Save the current state (e.g. everything you see, camera, color, visiblity, etc..) in a snapshot with optional title, key and description', side: 'right', align: 'start' } },
            { element: '#seleinfo', popover: { title: 'Seletion List', description: 'Save the current state (e.g. everything you see, camera, color, visiblity, etc..) in a snapshot with optional title, key and description', side: 'right', align: 'start' } },
            // { element: '#measurement', popover: { title: 'Measurements', description: 'Theses options olso get saved in the snapshot, set them prior to add a snapshot to see their effect when playing the animation', side: 'right', align: 'start' } },
            { element: '#graphicsquality', popover: { title: 'Graphics quality', description: 'Theses options olso get saved in the snapshot, set them prior to add a snapshot to see their effect when playing the animation', side: 'right', align: 'start' } },
            { element: '#searchtree', popover: { title: 'Search', description: 'Theses options olso get saved in the snapshot, set them prior to add a snapshot to see their effect when playing the animation', side: 'right', align: 'start' } },
            { element: '#grouptree', popover: { title: 'Group', description: 'Theses options olso get saved in the snapshot, set them prior to add a snapshot to see their effect when playing the animation', side: 'right', align: 'start' } },
            { element: '#tree', popover: { title: 'Tree hierarchy', description: 'Theses options olso get saved in the snapshot, set them prior to add a snapshot to see their effect when playing the animation', side: 'right', align: 'start' } },
            { popover: { title: 'Happy Exploring', description: 'And that is all, go ahead and start exploring or creating mesoscale tours.' } }
        ]
        );
        driverObj.refresh();
    };

    openHelp = () => {
        // open a new page with the documentation
        window.open('https://mesoscope.scripps.edu/explorer/docs/', '_blank');
    };

    toggleHelp = () => {
        if (!driverObj || !driverObj.hasNextStep()) {
            this.setupDriver();
        }
        this.setState({ showHelp: !this.state.showHelp }, () => {
            if (this.state.showHelp && driverObj) {
                driverObj.drive(); // start at 0
            }
        });
    };

    render() {
        const legend = `## Welcome to Mol* Mesoscale Explorer

Click (?) below to open the [Documentation](https://mesoscope.scripps.edu/explorer/docs/) in a new tab. 

Click &#9873; below to start a tour of the interface with contextual help.

v0.1`;
        const help = <IconButton svg={HelpOutlineSvg} toggleState={false} small onClick={this.openHelp} />;
        const tour = <IconButton svg={TourSvg} toggleState={false} small onClick={this.toggleHelp} />;
        return <>
            <div id='explorerinfo' className='msp-help-text'>
                <Markdown>{legend}</Markdown>
                {tour}{help}
            </div>
        </>;
    }
}


