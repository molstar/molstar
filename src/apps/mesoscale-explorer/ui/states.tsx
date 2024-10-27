/**
 * Copyright (c) 2022-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';
import { MmcifProvider } from '../../../mol-plugin-state/formats/trajectory';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { Button, ExpandGroup, IconButton } from '../../../mol-plugin-ui/controls/common';
import { GetAppSvg, HelpOutlineSvg, MagicWandSvg, TourSvg, Icon, OpenInBrowserSvg } from '../../../mol-plugin-ui/controls/icons';
import { CollapsableControls, PluginUIComponent } from '../../../mol-plugin-ui/base';
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
import { MesoscaleState, MesoscaleStateObject, setGraphicsCanvas3DProps, updateStyle } from '../data/state';
import { isTimingMode } from '../../../mol-util/debug';
import { now } from '../../../mol-util/now';

function adjustPluginProps(ctx: PluginContext) {
    const customState = ctx.customState as MesoscaleExplorerState;

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
                    transparentThreshold: 0.4,
                }
            },
            shadow: {
                name: 'on',
                params: {
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
            },
        },
        illumination: {
            enabled: customState.illumination,
            firstStepSize: 0.1,
            rayDistance: 1024,
        },
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
    let startTime = 0;
    if (isTimingMode) {
        startTime = now();
    }
    if (type === 'molx' || type === 'molj') {
        const customState = ctx.customState as MesoscaleExplorerState;
        delete customState.stateRef;
        customState.stateCache = {};
        ctx.managers.asset.clear();

        await PluginCommands.State.Snapshots.Clear(ctx);
        await PluginCommands.State.Snapshots.OpenUrl(ctx, { url, type });

        const cell = ctx.state.data.selectQ(q => q.ofType(MesoscaleStateObject))[0];
        if (!cell) throw new Error('Missing MesoscaleState');

        customState.stateRef = cell.transform.ref;
        customState.graphicsMode = cell.obj?.data.graphics || customState.graphicsMode;
    } else {
        await reset(ctx);
        const isBinary = type === 'bcif';
        const data = await ctx.builders.data.download({ url, isBinary });
        await createHierarchy(ctx, data.ref);
    }
    if (isTimingMode) {
        const endTime = now();
        // Calculate the elapsed time
        const timeTaken = endTime - startTime;
        console.log(`Model loaded in ${timeTaken} milliseconds`);
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
    let url: string;
    // 4 character PDB id, TODO: support extended PDB ID
    if (id.match(/^[1-9][A-Z0-9]{3}$/i) !== null) {
        url = `https://pdb-dev.wwpdb.org/bcif/${id.toLowerCase()}.bcif`;
    } else {
        const nId = id.toUpperCase().startsWith('PDBDEV_') ? id : `PDBDEV_${id.padStart(8, '0')}`;
        url = `https://pdb-dev.wwpdb.org/bcif/${nId.toUpperCase()}.bcif`;
    }
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
        return <div id='database' style={{ margin: '5px' }}>
            <ApplyActionControl state={this.plugin.state.data} action={LoadDatabase} nodeRef={this.plugin.state.data.tree.root.ref} applyLabel={'Load'} hideHeader />
        </div>;
    }
}

export class LoaderControls extends PluginUIComponent {
    componentDidMount() {

    }

    render() {
        return <div id='loader' style={{ margin: '5px' }}>
            <ApplyActionControl state={this.plugin.state.data} action={LoadModel} nodeRef={this.plugin.state.data.tree.root.ref} applyLabel={'Load'} hideHeader />
        </div>;
    }
}

export class ExampleControls extends PluginUIComponent {
    componentDidMount() {

    }

    render() {
        return <div id='example' style={{ margin: '5px' }}>
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
        return <div id='session' style={{ margin: '5px' }}>
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
            <div id='snaplist' style={{ marginBottom: '10px' }}>
                <LocalStateSnapshotList />
            </div>
            <div id='snap' style={{ marginBottom: '10px' }}>
                <LocalStateSnapshots />
            </div>

            <div id='snapoption' style={{ marginBottom: '10px' }}>
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
        // setup the tour of the interface
        const driver = (this.plugin.customState as MesoscaleExplorerState).driver;
        if (!driver) return;

        driver.setSteps([
            // Left panel
            { element: '#explorerinfo', popover: { title: 'Explorer Header Info', description: 'This section displays the explorer header with version information, documentation access, and tour navigation. Use the right and left arrow keys to navigate the tour.', side: 'left', align: 'start' } },
            { element: '#database', popover: { title: 'Import from PDB', description: 'Load structures directly from PDB and PDB-DEV databases.', side: 'bottom', align: 'start' } },
            { element: '#loader', popover: { title: 'Import from File', description: 'Load local files (.molx, .molj, .zip, .cif, .bcif) using this option.', side: 'bottom', align: 'start' } },
            { element: '#example', popover: { title: 'Example Models and Tours', description: 'Select from a range of example models and tours provided.', side: 'left', align: 'start' } },
            { element: '#session', popover: { title: 'Session Management', description: 'Download the current session in .molx format.', side: 'top', align: 'start' } },
            { element: '#snaplist', popover: { title: 'Snapshot List', description: 'View and manage the list of snapshots. You can reorder them and edit their titles, keys, and descriptions. Snapshot states cannot be edited.', side: 'right', align: 'start' } },
            { element: '#snap', popover: { title: 'Add Snapshot', description: 'Save the current state (e.g., camera position, color, visibility, etc.) in a snapshot with an optional title, key, and description.', side: 'right', align: 'start' } },
            { element: '#snapoption', popover: { title: 'Snapshot Options', description: 'These options are saved in the snapshot. Set them before adding a snapshot to see their effect during animation playback.', side: 'right', align: 'start' } },
            { element: '#exportanimation', popover: { title: 'Export Animation', description: 'Create movies or scenes with rocking, rotating, or snapshots animations.', side: 'right', align: 'start' } },
            { element: '#viewportsettings', popover: { title: 'Viewport Settings', description: 'Advanced settings for the renderer and trackball.', side: 'right', align: 'start' } },
            // Viewport
            { element: '#snapinfo', popover: { title: 'Snapshot Description', description: 'Save the current state (e.g., camera position, color, visibility, etc.) in a snapshot with an optional title, key, and description.', side: 'right', align: 'start' } },
            { element: '#snapinfoctrl', popover: { title: 'Snapshot Description Control', description: 'Control the visibility and text size of the snapshot description widget.', side: 'right', align: 'start' } },
            // Right panel
            { element: '#modelinfo', popover: { title: 'Model Information', description: 'Summary information about the model, if available.', side: 'right', align: 'start' } },
            { element: '#selestyle', popover: { title: 'Selection Style', description: 'Choose the rendering style for entity selection accessed via Shift/Ctrl mouse. Options include: Color & Outline, Color, Outline.', side: 'right', align: 'start' } },
            { element: '#seleinfo', popover: { title: 'Selection List', description: 'View the current list of selected entities.', side: 'right', align: 'start' } },
            { element: '#measurements', popover: { title: 'Measurements', description: 'Use this widget to create labels, measure distances, angles, dihedral orientations, and planes for the selected entities.', side: 'right', align: 'start' } },
            { element: '#quickstyles', popover: { title: 'Quick Styles', description: 'Change between a selection of style presets.', side: 'right', align: 'start' } },
            { element: '#graphicsquality', popover: { title: 'Graphics Quality', description: 'Adjust the overall graphics quality. Lower quality improves performance. Options are: Ultra, Quality (Default), Balanced, Performance, Custom. Custom settings use the Culling & LOD values set in the Tree.', side: 'right', align: 'start' } },
            { element: '#searchtree', popover: { title: 'Search', description: 'Filter the entity tree based on your queries.', side: 'right', align: 'start' } },
            { element: '#grouptree', popover: { title: 'Group By', description: 'Change the grouping of the hierarchy tree, e.g., group by instance or by compartment.', side: 'right', align: 'start' } },
            { element: '#tree', popover: { title: 'Tree Hierarchy', description: 'View the hierarchical tree of entity types in the model.', side: 'right', align: 'start' } },
            { element: '#focusinfo', popover: { title: 'Selection Description', description: 'Detailed information about the current selection, if present in the loaded file.', side: 'right', align: 'start' } },
            { popover: { title: 'Happy Exploring!', description: 'That’s all! Go ahead and start exploring or creating mesoscale tours.' } }
        ]);
        driver.refresh();
    };

    openHelp = () => {
        // open a new page with the documentation
        window.open('https://molstar.org/me-docs/', '_blank');
    };

    toggleHelp = () => {
        const driver = (this.plugin.customState as MesoscaleExplorerState).driver;
        if (!driver || !driver.hasNextStep()) {
            this.setupDriver();
        }
        this.setState({ showHelp: !this.state.showHelp }, () => {
            if (this.state.showHelp && driver) {
                driver.drive(); // start at 0
            }
        });
    };

    render() {
        const driver = (this.plugin.customState as MesoscaleExplorerState).driver;
        if (!driver) return;

        const help = <IconButton svg={HelpOutlineSvg} toggleState={false} small onClick={this.openHelp} title='Open the Documentation' />;
        const tour = <IconButton svg={TourSvg} toggleState={false} small onClick={this.toggleHelp} title='Start the interactive tour' />;
        return <>
            <div id='explorerinfo' style={{ display: 'flex', alignItems: 'center', padding: '4px 0 4px 8px' }} className='msp-help-text'>
                <h2 style={{ flexGrow: 1 }}>Mol* Mesoscale Explorer</h2>
                {tour}{help}
            </div>
        </>;
    }
}


export class MesoQuickStylesControls extends CollapsableControls {
    defaultState() {
        return {
            isCollapsed: true,
            header: 'Quick Styles',
            brand: { accent: 'gray' as const, svg: MagicWandSvg }
        };
    }

    renderControls() {
        return <>
            <MesoQuickStyles />
        </>;
    }
}

export class MesoQuickStyles extends PluginUIComponent {
    async default() {
        if (!this.plugin.canvas3d) return;
        const p = this.plugin.canvas3d.props;
        this.plugin.canvas3d.setProps({
            renderer: {
                exposure: 1.1,
            },
            postprocessing: {
                ...p.postprocessing,
                shadow: {
                    name: 'on',
                    params: {
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
                },
                dof: { name: 'off', params: {} },
            }
        });
        await updateStyle(this.plugin, {
            ignoreLight: true,
            material: { metalness: 0, roughness: 1.0, bumpiness: 0 },
            celShaded: false,
            illustrative: false,
        });
    }

    async celshading() {
        if (!this.plugin.canvas3d) return;
        const p = this.plugin.canvas3d.props;
        this.plugin.canvas3d.setProps({
            renderer: {
                exposure: 1.5,
            },
            postprocessing: {
                ...p.postprocessing,
                shadow: {
                    name: 'on',
                    params: {
                        maxDistance: 256,
                        steps: 64,
                        tolerance: 1.0,
                    }
                },
                outline: { name: 'off', params: {} },
                dof: { name: 'off', params: {} },
            }
        });
        await updateStyle(this.plugin, {
            ignoreLight: false,
            material: { metalness: 0, roughness: 1.0, bumpiness: 0 },
            celShaded: true,
            illustrative: false,
        });
    }

    async shinyDof() {
        if (!this.plugin.canvas3d) return;
        const p = this.plugin.canvas3d.props;
        this.plugin.canvas3d.setProps({
            renderer: {
                exposure: 1.1,
            },
            postprocessing: {
                ...p.postprocessing,
                shadow: {
                    name: 'on',
                    params: {
                        maxDistance: 256,
                        steps: 64,
                        tolerance: 1.0,
                    }
                },
                outline: { name: 'off', params: {} },
                dof: {
                    name: 'on',
                    params: {
                        blurSize: 9,
                        blurSpread: 1.0,
                        inFocus: 0.0,
                        PPM: 200.0,
                        center: 'camera-target',
                        mode: 'sphere',
                    }
                }
            }
        });
        await updateStyle(this.plugin, {
            ignoreLight: false,
            material: { metalness: 0, roughness: 0.2, bumpiness: 0 },
            celShaded: false,
            illustrative: false,
        });
    }

    async illustrative() {
        if (!this.plugin.canvas3d) return;
        const p = this.plugin.canvas3d.props;
        this.plugin.canvas3d.setProps({
            renderer: {
                exposure: 1.5,
            },
            postprocessing: {
                ...p.postprocessing,
                shadow: {
                    name: 'on',
                    params: {
                        maxDistance: 256,
                        steps: 64,
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
                },
                dof: { name: 'off', params: {} },
            }
        });
        await updateStyle(this.plugin, {
            ignoreLight: true,
            material: { metalness: 0, roughness: 1.0, bumpiness: 0 },
            celShaded: false,
            illustrative: true,
        });
    }

    async shiny() {
        if (!this.plugin.canvas3d) return;
        const p = this.plugin.canvas3d.props;
        this.plugin.canvas3d.setProps({
            renderer: {
                exposure: 1.5,
            },
            postprocessing: {
                ...p.postprocessing,
                shadow: { name: 'off', params: {} },
                outline: { name: 'off', params: {} },
                dof: { name: 'off', params: {} },
            }
        });
        await updateStyle(this.plugin, {
            ignoreLight: false,
            material: { metalness: 0, roughness: 0.2, bumpiness: 0 },
            celShaded: false,
            illustrative: false,
        });
    }

    async stylized() {
        if (!this.plugin.canvas3d) return;
        const p = this.plugin.canvas3d.props;
        this.plugin.canvas3d.setProps({
            renderer: {
                exposure: 1.1,
            },
            postprocessing: {
                ...p.postprocessing,
                shadow: {
                    name: 'on',
                    params: {
                        maxDistance: 256,
                        steps: 64,
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
                },
                dof: { name: 'off', params: {} },
            }
        });
        await updateStyle(this.plugin, {
            ignoreLight: false,
            material: { metalness: 0, roughness: 0.2, bumpiness: 0 },
            celShaded: false,
            illustrative: true,
        });
    }

    render() {
        return <>
            <div className='msp-flex-row'>
                <Button noOverflow title='Applies default representation preset and sets outline and occlusion effects to default' onClick={() => this.default()} style={{ width: 'auto' }}>
                    Default
                </Button>
                <Button noOverflow title='Applies celShading' onClick={() => this.celshading()} style={{ width: 'auto' }}>
                    Cel-shaded
                </Button>
                <Button noOverflow title='Applies illustrative colors preset' onClick={() => this.illustrative()} style={{ width: 'auto' }}>
                    Illustrative
                </Button>
            </div>
            <div className='msp-flex-row'>
                <Button noOverflow title='Apply shiny material to default' onClick={() => this.shiny()} style={{ width: 'auto' }}>
                    Shiny
                </Button>
                <Button noOverflow title='Enable shiny material, outline, and illustrative colors' onClick={() => this.stylized()} style={{ width: 'auto' }}>
                    Shiny-Illustrative
                </Button>
                <Button noOverflow title='Enable DOF and shiny material' onClick={() => this.shinyDof()} style={{ width: 'auto' }}>
                    Shiny-DOF
                </Button>
            </div>
        </>;
    }
}
