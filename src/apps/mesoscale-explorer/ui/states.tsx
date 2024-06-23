/**
 * Copyright (c) 2022-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';
import { MmcifProvider } from '../../../mol-plugin-state/formats/trajectory';
import { StructureComponentManager } from '../../../mol-plugin-state/manager/structure/component';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { CollapsableControls, PluginUIComponent } from '../../../mol-plugin-ui/base';
import { Button, ExpandGroup } from '../../../mol-plugin-ui/controls/common';
import { GetAppSvg, Icon, MagicWandSvg, OpenInBrowserSvg } from '../../../mol-plugin-ui/controls/icons';
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
import { MesoscaleState, MesoscaleStateObject, setGraphicsCanvas3DProps, updateColors } from '../data/state';

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

export class MesoQuickStylesControls extends CollapsableControls {
    defaultState() {
        return {
            isCollapsed: true,
            header: 'Styles',
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
    state = {
        celShaded: false,
    };
    default_color_values = {
        type: 'group-generate',
        illustrative: false,
        value: [1, 1, 1, 1],
        variability: 20,
        shift: 0,
        lightness: 0,
        alpha: 1,
        emissive: 0
    };
    illustrative_color_values = {
        type: 'group-generate',
        illustrative: true,
        value: [1, 1, 1, 1],
        variability: 20,
        shift: 0,
        lightness: 0,
        alpha: 1,
        emissive: 0
    };
    async default() {
        if (!this.plugin.canvas3d) return;
        this.plugin.canvas3d.setProps({
            renderer: {
                exposure: 1.1,
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
                },
                dof: { name: 'off', params: {} },
            }
        });

        const loptions = { ignoreLight: true, materialStyle: { metalness: 0, roughness: 1.0, bumpiness: 0 } };
        const options = { ...loptions, celShaded: false, };
        await this.plugin.managers.structure.component.setOptions(loptions as StructureComponentManager.Options);
        await updateColors(this.plugin, this.default_color_values, options);
    }

    async cellshading() {
        if (!this.plugin.canvas3d) return;
        this.plugin.canvas3d.setProps({
            renderer: {
                exposure: 1.5,
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
                        bias: 2.0,
                        blurKernelSize: 11,
                        blurDepthBias: 0.5,
                        resolutionScale: 1,
                        color: Color(0x000000),
                    }
                },
                shadow: {
                    name: 'on',
                    params: {
                        bias: 0.4,
                        maxDistance: 256,
                        steps: 64,
                        tolerance: 1.0,
                    }
                },
                outline: { name: 'off', params: {} },
                dof: { name: 'off', params: {} },
            }
        });
        // ignore Light
        const loptions = { ignoreLight: false, materialStyle: { metalness: 0, roughness: 1.0, bumpiness: 0 } };
        const options = { ...loptions, celShaded: true, };
        await this.plugin.managers.structure.component.setOptions(loptions as StructureComponentManager.Options);
        await updateColors(this.plugin, this.illustrative_color_values, options);
    }

    async stylizeddof() {
        if (!this.plugin.canvas3d) return;
        this.plugin.canvas3d.setProps({
            renderer: {
                exposure: 1.1,
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
                        bias: 1.3,
                        blurKernelSize: 11,
                        blurDepthBias: 0.5,
                        resolutionScale: 1,
                        color: Color(0x000000),
                    }
                },
                shadow: {
                    name: 'on',
                    params: {
                        bias: 0.4,
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
        // ignore Light
        const loptions = { ignoreLight: false, materialStyle: { metalness: 0, roughness: 0.2, bumpiness: 0 } };
        const options = { ...loptions, celShaded: false };
        await this.plugin.managers.structure.component.setOptions(loptions as StructureComponentManager.Options);
        await updateColors(this.plugin, this.default_color_values, options);
    }

    async illustrative() {
        if (!this.plugin.canvas3d) return;
        this.plugin.canvas3d.setProps({
            renderer: {
                exposure: 1.5,
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
                        bias: 2.0,
                        blurKernelSize: 11,
                        blurDepthBias: 0.5,
                        resolutionScale: 1,
                        color: Color(0x000000),
                    }
                },
                shadow: {
                    name: 'on',
                    params: {
                        bias: 0.4,
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
        // ignore Light
        const loptions = { ignoreLight: true, materialStyle: { metalness: 0, roughness: 1.0, bumpiness: 0 } };
        const options = { ...loptions, celShaded: false, };
        await this.plugin.managers.structure.component.setOptions(loptions as StructureComponentManager.Options);
        await updateColors(this.plugin, this.illustrative_color_values, options);
    }

    async metal() {
        if (!this.plugin.canvas3d) return;
        this.plugin.canvas3d.setProps({
            renderer: {
                exposure: 1.5,
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
                        bias: 1.3,
                        blurKernelSize: 11,
                        blurDepthBias: 0.5,
                        resolutionScale: 1,
                        color: Color(0x000000),
                    }
                },
                shadow: { name: 'off', params: {} },
                outline: { name: 'off', params: {} },
                dof: { name: 'off', params: {} },
            }
        });
        // ignore Light
        const loptions = { ignoreLight: false, materialStyle: { metalness: 0, roughness: 0.2, bumpiness: 0 } };
        const options = { ...loptions, celShaded: false };
        await this.plugin.managers.structure.component.setOptions(loptions as StructureComponentManager.Options);
        await updateColors(this.plugin, this.default_color_values, options);
    }

    async stylized() {
        if (!this.plugin.canvas3d) return;
        this.plugin.canvas3d.setProps({
            renderer: {
                exposure: 1.1,
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
                        bias: 1.3,
                        blurKernelSize: 11,
                        blurDepthBias: 0.5,
                        resolutionScale: 1,
                        color: Color(0x000000),
                    }
                },
                shadow: {
                    name: 'on',
                    params: {
                        bias: 0.4,
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
        // ignore Light
        const loptions = { ignoreLight: false, materialStyle: { metalness: 0, roughness: 0.2, bumpiness: 0 } };
        const options = { ...loptions, celShaded: false };
        await this.plugin.managers.structure.component.setOptions(loptions as StructureComponentManager.Options);
        await updateColors(this.plugin, this.illustrative_color_values, options);
    }

    render() {
        return <>
            <div className='msp-flex-row'>
                <Button noOverflow title='Applies default representation preset and sets outline and occlusion effects to default' onClick={() => this.default()} style={{ width: 'auto' }}>
                    Default
                </Button>
                <Button noOverflow title='Applies illustrative representation preset and Stylize it with cellShading' onClick={() => this.cellshading()} style={{ width: 'auto' }}>
                    Cell-shaded
                </Button>
                <Button noOverflow title='Enable DOF and plastic material' onClick={() => this.stylizeddof()} style={{ width: 'auto' }}>
                    Stylized-DOF
                </Button>
            </div>
            <div className='msp-flex-row'>
                <Button noOverflow title='Applies default representation preset and sets outline and occlusion effects to default' onClick={() => this.metal()} style={{ width: 'auto' }}>
                    Metal
                </Button>
                <Button noOverflow title='Applies illustrative representation preset and Stylize it with cellShading' onClick={() => this.illustrative()} style={{ width: 'auto' }}>
                    Illustrative
                </Button>
                <Button noOverflow title='Enable DOF and plastic material' onClick={() => this.stylized()} style={{ width: 'auto' }}>
                    Stylized
                </Button>
            </div>
        </>;
    }
}