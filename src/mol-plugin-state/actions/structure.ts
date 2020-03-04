/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../mol-plugin/context';
import { StateAction, StateBuilder, StateSelection, StateTransformer, State } from '../../mol-state';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { Download, ParsePsf } from '../transforms/data';
import { CustomModelProperties, StructureSelectionFromExpression, CustomStructureProperties, CoordinatesFromDcd, TrajectoryFromModelAndCoordinates, TopologyFromPsf } from '../transforms/model';
import { DataFormatProvider, guessCifVariant, DataFormatBuilderOptions } from './data-format';
import { FileInfo } from '../../mol-util/file-info';
import { Task } from '../../mol-task';
import { StructureElement } from '../../mol-model/structure';
import { createDefaultStructureComplex } from '../../mol-plugin/util/structure-complex-helper';
import { RootStructureDefinition } from '../helpers/root-structure';

export const MmcifProvider: DataFormatProvider<any> = {
    label: 'mmCIF',
    description: 'mmCIF',
    stringExtensions: ['cif', 'mmcif', 'mcif'],
    binaryExtensions: ['bcif'],
    isApplicable: (info: FileInfo, data: Uint8Array | string) => {
        if (info.ext === 'mmcif' || info.ext === 'mcif') return true
        // assume cif/bcif files that are not DensityServer CIF are mmCIF
        if (info.ext === 'cif' || info.ext === 'bcif') return guessCifVariant(info, data) !== 'dscif'
        return false
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.Binary | PluginStateObject.Data.String>, options: DataFormatBuilderOptions, state: State) => {
        return Task.create('mmCIF default builder', async taskCtx => {
            const traj = createModelTree(data, 'cif');
            await state.updateTree(options.visuals ? createStructureAndVisuals(ctx, traj, false) : traj).runInContext(taskCtx)
        })
    }
}

export const PdbProvider: DataFormatProvider<any> = {
    label: 'PDB',
    description: 'PDB',
    stringExtensions: ['pdb', 'ent'],
    binaryExtensions: [],
    isApplicable: (info: FileInfo, data: string) => {
        return info.ext === 'pdb' || info.ext === 'ent'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.String>, options: DataFormatBuilderOptions, state: State) => {
        return Task.create('PDB default builder', async taskCtx => {
            const traj = createModelTree(data, 'pdb');
            await state.updateTree(options.visuals ? createStructureAndVisuals(ctx, traj, false) : traj).runInContext(taskCtx)
        })
    }
}

export const GroProvider: DataFormatProvider<any> = {
    label: 'GRO',
    description: 'GRO',
    stringExtensions: ['gro'],
    binaryExtensions: [],
    isApplicable: (info: FileInfo, data: string) => {
        return info.ext === 'gro'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.String>, options: DataFormatBuilderOptions, state: State) => {
        return Task.create('GRO default builder', async taskCtx => {
            const traj = createModelTree(data, 'gro');
            await state.updateTree(options.visuals ? createStructureAndVisuals(ctx, traj, false) : traj).runInContext(taskCtx)
        })
    }
}

export const Provider3dg: DataFormatProvider<any> = {
    label: '3DG',
    description: '3DG',
    stringExtensions: ['3dg'],
    binaryExtensions: [],
    isApplicable: (info: FileInfo, data: string) => {
        return info.ext === '3dg'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.String>, options: DataFormatBuilderOptions, state: State) => {
        return Task.create('3DG default builder', async taskCtx => {
            const traj = createModelTree(data, '3dg');
            await state.updateTree(options.visuals ? createStructureAndVisuals(ctx, traj, false) : traj).runInContext(taskCtx)
        })
    }
}

export const PsfProvider: DataFormatProvider<any> = {
    label: 'PSF',
    description: 'PSF',
    stringExtensions: ['psf'],
    binaryExtensions: [],
    isApplicable: (info: FileInfo, data: string) => {
        return info.ext === 'psf'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.String>, options: DataFormatBuilderOptions, state: State) => {
        return Task.create('PSF default builder', async taskCtx => {
            await state.updateTree(data.apply(ParsePsf, {}, { state: { isGhost: true } }).apply(TopologyFromPsf)).runInContext(taskCtx)
        })
    }
}

export const DcdProvider: DataFormatProvider<any> = {
    label: 'DCD',
    description: 'DCD',
    stringExtensions: [],
    binaryExtensions: ['dcd'],
    isApplicable: (info: FileInfo, data: string) => {
        return info.ext === 'dcd'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.Binary>, options: DataFormatBuilderOptions, state: State) => {
        return Task.create('DCD default builder', async taskCtx => {
            await state.updateTree(data.apply(CoordinatesFromDcd)).runInContext(taskCtx)
        })
    }
}

type StructureFormat = 'pdb' | 'cif' | 'gro' | '3dg'

//

const DownloadModelRepresentationOptions = PD.Group({
    type: RootStructureDefinition.getParams(void 0, 'assembly').type,
    noRepresentation: PD.Optional(PD.Boolean(false, { description: 'Omit creating default representation.' }))
}, { isExpanded: false });

const DownloadStructurePdbIdSourceOptions = PD.Group({
    supportProps: PD.Optional(PD.Boolean(false)),
    asTrajectory: PD.Optional(PD.Boolean(false, { description: 'Load all entries into a single trajectory.' }))
});

export { DownloadStructure };
type DownloadStructure = typeof DownloadStructure
const DownloadStructure = StateAction.build({
    from: PluginStateObject.Root,
    display: { name: 'Download Structure', description: 'Load a structure from the provided source and create its representation.' },
    params: {
        source: PD.MappedStatic('bcif-static', {
            'pdbe-updated': PD.Group({
                id: PD.Text('1cbs', { label: 'PDB Id(s)', description: 'One or more comma separated PDB ids.' }),
                structure: DownloadModelRepresentationOptions,
                options: DownloadStructurePdbIdSourceOptions
            }, { isFlat: true, label: 'PDBe Updated' }),
            'rcsb': PD.Group({
                id: PD.Text('1tqn', { label: 'PDB Id(s)', description: 'One or more comma separated PDB ids.' }),
                structure: DownloadModelRepresentationOptions,
                options: DownloadStructurePdbIdSourceOptions
            }, { isFlat: true, label: 'RCSB' }),
            'pdb-dev': PD.Group({
                id: PD.Text('PDBDEV_00000001', { label: 'PDBDev Id(s)', description: 'One or more comma separated ids.' }),
                structure: DownloadModelRepresentationOptions,
                options: DownloadStructurePdbIdSourceOptions
            }, { isFlat: true, label: 'PDBDEV' }),
            'bcif-static': PD.Group({
                id: PD.Text('1tqn', { label: 'PDB Id(s)', description: 'One or more comma separated PDB ids.' }),
                structure: DownloadModelRepresentationOptions,
                options: DownloadStructurePdbIdSourceOptions
            }, { isFlat: true, label: 'BinaryCIF (static PDBe Updated)' }),
            'swissmodel': PD.Group({
                id: PD.Text('Q9Y2I8', { label: 'UniProtKB AC(s)', description: 'One or more comma separated ACs.' }),
                structure: DownloadModelRepresentationOptions,
                options: DownloadStructurePdbIdSourceOptions
            }, { isFlat: true, label: 'SWISS-MODEL', description: 'Loads the best homology model or experimental structure' }),
            'url': PD.Group({
                url: PD.Text(''),
                format: PD.Select('cif', [['cif', 'CIF'], ['pdb', 'PDB']] as ['cif' | 'pdb', string][]),
                isBinary: PD.Boolean(false),
                structure: DownloadModelRepresentationOptions,
                options: PD.Group({
                    supportProps: PD.Optional(PD.Boolean(false)),
                    createRepresentation: PD.Optional(PD.Boolean(true))
                })
            }, { isFlat: true, label: 'URL' })
        })
    }
})(({ params, state }, plugin: PluginContext) => Task.create('Download Structure', async ctx => {
    plugin.behaviors.layout.leftPanelTabName.next('data');

    const src = params.source;
    let downloadParams: StateTransformer.Params<Download>[];
    let supportProps = false, asTrajectory = false, format: StructureFormat = 'cif';

    switch (src.name) {
        case 'url':
            downloadParams = [{ url: src.params.url, isBinary: src.params.isBinary }];
            supportProps = !!src.params.options.supportProps;
            format = src.params.format
            break;
        case 'pdbe-updated':
            downloadParams = getDownloadParams(src.params.id, id => `https://www.ebi.ac.uk/pdbe/static/entry/${id.toLowerCase()}_updated.cif`, id => `PDBe: ${id}`, false);
            supportProps = !!src.params.options.supportProps;
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'rcsb':
            downloadParams = getDownloadParams(src.params.id, id => `https://files.rcsb.org/download/${id.toUpperCase()}.cif`, id => `RCSB: ${id}`, false);
            supportProps = !!src.params.options.supportProps;
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'pdb-dev':
            downloadParams = getDownloadParams(src.params.id,
                id => {
                    const nId = id.toUpperCase().startsWith('PDBDEV_') ? id : `PDBDEV_${id.padStart(8, '0')}`
                    return `https://pdb-dev.wwpdb.org/static/cif/${nId.toUpperCase()}.cif`
                },
                id => id.toUpperCase().startsWith('PDBDEV_') ? id : `PDBDEV_${id.padStart(8, '0')}`,
                false
            );
            supportProps = !!src.params.options.supportProps;
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'bcif-static':
            downloadParams = getDownloadParams(src.params.id, id => `https://webchem.ncbr.muni.cz/ModelServer/static/bcif/${id.toLowerCase()}`, id => `BinaryCIF: ${id}`, true);
            supportProps = !!src.params.options.supportProps;
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'swissmodel':
            downloadParams = getDownloadParams(src.params.id, id => `https://swissmodel.expasy.org/repository/uniprot/${id.toUpperCase()}.pdb`, id => `SWISS-MODEL: ${id}`, false);
            supportProps = !!src.params.options.supportProps;
            asTrajectory = !!src.params.options.asTrajectory;
            format = 'pdb'
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    const createRepr = !params.source.params.structure.noRepresentation;

    state.transaction(async () => {
        if (downloadParams.length > 0 && asTrajectory) {
            const traj = await createSingleTrajectoryModel(plugin, state, downloadParams);
            const struct = createStructure(state.build().to(traj), supportProps, src.params.structure.type);
            await state.updateTree(struct, { revertIfAborted: true }).runInContext(ctx);
            if (createRepr) {
                await plugin.builders.structureRepresentation.apply(struct.ref, 'auto');
            }
        } else {
            for (const download of downloadParams) {
                const data = await plugin.builders.data.download(download, { state: { isGhost: true } });
                const traj = createModelTree(state.build().to(data), format);

                const struct = createStructure(traj, supportProps, src.params.structure.type);
                await state.updateTree(struct, { revertIfAborted: true }).runInContext(ctx);
                if (createRepr) {
                    await plugin.builders.structureRepresentation.apply(struct.ref, 'auto');
                }
            }
        }
    }).runInContext(ctx);
}));

function getDownloadParams(src: string, url: (id: string) => string, label: (id: string) => string, isBinary: boolean): StateTransformer.Params<Download>[] {
    const ids = src.split(',').map(id => id.trim()).filter(id => !!id && (id.length >= 4 || /^[1-9][0-9]*$/.test(id)));
    const ret: StateTransformer.Params<Download>[] = [];
    for (const id of ids) {
        ret.push({ url: url(id), isBinary, label: label(id) })
    }
    return ret;
}

async function createSingleTrajectoryModel(plugin: PluginContext, state: State, sources: StateTransformer.Params<Download>[]) {
    const data = await plugin.builders.data.downloadBlob({
        sources: sources.map((src, i) => ({ id: '' + i, url: src.url, isBinary: src.isBinary })),
        maxConcurrency: 6
    }, { state: { isGhost: true } });

    const trajectory = state.build().to(data)
        .apply(StateTransforms.Data.ParseBlob, {
            formats: sources.map((_, i) => ({ id: '' + i, format: 'cif' as 'cif' }))
        }, { state: { isGhost: true } })
        .apply(StateTransforms.Model.TrajectoryFromBlob)
        .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 });

    await plugin.runTask(state.updateTree(trajectory, { revertIfAborted: true }));
    return trajectory.selector;
}

export function createModelTree(b: StateBuilder.To<PluginStateObject.Data.Binary | PluginStateObject.Data.String>, format: StructureFormat = 'cif') {
    let parsed: StateBuilder.To<PluginStateObject.Molecule.Trajectory>
    switch (format) {
        case 'cif':
            parsed = b.apply(StateTransforms.Data.ParseCif, void 0, { state: { isGhost: true } })
                .apply(StateTransforms.Model.TrajectoryFromMmCif)
            break
        case 'pdb':
            parsed = b.apply(StateTransforms.Model.TrajectoryFromPDB);
            break
        case 'gro':
            parsed = b.apply(StateTransforms.Model.TrajectoryFromGRO);
            break
        case '3dg':
            parsed = b.apply(StateTransforms.Model.TrajectoryFrom3DG);
            break
        default:
            throw new Error('unsupported format')
    }

    return parsed.apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 });
}

function createStructure(b: StateBuilder.To<PluginStateObject.Molecule.Model>, supportProps: boolean, params?: RootStructureDefinition.Params) {
    let root = b;
    if (supportProps) {
        root = root.apply(StateTransforms.Model.CustomModelProperties);
    }
    return root.apply(StateTransforms.Model.StructureFromModel, { type: params || { name: 'assembly', params: { } } });
}

function createStructureAndVisuals(ctx: PluginContext, b: StateBuilder.To<PluginStateObject.Molecule.Model>, supportProps: boolean, params?: RootStructureDefinition.Params) {
    const structure = createStructure(b, supportProps, params);
    createDefaultStructureComplex(ctx, structure);
    return b;
}

export const Create3DRepresentationPreset = StateAction.build({
    display: { name: '3D Representation Preset', description: 'Create one of preset 3D representations.' },
    from: PluginStateObject.Molecule.Structure,
    isApplicable(a, _, plugin: PluginContext) { return plugin.builders.structureRepresentation.hasProvider(a.data); },
    params(a, plugin: PluginContext) {
        return {
            type: plugin.builders.structureRepresentation.getOptions(a.data)
        };
    }
})(({ ref, params }, plugin: PluginContext) => {
    plugin.builders.structureRepresentation.apply(ref, params.type.name, params.type.params);
});

export const Remove3DRepresentationPreset = StateAction.build({
    display: { name: 'Remove 3D Representation Preset', description: 'Remove 3D representations.' },
    from: PluginStateObject.Molecule.Structure,
    isApplicable(_, t, plugin: PluginContext) { return plugin.builders.structureRepresentation.hasManagedRepresentation(t.ref); },
    params(a, plugin: PluginContext) {
        return {
            type: plugin.builders.structureRepresentation.getOptions(a.data).select
        };
    }
})(({ ref, params }, plugin: PluginContext) => {
    // TODO: this will be completely handled by the managed and is just for testing purposes
    plugin.builders.structureRepresentation.remove(params.type, ref);
});

export const UpdateTrajectory = StateAction.build({
    display: { name: 'Update Trajectory' },
    params: {
        action: PD.Select<'advance' | 'reset'>('advance', [['advance', 'Advance'], ['reset', 'Reset']]),
        by: PD.Optional(PD.Numeric(1, { min: -1, max: 1, step: 1 }))
    }
})(({ params, state }) => {
    const models = state.selectQ(q => q.ofTransformer(StateTransforms.Model.ModelFromTrajectory));

    const update = state.build();

    if (params.action === 'reset') {
        for (const m of models) {
            update.to(m).update({ modelIndex: 0 });
        }
    } else {
        for (const m of models) {
            const parent = StateSelection.findAncestorOfType(state.tree, state.cells, m.transform.ref, [PluginStateObject.Molecule.Trajectory]);
            if (!parent || !parent.obj) continue;
            const traj = parent.obj;
            update.to(m).update(old => {
                let modelIndex = (old.modelIndex + params.by!) % traj.data.length;
                if (modelIndex < 0) modelIndex += traj.data.length;
                return { modelIndex };
            });
        }
    }

    return state.updateTree(update);
});

export const EnableModelCustomProps = StateAction.build({
    display: { name: 'Custom Model Properties', description: 'Enable parameters for custom properties of the model.' },
    from: PluginStateObject.Molecule.Model,
    params(a, ctx: PluginContext) {
        return ctx.customModelProperties.getParams(a?.data)
    },
    isApplicable(a, t, ctx: PluginContext) {
        return t.transformer !== CustomModelProperties;
    }
})(({ ref, params, state }, ctx: PluginContext) => {
    const root = state.build().to(ref).insert(CustomModelProperties, params);
    return state.updateTree(root);
});

export const EnableStructureCustomProps = StateAction.build({
    display: { name: 'Custom Structure Properties', description: 'Enable parameters for custom properties of the structure.' },
    from: PluginStateObject.Molecule.Structure,
    params(a, ctx: PluginContext) {
        return ctx.customStructureProperties.getParams(a?.data)
    },
    isApplicable(a, t, ctx: PluginContext) {
        return t.transformer !== CustomStructureProperties;
    }
})(({ ref, params, state }, ctx: PluginContext) => {
    const root = state.build().to(ref).insert(CustomStructureProperties, params);
    return state.updateTree(root);
});

export const TransformStructureConformation = StateAction.build({
    display: { name: 'Transform Conformation' },
    from: PluginStateObject.Molecule.Structure,
    params: StateTransforms.Model.TransformStructureConformation.definition.params,
})(({ ref, params, state }) => {
    const root = state.build().to(ref).insert(StateTransforms.Model.TransformStructureConformation, params as any);
    return state.updateTree(root);
});

export const StructureFromSelection = StateAction.build({
    display: { name: 'Structure from Current Selection', description: 'Create a new Structure from the current selection.' },
    from: PluginStateObject.Molecule.Structure,
    params: {
        label: PD.Text()
    }
    // isApplicable(a, t, ctx: PluginContext) {
    //     return t.transformer !== CustomModelProperties;
    // }
})(({ a, ref, params, state }, plugin: PluginContext) => {
    const sel = plugin.helpers.structureSelectionManager.get(a.data);
    if (sel.kind === 'empty-loci') return Task.constant('', void 0);

    const expression = StructureElement.Loci.toExpression(sel);
    const root = state.build().to(ref).apply(StructureSelectionFromExpression, { expression, label: params.label });
    return state.updateTree(root);
});

export const AddTrajectory = StateAction.build({
    display: { name: 'Add Trajectory', description: 'Add trajectory from existing model/topology and coordinates.' },
    from: PluginStateObject.Root,
    params(a, ctx: PluginContext) {
        const state = ctx.state.dataState
        const models = [
            ...state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Model)),
            ...state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Topology)),
        ]
        const modelOptions = models.map(t => [t.transform.ref, t.obj!.label]) as [string, string][]
        const coords = state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Coordinates))
        const coordOptions = coords.map(c => [c.transform.ref, c.obj!.label]) as [string, string][]
        return {
            model: PD.Select(modelOptions.length ? modelOptions[0][0] : '', modelOptions),
            coordinates: PD.Select(coordOptions.length ? coordOptions[0][0] : '', coordOptions)
        }
    }
})(({ ref, params, state }, ctx: PluginContext) => {
    const dependsOn = [params.model, params.coordinates];
    const root = state.build().toRoot()
        .apply(TrajectoryFromModelAndCoordinates, {
            modelRef: params.model,
            coordinatesRef: params.coordinates
        }, { dependsOn })
        .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 })
    return state.updateTree(createStructureAndVisuals(ctx, root, false));
});