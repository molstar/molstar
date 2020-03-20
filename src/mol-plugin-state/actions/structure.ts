/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../mol-plugin/context';
import { StateAction, StateBuilder, StateSelection, StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { FileInfo } from '../../mol-util/file-info';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { TrajectoryFormat } from '../builder/structure';
import { BuiltInTrajectoryFormat } from '../formats/trajectory';
import { RootStructureDefinition } from '../helpers/root-structure';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { Download, ParsePsf } from '../transforms/data';
import { CoordinatesFromDcd, CustomModelProperties, CustomStructureProperties, TopologyFromPsf, TrajectoryFromModelAndCoordinates } from '../transforms/model';
import { DataFormatProvider, guessCifVariant } from './data-format';
import { applyTrajectoryHierarchyPreset } from '../builder/structure/hierarchy-preset';
import { PresetStructureReprentations } from '../builder/structure/representation-preset';

// TODO make unitcell creation part of preset

export const MmcifProvider: DataFormatProvider<PluginStateObject.Data.String | PluginStateObject.Data.Binary> = {
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
    getDefaultBuilder: (ctx: PluginContext, data, options) => {
        return Task.create('mmCIF default builder', async () => {
            const trajectory = await ctx.builders.structure.parseTrajectory(data, 'mmcif');
            const representationPreset = options.visuals ? 'auto' : 'empty';
            await applyTrajectoryHierarchyPreset(ctx, trajectory, 'first-model', { showUnitcell: options.visuals, representationPreset });
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
    getDefaultBuilder: (ctx: PluginContext, data, options) => {
        return Task.create('PDB default builder', async () => {
            const trajectory = await ctx.builders.structure.parseTrajectory(data, 'pdb');
            const representationPreset = options.visuals ? 'auto' : 'empty';
            await applyTrajectoryHierarchyPreset(ctx, trajectory, 'first-model', { showUnitcell: options.visuals, representationPreset });
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
    getDefaultBuilder: (ctx: PluginContext, data, options) => {
        return Task.create('GRO default builder', async () => {
            const trajectory = await ctx.builders.structure.parseTrajectory(data, 'gro');
            const representationPreset = options.visuals ? 'auto' : 'empty';
            await applyTrajectoryHierarchyPreset(ctx, trajectory, 'first-model', { showUnitcell: options.visuals, representationPreset });
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
    getDefaultBuilder: (ctx: PluginContext, data, options) => {
        return Task.create('3DG default builder', async () => {
            const trajectory = await ctx.builders.structure.parseTrajectory(data, '3dg');
            const representationPreset = options.visuals ? 'auto' : 'empty';
            await applyTrajectoryHierarchyPreset(ctx, trajectory, 'first-model', { showUnitcell: options.visuals, representationPreset });
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
    getDefaultBuilder: (ctx: PluginContext, data, options, state) => {
        return Task.create('PSF default builder', async taskCtx => {
            const build = state.build().to(data).apply(ParsePsf, {}, { state: { isGhost: true } }).apply(TopologyFromPsf)
            await state.updateTree(build).runInContext(taskCtx)
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
    getDefaultBuilder: (ctx: PluginContext, data, options, state) => {
        return Task.create('DCD default builder', async taskCtx => {
            const build = state.build().to(data).apply(CoordinatesFromDcd);
            await state.updateTree(build).runInContext(taskCtx)
        })
    }
}

//

const DownloadModelRepresentationOptions = (plugin: PluginContext) => PD.Group({
    type: RootStructureDefinition.getParams(void 0, 'assembly').type,
    representation: PD.Select(PresetStructureReprentations.auto.id,
        plugin.builders.structure.representation.getPresets().map(p => [p.id, p.display.name] as any),
        { description: 'Which representation preset to use.' }),
    asTrajectory: PD.Optional(PD.Boolean(false, { description: 'Load all entries into a single trajectory.' }))
}, { isExpanded: false });

export { DownloadStructure };
type DownloadStructure = typeof DownloadStructure
const DownloadStructure = StateAction.build({
    from: PluginStateObject.Root,
    display: { name: 'Download Structure', description: 'Load a structure from the provided source and create its representation.' },
    params: (_, plugin: PluginContext) => {
        const options = DownloadModelRepresentationOptions(plugin);
        return {
            source: PD.MappedStatic('pdb', {
                'pdb': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('1tqn', { label: 'PDB Id(s)', description: 'One or more comma separated PDB ids.' }),
                        server: PD.MappedStatic('rcsb', {
                            'rcsb': PD.Group({
                                encoding: PD.Select('bcif', [['cif', 'cif'], ['bcif', 'bcif']] as ['cif' | 'bcif', string][]),
                            }, { label: 'RCSB PDB', isFlat: true }),
                            'pdbe': PD.Group({
                                variant: PD.Select('updated', [['updated', 'Updated'], ['archival', 'Archival']] as ['updated' | 'archival', string][]),
                            }, { label: 'PDBe', isFlat: true }),
                        }),
                    }, { pivot: 'id' }),
                    options
                }, { isFlat: true, label: 'PDB' }),
                'pdb-dev': PD.Group({
                    id: PD.Text('PDBDEV_00000001', { label: 'PDBDev Id(s)', description: 'One or more comma separated ids.' }),
                    options
                }, { isFlat: true, label: 'PDBDEV' }),
                'bcif-static': PD.Group({
                    id: PD.Text('1tqn', { label: 'PDB Id(s)', description: 'One or more comma separated PDB ids.' }),
                    options
                }, { isFlat: true, label: 'BinaryCIF (static PDBe Updated)' }),
                'swissmodel': PD.Group({
                    id: PD.Text('Q9Y2I8', { label: 'UniProtKB AC(s)', description: 'One or more comma separated ACs.' }),
                    options
                }, { isFlat: true, label: 'SWISS-MODEL', description: 'Loads the best homology model or experimental structure' }),
                'url': PD.Group({
                    url: PD.Text(''),
                    format: PD.Select('mmcif', [['mmcif', 'CIF'], ['pdb', 'PDB']] as ['mmcif' | 'pdb', string][]),
                    isBinary: PD.Boolean(false),
                    options
                }, { isFlat: true, label: 'URL' })
            })
        }
    }
})(({ params, state }, plugin: PluginContext) => Task.create('Download Structure', async ctx => {
    plugin.behaviors.layout.leftPanelTabName.next('data');

    const src = params.source;
    let downloadParams: StateTransformer.Params<Download>[];
    let asTrajectory = false, format: BuiltInTrajectoryFormat = 'mmcif';

    switch (src.name) {
        case 'url':
            downloadParams = [{ url: src.params.url, isBinary: src.params.isBinary }];
            format = src.params.format
            break;
        case 'pdb':
            downloadParams = src.params.provider.server.name === 'pdbe'
                ? src.params.provider.server.params.variant === 'updated'
                    ? getDownloadParams(src.params.provider.id, id => `https://www.ebi.ac.uk/pdbe/static/entry/${id.toLowerCase()}_updated.cif`, id => `PDBe: ${id} (updated cif)`, false)
                    : getDownloadParams(src.params.provider.id, id => `https://www.ebi.ac.uk/pdbe/static/entry/${id.toLowerCase()}.cif`, id => `PDBe: ${id} (cif)`, false)
                : src.params.provider.server.params.encoding === 'cif'
                    ? getDownloadParams(src.params.provider.id, id => `https://files.rcsb.org/download/${id.toUpperCase()}.cif`, id => `RCSB: ${id} (cif)`, false)
                    : getDownloadParams(src.params.provider.id, id => `https://models.rcsb.org/${id.toUpperCase()}.bcif`, id => `RCSB: ${id} (bcif)`, true);
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'pdb-dev':
            downloadParams = getDownloadParams(src.params.id,
                id => {
                    const nId = id.toUpperCase().startsWith('PDBDEV_') ? id : `PDBDEV_${id.padStart(8, '0')}`
                    return `https://pdb-dev.wwpdb.org/cif/${nId.toUpperCase()}.cif`
                },
                id => id.toUpperCase().startsWith('PDBDEV_') ? id : `PDBDEV_${id.padStart(8, '0')}`,
                false
            );
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'bcif-static':
            downloadParams = getDownloadParams(src.params.id, id => `https://webchem.ncbr.muni.cz/ModelServer/static/bcif/${id.toLowerCase()}`, id => `BinaryCIF: ${id}`, true);
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'swissmodel':
            downloadParams = getDownloadParams(src.params.id, id => `https://swissmodel.expasy.org/repository/uniprot/${id.toUpperCase()}.pdb`, id => `SWISS-MODEL: ${id}`, false);
            asTrajectory = !!src.params.options.asTrajectory;
            format = 'pdb'
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    const representationPreset: any = params.source.params.options.representation || PresetStructureReprentations.auto.id;
    const showUnitcell = representationPreset !== PresetStructureReprentations.empty.id;

    await state.transaction(async () => {
        if (downloadParams.length > 0 && asTrajectory) {
            const blob = await plugin.builders.data.downloadBlob({
                sources: downloadParams.map((src, i) => ({ id: '' + i, url: src.url, isBinary: src.isBinary })),
                maxConcurrency: 6
            }, { state: { isGhost: true } });
            const trajectory = await plugin.builders.structure.parseTrajectory(blob, { formats: downloadParams.map((_, i) => ({ id: '' + i, format: 'cif' as 'cif' })) });

            await applyTrajectoryHierarchyPreset(plugin, trajectory, 'first-model', {
                structure: src.params.options.type,
                showUnitcell,
                representationPreset
            });
        } else {
            for (const download of downloadParams) {
                const data = await plugin.builders.data.download(download, { state: { isGhost: true } });
                const trajectory = await plugin.builders.structure.parseTrajectory(data, format);

                await applyTrajectoryHierarchyPreset(plugin, trajectory, 'first-model', {
                    structure: src.params.options.type,
                    showUnitcell,
                    representationPreset
                });
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

export function createModelTree(b: StateBuilder.To<PluginStateObject.Data.Binary | PluginStateObject.Data.String>, format: TrajectoryFormat = 'cif') {
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

export const Create3DRepresentationPreset = StateAction.build({
    display: { name: '3D Representation Preset', description: 'Create one of preset 3D representations.' },
    from: PluginStateObject.Molecule.Structure,
    isApplicable(a, _, plugin: PluginContext) { return plugin.builders.structure.representation.hasPreset(a); },
    params(a, plugin: PluginContext) {
        return {
            type: plugin.builders.structure.representation.getPresetsWithOptions(a)
        };
    }
})(({ ref, params }, plugin: PluginContext) => {
    return plugin.builders.structure.representation.applyPreset(ref, params.type.name, params.type.params);
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
})(({ ref, params }, ctx: PluginContext) => ctx.builders.structure.insertModelProperties(ref, params));

export const EnableStructureCustomProps = StateAction.build({
    display: { name: 'Custom Structure Properties', description: 'Enable parameters for custom properties of the structure.' },
    from: PluginStateObject.Molecule.Structure,
    params(a, ctx: PluginContext) {
        return ctx.customStructureProperties.getParams(a?.data)
    },
    isApplicable(a, t, ctx: PluginContext) {
        return t.transformer !== CustomStructureProperties;
    }
})(({ ref, params }, ctx: PluginContext) => ctx.builders.structure.insertStructureProperties(ref, params));

export const TransformStructureConformation = StateAction.build({
    display: { name: 'Transform Conformation' },
    from: PluginStateObject.Molecule.Structure,
    params: StateTransforms.Model.TransformStructureConformation.definition.params,
})(({ ref, params, state }) => {
    const root = state.build().to(ref).insert(StateTransforms.Model.TransformStructureConformation, params as any);
    return state.updateTree(root);
});

export const AddTrajectory = StateAction.build({
    display: { name: 'Add Trajectory', description: 'Add trajectory from existing model/topology and coordinates.' },
    from: PluginStateObject.Root,
    params(a, ctx: PluginContext) {
        const state = ctx.state.data
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
})(({ params, state }, ctx: PluginContext) => Task.create('Add Trajectory', taskCtx => {
    return state.transaction(async () => {
        const dependsOn = [params.model, params.coordinates];
        const model = state.build().toRoot()
            .apply(TrajectoryFromModelAndCoordinates, {
                modelRef: params.model,
                coordinatesRef: params.coordinates
            }, { dependsOn })
            .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 });

        await state.updateTree(model).runInContext(taskCtx);
        const structure = await ctx.builders.structure.createStructure(model.selector);
        await ctx.builders.structure.representation.applyPreset(structure, 'auto');
    }).runInContext(taskCtx)
}));