/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../mol-plugin/context';
import { StateAction, StateSelection, StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PresetStructureRepresentations } from '../builder/structure/representation-preset';
import { BuiltInTrajectoryFormat, BuiltInTrajectoryFormats } from '../formats/trajectory';
import { RootStructureDefinition } from '../helpers/root-structure';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { Download } from '../transforms/data';
import { CustomModelProperties, CustomStructureProperties, TrajectoryFromModelAndCoordinates } from '../transforms/model';
import { Asset } from '../../mol-util/assets';
import { PluginConfig } from '../../mol-plugin/config';

const DownloadModelRepresentationOptions = (plugin: PluginContext) => PD.Group({
    type: RootStructureDefinition.getParams(void 0, 'auto').type,
    representation: PD.Select(PresetStructureRepresentations.auto.id,
        plugin.builders.structure.representation.getPresets().map(p => [p.id, p.display.name, p.display.group] as any),
        { description: 'Which representation preset to use.' }),
    asTrajectory: PD.Optional(PD.Boolean(false, { description: 'Load all entries into a single trajectory.' }))
}, { isExpanded: false });

export const PdbDownloadProvider = {
    'rcsb': PD.Group({
        encoding: PD.Select('bcif', [['cif', 'cif'], ['bcif', 'bcif']] as ['cif' | 'bcif', string][]),
    }, { label: 'RCSB PDB', isFlat: true }),
    'pdbe': PD.Group({
        variant: PD.Select('updated-bcif', [['updated-bcif', 'Updated (bcif)'], ['updated', 'Updated'], ['archival', 'Archival']] as ['updated' | 'updtaed-bcif' | 'archival', string][]),
    }, { label: 'PDBe', isFlat: true }),
};
export type PdbDownloadProvider = keyof typeof PdbDownloadProvider;

export { DownloadStructure };
type DownloadStructure = typeof DownloadStructure
const DownloadStructure = StateAction.build({
    from: PluginStateObject.Root,
    display: { name: 'Download Structure', description: 'Load a structure from the provided source and create its representation.' },
    params: (_, plugin: PluginContext) => {
        const options = DownloadModelRepresentationOptions(plugin);
        const defaultPdbProvider = plugin.config.get(PluginConfig.Download.DefaultPdbProvider) || 'pdbe';
        return {
            source: PD.MappedStatic('pdb', {
                'pdb': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('1tqn', { label: 'PDB Id(s)', description: 'One or more comma/space separated PDB ids.' }),
                        server: PD.MappedStatic(defaultPdbProvider, PdbDownloadProvider),
                    }, { pivot: 'id' }),
                    options
                }, { isFlat: true, label: 'PDB' }),
                'pdb-dev': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('PDBDEV_00000001', { label: 'PDBDev Id(s)', description: 'One or more comma/space separated ids.' }),
                        encoding: PD.Select('bcif', [['cif', 'cif'], ['bcif', 'bcif']] as ['cif' | 'bcif', string][]),
                    }, { pivot: 'id' }),
                    options
                }, { isFlat: true, label: 'PDBDEV' }),
                'swissmodel': PD.Group({
                    id: PD.Text('Q9Y2I8', { label: 'UniProtKB AC(s)', description: 'One or more comma/space separated ACs.' }),
                    options
                }, { isFlat: true, label: 'SWISS-MODEL', description: 'Loads the best homology model or experimental structure' }),
                'pubchem': PD.Group({
                    id: PD.Text('2244,2245', { label: 'PubChem ID', description: 'One or more comma/space separated IDs.' }),
                    options
                }, { isFlat: true, label: 'PubChem', description: 'Loads 3D conformer from PubChem.' }),
                'url': PD.Group({
                    url: PD.Url(''),
                    format: PD.Select<BuiltInTrajectoryFormat>('mmcif', PD.arrayToOptions(BuiltInTrajectoryFormats.map(f => f[0]), f => f)),
                    isBinary: PD.Boolean(false),
                    options
                }, { isFlat: true, label: 'URL' })
            })
        };
    }
})(({ params, state }, plugin: PluginContext) => Task.create('Download Structure', async ctx => {
    plugin.behaviors.layout.leftPanelTabName.next('data');

    const src = params.source;
    let downloadParams: StateTransformer.Params<Download>[];
    let asTrajectory = false, format: BuiltInTrajectoryFormat = 'mmcif';

    switch (src.name) {
        case 'url':
            downloadParams = [{ url: src.params.url, isBinary: src.params.isBinary }];
            format = src.params.format;
            break;
        case 'pdb':
            downloadParams = src.params.provider.server.name === 'pdbe'
                ? src.params.provider.server.params.variant === 'updated'
                    ? getDownloadParams(src.params.provider.id, id => `https://www.ebi.ac.uk/pdbe/static/entry/${id.toLowerCase()}_updated.cif`, id => `PDBe: ${id} (updated cif)`, false)
                    : src.params.provider.server.params.variant === 'updated-bcif'
                        ? getDownloadParams(src.params.provider.id, id => `https://www.ebi.ac.uk/pdbe/entry-files/download/${id.toLowerCase()}.bcif`, id => `PDBe: ${id} (updated cif)`, true)
                        : getDownloadParams(src.params.provider.id, id => `https://www.ebi.ac.uk/pdbe/static/entry/${id.toLowerCase()}.cif`, id => `PDBe: ${id} (cif)`, false)
                : src.params.provider.server.params.encoding === 'cif'
                    ? getDownloadParams(src.params.provider.id, id => `https://files.rcsb.org/download/${id.toUpperCase()}.cif`, id => `RCSB: ${id} (cif)`, false)
                    : getDownloadParams(src.params.provider.id, id => `https://models.rcsb.org/${id.toUpperCase()}.bcif`, id => `RCSB: ${id} (bcif)`, true);
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'pdb-dev':
            downloadParams = getDownloadParams(src.params.provider.id,
                id => {
                    const nId = id.toUpperCase().startsWith('PDBDEV_') ? id : `PDBDEV_${id.padStart(8, '0')}`;
                    return src.params.provider.encoding === 'bcif'
                        ? `https://pdb-dev.wwpdb.org/bcif/${nId.toUpperCase()}.bcif`
                        : `https://pdb-dev.wwpdb.org/cif/${nId.toUpperCase()}.cif`;
                },
                id => id.toUpperCase().startsWith('PDBDEV_') ? id : `PDBDEV_${id.padStart(8, '0')}`,
                src.params.provider.encoding === 'bcif'
            );
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'swissmodel':
            downloadParams = getDownloadParams(src.params.id, id => `https://swissmodel.expasy.org/repository/uniprot/${id.toUpperCase()}.pdb`, id => `SWISS-MODEL: ${id}`, false);
            asTrajectory = !!src.params.options.asTrajectory;
            format = 'pdb';
            break;
        case 'pubchem':
            downloadParams = getDownloadParams(src.params.id, id => `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/${id.trim()}/record/SDF/?record_type=3d`, id => `PubChem: ${id}`, false);
            asTrajectory = !!src.params.options.asTrajectory;
            format = 'mol';
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    const representationPreset: any = params.source.params.options.representation || PresetStructureRepresentations.auto.id;
    const showUnitcell = representationPreset !== PresetStructureRepresentations.empty.id;

    const structure = src.params.options.type.name === 'auto' ? void 0 : src.params.options.type;

    await state.transaction(async () => {
        if (downloadParams.length > 0 && asTrajectory) {
            const blob = await plugin.builders.data.downloadBlob({
                sources: downloadParams.map((src, i) => ({ id: '' + i, url: src.url, isBinary: src.isBinary })),
                maxConcurrency: 6
            }, { state: { isGhost: true } });
            const trajectory = await plugin.builders.structure.parseTrajectory(blob, { formats: downloadParams.map((_, i) => ({ id: '' + i, format: 'cif' as 'cif' })) });

            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
                structure,
                showUnitcell,
                representationPreset
            });
        } else {
            for (const download of downloadParams) {
                const data = await plugin.builders.data.download(download, { state: { isGhost: true } });
                const trajectory = await plugin.builders.structure.parseTrajectory(data, format);

                await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
                    structure,
                    showUnitcell,
                    representationPreset
                });
            }
        }
    }).runInContext(ctx);
}));

function getDownloadParams(src: string, url: (id: string) => string, label: (id: string) => string, isBinary: boolean): StateTransformer.Params<Download>[] {
    const ids = src.split(/[,\s]/).map(id => id.trim()).filter(id => !!id && (id.length >= 4 || /^[1-9][0-9]*$/.test(id)));
    const ret: StateTransformer.Params<Download>[] = [];
    for (const id of ids) {
        ret.push({ url: Asset.Url(url(id)), isBinary, label: label(id) });
    }
    return ret;
}

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
        return ctx.customModelProperties.getParams(a?.data);
    },
    isApplicable(a, t, ctx: PluginContext) {
        return t.transformer !== CustomModelProperties;
    }
})(({ ref, params }, ctx: PluginContext) => ctx.builders.structure.insertModelProperties(ref, params));

export const EnableStructureCustomProps = StateAction.build({
    display: { name: 'Custom Structure Properties', description: 'Enable parameters for custom properties of the structure.' },
    from: PluginStateObject.Molecule.Structure,
    params(a, ctx: PluginContext) {
        return ctx.customStructureProperties.getParams(a?.data);
    },
    isApplicable(a, t, ctx: PluginContext) {
        return t.transformer !== CustomStructureProperties;
    }
})(({ ref, params }, ctx: PluginContext) => ctx.builders.structure.insertStructureProperties(ref, params));

export const AddTrajectory = StateAction.build({
    display: { name: 'Add Trajectory', description: 'Add trajectory from existing model/topology and coordinates.' },
    from: PluginStateObject.Root,
    params(a, ctx: PluginContext) {
        const state = ctx.state.data;
        const models = [
            ...state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Model)),
            ...state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Topology)),
        ];
        const modelOptions = models.map(t => [t.transform.ref, t.obj!.label]) as [string, string][];
        const coords = state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Coordinates));
        const coordOptions = coords.map(c => [c.transform.ref, c.obj!.label]) as [string, string][];
        return {
            model: PD.Select(modelOptions.length ? modelOptions[0][0] : '', modelOptions),
            coordinates: PD.Select(coordOptions.length ? coordOptions[0][0] : '', coordOptions)
        };
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
    }).runInContext(taskCtx);
}));