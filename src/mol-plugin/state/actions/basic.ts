/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from 'mol-plugin/context';
import { StateTree, Transformer } from 'mol-state';
import { StateAction } from 'mol-state/action';
import { StateSelection } from 'mol-state/state/selection';
import { StateTreeBuilder } from 'mol-state/tree/builder';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { Download } from '../transforms/data';
import { StructureRepresentation3DHelpers, VolumeRepresentation3DHelpers } from '../transforms/representation';
import { getFileInfo } from 'mol-util/file-info';

// TODO: "structure parser provider"

export { DownloadStructure };
type DownloadStructure = typeof DownloadStructure
const DownloadStructure = StateAction.build({
    from: PluginStateObject.Root,
    display: { name: 'Download Structure', description: 'Load a structure from the provided source and create its default Assembly and visual.' },
    params: {
        source: PD.MappedStatic('bcif-static', {
            'pdbe-updated': PD.Group({
                id: PD.Text('1cbs', { label: 'Id' }),
                supportProps: PD.Boolean(false)
            }, { isFlat: true }),
            'rcsb': PD.Group({
                id: PD.Text('1tqn', { label: 'Id' }),
                supportProps: PD.Boolean(false)
            }, { isFlat: true }),
            'bcif-static': PD.Group({
                id: PD.Text('1tqn', { label: 'Id' }),
                supportProps: PD.Boolean(false)
            }, { isFlat: true }),
            'url': PD.Group({
                url: PD.Text(''),
                isBinary: PD.Boolean(false),
                supportProps: PD.Boolean(false)
            }, { isFlat: true })
        }, {
            options: [
                ['pdbe-updated', 'PDBe Updated'],
                ['rcsb', 'RCSB'],
                ['bcif-static', 'BinaryCIF (static PDBe Updated)'],
                ['url', 'URL']
            ]
        })
    }
})(({ params, state }, ctx: PluginContext) => {
    const b = state.build();
    const src = params.source;
    let url: Transformer.Params<Download>;

    switch (src.name) {
        case 'url':
            url = src.params;
            break;
        case 'pdbe-updated':
            url = { url: `https://www.ebi.ac.uk/pdbe/static/entry/${src.params.id.toLowerCase()}_updated.cif`, isBinary: false, label: `PDBe: ${src.params.id}` };
            break;
        case 'rcsb':
            url = { url: `https://files.rcsb.org/download/${src.params.id.toUpperCase()}.cif`, isBinary: false, label: `RCSB: ${src.params.id}` };
            break;
        case 'bcif-static':
            url = { url: `https://webchem.ncbr.muni.cz/ModelServer/static/bcif/${src.params.id.toLowerCase()}`, isBinary: true, label: `BinaryCIF: ${src.params.id}` };
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    const data = b.toRoot().apply(StateTransforms.Data.Download, url);
    return state.update(createStructureTree(ctx, data, params.source.params.supportProps));
});

export const OpenStructure = StateAction.build({
    display: { name: 'Open Structure', description: 'Load a structure from file and create its default Assembly and visual' },
    from: PluginStateObject.Root,
    params: { file: PD.File({ accept: '.cif,.bcif' }) }
})(({ params, state }, ctx: PluginContext) => {
    const b = state.build();
    const data = b.toRoot().apply(StateTransforms.Data.ReadFile, { file: params.file, isBinary: /\.bcif$/i.test(params.file.name) });
    return state.update(createStructureTree(ctx, data, false));
});

function createStructureTree(ctx: PluginContext, b: StateTreeBuilder.To<PluginStateObject.Data.Binary | PluginStateObject.Data.String>, supportProps: boolean): StateTree {
    let root = b
        .apply(StateTransforms.Data.ParseCif)
        .apply(StateTransforms.Model.TrajectoryFromMmCif)
        .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 });

    if (supportProps) {
        root = root.apply(StateTransforms.Model.CustomModelProperties);
    }
    const structure = root.apply(StateTransforms.Model.StructureAssemblyFromModel);
    complexRepresentation(ctx, structure);

    return root.getTree();
}

function complexRepresentation(ctx: PluginContext, root: StateTreeBuilder.To<PluginStateObject.Molecule.Structure>) {
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-sequence' })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'cartoon'));
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-het' })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'ball-and-stick'));
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'water' })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'ball-and-stick', { alpha: 0.51 }));
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'spheres' })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'spacefill'));
}

export const CreateComplexRepresentation = StateAction.build({
    display: { name: 'Create Complex', description: 'Split the structure into Sequence/Water/Ligands/... ' },
    from: PluginStateObject.Molecule.Structure
})(({ ref, state }, ctx: PluginContext) => {
    const root = state.build().to(ref);
    complexRepresentation(ctx, root);
    return state.update(root.getTree());
});

export const UpdateTrajectory = StateAction.build({
    display: { name: 'Update Trajectory' },
    params: {
        action: PD.Select<'advance' | 'reset'>('advance', [['advance', 'Advance'], ['reset', 'Reset']]),
        by: PD.makeOptional(PD.Numeric(1, { min: -1, max: 1, step: 1 }))
    }
})(({ params, state }) => {
    const models = state.select(q => q.rootsOfType(PluginStateObject.Molecule.Model)
        .filter(c => c.transform.transformer === StateTransforms.Model.ModelFromTrajectory));

    const update = state.build();

    if (params.action === 'reset') {
        for (const m of models) {
            update.to(m.transform.ref).update(StateTransforms.Model.ModelFromTrajectory,
                () => ({ modelIndex: 0 }));
        }
    } else {
        for (const m of models) {
            const parent = StateSelection.findAncestorOfType(state.tree, state.cells, m.transform.ref, [PluginStateObject.Molecule.Trajectory]);
            if (!parent || !parent.obj) continue;
            const traj = parent.obj as PluginStateObject.Molecule.Trajectory;
            update.to(m.transform.ref).update(StateTransforms.Model.ModelFromTrajectory,
                old => {
                    let modelIndex = (old.modelIndex + params.by!) % traj.data.length;
                    if (modelIndex < 0) modelIndex += traj.data.length;
                    return { modelIndex };
                });
        }
    }

    return state.update(update);
});

//

const VolumeFormats = { 'ccp4': '', 'mrc': '', 'dsn6': '', 'brix': '' }
type VolumeFormat = keyof typeof VolumeFormats

function getVolumeData(format: VolumeFormat, b: StateTreeBuilder.To<PluginStateObject.Data.Binary | PluginStateObject.Data.String>) {
    switch (format) {
        case 'ccp4': case 'mrc':
            return b.apply(StateTransforms.Data.ParseCcp4).apply(StateTransforms.Model.VolumeFromCcp4);
        case 'dsn6': case 'brix':
            return b.apply(StateTransforms.Data.ParseDsn6).apply(StateTransforms.Model.VolumeFromDsn6);
    }
}

function createVolumeTree(format: VolumeFormat, ctx: PluginContext, b: StateTreeBuilder.To<PluginStateObject.Data.Binary | PluginStateObject.Data.String>): StateTree {

    const root = getVolumeData(format, b)
        .apply(StateTransforms.Representation.VolumeRepresentation3D,
            VolumeRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'isosurface'));

    return root.getTree();
}

function getFileFormat(format: VolumeFormat | 'auto', file: File): VolumeFormat {
    if (format === 'auto') {
        const fileFormat = getFileInfo(file).ext
        if (fileFormat in VolumeFormats) {
            return fileFormat as VolumeFormat
        } else {
            throw new Error('unsupported format')
        }
    } else {
        return format
    }
}

export const OpenVolume = StateAction.build({
    display: { name: 'Open Volume', description: 'Load a volume from file and create its default visual' },
    from: PluginStateObject.Root,
    params: {
        file: PD.File({ accept: '.ccp4,.mrc,.dsn6,.brix'}),
        format: PD.Select('auto', [
            ['auto', 'Automatic'], ['ccp4', 'CCP4'], ['mrc', 'MRC'], ['dsn6', 'DSN6'], ['brix', 'BRIX']
        ]),
    }
})(({ params, state }, ctx: PluginContext) => {
    const b = state.build();
    const data = b.toRoot().apply(StateTransforms.Data.ReadFile, { file: params.file, isBinary: true });
    const format = getFileFormat(params.format, params.file)
    return state.update(createVolumeTree(format, ctx, data));
});