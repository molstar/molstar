/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeIsoValue } from '../../mol-model/volume';
import { PluginContext } from '../../mol-plugin/context';
import { State, StateAction, StateBuilder, StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { ColorNames } from '../../mol-util/color/names';
import { FileInfo, getFileInfo } from '../../mol-util/file-info';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { Download } from '../transforms/data';
import { VolumeRepresentation3DHelpers } from '../transforms/representation';
import { DataFormatProvider, guessCifVariant, DataFormatBuilderOptions } from './data-format';

export const Ccp4Provider: DataFormatProvider<any> = {
    label: 'CCP4/MRC/BRIX',
    description: 'CCP4/MRC/BRIX',
    stringExtensions: [],
    binaryExtensions: ['ccp4', 'mrc', 'map'],
    isApplicable: (info: FileInfo, data: Uint8Array) => {
        return info.ext === 'ccp4' || info.ext === 'mrc' || info.ext === 'map'
    },
    getDefaultBuilder: (ctx: PluginContext, data, options: DataFormatBuilderOptions, state: State) => {
        return Task.create('CCP4/MRC/BRIX default builder', async taskCtx => {
            let tree: StateBuilder.To<any> = state.build().to(data)
                .apply(StateTransforms.Data.ParseCcp4, {}, { state: { isGhost: true } })
                .apply(StateTransforms.Volume.VolumeFromCcp4)
            if (options.visuals) {
                tree = tree.apply(StateTransforms.Representation.VolumeRepresentation3D)
            }
            await state.updateTree(tree).runInContext(taskCtx)
        })
    }
}

export const Dsn6Provider: DataFormatProvider<any> = {
    label: 'DSN6/BRIX',
    description: 'DSN6/BRIX',
    stringExtensions: [],
    binaryExtensions: ['dsn6', 'brix'],
    isApplicable: (info: FileInfo, data: Uint8Array) => {
        return info.ext === 'dsn6' || info.ext === 'brix'
    },
    getDefaultBuilder: (ctx: PluginContext, data, options: DataFormatBuilderOptions, state: State) => {
        return Task.create('DSN6/BRIX default builder', async taskCtx => {
            let tree: StateBuilder.To<any> = state.build().to(data)
                .apply(StateTransforms.Data.ParseDsn6, {}, { state: { isGhost: true } })
                .apply(StateTransforms.Volume.VolumeFromDsn6)
            if (options.visuals) {
                tree = tree.apply(StateTransforms.Representation.VolumeRepresentation3D)
            }
            await state.updateTree(tree).runInContext(taskCtx)
        })
    }
}

export const DscifProvider: DataFormatProvider<any> = {
    label: 'DensityServer CIF',
    description: 'DensityServer CIF',
    stringExtensions: ['cif'],
    binaryExtensions: ['bcif'],
    isApplicable: (info: FileInfo, data: Uint8Array | string) => {
        return guessCifVariant(info, data) === 'dscif' ? true : false
    },
    getDefaultBuilder: (ctx: PluginContext, data, options: DataFormatBuilderOptions, state: State) => {
        return Task.create('DensityServer CIF default builder', async taskCtx => {
            const cifBuilder = state.build().to(data).apply(StateTransforms.Data.ParseCif)
            const cifStateObject = await state.updateTree(cifBuilder).runInContext(taskCtx)
            const b = state.build().to(cifBuilder.ref);
            const blocks = cifStateObject.data.blocks.slice(1); // zero block contains query meta-data
            let tree: StateBuilder.To<any, any>
            if (blocks.length === 1) {
                tree = b
                    .apply(StateTransforms.Volume.VolumeFromDensityServerCif, { blockHeader: blocks[0].header })
                if (options.visuals) {
                    tree = tree.apply(StateTransforms.Representation.VolumeRepresentation3D, VolumeRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'isosurface', { isoValue: VolumeIsoValue.relative(1.5), alpha: 0.3 }, 'uniform', { value: ColorNames.teal }))
                }
            } else if (blocks.length === 2) {
                tree = b
                    .apply(StateTransforms.Volume.VolumeFromDensityServerCif, { blockHeader: blocks[0].header })
                if (options.visuals) {
                    tree = tree.apply(StateTransforms.Representation.VolumeRepresentation3D, VolumeRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'isosurface', { isoValue: VolumeIsoValue.relative(1.5), alpha: 0.3 }, 'uniform', { value: ColorNames.blue }))
                }
                const vol = tree.to(cifBuilder.ref)
                    .apply(StateTransforms.Volume.VolumeFromDensityServerCif, { blockHeader: blocks[1].header })
                const posParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'isosurface', { isoValue: VolumeIsoValue.relative(3), alpha: 0.3 }, 'uniform', { value: ColorNames.green })
                tree = vol.apply(StateTransforms.Representation.VolumeRepresentation3D, posParams)
                const negParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'isosurface', { isoValue: VolumeIsoValue.relative(-3), alpha: 0.3 }, 'uniform', { value: ColorNames.red })
                tree = tree.to(vol.ref).apply(StateTransforms.Representation.VolumeRepresentation3D, negParams)
            } else {
                throw new Error('unknown number of blocks')
            }

            await state.updateTree(tree).runInContext(taskCtx);
        })
    }
}

export { DownloadDensity };
type DownloadDensity = typeof DownloadDensity
const DownloadDensity = StateAction.build({
    from: PluginStateObject.Root,
    display: { name: 'Download Density', description: 'Load a density from the provided source and create its default visual.' },
    params: (a, ctx: PluginContext) => {
        const { options } = ctx.dataFormat.registry
        return {
            source: PD.MappedStatic('pdb-xray', {
                'pdb-xray': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('1tqn', { label: 'Id' }),
                        server: PD.Select('rcsb', [['pdbe', 'PDBe'], ['rcsb', 'RCSB PDB']]),
                    }, { pivot: 'id' }),
                    type: PD.Select('2fofc', [['2fofc', '2Fo-Fc'], ['fofc', 'Fo-Fc']]),
                }, { isFlat: true }),
                'pdb-xray-ds': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('1tqn', { label: 'Id' }),
                        server: PD.Select('pdbe', [['pdbe', 'PDBe'], ['rcsb', 'RCSB PDB']]),
                    }, { pivot: 'id' }),
                    detail: PD.Numeric(3, { min: 0, max: 10, step: 1 }, { label: 'Detail' }),
                }, { isFlat: true }),
                'pdb-emd-ds': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('emd-8004', { label: 'Id' }),
                        server: PD.Select('pdbe', [['pdbe', 'PDBe'], ['rcsb', 'RCSB PDB']]),
                    }, { pivot: 'id' }),
                    detail: PD.Numeric(3, { min: 0, max: 10, step: 1 }, { label: 'Detail' }),
                }, { isFlat: true }),
                'url': PD.Group({
                    url: PD.Text(''),
                    isBinary: PD.Boolean(false),
                    format: PD.Select('auto', options),
                }, { isFlat: true })
            }, {
                options: [
                    ['pdb-xray', 'PDB X-ray maps'],
                    ['pdb-emd-ds', 'PDB EMD Density Server'],
                    ['pdb-xray-ds', 'PDB X-ray Density Server'],
                    ['url', 'URL']
                ]
            })
        }
    }
})(({ params, state }, ctx: PluginContext) => Task.create('Download Density', async taskCtx => {
    const src = params.source;
    let downloadParams: StateTransformer.Params<Download>;
    let provider: DataFormatProvider<any>

    switch (src.name) {
        case 'url':
            downloadParams = src.params;
            break;
        case 'pdb-xray':
            downloadParams = src.params.provider.server === 'pdbe' ? {
                url: src.params.type === '2fofc'
                    ? `http://www.ebi.ac.uk/pdbe/coordinates/files/${src.params.provider.id.toLowerCase()}.ccp4`
                    : `http://www.ebi.ac.uk/pdbe/coordinates/files/${src.params.provider.id.toLowerCase()}_diff.ccp4`,
                isBinary: true,
                label: `PDBe X-ray map: ${src.params.provider.id}`
            } : {
                url: src.params.type === '2fofc'
                    ? `https://edmaps.rcsb.org/maps/${src.params.provider.id.toLowerCase()}_2fofc.dsn6`
                    : `https://edmaps.rcsb.org/maps/${src.params.provider.id.toLowerCase()}_fofc.dsn6`,
                isBinary: true,
                label: `RCSB X-ray map: ${src.params.provider.id}`
            };
            break;
        case 'pdb-emd-ds':
            downloadParams = src.params.provider.server === 'pdbe' ? {
                url: `https://www.ebi.ac.uk/pdbe/densities/emd/${src.params.provider.id.toLowerCase()}/cell?detail=${src.params.detail}`,
                isBinary: true,
                label: `PDBe EMD Density Server: ${src.params.provider.id}`
            } : {
                url: `https://maps.rcsb.org/em/${src.params.provider.id.toLowerCase()}/cell?detail=${src.params.detail}`,
                isBinary: true,
                label: `RCSB PDB EMD Density Server: ${src.params.provider.id}`
            };
            break;
        case 'pdb-xray-ds':
            downloadParams = src.params.provider.server === 'pdbe' ? {
                url: `https://www.ebi.ac.uk/pdbe/densities/x-ray/${src.params.provider.id.toLowerCase()}/cell?detail=${src.params.detail}`,
                isBinary: true,
                label: `PDBe X-ray Density Server: ${src.params.provider.id}`
            } : {
                url: `https://maps.rcsb.org/x-ray/${src.params.provider.id.toLowerCase()}/cell?detail=${src.params.detail}`,
                isBinary: true,
                label: `RCSB PDB X-ray Density Server: ${src.params.provider.id}`
            };
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    const data = await ctx.builders.data.download(downloadParams);

    switch (src.name) {
        case 'url':
            downloadParams = src.params;
            provider = src.params.format === 'auto' ? ctx.dataFormat.registry.auto(getFileInfo(downloadParams.url), data.cell?.obj!) : ctx.dataFormat.registry.get(src.params.format)
            break;
        case 'pdb-xray':
            provider = src.params.provider.server === 'pdbe'
                ? ctx.dataFormat.registry.get('ccp4')
                : ctx.dataFormat.registry.get('dsn6')
            break;
        case 'pdb-emd-ds':
        case 'pdb-xray-ds':
            provider = ctx.dataFormat.registry.get('dscif')
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    await provider.getDefaultBuilder(ctx, data, { visuals: true }, state).runInContext(taskCtx)
}));