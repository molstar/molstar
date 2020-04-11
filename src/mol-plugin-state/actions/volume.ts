/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../mol-plugin/context';
import { StateAction, StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { getFileInfo } from '../../mol-util/file-info';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject } from '../objects';
import { Download } from '../transforms/data';
import { DataFormatProvider } from '../formats/provider';

export { DownloadDensity };
type DownloadDensity = typeof DownloadDensity
const DownloadDensity = StateAction.build({
    from: PluginStateObject.Root,
    display: { name: 'Download Density', description: 'Load a density from the provided source and create its default visual.' },
    params: (a, ctx: PluginContext) => {
        const { options } = ctx.dataFormats
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
})(({ params }, plugin: PluginContext) => Task.create('Download Density', async taskCtx => {
    const src = params.source;
    let downloadParams: StateTransformer.Params<Download>;
    let provider: DataFormatProvider | undefined;

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

    const data = await plugin.builders.data.download(downloadParams);

    switch (src.name) {
        case 'url':
            downloadParams = src.params;
            provider = src.params.format === 'auto' ? plugin.dataFormats.auto(getFileInfo(downloadParams.url), data.cell?.obj!) : plugin.dataFormats.get(src.params.format)
            break;
        case 'pdb-xray':
            provider = src.params.provider.server === 'pdbe'
                ? plugin.dataFormats.get('ccp4')
                : plugin.dataFormats.get('dsn6')
            break;
        case 'pdb-emd-ds':
        case 'pdb-xray-ds':
            provider = plugin.dataFormats.get('dscif')
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    if (!provider) {
        plugin.log.warn('DownloadDensity: Format provider not found.');
        return;
    }

    const volumes = await provider.parse(plugin, data);
    await provider.visuals?.(plugin, volumes);
}));