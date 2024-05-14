/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Dušan Veľký <dvelky@mail.muni.cz>
 */

import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { StateTransforms } from '../../../mol-plugin-state/transforms';
import { PluginContext } from '../../../mol-plugin/context';
import { StateAction } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { Asset } from '../../../mol-util/assets';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { assertUnreachable } from '../../../mol-util/type-helpers';
import { ChannelsDBdata, Tunnel, TunnelDB } from './data-model';
import { TunnelsServerConfig } from './props';
import { TunnelsFromRawData, SelectTunnel, TunnelShapeProvider } from './representation';

export const TunnelDownloadServer = {
    'channelsdb': PD.EmptyGroup({ label: 'ChannelsDB' })
};

type DownloadTunnels = typeof DownloadTunnels
export const DownloadTunnels = StateAction.build({
    display: { name: 'Download Tunnels', description: 'Load a tunnels from the provided source and create theirs representations' },
    from: PluginStateObject.Root,
    params: (_, plugin: PluginContext) => {
        return {
            source: PD.MappedStatic('pdb', {
                'pdb': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('1tqn', { label: 'PDB Id(s)', description: 'One or more comma/space separated PDB ids.' }),
                        server: PD.MappedStatic('channelsdb', TunnelDownloadServer),
                    }, { pivot: 'id' }),
                }, { isFlat: true, label: 'PDB' }),
                'alphafolddb': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('Q8W3K0', { label: 'UniProtKB AC(s)', description: 'One or more comma/space separated ACs.' }),
                        server: PD.MappedStatic('channelsdb', TunnelDownloadServer),
                    }, { pivot: 'id' })
                }, { isFlat: true, label: 'AlphaFold DB', description: 'Loads the predicted model if available' }),
                'url': PD.Group({
                    url: PD.Url(''),
                }, { isFlat: true, label: 'URL' })
            })
        };
    }
})(({ params, state }, plugin: PluginContext) => Task.create('Download Tunnels', async ctx => {
    plugin.behaviors.layout.leftPanelTabName.next('data');

    const src = params.source;
    let downloadParams: PD.Normalize<{ url: string | Asset.Url }>[];

    switch (src.name) {
        case 'url':
            downloadParams = [{ url: src.params.url }];
            break;
        case 'pdb':
            downloadParams = src.params.provider.server.name === 'channelsdb'
                ? [{ url: `${plugin?.config.get(TunnelsServerConfig.DefaultServerUrl)}/channels/pdb/${src.params.provider.id}` }]
                : assertUnreachable(src as never);
            break;
        case 'alphafolddb':
            downloadParams = src.params.provider.server.name === 'channelsdb'
                ? [{ url: `${plugin?.config.get(TunnelsServerConfig.DefaultServerUrl)}/channels/alphafill/${src.params.provider.id.toLowerCase()}` }]
                : assertUnreachable(src as never);
            break;
        default: assertUnreachable(src);
    }

    await state.transaction(async () => {
        const update = plugin.build();
        const webgl = plugin.canvas3dContext?.webgl;

        for (const download of downloadParams) {
            const response = await (await fetch(download.url.toString())).json();
            const tunnels: Tunnel[] = [];

            Object.entries(response.Channels as ChannelsDBdata).forEach(([key, values]) => {
                if (values.length > 0) {
                    values.forEach((item: TunnelDB) => {
                        tunnels.push({ data: item.Profile, props: { id: item.Id, type: item.Type } });
                    });
                }
            });

            update
                .toRoot()
                .apply(TunnelsFromRawData, { data: tunnels })
                .apply(SelectTunnel)
                .apply(TunnelShapeProvider, {
                    webgl,
                })
                .apply(StateTransforms.Representation.ShapeRepresentation3D);

            await update.commit();
        }
    }).runInContext(ctx);
}));
