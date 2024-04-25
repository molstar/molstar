/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Dušan Veľký <dvelky@mail.muni.cz>
 */

import { StateTransforms } from '../../../mol-plugin-state/transforms';
import { PluginContext } from '../../../mol-plugin/context';
import { Channels, Tunnel, TunnelDB } from './props';
import { tunnel_Path55 } from './1ymg_Path55';
import { TunnelsDataTransformer, TunnelsToTunnelTransformer, TunnelShapeProvider, TunnelDataTransformer } from './representation';


export const DB_URL = 'https://channelsdb2.biodata.ceitec.cz/api/channels/';
export const SUB_DB = 'pdb';
export const CHANNEL = '1ymg';

export const URL = `${DB_URL}${SUB_DB}/${CHANNEL}`;

export async function runVisualizeTunnels(plugin: PluginContext, url: string = URL) {
    const update = plugin.build();
    const webgl = plugin.canvas3dContext?.webgl;

    const response = await (await fetch(url)).json();

    const tunnels: Tunnel[] = [];
    Object.entries(response.Channels as Channels).forEach(([key, values]) => {
        if (values.length > 0) {
            values.forEach((item: TunnelDB) => {
                tunnels.push({ data: item.Profile, props: { id: item.Id, type: item.Type } });
            });
        }
    });

    update
        .toRoot()
        .apply(TunnelsDataTransformer, { data: tunnels })
        .apply(TunnelsToTunnelTransformer)
        .apply(TunnelShapeProvider, {
            webgl,
        })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);

    await update.commit();
}

export async function runVisualizeTunnel(plugin: PluginContext) {
    const update = plugin.build();
    const webgl = plugin.canvas3dContext?.webgl;

    const tunnel = { data: tunnel_Path55.Profile, props: { id: tunnel_Path55.Id, type: tunnel_Path55.Type } };

    update
        .toRoot()
        .apply(TunnelDataTransformer, { data: tunnel })
        .apply(TunnelShapeProvider, {
            webgl,
        })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);

    await update.commit();
}

