import { StateTransforms } from '../../../mol-plugin-state/transforms';
import { PluginContext } from '../../../mol-plugin/context';
import { TunnelDataProvider, TunnelsDataProvider, TunnelShapeProvider, TunnelShapeSimpleProvider } from './representation';
import { tunnelJson } from './1ymg_Path55_09';
import { Tunnel } from './props';


export const DB_URL = 'https://channelsdb2.biodata.ceitec.cz/api/channels/';
export const SUB_DB = 'pdb';
export const CHANNEL = '1ymg';

export const URL = `${DB_URL}${SUB_DB}/${CHANNEL}`;

export async function runVisualizeTunnelsFromUrl(plugin: PluginContext, url: string = URL) {
    const update = plugin.build();
    const webgl = plugin.canvas3dContext?.webgl;

    update
        .toRoot()
        .apply(StateTransforms.Data.Download, { url })
        .apply(TunnelsDataProvider, { kind: { name: 'raw', params: {} } })
        .apply(TunnelShapeProvider, {
            visual: { name: 'mesh', params: { resolution: 3 } },
            webgl
        })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);

    await update.commit();
}

export async function runVisualizeTunnelFromUrl(plugin: PluginContext, url: string = URL) {
    const update = plugin.build();
    const webgl = plugin.canvas3dContext?.webgl;

    update
        .toRoot()
        .apply(StateTransforms.Data.Download, { url })
        .apply(TunnelDataProvider, { kind: { name: 'raw', params: {} } })
        .apply(TunnelShapeProvider, {
            visual: { name: 'mesh', params: { resolution: 3 } },
            webgl
        })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);

    await update.commit();
}

export async function runVisualizeTunnel(plugin: PluginContext) {
    const tunnel: Tunnel = { data: tunnelJson.Profile, type: tunnelJson.Type, id: tunnelJson.Id };

    const update = plugin.build();
    const webgl = plugin.canvas3dContext?.webgl;

    update
        .toRoot()
        .apply(TunnelShapeSimpleProvider, {
            data: tunnel,
            webgl
        })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);

    await update.commit();
}

