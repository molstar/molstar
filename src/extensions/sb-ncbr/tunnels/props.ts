import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

export interface Profile {
    Charge: number,
    Radius: number,
    FreeRadius: number,
    T: number,
    Distance: number,
    X: number,
    Y: number,
    Z: number,
}

export interface Channels{
    'CSATunnels_MOLE': [],
    'CSATunnels_Caver': [],
    'ReviewedChannels_MOLE': [],
    'ReviewedChannels_Caver': [],
    'CofactorTunnels_MOLE': [],
    'CofactorTunnels_Caver': [],
    'TransmembranePores_MOLE': [],
    'TransmembranePores_Caver': [],
    'ProcognateTunnels_MOLE': [],
    'ProcognateTunnels_Caver': [],
    'AlphaFillTunnels_MOLE': [],
    'AlphaFillTunnels_Caver': []
}

export interface ChannelsCache {
    Channels: Channels,
}

export interface Tunnel {
    data: Profile[],
    type: string,
    id: number | string,
}

export const TunnelParams = {
    data: PD.Value<any[]>([], { isHidden: true }),
    visual: PD.MappedStatic(
        'mesh',
        {
            mesh: PD.Group({
                resolution: PD.Numeric(2),
            }),
            sphere: PD.Group({
                resolution: PD.Numeric(2)
            })
        }
    )
};

export class TunnelStateObject extends PluginStateObject.Create<{ tunnel: Tunnel }>({ name: 'Tunnel Entry', typeClass: 'Data' }) { }