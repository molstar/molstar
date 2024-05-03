/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Dušan Veľký <dvelky@mail.muni.cz>
 */

import { WebGLContext } from '../../../mol-gl/webgl/context';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

export interface Profile {
    Charge: number,
    Radius: number,
    FreeRadius: number,
    T: number,
    Distance: number,
    X: number,
    Y: number,
    Z: number
}

export interface Layerweightedproperties {
    Hydrophobicity: number,
    Hydropathy: number,
    Polarity: number,
    Mutability: number
}

export interface LayerGeometry {
    MinRadius: number,
    MinFreeRadius: number,
    StartDistance: number,
    EndDistance: number,
    LocalMinimum: boolean,
    Bottleneck: boolean,
    bottleneck: boolean
}

export interface Properties {
    Charge: number,
    NumPositives: number,
    NumNegatives: number,
    Hydrophobicity: number,
    Hydropathy: number,
    Polarity: number,
    Mutability: number
}

export interface LayersInfo {
    LayerGeometry: LayerGeometry,
    Residues: string[],
    FlowIndices: string[],
    Properties: Properties
}

export interface Layers {
    ResidueFlow: string[],
    HetResidues: any[],
    LayerWeightedProperties: Layerweightedproperties
    LayersInfo: LayersInfo[]
}

export interface TunnelDB {
    Type: string,
    Id: string,
    Cavity: string,
    Auto: boolean,
    Properties: Properties,
    Profile: Profile[],
    Layers: Layers
};

export interface ChannelsDBdata {
    'CSATunnels_MOLE': TunnelDB[],
    'CSATunnels_Caver': TunnelDB[],
    'ReviewedChannels_MOLE': TunnelDB[],
    'ReviewedChannels_Caver': TunnelDB[],
    'CofactorTunnels_MOLE': TunnelDB[],
    'CofactorTunnels_Caver': TunnelDB[],
    'TransmembranePores_MOLE': TunnelDB[],
    'TransmembranePores_Caver': TunnelDB[],
    'ProcognateTunnels_MOLE': TunnelDB[],
    'ProcognateTunnels_Caver': TunnelDB[],
    'AlphaFillTunnels_MOLE': TunnelDB[],
    'AlphaFillTunnels_Caver': TunnelDB[]
}

export interface ChannelsCache {
    Channels: ChannelsDBdata
}

export interface Tunnel {
    data: Profile[],
    props: {
        highlight_label?: string,
        type?: string,
        id?: string,
        label?: string,
        description?: string
    }
}

export const TunnelShapeParams = {
    webgl: PD.Value<WebGLContext | null>(null),
    colorTheme: PD.Color(Color(0xff0000)),
    visual: PD.MappedStatic(
        'mesh',
        {
            mesh: PD.Group({ resolution: PD.Numeric(2) }),
            spheres: PD.Group({ resolution: PD.Numeric(2) })
        }
    ),
    samplingRate: PD.Numeric(1, { min: 0.05, max: 1, step: 0.05 }),
    showRadii: PD.Boolean(false),
};

export class TunnelStateObject extends PluginStateObject.Create<{ tunnel: Tunnel }>({ name: 'Tunnel Entry', typeClass: 'Data' }) { }
export class TunnelsStateObject extends PluginStateObject.Create<{ tunnels: Tunnel[] }>({ name: 'Tunnels', typeClass: 'Data' }) { }
