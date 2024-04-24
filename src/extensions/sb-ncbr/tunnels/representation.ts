import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { StateTransformer } from '../../../mol-state';
import { TunnelStateObject, Tunnel, TunnelShapeParams, TunnelsStateObject } from './props';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { Task } from '../../../mol-task';
import { Color } from '../../../mol-util/color';
import { createTunnelShape, createSpheresShape } from './algorithm';

const Transform = StateTransformer.builderFactory('sb-ncbr-tunnels');

export const TunnelsDataProvider = Transform({
    name: 'tunnels-from-data',
    display: { name: 'Channels' },
    from: PluginStateObject.Data.String,
    to: TunnelStateObject,
    params: {
        kind: PD.MappedStatic('raw', {
            raw: PD.EmptyGroup(),
            trajectory: PD.Group({
                index: PD.Numeric(0)
            })
        }),
        channel_index: PD.Numeric(0),
        channel: PD.Select('CSATunnels_MOLE', PD.arrayToOptions([
            'CSATunnels_MOLE',
            'CSATunnels_Caver',
            'ReviewedChannels_MOLE',
            'ReviewedChannels_Caver',
            'CofactorTunnels_MOLE',
            'CofactorTunnels_Caver',
            'TransmembranePores_MOLE',
            'TransmembranePores_Caver',
            'ProcognateTunnels_MOLE',
            'ProcognateTunnels_Caver',
            'AlphaFillTunnels_MOLE',
            'AlphaFillTunnels_Caver'
        ]))
    } })({
    apply({ a, params, cache }) {
        const data = JSON.parse(a.data);
        (cache as any).Channels = data.Channels;
        const channel = data.Channels[params.channel][params.channel_index];
        if (!channel) {
            return new TunnelStateObject({
                tunnel: { data: [], type: 'Empty', id: '' } as Tunnel
            }, { label: 'Channel settings' });
        }
        return new TunnelStateObject({
            tunnel: { data: channel.Profile, type: channel.Type, id: channel.Id } as Tunnel
        }, { label: 'Channel settings' });
    },
    update({ a, b, newParams, cache }) {
        const channel = (cache as any).Channels[newParams.channel][newParams.channel_index];
        if (!channel) {
            b.data.tunnel.type = 'Empty';
            b.data.tunnel.id = '';
            b.data.tunnel.data = [];
            return StateTransformer.UpdateResult.Updated;
        }
        b.data.tunnel.type = channel.Type;
        b.data.tunnel.id = channel.Id;
        b.data.tunnel.data = channel.Profile;
        return StateTransformer.UpdateResult.Updated;
    }
});

export const TunnelDataProvider = Transform({
    name: 'tunnel-from-data',
    display: { name: 'Channel' },
    from: PluginStateObject.Data.String,
    to: TunnelStateObject,
    params: {
        kind: PD.MappedStatic('raw', {
            raw: PD.EmptyGroup(),
            trajectory: PD.Group({
                index: PD.Numeric(0)
            })
        })
    } })({
    apply({ a, params }) {
        const data = JSON.parse(a.data);
        return new TunnelStateObject({
            tunnel: { data: data.Profile, type: data.Type, id: data.Id } as Tunnel
        }, { label: 'Channel settings' });
    }
});

export const TunnelShapeProvider = Transform({
    name: 'tunnel-shape',
    display: { name: 'Channel' },
    from: TunnelStateObject,
    to: PluginStateObject.Shape.Provider,
    params: {
        webgl: PD.Value<WebGLContext | null>(null),
        colorTheme: PD.Color(Color(0xff0000)),
        visual: PD.MappedStatic(
            'mesh',
            {
                mesh: PD.Group({ resolution: PD.Numeric(2) }),
                spheres: PD.Group({ resolution: PD.Numeric(2) })
            }
        ),
    },
})({
    apply({ a, params }) {
        return Task.create('Tunnel Shape Representation', async ctx => {
            return new PluginStateObject.Shape.Provider({
                label: 'Channel',
                data: { tunnel: a.data.tunnel, params },
                params: Mesh.Params,
                geometryUtils: Mesh.Utils,
                getShape: (_, data) => data.params.visual.name === 'mesh'
                    ? createTunnelShape(data.tunnel, data.params.webgl, data.params.visual.params.resolution, data.params.colorTheme)
                    : createSpheresShape(data.tunnel, data.params.webgl, data.params.visual.params.resolution, data.params.colorTheme)
            }, { label: `${a.data.tunnel.type} ${a.data.tunnel.id}` });
        });
    },
});

export const TunnelShapeSimpleProvider = Transform({
    name: 'tunnel-shape',
    display: { name: 'Channel' },
    from: PluginStateObject.Root,
    to: PluginStateObject.Shape.Provider,
    params: {
        data: PD.Value<Tunnel>({ data: [], type: '', id: '' }),
        webgl: PD.Value<WebGLContext | null>(null),
        colorTheme: PD.Color(Color(0xff0000)),
        visual: PD.MappedStatic(
            'mesh',
            {
                mesh: PD.Group({ resolution: PD.Numeric(2) }),
                spheres: PD.Group({ resolution: PD.Numeric(2) })
            }
        ),
    },
})({
    apply({ params }) {
        return Task.create('Tunnel Shape Representation', async ctx => {
            return new PluginStateObject.Shape.Provider({
                label: 'Channel',
                data: { params },
                params: Mesh.Params,
                geometryUtils: Mesh.Utils,
                getShape: (_, data) => data.params.visual.name === 'mesh'
                    ? createTunnelShape(data.data, data.colorTheme, data.params.visual.params.resolution, data.webgl)
                    : createSpheresShape(data.data, data.colorTheme, data.params.visual.params.resolution, data.webgl)
            }, { label: `${params.data.type} ${params.data.id}` });
        });
    },
});