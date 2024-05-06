/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Dušan Veľký <dvelky@mail.muni.cz>
 */

import { PluginBehavior } from '../../../mol-plugin/behavior';
import { DownloadTunnels } from './actions';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PresetStructureRepresentations, StructureRepresentationPresetProvider } from '../../../mol-plugin-state/builder/structure/representation-preset';
import { Model, Structure } from '../../../mol-model/structure';
import { PluginContext } from '../../../mol-plugin/context';
import { StateObjectRef } from '../../../mol-state';
import { getTunnelsConfig, TunnelsDataParams } from './props';
import { StateTransforms } from '../../../mol-plugin-state/transforms';
import { Tunnel, ChannelsDBdata, TunnelDB } from './data-model';
import { TunnelShapeProvider, TunnelFromRawData } from './representation';
import { ColorGenerator } from '../../meshes/mesh-utils';

export const SbNcbrTunnels = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'sb-ncbr-tunnels',
    category: 'misc',
    display: {
        name: 'SB NCBR Tunnels',
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        register(): void {
            this.ctx.state.data.actions.add(DownloadTunnels);
            this.ctx.builders.structure.representation.registerPreset(TunnelsPreset);
        }
        unregister() {
            this.ctx.state.data.actions.remove(DownloadTunnels);
            this.ctx.builders.structure.representation.unregisterPreset(TunnelsPreset);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(true),
    })
});


export function isApplicable(structure?: Structure): boolean {
    return (
        !!structure && structure.models.length === 1 &&
        Model.hasPdbId(structure.models[0])
    );
}

export const TunnelsPreset = StructureRepresentationPresetProvider({
    id: 'sb-ncbr-preset-structure-tunnels',
    display: {
        name: 'Tunnels', group: 'Annotation',
        description: 'Shows Tunnels from ChannelsDB contained in the structure.'
    },
    isApplicable(a) {
        return isApplicable(a.data);
    },
    params: (a, plugin) => {
        return {
            ...StructureRepresentationPresetProvider.CommonParams,
            ...getConfiguredDefaultParams(plugin)
        };
    },
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        const structure = structureCell?.obj?.data;
        if (!structureCell || !structure) return {};

        const update = plugin.build();
        const webgl = plugin.canvas3dContext?.webgl;
        const response = await (await fetch(`${params.serverUrl}/channels/${params.serverType}/${structure.model.entryId.toLowerCase()}`)).json();
        const tunnels: Tunnel[] = [];

        Object.entries(response.Channels as ChannelsDBdata).forEach(([key, values]) => {
            if (values.length > 0) {
                values.forEach((item: TunnelDB) => {
                    tunnels.push({ data: item.Profile, props: { id: item.Id, type: item.Type } });
                });
            }
        });

        await tunnels.forEach(async (tunnel) => {
            await update
                .toRoot()
                .apply(TunnelFromRawData, { data: tunnel })
                .apply(TunnelShapeProvider, {
                    webgl,
                    colorTheme: ColorGenerator.next().value,
                })
                .apply(StateTransforms.Representation.ShapeRepresentation3D);
            await update.commit();
        });

        const preset = await PresetStructureRepresentations.auto.apply(ref, { ...params }, plugin);

        return { components: preset.components, representations: { ...preset.representations, } };
    }
});

function getConfiguredDefaultParams(plugin: PluginContext) {
    const config = getTunnelsConfig(plugin);
    const params = PD.clone(TunnelsDataParams);
    PD.setDefaultValues(params, { serverType: config.DefaultServerType, serverUrl: config.DefaultServerUrl });
    return params;
}
