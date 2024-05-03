import { PluginConfigItem } from '../../../mol-plugin/config';
import { PluginContext } from '../../../mol-plugin/context';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

export namespace TunnelsData {
    export const DefaultServerUrl = 'https://channelsdb2.biodata.ceitec.cz/api';
    export const DefaultServerType = 'pdb';
}

export const TunnelsDataParams = {
    serverType: PD.Select('pdb', [['pdb', 'pdb']] as const),
    serverUrl: PD.Text(TunnelsData.DefaultServerUrl)
};
export type TunnelsDataParams = typeof TunnelsDataParams

export const TunnelsServerConfig = {
    DefaultServerUrl: new PluginConfigItem('channelsdb-server', 'https://channelsdb2.biodata.ceitec.cz/api'),
    DefaultServerType: new PluginConfigItem<'pdb'>('serverType', 'pdb')
};

export function getTunnelsConfig(plugin: PluginContext): { [key in keyof typeof TunnelsServerConfig]: NonNullable<typeof TunnelsServerConfig[key]['defaultValue']> } {
    return {
        DefaultServerUrl: plugin.config.get(TunnelsServerConfig.DefaultServerUrl) ?? TunnelsServerConfig.DefaultServerUrl.defaultValue ?? TunnelsData.DefaultServerUrl,
        DefaultServerType: plugin.config.get(TunnelsServerConfig.DefaultServerType) ?? TunnelsServerConfig.DefaultServerType.defaultValue ?? TunnelsData.DefaultServerType,
    };
}

