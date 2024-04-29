import { PluginBehavior } from '../../../mol-plugin/behavior';
import { DownloadTunnels } from './actions';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

export const SbNcbrTunnels = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'sb-ncbr-tunnels',
    category: 'misc',
    display: {
        name: 'SB NCBR Tunnels',
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        register(): void {
            this.ctx.state.data.actions.add(DownloadTunnels);
        }
        unregister() {
            this.ctx.state.data.actions.remove(DownloadTunnels);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(true),
    })
});