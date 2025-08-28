/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PluginBehavior } from '../../mol-plugin/behavior';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { SequenceColoringProvider } from './coloring-provider';
import { SequenceColorProperty } from './prop';


/** Allows coloring residues in sequence panel */
export const SequenceColor = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'sequence-color',
    category: 'misc',
    display: {
        name: 'Sequence Color',
        description: 'Sequence Color extension, allows assigning custom residue colors to be shown in the sequence panel, based on a custom structure property',
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        register(): void {
            this.ctx.customStructureProperties.register(SequenceColorProperty.Provider, this.params.autoAttach);
            this.ctx.customSequenceColoringRegistry.register(SequenceColoringProvider);
        }
        update(p: { autoAttach: boolean }) {
            const updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            this.ctx.customStructureProperties.setDefaultAutoAttach(SequenceColorProperty.Provider.descriptor.name, this.params.autoAttach);
            return updated;
        }
        unregister() {
            this.ctx.customSequenceColoringRegistry.unregister(SequenceColoringProvider.name);
            this.ctx.customStructureProperties.unregister(SequenceColorProperty.Provider.descriptor.name);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(true),
    })
});
