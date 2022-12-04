/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../mol-plugin/behavior';
import { PetworldColorThemeProvider } from './color';
import { LoadPetworldModel, StructureFromPetworld } from './model';
import { PetworldPreset } from './preset';

export const PetWorld = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'petworld',
    category: 'custom-props',
    display: {
        name: 'PetWorld',
        description: 'PetWorld Model Loading and Viewing.'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        register(): void {
            this.ctx.state.data.actions.add(LoadPetworldModel);
            this.ctx.state.data.actions.add(StructureFromPetworld);
            this.ctx.builders.structure.hierarchy.registerPreset(PetworldPreset);
            this.ctx.representation.structure.themes.colorThemeRegistry.add(PetworldColorThemeProvider);
        }

        unregister() {
            this.ctx.state.data.actions.remove(LoadPetworldModel);
            this.ctx.state.data.actions.remove(StructureFromPetworld);
            this.ctx.builders.structure.hierarchy.unregisterPreset(PetworldPreset);
            this.ctx.representation.structure.themes.colorThemeRegistry.remove(PetworldColorThemeProvider);
        }
    }
});
