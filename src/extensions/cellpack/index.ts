/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../mol-plugin/behavior';
import { LoadCellPackModel } from './model';
import { CellPackGenerateColorThemeProvider } from './color/generate';
import { CellPackProvidedColorThemeProvider } from './color/provided';
import { CellPackColorThemeProvider } from './color/basic';
import { CellpackUniformColorThemeProvider } from './color/uniform';

export const CellPack = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'cellpack',
    category: 'custom-props',
    display: {
        name: 'CellPack',
        description: 'CellPack Model Loading and Viewing.'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        register(): void {
            this.ctx.state.data.actions.add(LoadCellPackModel);
            this.ctx.representation.structure.themes.colorThemeRegistry.add(CellPackGenerateColorThemeProvider);
            this.ctx.representation.structure.themes.colorThemeRegistry.add(CellPackProvidedColorThemeProvider);
            this.ctx.representation.structure.themes.colorThemeRegistry.add(CellPackColorThemeProvider);
            this.ctx.representation.structure.themes.colorThemeRegistry.add(CellpackUniformColorThemeProvider);
        }

        unregister() {
            this.ctx.state.data.actions.remove(LoadCellPackModel);
            this.ctx.representation.structure.themes.colorThemeRegistry.remove(CellPackGenerateColorThemeProvider);
            this.ctx.representation.structure.themes.colorThemeRegistry.remove(CellPackProvidedColorThemeProvider);
            this.ctx.representation.structure.themes.colorThemeRegistry.remove(CellPackColorThemeProvider);
            this.ctx.representation.structure.themes.colorThemeRegistry.remove(CellpackUniformColorThemeProvider);
        }
    }
});