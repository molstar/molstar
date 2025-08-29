/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { BehaviorSubject, Unsubscribable } from 'rxjs';
import { ColorTypeLocation } from '../../mol-geo/geometry/color-data';
import { AccessibleSurfaceAreaColorThemeProvider } from '../../mol-model-props/computed/themes/accessible-surface-area';
import { CustomStructureProperties } from '../../mol-plugin-state/transforms/model';
import { PluginUIContext } from '../../mol-plugin-ui/context';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { ColorTheme } from '../../mol-theme/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { CustomSequenceColorThemeProvider } from './coloring-provider';
import { SequenceColorProperty } from './prop';


type ColorThemeProvider = ColorTheme.Provider<any, string, ColorTypeLocation> | undefined;

/** Allows coloring residues in sequence panel */
export const SequenceColor = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'sequence-color',
    category: 'misc',
    display: {
        name: 'Sequence Color',
        description: 'Sequence Color extension, allows assigning custom residue colors to be shown in the sequence panel, based on a custom structure property',
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        sub?: Unsubscribable;

        register(): void {
            this.ctx.customStructureProperties.register(SequenceColorProperty.Provider, this.params.autoAttach);
            if (this.ctx instanceof PluginUIContext) {
                console.log('register with PluginUIContext')
                const theme: BehaviorSubject<ColorThemeProvider> = this.ctx.customUIState.experimentalSequenceColorTheme ??= new BehaviorSubject<ColorThemeProvider>(undefined);
                this.sub = this.ctx.state.events.cell.stateUpdated.subscribe(s => {
                    if (s.cell.transform.transformer === CustomStructureProperties) {
                        theme.next(CustomSequenceColorThemeProvider);
                        // theme.next(AccessibleSurfaceAreaColorThemeProvider);
                    }
                });
            }
        }
        update(p: { autoAttach: boolean }) {
            const updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            this.ctx.customStructureProperties.setDefaultAutoAttach(SequenceColorProperty.Provider.descriptor.name, this.params.autoAttach);
            return updated;
        }
        unregister() {
            this.ctx.customStructureProperties.unregister(SequenceColorProperty.Provider.descriptor.name);
            this.sub?.unsubscribe();
            this.sub = undefined;
            if (this.ctx instanceof PluginUIContext) {
                console.log('unregister with PluginUIContext')
                const theme: BehaviorSubject<ColorThemeProvider> | undefined = this.ctx.customUIState.experimentalSequenceColorTheme;
                theme?.next(undefined);
            }
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(true),
    })
});
