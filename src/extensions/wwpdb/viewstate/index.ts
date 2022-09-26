/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '../../../mol-plugin/behavior';
import { StateObjectCell } from '../../../mol-state';
import { LoadWwPDBViewState, wwPDBViewStateObject } from './model';
import { wwPDBViewStateEditUI } from './ui';

export const wwPDBViewState = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'wwpdb-viewstate',
    category: 'misc',
    display: {
        name: 'wwPDB View State',
        description: 'wwPDB View State Handling'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        register(): void {
            this.ctx.state.data.actions.add(LoadWwPDBViewState);
            this.ctx.customStructureControls.set('wwpdb-viwestate-edit', wwPDBViewStateEditUI as any);

            this.subscribeObservable(this.ctx.state.data.events.changed, e => {
                if (e.inTransaction || this.ctx.behaviors.state.isAnimating.value) return;
                this.syncVersions();
            });
        }

        unregister() {
            this.ctx.state.data.actions.remove(LoadWwPDBViewState);
            this.ctx.customStructureControls.delete('wwpdb-viwestate-edit');
        }

        private versions: Record<string, string> = {};
        private async syncVersions() {
            const cells = this.ctx.state.data.selectQ(q => q.ofType(wwPDBViewStateObject));

            // TODO: wrap in transation?
            for (const c of cells) {
                if (this.versions[c.transform.ref] !== c.transform.version) {
                    try {
                        this.versions[c.transform.ref] = c.transform.version;
                        await this.syncCell(c);
                    } catch (err) {
                        console.error(err);
                        this.ctx.log.error(`wwPDB View State update: ${err}`);
                    }
                }
            }

        }

        private async syncCell(root: StateObjectCell<wwPDBViewStateObject>) {
            const plugin = this.ctx;
            const deleteAction = plugin.build();
            for (const c of plugin.state.data.selectQ(q => q.byValue(root).children())) {
                deleteAction.delete(c);
            }
            await deleteAction.commit();

            const state = root.obj?.data.state;
            if (!state) return;

            const data = await plugin.builders.data.download({ url: state.url, isBinary: state.format.isBinary }, { parent: root });
            const trajectory = await plugin.builders.structure.parseTrajectory(data, state.format.name as any);

            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
                representationPreset: state.presetName as any ?? 'auto',
            });
        }
    }
});