/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../../behavior';
import { Structure, Model } from '../../../../mol-model/structure';
import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { StateSelection, StateObject } from '../../../../mol-state';

export const StructureInfo = PluginBehavior.create({
    name: 'structure-info-prop',
    category: 'custom-props',
    display: { name: 'Structure Info' },
    ctor: class extends PluginBehavior.Handler {
        private get maxModelIndex() {
            let maxIndex = -1;
            const cells = this.ctx.state.data.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Model));
            for (const c of cells) {
                const index = c.obj?.data && Model.Index.get(c.obj?.data).value;
                if (index !== undefined && index > maxIndex) maxIndex = index;
            }
            return maxIndex;
        }

        private get maxStructureIndex() {
            let maxIndex = -1;
            const cells = this.ctx.state.data.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure));
            for (const c of cells) {
                const index = c.obj?.data && Structure.Index.get(c.obj?.data).value;
                if (index !== undefined && index > maxIndex) maxIndex = index;
            }
            return maxIndex;
        }

        private get asymIdOffset() {
            let auth = 0;
            let label = 0;
            const cells = this.ctx.state.data.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Model));
            for (const c of cells) {
                const m = c.obj?.data;
                if (m) {
                    const count = Model.AsymIdCount.get(m);
                    const offset = Model.AsymIdOffset.get(m).value;
                    if (count !== undefined && offset !== undefined) {
                        auth = Math.max(auth, offset.auth + count.auth);
                        label = Math.max(label, offset.label + count.label);
                    }
                }
            }
            return { auth, label };
        }

        private handleModel(model: Model, oldModel?: Model) {
            if (Model.Index.get(model).value === undefined) {
                const oldIndex = oldModel && Model.Index.get(oldModel).value;
                const value = oldIndex ?? (this.maxModelIndex + 1);
                Model.Index.set(model, { value });
            }

            if (Model.AsymIdOffset.get(model).value === undefined) {
                const oldOffset = oldModel && Model.AsymIdOffset.get(oldModel).value;
                const value = oldOffset ?? { ...this.asymIdOffset };
                Model.AsymIdOffset.set(model, { value });
            }
        }

        private handleStructure(structure: Structure, oldStructure?: Structure) {
            if (structure.parent !== undefined) return;
            if (Structure.Index.get(structure).value !== undefined) return;

            const oldIndex = oldStructure && Structure.Index.get(oldStructure).value;
            const value = oldIndex ?? (this.maxStructureIndex + 1);
            Structure.Index.set(structure, { value });
        }

        private handle(ref: string, obj: StateObject<any, StateObject.Type<any>>, oldObj?: StateObject<any, StateObject.Type<any>>) {
            if (PluginStateObject.Molecule.Structure.is(obj)) {
                const transform = this.ctx.state.data.tree.transforms.get(ref);
                if (!transform.transformer.definition.isDecorator && obj.data.parent === undefined) {
                    this.handleStructure(obj.data, oldObj?.data);
                }
            } else if (PluginStateObject.Molecule.Model.is(obj)) {
                const transform = this.ctx.state.data.tree.transforms.get(ref);
                if (!transform.transformer.definition.isDecorator) {
                    this.handleModel(obj.data, oldObj?.data);
                }
            }
        }

        register(): void {
            this.ctx.customModelProperties.register(Model.AsymIdOffset, true);
            this.ctx.customModelProperties.register(Model.Index, true);
            this.ctx.customStructureProperties.register(Structure.Index, true);

            this.subscribeObservable(this.ctx.state.data.events.object.created, o => {
                this.handle(o.ref, o.obj);
            });

            this.subscribeObservable(this.ctx.state.data.events.object.updated, o => {
                this.handle(o.ref, o.obj, o.oldObj);
            });
        }

        unregister() {
            this.ctx.customModelProperties.unregister(Model.AsymIdOffset.descriptor.name);
            this.ctx.customModelProperties.unregister(Model.Index.descriptor.name);
            this.ctx.customStructureProperties.unregister(Structure.Index.descriptor.name);
        }
    }
});