/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '../../behavior.js';
import { Structure, Model } from '@molstar/model/model/structure';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { StateSelection, StateObject } from '@molstar/core/state';

const isModelOrStructure = (obj: StateObject | undefined) =>
  !!obj && (PluginStateObject.Molecule.Model.is(obj) || PluginStateObject.Molecule.Structure.is(obj));

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
      const cells = this.ctx.state.data.select(
        StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure),
      );
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

    /** The max model index the last walk found: what a model or structure replaced in place (`updated`) is given. */
    private maxIndex = -1;

    private setMaxIndices() {
      const value = this.maxModelIndex;
      this.maxIndex = value;
      this.setModelMaxIndex(value);
      this.setStructureMaxIndex(value);
    }

    private setModelMaxIndex(value: number) {
      const cells = this.ctx.state.data.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Model));
      for (const c of cells) {
        const m = c.obj?.data;
        if (m) {
          if (Model.MaxIndex.get(m).value !== value) {
            Model.MaxIndex.set(m, { value }, value);
          }
        }
      }
    }

    private setStructureMaxIndex(value: number) {
      const cells = this.ctx.state.data.select(
        StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure),
      );
      for (const c of cells) {
        const s = c.obj?.data;
        if (s) {
          if (Structure.MaxIndex.get(s).value !== value) {
            Structure.MaxIndex.set(s, { value }, value);
          }
        }
      }
    }

    private setMaxIndexOf(obj: StateObject) {
      const value = this.maxIndex;
      if (value < 0) return;
      if (PluginStateObject.Molecule.Model.is(obj)) {
        if (Model.MaxIndex.get(obj.data).value !== value) Model.MaxIndex.set(obj.data, { value }, value);
      } else if (PluginStateObject.Molecule.Structure.is(obj) && obj.data.parent === undefined) {
        if (Structure.MaxIndex.get(obj.data).value !== value) Structure.MaxIndex.set(obj.data, { value }, value);
      }
    }

    private handleModel(model: Model, oldModel?: Model) {
      if (Model.Index.get(model).value === undefined) {
        const oldIndex = oldModel && Model.Index.get(oldModel).value;
        const value = oldIndex ?? this.maxModelIndex + 1;
        Model.Index.set(model, { value }, value);
      }

      if (Model.AsymIdOffset.get(model).value === undefined) {
        const oldOffset = oldModel && Model.AsymIdOffset.get(oldModel).value;
        const value = oldOffset ?? { ...this.asymIdOffset };
        Model.AsymIdOffset.set(model, { value }, value);
      }
    }

    private handleStructure(structure: Structure, oldStructure?: Structure) {
      if (structure.parent !== undefined) return;
      if (Structure.Index.get(structure).value !== undefined) return;

      const oldIndex = oldStructure && Structure.Index.get(oldStructure).value;
      const value = oldIndex ?? this.maxStructureIndex + 1;
      Structure.Index.set(structure, { value }, value);
    }

    private handle(
      ref: string,
      obj: StateObject<any, StateObject.Type<any>>,
      oldObj?: StateObject<any, StateObject.Type<any>>,
    ) {
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
      this.ctx.customModelProperties.register(Model.MaxIndex, true);
      this.ctx.customStructureProperties.register(Structure.Index, true);
      this.ctx.customStructureProperties.register(Structure.MaxIndex, true);

      this.subscribeObservable(this.ctx.state.data.events.object.created, (o) => {
        this.handle(o.ref, o.obj);
        // Only a model or a structure can move the max indices or need them set: anything else (a representation, a
        // component, a layer) skips the walks, which made an op over N structures O(N²).
        if (isModelOrStructure(o.obj)) this.setMaxIndices();
      });

      // A removed model lowers the max at once (before, the next created object of any kind did it).
      this.subscribeObservable(this.ctx.state.data.events.object.removed, (o) => {
        if (isModelOrStructure(o.obj)) this.setMaxIndices();
      });

      this.subscribeObservable(this.ctx.state.data.events.object.updated, (o) => {
        this.handle(o.ref, o.obj, o.oldObj);
        // A model or a structure made anew (a trajectory frame, other params) has the max at once, without a walk (before,
        // the next created object of any kind set it).
        if (o.action === 'recreate') this.setMaxIndexOf(o.obj);
      });
    }

    unregister() {
      this.ctx.customModelProperties.unregister(Model.AsymIdOffset.descriptor.name);
      this.ctx.customModelProperties.unregister(Model.Index.descriptor.name);
      this.ctx.customModelProperties.unregister(Model.MaxIndex.descriptor.name);
      this.ctx.customStructureProperties.unregister(Structure.Index.descriptor.name);
      this.ctx.customStructureProperties.unregister(Structure.MaxIndex.descriptor.name);
    }
  },
});
