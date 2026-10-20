/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '../../behavior.js';
import { Structure, Model } from '@molstar/model/model/structure';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { StateSelection, StateObject, StateTransform } from '@molstar/core/state';

export const StructureInfo = PluginBehavior.create({
  name: 'structure-info-prop',
  category: 'custom-props',
  display: { name: 'Structure Info' },
  ctor: class StructureInfoHandler extends PluginBehavior.Handler {
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

    /**
     * The max indices depend only on the root model and structure objects (`rootsOfType`): they are walked for again
     * only after one of those was created or updated, or a model or structure was removed, not on every created object
     * (components, representations: O(objects × cells) otherwise). The result is the same, at the same moment.
     */
    private maxIndexDirty = true;

    private setMaxIndex() {
      if (!this.maxIndexDirty) return;
      this.maxIndexDirty = false;
      const value = this.maxModelIndex;
      this.setModelMaxIndex(value);
      this.setStructureMaxIndex(value);
    }

    private static typeOf(obj: StateObject | undefined) {
      if (PluginStateObject.Molecule.Model.is(obj)) return PluginStateObject.Molecule.Model.type;
      if (PluginStateObject.Molecule.Structure.is(obj)) return PluginStateObject.Molecule.Structure.type;
      return undefined;
    }

    /** a created or updated object below another of its type (a component, a decorator) is not one of the roots */
    private markIfRoot(ref: string, obj: StateObject | undefined) {
      const type = StructureInfoHandler.typeOf(obj);
      if (!type || this.maxIndexDirty) return;
      const { tree, cells } = this.ctx.state.data;
      for (let t = tree.transforms.get(ref); t && t.ref !== StateTransform.RootRef; ) {
        t = tree.transforms.get(t.parent);
        if (t && cells.get(t.ref)?.obj?.type === type) return;
      }
      this.maxIndexDirty = true;
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
        this.markIfRoot(o.ref, o.obj);
        this.setMaxIndex();
      });

      this.subscribeObservable(this.ctx.state.data.events.object.updated, (o) => {
        this.handle(o.ref, o.obj, o.oldObj);
        this.markIfRoot(o.ref, o.obj);
        this.markIfRoot(o.ref, o.oldObj);
      });

      // the cell is gone from the tree: any model or structure counts
      this.subscribeObservable(this.ctx.state.data.events.object.removed, (o) => {
        if (StructureInfoHandler.typeOf(o.obj)) this.maxIndexDirty = true;
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
