/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '../../behavior.js';
import { Structure, Model } from '@molstar/model/model/structure';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { StateSelection, StateObject, StateTransform } from '@molstar/core/state';

const isModelOrStructure = (obj: StateObject | undefined) =>
  !!obj && (PluginStateObject.Molecule.Model.is(obj) || PluginStateObject.Molecule.Structure.is(obj));

export const StructureInfo = PluginBehavior.create({
  name: 'structure-info-prop',
  category: 'custom-props',
  display: { name: 'Structure Info' },
  ctor: class StructureInfoHandler extends PluginBehavior.Handler {
    /**
     * The root model and structure cells (what `rootsOfType` selects: no ancestor of the same type), kept from the
     * object events instead of walking the state for each created model and structure (O(structures²) otherwise), plus
     * the cells new in the current update whose events have not come yet (`pending`: the walk would see them too). A
     * cell's object is checked when it is read. The values below are maxima over them, so their order does not matter.
     */
    private rootModels = new Set<string>();
    private rootStructures = new Set<string>();
    private pending = new Set<string>();

    private isRoot(ref: string, type: StateObject.Type) {
      const { tree, cells } = this.ctx.state.data;
      for (let t = tree.transforms.get(ref); t && t.ref !== StateTransform.RootRef; ) {
        t = tree.transforms.get(t.parent);
        if (t && cells.get(t.ref)?.obj?.type === type) return false;
      }
      return true;
    }

    private *rootObjects<T>(roots: Set<string>, type: StateObject.Type): Generator<T> {
      const { cells } = this.ctx.state.data;
      for (const ref of roots) {
        const obj = cells.get(ref)?.obj;
        if (obj?.type === type) yield obj.data as T;
      }
      for (const ref of this.pending) {
        const obj = cells.get(ref)?.obj;
        if (obj?.type === type && !roots.has(ref) && this.isRoot(ref, type)) yield obj.data as T;
      }
    }

    private models() {
      return this.rootObjects<Model>(this.rootModels, PluginStateObject.Molecule.Model.type);
    }

    private structures() {
      return this.rootObjects<Structure>(this.rootStructures, PluginStateObject.Molecule.Structure.type);
    }

    private get maxModelIndex() {
      let maxIndex = -1;
      for (const m of this.models()) {
        const index = Model.Index.get(m).value;
        if (index !== undefined && index > maxIndex) maxIndex = index;
      }
      return maxIndex;
    }

    private get maxStructureIndex() {
      let maxIndex = -1;
      for (const s of this.structures()) {
        const index = Structure.Index.get(s).value;
        if (index !== undefined && index > maxIndex) maxIndex = index;
      }
      return maxIndex;
    }

    private get asymIdOffset() {
      let auth = 0;
      let label = 0;
      for (const m of this.models()) {
        const count = Model.AsymIdCount.get(m);
        const offset = Model.AsymIdOffset.get(m).value;
        if (count !== undefined && offset !== undefined) {
          auth = Math.max(auth, offset.auth + count.auth);
          label = Math.max(label, offset.label + count.label);
        }
      }
      return { auth, label };
    }

    /** The max model index last set: what a model or structure replaced in place (`updated`) is given. */
    private maxIndex = -1;

    /** After a model or structure was created or removed: over the roots, not a walk of the state */
    private setMaxIndices() {
      const value = this.maxModelIndex;
      this.maxIndex = value;
      this.setModelMaxIndex(value);
      this.setStructureMaxIndex(value);
    }

    private rootsOf(obj: StateObject | undefined) {
      if (PluginStateObject.Molecule.Model.is(obj)) return this.rootModels;
      if (PluginStateObject.Molecule.Structure.is(obj)) return this.rootStructures;
      return undefined;
    }

    /** a created or updated object below another of its type (a component, a decorator) is not one of the roots */
    private trackIfRoot(ref: string, obj: StateObject | undefined) {
      const roots = this.rootsOf(obj);
      if (!roots || !this.isRoot(ref, obj!.type)) return;
      roots.add(ref);
    }

    /** whatever the object was (an update in place to a Null object carries no old one) */
    private untrack(ref: string) {
      this.rootModels.delete(ref);
      this.rootStructures.delete(ref);
    }

    private setModelMaxIndex(value: number) {
      for (const m of this.models()) {
        if (Model.MaxIndex.get(m).value !== value) {
          Model.MaxIndex.set(m, { value }, value);
        }
      }
    }

    private setStructureMaxIndex(value: number) {
      for (const s of this.structures()) {
        if (Structure.MaxIndex.get(s).value !== value) {
          Structure.MaxIndex.set(s, { value }, value);
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

      // the roots already in the state (the behavior may be registered after loading)
      const { data } = this.ctx.state;
      for (const c of data.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Model))) {
        this.rootModels.add(c.transform.ref);
      }
      for (const c of data.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure))) {
        this.rootStructures.add(c.transform.ref);
      }

      this.subscribeObservable(data.events.cell.created, (o) => this.pending.add(o.ref));
      // after the object events of an update (a cell that failed has none)
      this.subscribeObservable(data.events.changed, () => this.pending.clear());
      this.subscribeObservable(data.events.cell.removed, (o) => {
        this.pending.delete(o.ref);
        this.untrack(o.ref);
      });

      this.subscribeObservable(data.events.object.created, (o) => {
        this.pending.delete(o.ref);
        this.trackIfRoot(o.ref, o.obj);
        this.handle(o.ref, o.obj);
        // Only a model or a structure can move the max indices or need them set: anything else (a representation, a
        // layer) skips them, which made an op over N structures O(N²).
        if (isModelOrStructure(o.obj)) this.setMaxIndices();
      });

      // A removed model lowers the max at once (before, the next created object of any kind did it).
      this.subscribeObservable(data.events.object.removed, (o) => {
        this.untrack(o.ref);
        if (isModelOrStructure(o.obj)) this.setMaxIndices();
      });

      this.subscribeObservable(data.events.object.updated, (o) => {
        this.untrack(o.ref);
        this.trackIfRoot(o.ref, o.obj);
        // An in-place update gives the old data, not the old object (a structure remapped onto a trajectory's next
        // model): its index carries over all the same.
        const old = o.oldObj ?? (o.oldData !== undefined ? ({ data: o.oldData } as StateObject) : undefined);
        this.handle(o.ref, o.obj, old);
        // A model or a structure made anew or remapped in place (a trajectory frame, other params) has the max at once,
        // without a walk (before, the next created object of any kind set it).
        this.setMaxIndexOf(o.obj);
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
