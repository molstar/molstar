/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russ Taylor <russ@reliasolve.com>
 */

/** Based on the ../anvil extension. */

import { Vec3, Mat3 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { KinemageDataProvider, KinemageData } from './prop';
import { StateTransformer, StateBuilder } from '../../mol-state';
import { Task } from '../../mol-task';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginDragAndDropHandler } from '../../mol-plugin-state/manager/drag-and-drop';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginContext } from '../../mol-plugin/context';
import { DefaultQueryRuntimeTable } from '../../mol-script/runtime/query/compiler';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { shapePointsFromKin, shapeLinesFromKin, shapeMeshFromKin, shapeSpheresFromKin } from './kin';
import { Kinemage } from './reader/schema';
import { DataFormatProvider } from '../../mol-plugin-state/formats/provider';
import { Camera } from '../../mol-canvas3d/camera';
import { PluginCommands } from '../../mol-plugin/commands';
import { getPluginBoundingSphere } from '../../mol-plugin-state/manager/focus-camera/focus-object';
import { KinemageControls } from './ui';
import { StateObjectSelector } from '../../mol-state';

const Tag = KinemageData.Tag;

const Transform = StateTransformer.builderFactory('sb-kinemage');

/**
 * State object to hold parsed Kinemage data
 */
export class KinemageObject extends PluginStateObject.Create<KinemageData>({ name: 'Kinemage', typeClass: 'Object' }) { }

/**
 * Apply a saved snapshot object (from a view state node) to the plugin camera.
 * Use PluginCommands.Camera.SetSnapshot so transitions and canvas props are handled properly.
 */
export async function applyViewSnapshot(plugin: PluginContext, snapshot: Partial<Camera.Snapshot>) {
  if (!snapshot) return;
  // If the snapshot provides a target, adjust the canvas `sceneRadiusFactor` so the scene isn't clipped
  // when we switch camera.
  if (snapshot.target) {
    try {
      const boundingSphere = getPluginBoundingSphere(plugin);
      if (boundingSphere && boundingSphere.radius > 0) {
        const offset = Vec3.distance(snapshot.target as Vec3, boundingSphere.center);
        const sceneRadiusFactor = (boundingSphere.radius + offset) / boundingSphere.radius;
        plugin.canvas3d?.setProps({ sceneRadiusFactor });
      }
    } catch (e) {
      // fallback: ignore errors and continue to set the camera snapshot
      console.warn('Failed to adjust sceneRadiusFactor for view snapshot', e);
    }
  }
  await PluginCommands.Camera.SetSnapshot(plugin, { snapshot });
}

/**
 * Transform to parse Kinemage data from string/data input
 */
export const ParseKinemage = Transform({
  name: 'sb-kinemage-parse',
  display: { name: 'Parse Kinemage' },
  from: [PluginStateObject.Data.String],
  to: KinemageObject,
  params: {
    label: PD.Optional(PD.Text('', { description: 'Label for the Kinemage data' }))
  }
})({
  apply({ a, params }) {
    return Task.create('Parse Kinemage', async ctx => {
      const input = a.data;
      let data: KinemageData;

      if (typeof input === 'string') {
        // Parse from string content
        const file = new File([input], 'input.kin', { type: 'text/plain' });
        data = await KinemageData.open(file);
      } else {
        throw new Error('Unsupported input type for ParseKinemage');
      }

      // Precompute camera snapshots for all views in all kinemages
      for (const kinData of data.kinemages) {
        (kinData as any).viewSnapshots = (kinData as any).viewSnapshots || Object.create(null);
        for (const [viewKey, viewObj] of Object.entries(kinData.viewDict)) {
          const center = Vec3.create(0, 0, 0);
          if (viewObj.center) {
            Vec3.set(center, viewObj.center[0], viewObj.center[1], viewObj.center[2]);
          }

          const orientation: Mat3 = Mat3.identity();
          if (viewObj.matrix) {
            Mat3.fromArray(orientation, viewObj.matrix, 0);
            Mat3.transpose(orientation, orientation);
          }

          const zAxis = Vec3.create(0, 0, 1);
          Vec3.transformMat3(zAxis, zAxis, orientation);

          const yAxis = Vec3.create(0, 1, 0);
          Vec3.transformMat3(yAxis, yAxis, orientation);

          let distance = 100;
          if (viewObj.span) {
            distance = viewObj.span;
          }
          Vec3.scale(zAxis, zAxis, distance);
          const position = Vec3.create(0, 0, 100);
          Vec3.add(position, center, zAxis);

          let radius = 100;
          if (viewObj.zslab) {
            const scale = viewObj.zslab / 200;
            radius = 0.5 * distance * scale;
          }

          const snap: Camera.Snapshot = {
            mode: 'orthographic',
            fov: Math.PI / 4,
            position,
            up: yAxis,
            target: center,
            radius,
            radiusMax: 1e4,
            fog: 0,
            clipFar: true,
            minNear: 1,
            minFar: 1
          };

          (kinData as any).viewSnapshots[viewKey] = snap;
        }
      }

      const label = params.label || data.kinemages[0]?.caption || 'Kinemage';
      return new KinemageObject(data, { label, description: `Kinemage with ${data.kinemages.length} view(s)` });
    });
  }
});

/**
 * Transform to select a specific kinemage from parsed data
 */
export const SelectKinemage = Transform({
  name: 'sb-kinemage-select',
  display: { name: 'Select Kinemage' },
  from: KinemageObject,
  to: PluginStateObject.Format.Json,
  params: (a) => {
    const kinemages = a?.data?.kinemages || [];
    const options = kinemages.map((k: Kinemage, i: number) => [i, k.pdbfile || k.caption || `Kinemage ${i}`] as const);
    return {
      index: PD.Select(0, options, { description: 'Which kinemage to use' })
    };
  }
})({
  apply({ a, params }) {
    return Task.create('Select Kinemage', async ctx => {
      const kinData = a.data.kinemages[params.index];
      if (!kinData) {
        throw new Error(`No kinemage found at index ${params.index}`);
      }

      const label = kinData.pdbfile || kinData.caption || `Kinemage ${params.index}`;

      // Store the kinemage data in a Format.Json node so downstream transforms can access it
      return new PluginStateObject.Format.Json(
        { kinData },
        { label, description: kinData.text || '' }
      );
    });
  }
});

export const KinemageShapePointsProvider = Transform({
  name: 'sb-kinemage-shape-points-provider',
  display: { name: 'Kinemage Shape Points Provider' },
  from: PluginStateObject.Format.Json,
  to: PluginStateObject.Shape.Provider,
  params: {}
})({
  apply({ a }) {
    return Task.create('Kinemage Points Shape Provider', async ctx => {
      const kinData = (a.data as any).kinData as Kinemage;
      if (!kinData) {
        throw new Error('No kinData found in parent Format.Json node');
      }

      const provider = await shapePointsFromKin(kinData, { transforms: undefined }, 'Dots').runInContext(ctx);
      return new PluginStateObject.Shape.Provider(provider as any, {
        label: kinData.pdbfile || kinData.caption || 'Kinemage Points',
        description: kinData.text || ''
      });
    });
  }
});

export const KinemageShapeLinesProvider = Transform({
  name: 'sb-kinemage-shape-lines-provider',
  display: { name: 'Kinemage Shape Lines Provider' },
  from: PluginStateObject.Format.Json,
  to: PluginStateObject.Shape.Provider,
  params: {}
})({
  apply({ a }) {
    return Task.create('Kinemage Lines Shape Provider', async ctx => {
      const kinData = (a.data as any).kinData as Kinemage;
      if (!kinData) {
        throw new Error('No kinData found in parent Format.Json node');
      }

      const provider = await shapeLinesFromKin(kinData).runInContext(ctx);
      return new PluginStateObject.Shape.Provider(provider as any, {
        label: kinData.pdbfile || kinData.caption || 'Kinemage Lines',
        description: kinData.text || ''
      });
    });
  }
});

export const KinemageShapeMeshProvider = Transform({
  name: 'sb-kinemage-shape-mesh-provider',
  display: { name: 'Kinemage Shape Mesh Provider' },
  from: PluginStateObject.Format.Json,
  to: PluginStateObject.Shape.Provider,
  params: {}
})({
  apply({ a }) {
    return Task.create('Kinemage Mesh Shape Provider', async ctx => {
      const kinData = (a.data as any).kinData as Kinemage;
      if (!kinData) {
        throw new Error('No kinData found in parent Format.Json node');
      }

      const provider = await shapeMeshFromKin(kinData).runInContext(ctx);
      return new PluginStateObject.Shape.Provider(provider as any, {
        label: kinData.pdbfile || kinData.caption || 'Kinemage Meshes',
        description: kinData.text || ''
      });
    });
  }
});

export const KinemageShapeSpheresProvider = Transform({
  name: 'sb-kinemage-shape-spheres-provider',
  display: { name: 'Kinemage Shape Spheres Provider' },
  from: PluginStateObject.Format.Json,
  to: PluginStateObject.Shape.Provider,
  params: {}
})({
  apply({ a }) {
    return Task.create('Kinemage Spheres Shape Provider', async ctx => {
      const kinData = (a.data as any).kinData as Kinemage;
      if (!kinData) {
        throw new Error('No kinData found in parent Format.Json node');
      }

      const provider = await shapeSpheresFromKin(kinData).runInContext(ctx);
      return new PluginStateObject.Shape.Provider(provider as any, {
        label: kinData.pdbfile || kinData.caption || 'Kinemage Spheres',
        description: kinData.text || ''
      });
    });
  }
});

export const KinemageExtension = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'kinemage-data-prop',
    category: 'custom-props',
    display: {
        name: 'Kinemage data',
        description: 'Data loaded from Kinemage.'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private provider = KinemageDataProvider;

        register(): void {
            DefaultQueryRuntimeTable.addCustomProp(this.provider.descriptor);

            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);

            // Register right-panel controls for Kinemage (show in the right-hand inspector)
            this.ctx.customStructureControls.set(Tag.Representation, KinemageControls as any);
            // Some app hosts expose a global customControls registry; register there too so the card is visible
            // even when no structure is loaded. Use `any` guards to avoid type errors if customControls isn't present.
            if ((this.ctx as any).customControls && typeof (this.ctx as any).customControls.set === 'function') {
              (this.ctx as any).customControls.set('kinemage', KinemageControls as any);
            }

            this.ctx.managers.dragAndDrop.addHandler(KinemageDragAndDropHandler.name, KinemageDragAndDropHandler.handle);

            // Register .kin file handler so opening/dropping .kin is supported via the data formats system
            this.ctx.dataFormats.add('KIN', KINFormatProvider);
        }

        update(p: { autoAttach: boolean }) {
            const updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        unregister() {
            DefaultQueryRuntimeTable.removeCustomProp(this.provider.descriptor);

            this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);

            this.ctx.genericRepresentationControls.delete(Tag.Representation);

            this.ctx.managers.dragAndDrop.removeHandler(KinemageDragAndDropHandler.name);

            // Unregister the .kin data format provider
            this.ctx.dataFormats.remove('KIN');

            // Remove right-panel controls
            try { this.ctx.customStructureControls.delete(Tag.Representation); } catch {}
            if ((this.ctx as any).customControls && typeof (this.ctx as any).customControls.delete === 'function') {
              try { (this.ctx as any).customControls.delete('kinemage'); } catch {}
            }
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false)
    })
});

/** Registerable method for handling dragged-and-dropped files */
interface DragAndDropHandler {
  name: string,
  handle: PluginDragAndDropHandler,
}

/** Helper function to create all shapes for a kinemage via proper transform chain */
async function createShapesForKinemage(plugin: PluginContext, update: StateBuilder.Root, kinDataSelector: StateObjectSelector<PluginStateObject.Format.Json>) {
  const kinDataCell = plugin.state.data.cells.get(kinDataSelector.ref);
  if (!kinDataCell?.obj?.data) return;

  const kinData = (kinDataCell.obj.data as any).kinData as Kinemage;
  if (!kinData) return;

  // Generate all shape types that have data, each as child of the selected kinemage
  if (kinData.dotLists.length > 0) {
    await update
      .to(kinDataSelector)
      .apply(KinemageShapePointsProvider, {}, { state: { isGhost: true } })
      .apply(StateTransforms.Representation.ShapeRepresentation3D);
  }
  if (kinData.vectorLists.length > 0) {
    await update
      .to(kinDataSelector)
      .apply(KinemageShapeLinesProvider, {}, { state: { isGhost: true } })
      .apply(StateTransforms.Representation.ShapeRepresentation3D);
  }
  if (kinData.ribbonLists.length > 0) {
    await update
      .to(kinDataSelector)
      .apply(KinemageShapeMeshProvider, {}, { state: { isGhost: true } })
      .apply(StateTransforms.Representation.ShapeRepresentation3D, { doubleSided: true });
  }
  if (kinData.ballLists.length > 0) {
    await update
      .to(kinDataSelector)
      .apply(KinemageShapeSpheresProvider, {}, { state: { isGhost: true } })
      .apply(StateTransforms.Representation.ShapeRepresentation3D);
  }
}

/** Helper function to rebuild shapes for a kinemage (remove and recreate) */
export async function rebuildShapesForKinemage(plugin: PluginContext, kinDataSelector: StateObjectSelector<PluginStateObject.Format.Json>) {
  // Store current camera snapshot
  const curSnap = (plugin.canvas3d && (plugin.canvas3d as any).camera && (plugin.canvas3d as any).camera.getSnapshot)
    ? (plugin.canvas3d as any).camera.getSnapshot()
    : undefined;

  const update = plugin.state.data.build();

  // Remove all children of this kinemage node (shapes/representations)
  const children = plugin.state.data.tree.children.get(kinDataSelector.ref);
  if (children) {
    for (const childRef of children.values()) {
      update.delete(childRef);
    }
  }

  // Recreate shapes
  await createShapesForKinemage(plugin, update, kinDataSelector);
  await update.commit();

  // Restore camera
  if (curSnap) {
    try {
      await applyViewSnapshot(plugin, curSnap);
    } catch (e) {
      console.warn('Failed to restore camera snapshot after recreating shapes', e);
    }
  }
}

/** Centralized helper to apply kinemage content into plugin state */
async function applyKinemageToState(plugin: PluginContext, data: string, label?: string) {
  const update = plugin.state.data.build();

  // Create String data node
  const dataNode = update
    .toRoot()
    .apply(StateTransforms.Data.RawData, { data, label: label || 'Kinemage File' });

  // Parse into KinemageObject
  const parsedNode = dataNode
    .apply(ParseKinemage, { label });

  // Select first kinemage (default)
  const selectedNode = parsedNode
    .apply(SelectKinemage, { index: 0 });

  await update.commit();

  // Now create shapes from the selected kinemage
  const shapeUpdate = plugin.state.data.build();
  await createShapesForKinemage(plugin, shapeUpdate, selectedNode.selector);
  await shapeUpdate.commit();

  // Wait for bounding sphere and focus camera
  async function waitForNonEmptyBoundingSphere(plugin: PluginContext, timeoutMs = 2000, pollMs = 50) {
    const start = Date.now();
    while (Date.now() - start < timeoutMs) {
      try {
        const bs = getPluginBoundingSphere(plugin);
        if (bs && bs.radius > 0) return bs;
      } catch { /* ignore */ }
      await new Promise<void>(r => setTimeout(r, pollMs));
    }
    return null;
  }

  try {
    const bs = await waitForNonEmptyBoundingSphere(plugin);
    if (bs && bs.radius > 0 && plugin.canvas3d) {
      await PluginCommands.Camera.Focus(plugin, { center: bs.center, radius: bs.radius, durationMs: 250 });
      plugin.canvas3d?.commit();
    }
  } catch (e) {
    console.warn('Failed to apply initial kinemage view snapshot', e);
  }

  return selectedNode.selector;
}

/** Programmatic loader: load a single File (a .kin) into the plugin state.
 * Returns the ref to the selected kinemage node.
 */
export async function loadKinemageFile(plugin: PluginContext, file: File): Promise<StateObjectSelector<PluginStateObject.Format.Json> | undefined> {
  const content = await file.text();
  return await applyKinemageToState(plugin, content, file.name);
}

/** DragAndDropHandler handler for `.kin` files */
const KinemageDragAndDropHandler: DragAndDropHandler = {
  name: 'kin',
  async handle(files: File[], plugin: PluginContext): Promise<boolean> {
    let applied = false;
    for (const file of files) {
      if (file.name.toLowerCase().endsWith('.kin')) {
        const ref = await loadKinemageFile(plugin, file);
        applied = applied || !!ref;
      }
    }
    return applied;
  },
};

const KINFormatProvider: DataFormatProvider<{}, any, any> = DataFormatProvider({
  label: 'KIN',
  description: 'Kinemage',
  category: 'Miscellaneous',
  stringExtensions: ['kin', 'KIN'],
  parse: async (plugin, data) => {
    try {
      // data is already a StateObjectRef to the raw data in the tree
      // Build the transform chain from it
      const builder = plugin.state.data.build()
        .to(data)
        .apply(ParseKinemage, {});

      const selectedKin = builder
        .apply(SelectKinemage, { index: 0 });

      await builder.commit();

      // Return the selector for the selected kinemage so visuals can use it
      return { selectedKin: selectedKin.selector };
    } catch (e) {
      console.error('Failed to parse KIN file', e);
      throw e;
    }
  },
  visuals: async (plugin, data) => {
    if (!data?.selectedKin) {
      console.warn('[Kinemage] visuals: no selectedKin ref provided');
      return;
    }

    // Create shapes from the selected kinemage
    const shapeBuilder = plugin.state.data.build();
    await createShapesForKinemage(plugin, shapeBuilder, data.selectedKin);
    await shapeBuilder.commit();

    // Wait for bounding sphere and focus camera
    async function waitForNonEmptyBoundingSphere(plugin: PluginContext, timeoutMs = 2000, pollMs = 50) {
      const start = Date.now();
      while (Date.now() - start < timeoutMs) {
        try {
          const bs = getPluginBoundingSphere(plugin);
          if (bs && bs.radius > 0) return bs;
        } catch { /* ignore */ }
        await new Promise<void>(r => setTimeout(r, pollMs));
      }
      return null;
    }

    try {
      const bs = await waitForNonEmptyBoundingSphere(plugin);
      if (bs && bs.radius > 0 && plugin.canvas3d) {
        await PluginCommands.Camera.Focus(plugin, { center: bs.center, radius: bs.radius, durationMs: 250 });
        plugin.canvas3d?.commit();
      }
    } catch (e) {
      console.warn('Failed to focus camera on kinemage', e);
    }

    return undefined;
  }
});
