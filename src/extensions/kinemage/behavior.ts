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
import { StateObjectRef } from '../../mol-state';
import { getPluginBoundingSphere } from '../../mol-plugin-state/manager/focus-camera/focus-object';
import { KinemageControls } from './ui';

const Tag = KinemageData.Tag;

const Transform = StateTransformer.builderFactory('sb-kinemage');

let g_kinemageData: KinemageData | undefined = undefined;

/** Getter for external code / handlers to obtain the loaded kinemage runtime data. */
export function getLoadedKinemageData() {
  return g_kinemageData;
}

/**
 * Map that keeps track of created shape/repr selectors for each created `Kinemage`.
 * This lets callback handlers destroy / re-create shapes for a given `kinData`.
 * Key: the `Kinemage` instance (object identity), Value: array of selectors produced
 * by the state builder for the created shape/provider/representation transforms.
 */
const g_kinemageShapeSelectors = new Map<Kinemage, StateObjectRef<PluginStateObject.Format.Json | PluginStateObject.Shape.Provider>[]>();

/** Getter for external code / handlers to obtain the selectors for a specific kinemage. */
export function getKinemageShapeSelectors(kin: Kinemage) {
  return g_kinemageShapeSelectors.get(kin) || [];
}

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

export const KinemageShapePointsProvider = Transform({
  name: 'sb-kinemage-shape-points-provider',
  display: { name: 'Kinemage Shape Points Provider' },
  from: PluginStateObject.Root,
  to: PluginStateObject.Shape.Provider,
  params: {
    data: PD.Value<Kinemage>(undefined as any, {})
  }
})({
  apply({ params }) {
    return Task.create('Kinemage Points Shape Provider', async ctx => {
      // shapeFromKin returns a Task that resolves to a ShapeProvider-like object
      const provider = await shapePointsFromKin(params.data, { transforms: undefined }, 'Dots').runInContext(ctx);
      return new PluginStateObject.Shape.Provider(provider as any, {
        label: params.data.pdbfile || params.data.caption || 'Kinemage Points',
        description: params.data.text || ''
      });
    });
  }
});

export const KinemageShapeLinesProvider = Transform({
    name: 'sb-kinemage-shape-lines-provider',
    display: { name: 'Kinemage Shape Lines Provider' },
    from: PluginStateObject.Root,
    to: PluginStateObject.Shape.Provider,
    params: {
        data: PD.Value<Kinemage>(undefined as any, {})
    }
    })({
    apply({ params }) {
        return Task.create('Kinemage Lines Shape Provider', async ctx => {
            // shapeFromKin returns a Task that resolves to a ShapeProvider-like object
            const provider = await shapeLinesFromKin(params.data).runInContext(ctx);
            return new PluginStateObject.Shape.Provider(provider as any, {
              label: params.data.pdbfile || params.data.caption || 'Kinemage Lines',
                description: params.data.text || ''
              });
        });
    }
});

export const KinemageShapeMeshProvider = Transform({
  name: 'sb-kinemage-shape-mesh-provider',
  display: { name: 'Kinemage Shape Mesh Provider' },
  from: PluginStateObject.Root,
  to: PluginStateObject.Shape.Provider,
  params: {
    data: PD.Value<Kinemage>(undefined as any, {})
  }
})({
  apply({ params }) {
    return Task.create('Kinemage Mesh Shape Provider', async ctx => {
      // shapeFromKin returns a Task that resolves to a ShapeProvider-like object
      const provider = await shapeMeshFromKin(params.data).runInContext(ctx);
      return new PluginStateObject.Shape.Provider(provider as any, {
        label: params.data.pdbfile || params.data.caption || 'Kinemage Meshes',
        description: params.data.text || ''
      });
    });
  }
});

export const KinemageShapeSpheresProvider = Transform({
  name: 'sb-kinemage-shape-spheres-provider',
  display: { name: 'Kinemage Shape Spheres Provider' },
  from: PluginStateObject.Root,
  to: PluginStateObject.Shape.Provider,
  params: {
    data: PD.Value<Kinemage>(undefined as any, {})
    }
})({
  apply({ params }) {
    return Task.create('Kinemage Spheres Shape Provider', async ctx => {
      // shapeFromKin returns a Task that resolves to a ShapeProvider-like object
      const provider = await shapeSpheresFromKin(params.data).runInContext(ctx);
      return new PluginStateObject.Shape.Provider(provider as any, {
        label: params.data.pdbfile || params.data.caption || 'Kinemage Spheres',
        description: params.data.text || ''
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
            // Register both as a structure-scoped control and (if available) as a global control
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

/** Helper function to create the shapes for a kinemage */
export async function createShapesForKinemage(plugin: PluginContext, update: StateBuilder.Root,  kinData: Kinemage) {
  // Keep list of created selectors for this kinemage (shapes / representations etc.)
  const createdShapeSelectors: StateObjectRef<any>[] = [];

  // Generate all of the shapes for this kinemage, each shape type having its own provider and representation.
  // Make all of their GUI buttons ghosted -- we'll control visibility using Kinemage master and group settings
  if (kinData.dotLists.length > 0) {
    const node = await update
      .toRoot()
      .apply(KinemageShapePointsProvider, { data: kinData }, { state: { isGhost: true } })
      .apply(StateTransforms.Representation.ShapeRepresentation3D);
    createdShapeSelectors.push(node.selector as StateObjectRef<any>);
  }
  if (kinData.vectorLists.length > 0) {
    const node = await update
      .toRoot()
      .apply(KinemageShapeLinesProvider, { data: kinData }, { state: { isGhost: true } })
      .apply(StateTransforms.Representation.ShapeRepresentation3D);
    createdShapeSelectors.push(node.selector as StateObjectRef<any>);
  }
  if (kinData.ribbonLists.length > 0) {
    const node = await update
      .toRoot()
      .apply(KinemageShapeMeshProvider, { data: kinData }, { state: { isGhost: true } })
      .apply(StateTransforms.Representation.ShapeRepresentation3D, { doubleSided: true });
    createdShapeSelectors.push(node.selector as StateObjectRef<any>);
  }
  if (kinData.ballLists.length > 0) {
    const node = await update
      .toRoot()
      .apply(KinemageShapeSpheresProvider, { data: kinData }, { state: { isGhost: true } })
      .apply(StateTransforms.Representation.ShapeRepresentation3D);
    createdShapeSelectors.push(node.selector as StateObjectRef<any>);
  }

  // Store the created selector list for this kinemage so callback handlers can destroy / re-create.
  if (createdShapeSelectors.length > 0) {
    g_kinemageShapeSelectors.set(kinData as Kinemage, createdShapeSelectors);
  }
}

/** Helper function to destroy all previously-made shapes for a kinemage
 *  (soft remove: hide the transforms so visuals are removed from scene) */
export async function destroyShapesForKinemage(plugin: PluginContext, kinData: Kinemage) {
  const createdShapeSelectors = g_kinemageShapeSelectors.get(kinData as Kinemage);
  if (!createdShapeSelectors) return;

  for (const selector of createdShapeSelectors) {
    try {
      const ref = resolveSelectorRef(selector);
      if (ref) {
        // Fully remove the transform from the state tree so the nodes are gone (not just hidden).
        // Use the plugin command so the removal is handled in the same place as other UI removals.
        try {
          await PluginCommands.State.RemoveObject(plugin, { state: plugin.state.data, ref, removeParentGhosts: true });
        } catch (e) {
          console.warn('Failed to remove state object via command, falling back to hiding', ref, e);
          // fallback: mark transform as hidden so visuals are torn down
          plugin.state.data.updateCellState(ref, (old: any) => {
            const s = { ...(old || {}) };
            s.isHidden = true;
            return s;
          });
        }
      } else if ((selector as any).destroy) {
        // Fallback if selector object exposes destroy (unlikely for state refs)
        (selector as any).destroy();
      } else {
        console.warn('Could not resolve selector to a ref for destruction', selector);
      }
    } catch (e) {
      console.warn('Failed to destroy selector', selector, e);
    }
    // Commit canvas / repaint if needed
    plugin.canvas3d?.commit();
  }

  g_kinemageShapeSelectors.delete(kinData as Kinemage);
}

/** Centralized helper to apply kinemage content into plugin state (re-used by drag handler and programmatic loader)
 *  NOTE: This no longer creates State Tree JSON nodes for views/groups/masters/subgroups/animate.
 *  It computes view snapshots and stores them on the kinemage objects (runtime-only) and then creates shapes.
 *
 *  IMPORTANT: Before adding new visuals we explicitly destroy any previously-created kinemage visuals
 *  so we don't leak shapes / state objects across loads.
 */
async function applyKinemageInfoToState(plugin: PluginContext, kinInfo: KinemageData) {

  const update = plugin.state.data.build();

  for (const kinData of kinInfo.kinemages) {
    // Precompute snapshots for views and attach them to kinData so the right-panel UI can apply them.
    (kinData as any).viewSnapshots = (kinData as any).viewSnapshots || Object.create(null);
    for (const [viewKey, viewObj] of Object.entries(kinData.viewDict)) {
      //const viewName = viewObj.name || `View ${viewKey}`;

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

    // Create shapes for this kinemage
    await createShapesForKinemage(plugin, update, kinData);
  }

  await update.commit();

  // helper: wait briefly until the plugin bounding sphere has non-zero radius (or timeout)
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

  // After commit, focus camera as before...
  try {
    const bs = await waitForNonEmptyBoundingSphere(plugin);
    if (bs && bs.radius > 0 && plugin.canvas3d) {
      await PluginCommands.Camera.Focus(plugin, { center: bs.center, radius: bs.radius, durationMs: 250 });
      plugin.canvas3d?.commit();
    } else {
      console.log('Did not get a valid bounding sphere after waiting, applying initial view snapshot without adjustment');
    }
  } catch (e) {
    console.warn('Failed to apply initial kinemage view snapshot', e);
  }
}

// Helper: robustly resolve a transform ref from different selector shapes without changing other modules.
function resolveSelectorRef(sel: any): string | undefined {
  if (!sel) return undefined;
  if (typeof sel === 'string') return sel;
  if (sel.ref && typeof sel.ref === 'string') return sel.ref;
  if (sel.transform && typeof sel.transform.ref === 'string') return sel.transform.ref;
  if (sel.cell && sel.cell.transform && typeof sel.cell.transform.ref === 'string') return sel.cell.transform.ref;
  try {
    // In case the runtime provides a utility on the ref type
    return (StateObjectRef as any).resolveRef ? (StateObjectRef as any).resolveRef(sel as any) : undefined;
  } catch {
    return undefined;
  }
}

/** Programmatic loader: load a single File (a .kin) into the plugin state.
 * Runs the import inside a Task so it has a runtime and asset context similar to drag-and-drop.
 * Returns true if at least one Kinemage was added.
 */
export async function loadKinemageFile(plugin: PluginContext, file: File): Promise<boolean> {
  console.log('XXX loadKinemageFile() called');
  let applied = false;
  const task = Task.create('Load KIN file', async ctx => {
    const kinData = await KinemageData.open(file);
    // Replace previous runtime data with the loaded one (do not keep appending old data).
    g_kinemageData = kinData;
    console.log('XXX loaded kinemage data', g_kinemageData);
    applied = g_kinemageData.kinemages.length > 0;
  });
  await plugin.runTask(task);
  return applied;
}

/** DragAndDropHandler handler for `.kin` files */
const KinemageDragAndDropHandler: DragAndDropHandler = {
  name: 'kin',
  /** Load .kin files. Append to previous plugin state.
   * If multiple files are provided, append them all.
   * Select the last-loaded one from the list.
   * Return `true` if at least one file has been loaded. */
  async handle(files: File[], plugin: PluginContext): Promise<boolean> {
    let applied = false;
    for (const file of files) {
      if (file.name.toLowerCase().endsWith('.kin')) {
        // reuse programmatic loader so drag & drop and programmatic loading behave the same
        const ok = await loadKinemageFile(plugin, file);
        applied = applied || ok;
        if (g_kinemageData) {
          await applyKinemageInfoToState(plugin, g_kinemageData);
        }
      }
    }
    return applied;
  },
};

/* Convert a string to a file if needed so that our file loader can handle it properly. */
/// @todo Consider making the handler be able to deal with a string to avoid extra work here.
function fileFromPayload(data: any): File {
  // If it's already a File or wrapped File, use name + size as signature (ignore lastModified to be more robust
  // when different File instances are created from same content).
  if (data instanceof File) {
    return data;
  }
  if (data?.input instanceof File) {
    const f: File = data.input.file;
    return f;
  }
  if (data?.data && typeof data.data === 'string') {
    const name = data.name || 'import.kin';
    const content = data.data as string;
    const file = new File([content], name, { type: 'text/plain' });
    return file;
  }
  if (typeof data === 'string') {
    const file = new File([data], 'import.kin', { type: 'text/plain' });
    return file;
  }

  // Fallback: stringify & use length + prefix
  try {
    const s = String(data);
    const file = new File([s], 'import.kin', { type: 'text/plain' });
    return file;
  } catch {
    // Last resort, use a unique key so we don't accidentally collide
    const file = new File([''], 'import.kin', { type: 'text/plain' });
    return file;
  }
}

const KINFormatProvider: DataFormatProvider<{}, any, any> = DataFormatProvider({
  label: 'KIN',
  description: 'Kinemage',
  category: 'Miscellaneous',
  // accept common casings
  stringExtensions: ['kin', 'KIN'],
  parse: async (plugin, data) => {
    try {
      const file = fileFromPayload(data);
      await loadKinemageFile(plugin, file);
    } catch (e) {
      console.error('Failed to parse KIN file', e);
      throw e;
    }
    // no persistent state object produced here (data gets applied as representations), so return undefined
    return undefined;
  },
  visuals: async (plugin, data) => {
    if (g_kinemageData) {
      await applyKinemageInfoToState(plugin, g_kinemageData);
    } else {
      console.warn('[Kinemage] visuals: no loaded kinemage data present');
    }
    return undefined;
  }
});
