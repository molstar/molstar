/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russ Taylor <russ@reliasolve.com>
 */

/** Based on the ../anvil extension. */

import { Vec3, Mat3 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { KinemageDataProvider, KinemageData } from './prop';
import { StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginDragAndDropHandler } from '../../mol-plugin-state/manager/drag-and-drop';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginContext } from '../../mol-plugin/context';
import { DefaultQueryRuntimeTable } from '../../mol-script/runtime/query/compiler';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { shapePointsFromKin, shapeLinesFromKin, shapeMeshFromKin, shapeSpheresFromKin } from '../../mol-model-formats/shape/kin';
import { Kinemage } from '../../mol-io/reader/kin/schema';
import { DataFormatProvider } from '../../mol-plugin-state/formats/provider';
import { Camera } from '../../mol-canvas3d/camera';
import { PluginCommands } from '../../mol-plugin/commands';
import { StateObjectRef } from '../../mol-state';
import { getPluginBoundingSphere } from '../../mol-plugin-state/manager/focus-camera/focus-object';

const Tag = KinemageData.Tag;

const Transform = StateTransformer.builderFactory('sb-kinemage');

let g_kinemageData: KinemageData | undefined = undefined;

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
    /// @todo What is the point of setting this and others to isHidden?
    data: PD.Value<Kinemage>(undefined as any, { isHidden: true })
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
        data: PD.Value<Kinemage>(undefined as any, { isHidden: true })
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
    data: PD.Value<Kinemage>(undefined as any, { isHidden: true })
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
    data: PD.Value<Kinemage>(undefined as any, { isHidden: true })
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

export const KinemageViewProvider = Transform({
  name: 'sb-kinemage-view-provider',
  display: { name: 'Kinemage View Provider' },
  from: PluginStateObject.Root,
  to: PluginStateObject.Format.Json, // store view metadata as JSON data node
  params: {
    name: PD.Text(''),
    snapshot: PD.Value<Partial<Camera.Snapshot>>(undefined as any, { isHidden: true })
  }
})({
  apply({ params }) {
    return Task.create('Kinemage View Provider', async ctx => {
      // PluginStateObject.Format.Json holds arbitrary JSON-like data; create instance with the payload
      // Pass the view name as the node label so the State Tree shows the provided name instead of "JSON Data"
      const viewName = 'View ' + String(params.name || '');
      return new PluginStateObject.Format.Json(
        { name: viewName, snapshot: params.snapshot } as any,
        { label: viewName }
      );
    });
  }
});

export const KinemageMasterProvider = Transform({
  name: 'sb-kinemage-master-provider',
  display: { name: 'Kinemage Master Provider' },
  from: PluginStateObject.Root,
  to: PluginStateObject.Format.Json, // store view metadata as JSON data node
  params: {
    name: PD.Text(''),
    masterData: PD.Text('') /// @todo Fill this in with actual master data if needed, and parse it in the provider
  }
})({
  apply({ params }) {
    return Task.create('Kinemage Master Provider', async ctx => {
      // PluginStateObject.Format.Json holds arbitrary JSON-like data; create instance with the payload
      // Pass the view name as the node label so the State Tree shows the provided name instead of "JSON Data"
      const masterName = String(params.name);
      return new PluginStateObject.Format.Json(
        { name: masterName, masterData: params.masterData } as any,
        { label: masterName }
      );
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
        private applyViewAction?: any;
        private selectedSub?: any;
        private visibilitySub?: any;
        private visibilityMap = new Map<string, boolean>();

        register(): void {
            DefaultQueryRuntimeTable.addCustomProp(this.provider.descriptor);

            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);

            this.ctx.managers.dragAndDrop.addHandler(KinemageDragAndDropHandler.name, KinemageDragAndDropHandler.handle);

            // Register .kin file handler so opening/dropping .kin is supported via the data formats system
            this.ctx.dataFormats.add('KIN', KINFormatProvider);

            const applyViewAction = {
                // keep fields minimal and compatible with the actions system
                // adjust names if your repo expects different keys (use MVS LoadMvsData as example)
                label: 'Apply View',
                isApplicable: (a: any) => {
                    const d = a?.data;
                    return !!(d && d.data && (d.data as any).snapshot);
                },
                apply: async (a: any) => {
                    const node = a?.data;
                    const snapshot = node?.data?.snapshot as Partial<Camera.Snapshot> | undefined;
                    if (snapshot) await applyViewSnapshot(this.ctx, snapshot);
                }
            };

            // When one of the state objects is selected in the GUI, handle the appropriate behavior:
            // For a view object (which has a snapshot), update the camera to its snapshot.
            this.selectedSub = this.ctx.state.data.behaviors.currentObject.subscribe((e: any) => {
              const ref = e.ref;
              // state.select returns an array of cells; the first is the matching cell
              const cell = this.ctx.state.data.select(ref)[0];
              const obj = cell?.obj;
              const nodeData = obj?.data;
              if (nodeData && nodeData.snapshot) {
                applyViewSnapshot(this.ctx, nodeData.snapshot);
              }
            });

            // When one of the state objects is has its visibility changed in the GUI, handle the appropriate behavior:
            // For a master, turn on or off the value and regenerate appropriate geometry.
            this.visibilitySub = this.ctx.state.data.events.cell.stateUpdated.subscribe((e: any) => {
              const ref = e.ref;
              // state.select returns an array of cells; the first is the matching cell
              const cell = this.ctx.state.data.select(ref)[0];

              const obj = cell?.obj;
              const nodeData = obj?.data;
              if (nodeData && nodeData.masterData) {
                // transform.state usually holds the runtime state; inspect it to find the exact visibility flag used
                const st = (cell.transform && cell.transform.state) || cell.state || {};
                const nowHidden = !!st.isHidden;
                const prev = this.visibilityMap.get(ref);
                if (prev === undefined) {
                  // first time seeing this cell, store and return (or handle initial state)
                  this.visibilityMap.set(ref, nowHidden);
                  return;
                }
                if (prev !== nowHidden) {
                  this.visibilityMap.set(ref, nowHidden);
                  // visibility changed for the object referenced by `ref`
                  if (nowHidden) {
                    // Became hidden
                    console.log(`XXX Object ${nodeData.masterData} became hidden`);
                  } else {
                    // Became visible
                    console.log(`XXX Object ${nodeData.masterData} became visible`);
                  }
                  /// @todo Handle properly
                }
              }
            });

            // store and register the exact object reference
            this.applyViewAction = applyViewAction;
            this.ctx.state.data.actions.add(this.applyViewAction);
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

            // Remove the state action if we registered one
            if (this.applyViewAction) {
              this.ctx.state.data.actions.remove(this.applyViewAction);
              this.applyViewAction = undefined;
            }

            // Remove the action subscriptions if we created them
            if (this.selectedSub) {
              this.selectedSub.unsubscribe();
              this.selectedSub = undefined;
            }
            if (this.visibilitySub) {
              this.visibilitySub.unsubscribe();
              this.visibilitySub = undefined;
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

/** Centralized helper to apply kinemage content into plugin state (re-used by drag handler and programmatic loader) */
async function applyKinemageInfoToState(plugin: PluginContext, kinInfo: KinemageData) {
  const update = plugin.state.data.build();
  const createdMasterPairs: { selector: StateObjectRef<PluginStateObject.Format.Json>, visible: boolean }[] = [];

  for (const kinData of kinInfo.kinemages) {

    // Keep list of created selectors for this kinemage (shapes / representations etc.)
    const createdShapeSelectors: StateObjectRef<any>[] = [];

    // Iterate over all entries in the view dictionary. Do this before creating shapes so that the views show up
    // in the state tree first and don't change order when we update the masters.
    const createdViewRefs: StateObjectRef<PluginStateObject.Format.Json>[] = [];
    for (const [viewKey, viewObj] of Object.entries(kinData.viewDict)) {
      const viewName = viewObj.name || `View ${viewKey}`;

      // If center is specified, then we will use that as the camera target. Otherwise, we will use the origin.
      const center = Vec3.create(0, 0, 0);
      if (viewObj.center) {
        Vec3.set(center, viewObj.center[0], viewObj.center[1], viewObj.center[2]);
      }

      // Make an orientation matrix based on the matrix provided, otherwise make the identity matrix.
      const orientation: Mat3 = Mat3.identity();
      if (viewObj.matrix) {
        /// Transpose this so that it matches the order for Molstar's Mat3 (row-major vs column-major).
        Mat3.fromArray(orientation, viewObj.matrix, 0);
        Mat3.transpose(orientation, orientation);
      }

      // Rotate the +Z axis by the orientation to see which way points to the camera.
      const zAxis = Vec3.create(0, 0, 1);
      Vec3.transformMat3(zAxis, zAxis, orientation);

      // Rotate the +Y axis by the orientation to see which way points up.
      const yAxis = Vec3.create(0, 1, 0);
      Vec3.transformMat3(yAxis, yAxis, orientation);

      // If span is specified, then we go that distance along Z to find the camera position (90 degree FOV).
      // Otherwise, we go a default distance of 100 along Z.
      let distance = 100;
      if (viewObj.span) {
        distance = viewObj.span;
      } else if (viewObj.zoom) {
        /// @todo If zoom is specified, then we need to do more computations based on the bounds on the geometry.
      }
      Vec3.scale(zAxis, zAxis, distance);
      const position = Vec3.create(0, 0, 100);
      Vec3.add(position, center, zAxis);

      // If the zslab is specified, then we set the radius to it; otherwise, we use a default of 100.
      // When the zslab value is 200, it should match the same as the span (half is percent of half span).
      let radius = 100;
      if (viewObj.zslab) {
        const scale = viewObj.zslab / 200;
        // Scale by 0.5 here to match the behavior of KiNG
        radius = 0.5 * distance * scale;
      }

      // Fill in the camera shapshot
      let snap: Camera.Snapshot;
      snap = {
        mode: 'orthographic',   ///< Make this orthographic by default, to match Kinemage defaults
        fov: Math.PI / 4,       ///< 90-degree view by default for Molstar

        position: position,
        up: yAxis,
        target: center,

        radius: radius,
        radiusMax: 1e4,
        fog: 0,
        clipFar: true,
        minNear: 1,
        minFar: 1
      };

      // Create a state object for the view (visible in State Tree)
      const viewNode = update
        .toRoot()
        .apply(KinemageViewProvider, { name: viewName, snapshot: snap });

      // Store the selector for UI wiring
      createdViewRefs.push(viewNode.selector as StateObjectRef<PluginStateObject.Format.Json>);
    }

    // Iterate over all of the masterDict entries and create a state object for each master.
    // Add a callback handler for visibility changes on each and print whether it is visible or not.
    // Name each after the master dictionary key.
    // Set its visibility according to the visible entry in the masterDict.
    const createdMasterRefs: StateObjectRef<PluginStateObject.Format.Json>[] = [];
    for (const [masterKey, masterInfo] of Object.entries(kinData.masterDict)) {
      const masterName = masterKey;

      const masterNode = update
        .toRoot()
        .apply(KinemageMasterProvider, { name: masterName, masterData: masterKey });

      // Store the selector for UI wiring
      createdMasterRefs.push(masterNode.selector as StateObjectRef<PluginStateObject.Format.Json>);

      // capture desired visibility so we can set transform state after commit
      const visible = !!(masterInfo && (masterInfo as any).visible);
      createdMasterPairs.push({ selector: masterNode.selector as StateObjectRef<PluginStateObject.Format.Json>, visible });
    }

    // Generate all of the shapes for this kinemage, each shape type having its own provider and representation.
    if (kinData.dotLists.length > 0) {
      const node = await update
        .toRoot()
        .apply(KinemageShapePointsProvider, { data: kinData })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);
      createdShapeSelectors.push(node.selector as StateObjectRef<any>);
    }
    if (kinData.vectorLists.length > 0) {
      const node = await update
        .toRoot()
        .apply(KinemageShapeLinesProvider, { data: kinData })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);
      createdShapeSelectors.push(node.selector as StateObjectRef<any>);
    }
    if (kinData.ribbonLists.length > 0) {
      const node = await update
        .toRoot()
        .apply(KinemageShapeMeshProvider, { data: kinData })
        .apply(StateTransforms.Representation.ShapeRepresentation3D, { doubleSided: true });
      createdShapeSelectors.push(node.selector as StateObjectRef<any>);
    }
    if (kinData.ballLists.length > 0) {
      const node = await update
        .toRoot()
        .apply(KinemageShapeSpheresProvider, { data: kinData })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);
      createdShapeSelectors.push(node.selector as StateObjectRef<any>);
    }

    // Store the created selector list for this kinemage so callback handlers can destroy / re-create.
    if (createdShapeSelectors.length > 0) {
      g_kinemageShapeSelectors.set(kinData as Kinemage, createdShapeSelectors);
    }
  }
  await update.commit();

  // Helper: robustly resolve a transform ref from different selector shapes without changing other modules.
  function resolveSelectorRef(sel: any): string | undefined {
    if (!sel) return undefined;
    if (typeof sel === 'string') return sel;
    // StateObjectSelector has a .ref property containing the transform ref
    if (sel.ref && typeof sel.ref === 'string') return sel.ref;
    // StateObjectCell or StateObjectCell-like has transform.ref
    if (sel.transform && typeof sel.transform.ref === 'string') return sel.transform.ref;
    // Some callers may wrap selector as { cell: ... }
    if (sel.cell && sel.cell.transform && typeof sel.cell.transform.ref === 'string') return sel.cell.transform.ref;
    // fallback: try existing util (in case shape matches)
    try {
      return StateObjectRef.resolveRef(sel as any);
    } catch {
      return undefined;
    }
  }

  // Ensure the State Tree visibility matches the masters' initial 'visible' flags.
  // The UI commonly uses `isHidden` on transform state; set it here so the created
  // master nodes show the expected checked/unchecked visibility in the GUI.
  for (const pair of createdMasterPairs) {
    try {
      const ref = resolveSelectorRef(pair.selector);
      if (!ref) continue;

      // Set the isHidden state for this master based on the visibile flag stored above.
      plugin.state.data.updateCellState(ref, (old: any) => {
        const s = { ...(old || {}) };
        s.isHidden = !pair.visible;
        return s;
      });
    } catch (e) {
      // Non-fatal: log and continue
      console.warn('Failed to set master visibility for', pair, e);
    }
  }
}

/** Programmatic loader: load a single File (a .kin) into the plugin state.
 * Runs the import inside a Task so it has a runtime and asset context similar to drag-and-drop.
 * Returns true if at least one Kinemage was added.
 */
export async function loadKinemageFile(plugin: PluginContext, file: File): Promise<boolean> {
  let applied = false;
  const task = Task.create('Load KIN file', async ctx => {
    const kinData = await KinemageData.open(file);
    if (!g_kinemageData) {
      g_kinemageData = kinData;
    } else {
      // If we already have kinemage data loaded, append to the list of kinemages and make the last one active
      g_kinemageData.kinemages.push(...kinData.kinemages);
      g_kinemageData.activeKinemage = g_kinemageData.kinemages.length - 1;
    }
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
  stringExtensions: ['kin'],
  parse: async (plugin, data) => {
    try {
      /// @todo Consider allowing the handler to directly take a string and not to open and read the file.
      /// This avoids having to create a File object from the string content.
      const file = fileFromPayload(data);
      let p = loadKinemageFile(plugin, file);
      await p;
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
    }
    return undefined;
  }
});
