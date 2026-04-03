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

export const KinemageViewProvider = Transform({
  name: 'sb-kinemage-view-provider',
  display: { name: 'Kinemage View Provider' },
  from: PluginStateObject.Root,
  to: PluginStateObject.Format.Json, // store view metadata as JSON data node
  params: {
    name: PD.Text(''),
    snapshot: PD.Value<Partial<Camera.Snapshot>>(undefined as any, {})
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

export const KinemageGroupProvider = Transform({
  name: 'sb-kinemage-group-provider',
  display: { name: 'Kinemage Group Provider' },
  from: PluginStateObject.Root,
  to: PluginStateObject.Format.Json, // store view metadata as JSON data node
  params: {
    name: PD.Text(''),
    groupData: PD.Text(''), /// @todo Fill this in with actual group data if needed, and parse it in the provider
    data: PD.Value<Kinemage>(undefined as any, {}) // store kinData reference so visibility handlers can access it
  }
})({
  apply({ params }) {
    return Task.create('Kinemage Group Provider', async ctx => {
      // PluginStateObject.Format.Json holds arbitrary JSON-like data; create instance with the payload
      // Pass the view name as the node label so the State Tree shows the provided name instead of "JSON Data"
      const groupName = String(params.name);
      return new PluginStateObject.Format.Json(
        { name: groupName, groupData: params.groupData, kinData: params.data } as any,
        { label: groupName }
      );
    });
  }
});

export const KinemageSubgroupProvider = Transform({
  name: 'sb-kinemage-subgroup-provider',
  display: { name: 'Kinemage Subgroup Provider' },
  from: PluginStateObject.Root,
  to: PluginStateObject.Format.Json, // store view metadata as JSON data node
  params: {
    name: PD.Text(''),
    subgroupData: PD.Text(''), /// @todo Fill this in with actual subgroup data if needed, and parse it in the provider
    data: PD.Value<Kinemage>(undefined as any, {}) // store kinData reference so visibility handlers can access it
  }
})({
  apply({ params }) {
    return Task.create('Kinemage Subgroup Provider', async ctx => {
      // PluginStateObject.Format.Json holds arbitrary JSON-like data; create instance with the payload
      // Pass the view name as the node label so the State Tree shows the provided name instead of "JSON Data"
      const subgroupName = String(params.name);
      return new PluginStateObject.Format.Json(
        { name: subgroupName, subgroupData: params.subgroupData, kinData: params.data } as any,
        { label: subgroupName }
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
    masterData: PD.Text(''), /// @todo Fill this in with actual master data if needed, and parse it in the provider
    data: PD.Value<Kinemage>(undefined as any, {}) // store kinData reference so visibility handlers can access it
  }
})({
  apply({ params }) {
    return Task.create('Kinemage Master Provider', async ctx => {
      // PluginStateObject.Format.Json holds arbitrary JSON-like data; create instance with the payload
      // Pass the view name as the node label so the State Tree shows the provided name instead of "JSON Data"
      const masterName = String(params.name);
      return new PluginStateObject.Format.Json(
        { name: masterName, masterData: params.masterData, kinData: params.data } as any,
        { label: masterName }
      );
    });
  }
});

export const KinemageAnimateProvider = Transform({
  name: 'sb-kinemage-animate-provider',
  display: { name: 'Kinemage Animate Provider' },
  from: PluginStateObject.Root,
  to: PluginStateObject.Format.Json, // store view metadata as JSON data node
  params: {
    name: PD.Text(''),
    animateData: PD.Text(''), /// @todo Fill this in with actual animate data if needed, and parse it in the provider
    data: PD.Value<Kinemage>(undefined as any, {}) // store kinData reference so visibility handlers can access it
  }
})({
  apply({ params }) {
    return Task.create('Kinemage Animate Provider', async ctx => {
      // PluginStateObject.Format.Json holds arbitrary JSON-like data; create instance with the payload
      // Pass the view name as the node label so the State Tree shows the provided name instead of "JSON Data"
      const animateName = String(params.name);
      return new PluginStateObject.Format.Json(
        { name: animateName, animateData: params.animateData, kinData: params.data, firedOnce: false } as any,
        { label: animateName }
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
        private selectedSub?: any;
        private visibilitySub?: any;
        private visibilityMap = new Map<string, boolean>();

        register(): void {
            DefaultQueryRuntimeTable.addCustomProp(this.provider.descriptor);

            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);

            this.ctx.managers.dragAndDrop.addHandler(KinemageDragAndDropHandler.name, KinemageDragAndDropHandler.handle);

            // Register .kin file handler so opening/dropping .kin is supported via the data formats system
            this.ctx.dataFormats.add('KIN', KINFormatProvider);

            // When one of the state objects is selected in the GUI, handle the appropriate behavior:
            // For a view object (which has a snapshot), update the camera to its snapshot.
            // For an animation button, have it adjust the various visibilities and then regenerate shapes.
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
            // For an animation object, adjust which group is visible and then regenerate shapes.
            // For a master, group, or subgroup, turn on or off the value and regenerate appropriate geometry.
            this.visibilitySub = this.ctx.state.data.events.cell.stateUpdated.subscribe(async (e: any) => {
              const ref = e.ref;
              const cell = this.ctx.state.data.select(ref)[0];
              const obj = cell?.obj;
              const nodeData = obj?.data;
              let madeChanges = false;
              if (nodeData && nodeData.animateData) {
                // If we have not yet fired, ignore this event because it is just the creation of the node.
                if (!nodeData.firedOnce) {
                  nodeData.firedOnce = true;
                  return;
                }
                const kinData = nodeData.kinData as Kinemage;
                if (nodeData.animateData === 'animate') {
                  // Increment the activeAnimateGroup index and wrap around if needed,
                  // then make the selected group visible and the others not.
                  kinData.activeAnimateGroup = (kinData.activeAnimateGroup + 1) % kinData.groupsAnimate.length;
                  for (let i = 0; i < kinData.groupsAnimate.length; i++) {
                    const groupName = kinData.groupsAnimate[i];
                    const groupInfo = kinData.groupDict[groupName];
                    groupInfo.off = i !== kinData.activeAnimateGroup;
                    // Also set the GUI element visibility state to match the kinemage data,
                    // so that the GUI reflects which group is currently active.
                    /// @todo
                  }

                } else if (nodeData.animateData === '2animate') {
                  // Increment the activeAnimateGroup2 index and wrap around if needed,
                  // then make the selected group visible and the others not.
                  kinData.activeAnimateGroup2 = (kinData.activeAnimateGroup2 + 1) % kinData.groupsAnimate2.length;
                  for (let i = 0; i < kinData.groupsAnimate2.length; i++) {
                    const groupName = kinData.groupsAnimate2[i];
                    const groupInfo = kinData.groupDict[groupName];
                    groupInfo.off = i !== kinData.activeAnimateGroup2;
                    // Also set the GUI element visibility state to match the kinemage data,
                    // so that the GUI reflects which group is currently active.
                    /// @todo
                  }
                }

                // Indicate that we need to rebuild the shapes for this kinemage based on the animation change.
                madeChanges = true;
              }

              if (nodeData && (nodeData.masterData || nodeData.groupData || nodeData.subgroupData)) {
                const st = (cell.transform && cell.transform.state) || cell.state || {};
                const nowHidden = !!st.isHidden;

                // Change the record of visibility so we know whether we became visible or invisible.
                const prev = this.visibilityMap.get(ref);
                if (prev === undefined) {
                  this.visibilityMap.set(ref, nowHidden);
                  return;
                }
                if (prev !== nowHidden) {
                  this.visibilityMap.set(ref, nowHidden);
                  const kinRef: Kinemage | undefined = nodeData.kinData;
                  if (!kinRef) return;

                  // Set the Kinemage master visibility based on the isHidden state of the transform.
                  // When the transform is hidden, the master is invisible, and vice versa.
                  if (nodeData.groupData) kinRef.groupDict[nodeData.groupData].off = nowHidden;
                  if (nodeData.subgroupData) kinRef.subgroupDict[nodeData.subgroupData].off = nowHidden;
                  if (nodeData.masterData) kinRef.masterDict[nodeData.masterData].visible = !nowHidden;

                  // Indicate that we need to rebuild the shapes for this kinemage based on the animation change.
                  madeChanges = true;
                }
              }

              if (madeChanges) {
                // capture current camera snapshot so we can restore view after re-creating shapes
                const curSnap = (this.ctx.canvas3d && (this.ctx.canvas3d as any).camera && (this.ctx.canvas3d as any).camera.getSnapshot)
                  ? (this.ctx.canvas3d as any).camera.getSnapshot()
                  : undefined;

                // recreate: ensure old selectors are cleared, then build new ones with a fresh builder
                const kinRef: Kinemage | undefined = nodeData.kinData;
                if (!kinRef) return;
                await destroyShapesForKinemage(this.ctx, kinRef);
                const update = this.ctx.state.data.build();
                try {
                  await createShapesForKinemage(this.ctx, update, kinRef);
                  await update.commit();

                  // restore camera snapshot to avoid the temporary zoom-out caused by removing geometry
                  if (curSnap) {
                    try {
                      await applyViewSnapshot(this.ctx, curSnap);
                    } catch (e) {
                      console.warn('Failed to restore camera snapshot after recreating shapes', e);
                    }
                  }
                } catch (err) {
                  console.error('Failed to recreate kinemage shapes', err);
                }
              }

            });
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

/** Helper function to create the shapes for a kinemage */
async function createShapesForKinemage(plugin: PluginContext, update: StateBuilder.Root,  kinData: Kinemage) {
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
async function destroyShapesForKinemage(plugin: PluginContext, kinData: Kinemage) {
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

/** Centralized helper to apply kinemage content into plugin state (re-used by drag handler and programmatic loader) */
async function applyKinemageInfoToState(plugin: PluginContext, kinInfo: KinemageData) {
  const update = plugin.state.data.build();

  for (const kinData of kinInfo.kinemages) {

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

    // If there are any entries in the groupsAnimate list, create a state object for the animate provider so it shows up in the State Tree.
    if (kinData.groupsAnimate.length > 0) {
      update
        .toRoot()
        .apply(KinemageAnimateProvider, { name: 'Animate (change vis)', animateData: 'animate', data: kinData });
    }

    // If there are any entries in the groupsAnimate2 list, create a state object for the animate provider so it shows up in the State Tree.
    if (kinData.groupsAnimate2.length > 0) {
      update
        .toRoot()
        .apply(KinemageAnimateProvider, { name: 'Animate2 (change vis)', animateData: '2animate', data: kinData });
    }

    // Iterate over all of the groupDict entries and create a state object for each group.
    // Name each after the group dictionary key.
    for (const [groupKey, groupInfo] of Object.entries(kinData.groupDict)) {

      // Only make state object for this group if it did not have noButton set.
      if (!(groupInfo as any).nobutton) {

        // capture desired visibility and set it as the transform state at creation time
        const visible = !(groupInfo as any).off;
        update
          .toRoot()
          .apply(KinemageGroupProvider, { name: groupKey, groupData: groupKey, data: kinData }, { state: { isHidden: !visible } });
      }

      // Iterate over all of the subgroupDict entries under this group and create a state object for each subgroup.
      // Name each after the subgroup dictionary key, which is in the format "GroupName:SubgroupName" to preserve tree structure.
      for (const [subgroupKey, subgroupInfo] of Object.entries(kinData.subgroupDict)) {
        // Skip subgroups that don't belong to this group (based on the naming convention of "GroupName:SubgroupName")
        if (!subgroupKey.startsWith(groupKey + ':')) continue;

        // Skip this subgroup if its parent group is dominant.
        if ((groupInfo as any).dominant) continue;

        // Skip this subgroup if it has noButton set.
        if ((subgroupInfo as any).nobutton) continue;

        // capture desired visibility and set it as the transform state at creation time
        const visible = !(subgroupInfo as any).off;
        update
          .toRoot()
          .apply(KinemageSubgroupProvider, { name: subgroupKey, subgroupData: subgroupKey, data: kinData }, { state: { isHidden: !visible } });
      }
    }

    // Iterate over all subgroupDict entries that don't have a parent group and create state objects for them as well,
    // so they show up in the State Tree.
    for (const [subgroupKey, subgroupInfo] of Object.entries(kinData.subgroupDict)) {
      // Skip subgroups that belong to a group (characters before the ":" indicate the parent group)
      if (subgroupKey[0] !== ':' || subgroupKey.indexOf(':') === -1) continue;

      // capture desired visibility and set it as the transform state at creation time
      const visible = !(subgroupInfo as any).off;
      update
        .toRoot()
        .apply(KinemageSubgroupProvider, { name: subgroupKey, subgroupData: subgroupKey, data: kinData }, { state: { isHidden: !visible } });
    }

    // Iterate over all of the masterDict entries and create a state object for each master.
    // Name each after the master dictionary key.
    for (const [masterKey, masterInfo] of Object.entries(kinData.masterDict)) {
      const masterName = masterKey;

      // capture desired visibility and set it as the transform state at creation time
      const visible = !!(masterInfo && (masterInfo as any).visible);
      update
        .toRoot()
        .apply(KinemageMasterProvider, { name: masterName, masterData: masterKey, data: kinData }, { state: { isHidden: !visible } });
    }

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
      // Clear the kinemage data after applying so that if the user drags in another file, it doesn't get merged with the previous one.
      g_kinemageData = undefined;
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
    // Clear the kinemage data after applying so that if the user drags in another file, it doesn't get merged with the previous one.
    g_kinemageData = undefined;

    return undefined;
  }
});
