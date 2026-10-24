/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russ Taylor <russ@reliasolve.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Based on the ../anvil extension. */

// The kinemage data model and its plugin state (transforms, shapes, the .kin format), without UI: the v5 behavior
// (behavior.ts, which adds the right-panel controls) and the Prism extension (prism-viewer.ts) both stand on it.

import { StringLike } from '@molstar/core/util/string-like';
import { Vec3, Mat3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { KinemageData } from '@molstar/kinemage-extension/prop';
import { StateTransformer, StateBuilder } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import type { PluginDragAndDropEntry } from '@molstar/plugin/state/manager/drag-and-drop';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { PluginContext } from '@molstar/plugin/context';
import { ShapeRepresentation3D } from '@molstar/plugin/state/transforms/shape/representation';
import { RawData } from '@molstar/plugin/state/transforms/data/fetch';
import {
  shapePointsFromKin,
  shapeLinesFromKin,
  shapeMeshFromKin,
  shapeSpheresFromKin,
} from '@molstar/kinemage-extension/kin';
import type { Kinemage } from '@molstar/kinemage-extension/reader/schema';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { PluginCommands } from '@molstar/plugin/commands';
import { getPluginBoundingSphere } from '@molstar/plugin/state/manager/focus-camera/focus-object';
import { StateObjectSelector } from '@molstar/core/state';
import { Color } from '@molstar/core/util/color';

const Transform = StateTransformer.builderFactory('sb-kinemage');

/**
 * State object to hold parsed Kinemage data
 */
export class KinemageObject extends PluginStateObject.Create<KinemageData>({ name: 'Kinemage', typeClass: 'Object' }) {}

/**
 * Visibility state for kinemage elements - stores which items are VISIBLE (not hidden)
 */
export interface KinemageVisibilityState {
  /** Map of group name -> visibility (true = visible, false = hidden/off) */
  groupVisibility: Map<string, boolean>;
  /** Map of subgroup name -> visibility (true = visible, false = hidden/off) */
  subgroupVisibility: Map<string, boolean>;
  /** Map of master name -> visibility (true = visible, false = hidden) */
  masterVisibility: Map<string, boolean>;
  activeAnimateGroup: number;
  activeAnimateGroup2: number;
}

/**
 * Apply a saved snapshot object (from a view state node) to the plugin camera.
 * Use PluginCommands.Camera.SetSnapshot so transitions and canvas props are handled properly.
 */
export async function applyViewSnapshot(plugin: PluginContext, snapshot: Partial<Camera.Snapshot>) {
  if (!snapshot) return;

  // Set background color to black
  await setKinemageBackground(plugin);

  // If the snapshot provides a target, adjust the canvas `sceneRadiusFactor` so the scene isn't clipped
  // when we switch camera.
  if (snapshot.target) {
    try {
      const boundingSphere = getPluginBoundingSphere(plugin);
      if (boundingSphere && boundingSphere.radius > 0) {
        const offset = Vec3.distance(snapshot.target as Vec3, boundingSphere.center);
        const sceneRadiusFactor = (boundingSphere.radius + offset) / boundingSphere.radius;
        plugin.setCanvas3DProps({ sceneRadiusFactor });
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
  from: [PluginStateObject.Data.String, PluginStateObject.Data.Binary, PluginStateObject.Data.Response],
  to: KinemageObject,
  params: {
    label: PD.Optional(PD.Text('', { description: 'Label for the Kinemage data' })),
  },
})({
  apply({ a, params }) {
    return Task.create('Parse Kinemage', async (ctx) => {
      const input = StringLike.toString(await PluginStateObject.Data.readText(ctx, a));
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
            minFar: 1,
          };

          (kinData as any).viewSnapshots[viewKey] = snap;
        }
      }

      const label = params.label || data.kinemages[0]?.caption || 'Kinemage';
      return new KinemageObject(data, { label, description: `Kinemage with ${data.kinemages.length} kinemage(s)` });
    });
  },
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
      index: PD.Select(0, options, { description: 'Which kinemage to use' }),
    };
  },
})({
  apply({ a, params }) {
    return Task.create('Select Kinemage', async (ctx) => {
      const kinData = a.data.kinemages[params.index];
      if (!kinData) {
        throw new Error(`No kinemage found at index ${params.index}`);
      }

      const label = kinData.pdbfile || kinData.caption || `Kinemage ${params.index}`;

      // Store the kinemage data in a Format.Json node so downstream transforms can access it
      return new PluginStateObject.Format.Json({ kinData }, { label, description: kinData.text || '' });
    });
  },
});

/**
 * Visibility Controller Transform - centralizes visibility state for all shape types
 * Stores visibility as key-value pairs where key is the item name and value is boolean (true = visible)
 */
export const KinemageVisibilityController = Transform({
  name: 'sb-kinemage-visibility-controller',
  display: { name: 'Kinemage Visibility Controller' },
  from: PluginStateObject.Format.Json,
  to: PluginStateObject.Format.Json,
  params: (a) => {
    const kinData = (a?.data as any)?.kinData as Kinemage | undefined;
    if (!kinData) {
      return {
        groupVisibility: PD.Value<{ [key: string]: boolean }>({}),
        subgroupVisibility: PD.Value<{ [key: string]: boolean }>({}),
        masterVisibility: PD.Value<{ [key: string]: boolean }>({}),
        activeAnimateGroup: PD.Numeric(0, { min: 0, max: 0, step: 1 }, { description: 'Active animate group index' }),
        activeAnimateGroup2: PD.Numeric(0, { min: 0, max: 0, step: 1 }, { description: 'Active animate2 group index' }),
      };
    }

    // Build initial visibility from parsed data
    const groupVisibility: { [key: string]: boolean } = {};
    const subgroupVisibility: { [key: string]: boolean } = {};
    const masterVisibility: { [key: string]: boolean } = {};

    for (const [groupKey, groupInfo] of Object.entries(kinData.groupDict)) {
      groupVisibility[groupKey] = !(groupInfo as any).off;
    }

    for (const [subgroupKey, subgroupInfo] of Object.entries(kinData.subgroupDict)) {
      subgroupVisibility[subgroupKey] = !(subgroupInfo as any).off;
    }

    for (const [masterKey, masterInfo] of Object.entries(kinData.masterDict)) {
      masterVisibility[masterKey] = !!(masterInfo as any).visible;
    }

    return {
      groupVisibility: PD.Value(groupVisibility, { isHidden: true }),
      subgroupVisibility: PD.Value(subgroupVisibility, { isHidden: true }),
      masterVisibility: PD.Value(masterVisibility, { isHidden: true }),
      activeAnimateGroup: PD.Numeric(
        0,
        { min: 0, max: Math.max(0, kinData.groupsAnimate.length - 1), step: 1 },
        { description: 'Active animate group index', isHidden: true },
      ),
      activeAnimateGroup2: PD.Numeric(
        0,
        { min: 0, max: Math.max(0, kinData.groupsAnimate2.length - 1), step: 1 },
        { description: 'Active animate2 group index', isHidden: true },
      ),
    };
  },
})({
  apply({ a, params }) {
    return Task.create('Kinemage Visibility Controller', async (ctx) => {
      const kinData = (a.data as any).kinData as Kinemage;
      if (!kinData) {
        throw new Error('No kinData found in parent Format.Json node');
      }

      // Store visibility state alongside kinData
      const visibilityState: KinemageVisibilityState = {
        groupVisibility: new Map(Object.entries(params.groupVisibility)),
        subgroupVisibility: new Map(Object.entries(params.subgroupVisibility)),
        masterVisibility: new Map(Object.entries(params.masterVisibility)),
        activeAnimateGroup: params.activeAnimateGroup,
        activeAnimateGroup2: params.activeAnimateGroup2,
      };

      return new PluginStateObject.Format.Json(
        { kinData, visibilityState },
        { label: a.label, description: a.description },
      );
    });
  },
});

export const KinemageShapePointsProvider = Transform({
  name: 'sb-kinemage-shape-points-provider',
  display: { name: 'Kinemage Shape Points Provider' },
  from: PluginStateObject.Format.Json,
  to: PluginStateObject.Shape.Provider,
  params: {},
})({
  apply({ a }) {
    return Task.create('Kinemage Points Shape Provider', async (ctx) => {
      const kinData = (a.data as any).kinData as Kinemage;
      const visibilityState = (a.data as any).visibilityState as KinemageVisibilityState | undefined;

      if (!kinData) {
        throw new Error('No kinData found in parent Format.Json node');
      }

      const provider = await shapePointsFromKin(
        kinData,
        visibilityState,
        { transforms: undefined },
        'Dots',
      ).runInContext(ctx);
      return new PluginStateObject.Shape.Provider(provider as any, {
        label: kinData.pdbfile || kinData.caption || 'Kinemage Points',
        description: kinData.text || '',
      });
    });
  },
});

export const KinemageShapeLinesProvider = Transform({
  name: 'sb-kinemage-shape-lines-provider',
  display: { name: 'Kinemage Shape Lines Provider' },
  from: PluginStateObject.Format.Json,
  to: PluginStateObject.Shape.Provider,
  params: {},
})({
  apply({ a }) {
    return Task.create('Kinemage Lines Shape Provider', async (ctx) => {
      const kinData = (a.data as any).kinData as Kinemage;
      const visibilityState = (a.data as any).visibilityState as KinemageVisibilityState | undefined;

      if (!kinData) {
        throw new Error('No kinData found in parent Format.Json node');
      }

      const provider = await shapeLinesFromKin(kinData, visibilityState).runInContext(ctx);
      return new PluginStateObject.Shape.Provider(provider as any, {
        label: kinData.pdbfile || kinData.caption || 'Kinemage Lines',
        description: kinData.text || '',
      });
    });
  },
});

export const KinemageShapeMeshProvider = Transform({
  name: 'sb-kinemage-shape-mesh-provider',
  display: { name: 'Kinemage Shape Mesh Provider' },
  from: PluginStateObject.Format.Json,
  to: PluginStateObject.Shape.Provider,
  params: {},
})({
  apply({ a }) {
    return Task.create('Kinemage Mesh Shape Provider', async (ctx) => {
      const kinData = (a.data as any).kinData as Kinemage;
      const visibilityState = (a.data as any).visibilityState as KinemageVisibilityState | undefined;

      if (!kinData) {
        throw new Error('No kinData found in parent Format.Json node');
      }

      const provider = await shapeMeshFromKin(kinData, visibilityState).runInContext(ctx);
      return new PluginStateObject.Shape.Provider(provider as any, {
        label: kinData.pdbfile || kinData.caption || 'Kinemage Meshes',
        description: kinData.text || '',
      });
    });
  },
});

export const KinemageShapeSpheresProvider = Transform({
  name: 'sb-kinemage-shape-spheres-provider',
  display: { name: 'Kinemage Shape Spheres Provider' },
  from: PluginStateObject.Format.Json,
  to: PluginStateObject.Shape.Provider,
  params: {},
})({
  apply({ a }) {
    return Task.create('Kinemage Spheres Shape Provider', async (ctx) => {
      const kinData = (a.data as any).kinData as Kinemage;
      const visibilityState = (a.data as any).visibilityState as KinemageVisibilityState | undefined;

      if (!kinData) {
        throw new Error('No kinData found in parent Format.Json node');
      }

      const provider = await shapeSpheresFromKin(kinData, visibilityState).runInContext(ctx);
      return new PluginStateObject.Shape.Provider(provider as any, {
        label: kinData.pdbfile || kinData.caption || 'Kinemage Spheres',
        description: kinData.text || '',
      });
    });
  },
});

async function waitForNonEmptyBoundingSphere(plugin: PluginContext, timeoutMs = 2000, pollMs = 50) {
  // Without a canvas there is no scene to wait for (headless).
  if (!plugin.canvas3d) return null;
  const start = Date.now();
  while (Date.now() - start < timeoutMs) {
    try {
      const bs = getPluginBoundingSphere(plugin);
      if (bs && bs.radius > 0) return bs;
    } catch {
      /* ignore */
    }
    await new Promise<void>((r) => setTimeout(r, pollMs));
  }
  return null;
}

/** Helper function to create all shapes for a kinemage via proper transform chain */
export async function createShapesForKinemage(
  plugin: PluginContext,
  update: StateBuilder.Root,
  visControllerSelector: StateObjectSelector<PluginStateObject.Format.Json>,
) {
  const visControllerCell = plugin.state.data.cells.get(visControllerSelector.ref);
  if (!visControllerCell?.obj?.data) return;

  const kinData = (visControllerCell.obj.data as any).kinData as Kinemage;
  if (!kinData) return;

  // Generate all shape types that have data, each as child of the visibility controller
  if (kinData.dotLists.length > 0) {
    update
      .to(visControllerSelector.ref)
      .apply(KinemageShapePointsProvider, {}, { state: { isGhost: true } })
      .apply(ShapeRepresentation3D);
  }
  if (kinData.vectorLists.length > 0) {
    update
      .to(visControllerSelector.ref)
      .apply(KinemageShapeLinesProvider, {}, { state: { isGhost: true } })
      .apply(ShapeRepresentation3D);
  }
  if (kinData.ribbonLists.length > 0) {
    update
      .to(visControllerSelector.ref)
      .apply(KinemageShapeMeshProvider, {}, { state: { isGhost: true } })
      .apply(ShapeRepresentation3D, { doubleSided: true });
  }
  if (kinData.ballLists.length > 0) {
    update
      .to(visControllerSelector.ref)
      .apply(KinemageShapeSpheresProvider, {}, { state: { isGhost: true } })
      .apply(ShapeRepresentation3D);
  }
}

/** Centralized helper to apply kinemage content into plugin state */
export async function applyKinemageToState(plugin: PluginContext, data: string, label?: string) {
  const update = plugin.state.data.build();

  // Create String data node
  const dataNode = update.toRoot().apply(RawData, { data, label: label || 'Kinemage File' });

  // Parse into KinemageObject
  const parsedNode = dataNode.apply(ParseKinemage, { label });

  await update.commit();

  // Get the parsed kinemage object to see how many kinemages it contains
  const parsedCell = plugin.state.data.cells.get(parsedNode.ref);
  const kinemageData = parsedCell?.obj?.data as KinemageData | undefined;

  if (!kinemageData || !kinemageData.kinemages || kinemageData.kinemages.length === 0) {
    console.warn('No kinemages found in parsed data');
    return undefined;
  }

  // Kinemages are traditionally displayed on a black background.
  await setKinemageBackground(plugin);

  // Create a separate visibility controller and shapes for EACH kinemage
  const visControllerSelectors: StateObjectSelector<PluginStateObject.Format.Json>[] = [];

  for (let i = 0; i < kinemageData.kinemages.length; i++) {
    const kinUpdate = plugin.state.data.build();

    // Select this specific kinemage
    const selectedNode = kinUpdate.to(parsedNode.ref).apply(SelectKinemage, { index: i });

    // Add visibility controller for this kinemage
    const visControllerNode = selectedNode.apply(KinemageVisibilityController, {});

    await kinUpdate.commit();

    visControllerSelectors.push(visControllerNode.selector);
  }

  // Now create shapes for all kinemages
  const shapeUpdate = plugin.state.data.build();
  for (const visControllerSelector of visControllerSelectors) {
    await createShapesForKinemage(plugin, shapeUpdate, visControllerSelector);
  }
  await shapeUpdate.commit();

  // Wait for bounding sphere and focus camera
  try {
    const bs = await waitForNonEmptyBoundingSphere(plugin);
    if (bs && bs.radius > 0 && plugin.canvas3d) {
      await PluginCommands.Camera.Focus(plugin, { center: bs.center, radius: bs.radius, durationMs: 250 });
      plugin.canvas3d?.commit();
    }
  } catch (e) {
    console.warn('Failed to apply initial kinemage view snapshot', e);
  }

  return visControllerSelectors[0]; // Return first for backward compatibility
}

/** Programmatic loader: load a single File (a .kin) into the plugin state.
 * Returns the ref to the first visibility controller node.
 */
export async function loadKinemageFile(
  plugin: PluginContext,
  file: File,
): Promise<StateObjectSelector<PluginStateObject.Format.Json> | undefined> {
  const content = await file.text();
  return await applyKinemageToState(plugin, content, file.name);
}

/** DragAndDropHandler handler for `.kin` files */
export const KinemageDragAndDropHandler: PluginDragAndDropEntry = {
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

export const KINFormatProvider = DataFormatProvider({
  name: 'KIN',
  label: 'KIN',
  description: 'Kinemage',
  category: 'Miscellaneous',
  stringExtensions: ['kin', 'KIN'],
  parse: async (plugin, data) => {
    try {
      // data is already a StateObjectRef to the raw data in the tree
      // Build the transform chain from it
      const builder = plugin.state.data.build().to(data).apply(ParseKinemage, {});

      await builder.commit();

      // Get the parsed data to see how many kinemages
      const parsedRef = builder.selector.ref;
      const parsedCell = plugin.state.data.cells.get(parsedRef);
      const kinemageData = parsedCell?.obj?.data as KinemageData | undefined;

      if (!kinemageData || !kinemageData.kinemages || kinemageData.kinemages.length === 0) {
        console.warn('No kinemages found in parsed data');
        return {};
      }

      // Kinemages are traditionally displayed on a black background.
      await setKinemageBackground(plugin);

      // Create visibility controllers for all kinemages
      const visControllers: StateObjectSelector<PluginStateObject.Format.Json>[] = [];

      for (let i = 0; i < kinemageData.kinemages.length; i++) {
        const kinBuilder = plugin.state.data.build();

        const selectedKin = kinBuilder.to(parsedRef).apply(SelectKinemage, { index: i });

        const visController = selectedKin.apply(KinemageVisibilityController, {});

        await kinBuilder.commit();
        visControllers.push(visController.selector);
      }

      // Return all visibility controllers
      return { visControllers };
    } catch (e) {
      console.error('Failed to parse KIN file', e);
      throw e;
    }
  },
  visuals: async (plugin, data) => {
    if (!data?.visControllers || !Array.isArray(data.visControllers)) {
      console.warn('[Kinemage] visuals: no visControllers array provided');
      return;
    }

    // Create shapes for all kinemages
    const shapeBuilder = plugin.state.data.build();
    for (const visController of data.visControllers) {
      await createShapesForKinemage(plugin, shapeBuilder, visController);
    }
    await shapeBuilder.commit();

    // Wait for bounding sphere and focus camera
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
  },
});

// What the right-panel controls (ui.tsx) and the Prism operations (operations.ts) both do to the loaded kinemages.

export type KinemageVisibilityTarget = { type: 'group' | 'subgroup' | 'master'; key: string };

/** One loaded kinemage: its data and the visibility controller node that owns its shapes. */
export interface LoadedKinemage {
  kinData: Kinemage;
  /** The visibility controller's ref (the kinemage's node in the state tree). */
  ref: string;
}

/** The kinemages in the state tree: visibility controllers carry `kinData` and `visibilityState`. */
export function getLoadedKinemages(plugin: PluginContext): LoadedKinemage[] {
  const result: LoadedKinemage[] = [];
  try {
    for (const [ref, entry] of plugin.state.data.cells) {
      const data = entry.obj?.data as { kinData?: Kinemage; visibilityState?: unknown } | undefined;
      if (data?.kinData && data.visibilityState) result.push({ kinData: data.kinData, ref });
    }
  } catch (e) {
    console.warn('Failed to enumerate kinemage nodes', e);
  }
  return result;
}

function getAllDescendants(plugin: PluginContext, nodeRef: string): string[] {
  const result: string[] = [];
  const tree = plugin.state.data.tree;
  const queue = [nodeRef];
  while (queue.length > 0) {
    const current = queue.shift()!;
    const children = tree.children.get(current);
    if (children) {
      for (const childRef of children.values()) {
        result.push(childRef);
        queue.push(childRef);
      }
    }
  }
  return result;
}

/** Drops a kinemage's shapes and builds them again from its visibility controller (the controller's data changed). */
export async function rebuildKinemageShapes(plugin: PluginContext, visControllerRef: string, kinData: Kinemage) {
  const update = plugin.state.data.build();
  for (const nodeRef of getAllDescendants(plugin, visControllerRef)) update.delete(nodeRef);
  await update.commit();

  const rebuildUpdate = plugin.state.data.build();
  if (kinData.dotLists.length > 0) {
    rebuildUpdate
      .to(visControllerRef)
      .apply(KinemageShapePointsProvider, {}, { state: { isGhost: true } })
      .apply(ShapeRepresentation3D);
  }
  if (kinData.vectorLists.length > 0) {
    rebuildUpdate
      .to(visControllerRef)
      .apply(KinemageShapeLinesProvider, {}, { state: { isGhost: true } })
      .apply(ShapeRepresentation3D);
  }
  if (kinData.ribbonLists.length > 0) {
    rebuildUpdate
      .to(visControllerRef)
      .apply(KinemageShapeMeshProvider, {}, { state: { isGhost: true } })
      .apply(ShapeRepresentation3D, { doubleSided: true });
  }
  if (kinData.ballLists.length > 0) {
    rebuildUpdate
      .to(visControllerRef)
      .apply(KinemageShapeSpheresProvider, {}, { state: { isGhost: true } })
      .apply(ShapeRepresentation3D);
  }
  await rebuildUpdate.commit();
}

/** Whether a group, subgroup or master is shown (a missing one counts as shown, as in the controls). */
export function isKinemageVisible(plugin: PluginContext, visControllerRef: string, target: KinemageVisibilityTarget) {
  try {
    const params = plugin.state.data.cells.get(visControllerRef)?.transform?.params;
    if (!params) return true;
    const map =
      target.type === 'group'
        ? params.groupVisibility
        : target.type === 'subgroup'
          ? params.subgroupVisibility
          : params.masterVisibility;
    return map[target.key] !== false;
  } catch {
    return true;
  }
}

/** Sets a group, subgroup or master shown or hidden, or flips it (`visible` absent), and rebuilds the shapes. */
export async function setKinemageVisibility(
  plugin: PluginContext,
  visControllerRef: string,
  kinData: Kinemage,
  target: KinemageVisibilityTarget,
  visible?: boolean,
) {
  const cell = plugin.state.data.cells.get(visControllerRef);
  if (!cell?.transform?.params) return;
  const p = cell.transform.params;
  const groupVisibility = { ...p.groupVisibility };
  const subgroupVisibility = { ...p.subgroupVisibility };
  const masterVisibility = { ...p.masterVisibility };
  const map =
    target.type === 'group' ? groupVisibility : target.type === 'subgroup' ? subgroupVisibility : masterVisibility;
  map[target.key] = visible ?? !map[target.key];

  const update = plugin.state.data.build();
  update.to(visControllerRef).update({ groupVisibility, subgroupVisibility, masterVisibility });
  await update.commit();
  await rebuildKinemageShapes(plugin, visControllerRef, kinData);
}

/** Shows the next frame of the kinemage's animate (or 2animate) groups and hides the others; returns the frame's group. */
export async function animateKinemage(
  plugin: PluginContext,
  visControllerRef: string,
  kinData: Kinemage,
  mode: 'animate' | '2animate',
): Promise<string | undefined> {
  const cell = plugin.state.data.cells.get(visControllerRef);
  if (!cell?.transform?.params) return undefined;
  const currentParams = cell.transform.params;
  const animateGroups = mode === 'animate' ? kinData.groupsAnimate : kinData.groupsAnimate2;
  const currentActive = mode === 'animate' ? currentParams.activeAnimateGroup : currentParams.activeAnimateGroup2;
  const nextActive = (currentActive + 1) % Math.max(1, animateGroups.length);

  // The controller's current data, not its params: it keeps what interactions changed.
  const state = (cell.obj?.data as { visibilityState?: KinemageVisibilityState } | undefined)?.visibilityState;
  const groupVisibility = state ? Object.fromEntries(state.groupVisibility) : { ...currentParams.groupVisibility };
  for (let i = 0; i < animateGroups.length; i++) groupVisibility[animateGroups[i]] = i === nextActive;

  const updateParams: Record<string, unknown> = {
    groupVisibility,
    subgroupVisibility: state ? Object.fromEntries(state.subgroupVisibility) : currentParams.subgroupVisibility,
    masterVisibility: state ? Object.fromEntries(state.masterVisibility) : currentParams.masterVisibility,
  };
  if (mode === 'animate') updateParams.activeAnimateGroup = nextActive;
  else updateParams.activeAnimateGroup2 = nextActive;

  const update = plugin.state.data.build();
  update.to(visControllerRef).update(updateParams);
  await update.commit();
  await rebuildKinemageShapes(plugin, visControllerRef, kinData);
  return animateGroups[nextActive];
}

/** Set the viewer background to black, keeping the plugin's settings UI in sync. */
export async function setKinemageBackground(plugin: PluginContext) {
  if (!plugin.canvas3d) return;
  await PluginCommands.Canvas3D.SetSettings(plugin, {
    settings: (props) => {
      props.renderer.backgroundColor = Color(0x000000);
    },
  });
}
