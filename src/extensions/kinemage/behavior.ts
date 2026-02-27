/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russ Taylor <russ@reliasolve.com>
 */

/** Based on the ../anvil extension. */

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

const Tag = KinemageData.Tag;

const Transform = StateTransformer.builderFactory('sb-kinemage');

export const KinemageShapePointsProvider = Transform({
  name: 'sb-kinemage-shape-points-provider',
  display: { name: 'Kinemage Shape Points Provider' },
  from: PluginStateObject.Root,
  to: PluginStateObject.Shape.Provider,
  params: {
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
  for (const kinData of kinInfo.kinemages) {
    if (kinData.dotLists.length > 0) {
      await update
        .toRoot()
        .apply(KinemageShapePointsProvider, { data: kinData })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);
    }
    if (kinData.vectorLists.length > 0) {
      await update
        .toRoot()
        .apply(KinemageShapeLinesProvider, { data: kinData })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);
    }
    if (kinData.ribbonLists.length > 0) {
      await update
        .toRoot()
        .apply(KinemageShapeMeshProvider, { data: kinData })
        .apply(StateTransforms.Representation.ShapeRepresentation3D, { doubleSided: true });
    }
    if (kinData.ballLists.length > 0) {
      await update
        .toRoot()
        .apply(KinemageShapeSpheresProvider, { data: kinData })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);
    }
  }
  update.commit();
}

/** Programmatic loader: load a single File (a .kin) into the plugin state.
 * Runs the import inside a Task so it has a runtime and asset context similar to drag-and-drop.
 * Returns true if at least one Kinemage was added.
 */
export async function loadKinemageFile(plugin: PluginContext, file: File): Promise<boolean> {
  let applied = false;
  const task = Task.create('Load KIN file', async ctx => {
    const kinInfo = await KinemageData.open(file);
    await applyKinemageInfoToState(plugin, kinInfo);
    applied = kinInfo.kinemages.length > 0;
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
    /// @todo A more standard approach would be to have the parse call generate data and the visuals call create
    /// the representations. However, since the KIN loader is already implemented as a side-effecting function that
    /// applies directly to the plugin state, we can just call it from parse and have visuals be a no-op.  If we wanted
    // to split it up more cleanly, we would need to refactor the KIN loader to separate parsing from applying to state.
    //await (KINFormatProvider.parse as any)(plugin, data);
    return undefined;
  }
});
