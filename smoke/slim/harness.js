// Smoke harness for the slim plugin example (`node smoke/run.mjs slim`). It is bundled together with a copy of the
// example's entry (`examples/slim-plugin/src`), so the registered providers are exactly those of the example's spec.
// The assertions live in `smoke/run.mjs`; this module only exposes what they need.
import './slim-example/index.ts';
import { Script } from '@molstar/model/script/script';

const StructureRepresentation3D = 'ms-plugin.structure-representation-3d';
const warnings = [];

function structureReprs(plugin) {
  return [...plugin.state.data.cells.values()]
    .filter((cell) => cell.transform.transformer.id === StructureRepresentation3D)
    .map((cell) => {
      const repr = cell.obj?.data?.repr;
      return {
        ref: cell.transform.ref,
        type: cell.transform.params?.type?.name,
        status: cell.status,
        visible: repr ? repr.state.visible : false,
        renderObjects: repr ? repr.renderObjects.length : 0,
        renderObjectsWithGeometry: repr ? repr.renderObjects.filter((o) => o.values.drawCount.ref.value > 0).length : 0,
      };
    });
}

function describe(plugin) {
  return {
    registered: {
      structureRepresentations: plugin.representation.structure.registry.list.map((e) => e.name),
      formats: plugin.dataFormats.list.map((e) => e.name),
      hierarchyPresets: plugin.builders.structure.hierarchy.providers.map((p) => p.id),
      representationPresets: plugin.builders.structure.representation.providers.map((p) => p.id),
      scriptLanguages: Script.getAvailableLanguages(),
    },
    structures: plugin.managers.structure.hierarchy.current.structures.length,
    structureReprs: structureReprs(plugin),
    sceneRenderObjects: plugin.canvas3d?.getRenderObjects().length ?? 0,
    reprCount: plugin.canvas3d?.reprCount.value ?? 0,
  };
}

window.slimSmoke = {
  describe: () => describe(window.slimPlugin),
  /** Warnings logged while restoring the snapshot. */
  warnings,
  /** Replaces the type of every structure representation of the current state by `cartoon` in a snapshot and restores it. */
  async restoreCartoonSnapshot() {
    const plugin = window.slimPlugin;
    const snapshot = JSON.parse(JSON.stringify(plugin.state.getSnapshot({ data: true, behavior: true })));
    let changed = 0;
    for (const t of snapshot.data.tree.transforms) {
      if (t.transformer !== StructureRepresentation3D) continue;
      t.params = { ...t.params, type: { name: 'cartoon', params: {} } };
      changed++;
    }
    const sub = plugin.events.log.subscribe((e) => {
      if (e.type === 'warning') warnings.push(e.message);
    });
    try {
      await plugin.state.setSnapshot(snapshot);
    } finally {
      sub.unsubscribe();
    }
    plugin.canvas3d?.commit(true);
    return { changed };
  },
  /** Evaluates a PyMOL script. Returns the error message, or undefined when the evaluation did not fail. */
  evaluatePyMol() {
    try {
      Script.toExpression(Script('resn ALA', 'pymol'));
    } catch (e) {
      return e instanceof Error ? e.message : String(e);
    }
    return undefined;
  },
};

window.SlimPlugin.init('app');
