/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO, type PluginStateObject } from '@molstar/plugin/state/objects';
import { DistanceParams, DistanceRepresentation, type DistanceData } from '@molstar/graphics/repr/shape/loci/distance';
import type { PluginContext } from '@molstar/plugin/context';
import { Task } from '@molstar/core/task';
import { StateTransformer } from '@molstar/core/state';
import { AngleParams, AngleRepresentation, type AngleData } from '@molstar/graphics/repr/shape/loci/angle';
import { DihedralParams, DihedralRepresentation, type DihedralData } from '@molstar/graphics/repr/shape/loci/dihedral';
import { LabelParams, LabelRepresentation, type LabelData } from '@molstar/graphics/repr/shape/loci/label';
import { MarkerActions, MarkerAction } from '@molstar/core/util/marker-action';
import {
  OrientationParams,
  OrientationRepresentation,
  type OrientationData,
} from '@molstar/graphics/repr/shape/loci/orientation';
import { PlaneParams, PlaneRepresentation, type PlaneData } from '@molstar/graphics/repr/shape/loci/plane';

export { StructureSelectionsDistance3D };
type StructureSelectionsDistance3D = typeof StructureSelectionsDistance3D;
const StructureSelectionsDistance3D = PluginStateTransform.BuiltIn({
  name: 'structure-selections-distance-3d',
  display: '3D Distance',
  from: SO.Molecule.Structure.Selections,
  to: SO.Shape.Representation3D,
  params: () => ({
    ...DistanceParams,
  }),
})({
  canAutoUpdate({ oldParams, newParams }) {
    return true;
  },
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Structure Distance', async (ctx) => {
      const data = getDistanceDataFromStructureSelections(a.data);
      const repr = DistanceRepresentation(
        { webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes },
        () => DistanceParams,
      );
      await repr.createOrUpdate(params, data).runInContext(ctx);
      return new SO.Shape.Representation3D({ repr, sourceData: data }, { label: `Distance` });
    });
  },
  update({ a, b, oldParams, newParams }, plugin: PluginContext) {
    return Task.create('Structure Distance', async (ctx) => {
      const props = { ...b.data.repr.props, ...newParams };
      const data = getDistanceDataFromStructureSelections(a.data);
      await b.data.repr.createOrUpdate(props, data).runInContext(ctx);
      b.data.sourceData = data;
      return StateTransformer.UpdateResult.Updated;
    });
  },
});

export { StructureSelectionsAngle3D };
type StructureSelectionsAngle3D = typeof StructureSelectionsAngle3D;
const StructureSelectionsAngle3D = PluginStateTransform.BuiltIn({
  name: 'structure-selections-angle-3d',
  display: '3D Angle',
  from: SO.Molecule.Structure.Selections,
  to: SO.Shape.Representation3D,
  params: () => ({
    ...AngleParams,
  }),
})({
  canAutoUpdate({ oldParams, newParams }) {
    return true;
  },
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Structure Angle', async (ctx) => {
      const data = getAngleDataFromStructureSelections(a.data);
      const repr = AngleRepresentation(
        { webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes },
        () => AngleParams,
      );
      await repr.createOrUpdate(params, data).runInContext(ctx);
      return new SO.Shape.Representation3D({ repr, sourceData: data }, { label: `Angle` });
    });
  },
  update({ a, b, oldParams, newParams }, plugin: PluginContext) {
    return Task.create('Structure Angle', async (ctx) => {
      const props = { ...b.data.repr.props, ...newParams };
      const data = getAngleDataFromStructureSelections(a.data);
      await b.data.repr.createOrUpdate(props, data).runInContext(ctx);
      b.data.sourceData = data;
      return StateTransformer.UpdateResult.Updated;
    });
  },
});

export { StructureSelectionsDihedral3D };
type StructureSelectionsDihedral3D = typeof StructureSelectionsDihedral3D;
const StructureSelectionsDihedral3D = PluginStateTransform.BuiltIn({
  name: 'structure-selections-dihedral-3d',
  display: '3D Dihedral',
  from: SO.Molecule.Structure.Selections,
  to: SO.Shape.Representation3D,
  params: () => ({
    ...DihedralParams,
  }),
})({
  canAutoUpdate({ oldParams, newParams }) {
    return true;
  },
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Structure Dihedral', async (ctx) => {
      const data = getDihedralDataFromStructureSelections(a.data);
      const repr = DihedralRepresentation(
        { webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes },
        () => DihedralParams,
      );
      await repr.createOrUpdate(params, data).runInContext(ctx);
      return new SO.Shape.Representation3D({ repr, sourceData: data }, { label: `Dihedral` });
    });
  },
  update({ a, b, oldParams, newParams }, plugin: PluginContext) {
    return Task.create('Structure Dihedral', async (ctx) => {
      const props = { ...b.data.repr.props, ...newParams };
      const data = getDihedralDataFromStructureSelections(a.data);
      await b.data.repr.createOrUpdate(props, data).runInContext(ctx);
      b.data.sourceData = data;
      return StateTransformer.UpdateResult.Updated;
    });
  },
});

export { StructureSelectionsLabel3D };
type StructureSelectionsLabel3D = typeof StructureSelectionsLabel3D;
const StructureSelectionsLabel3D = PluginStateTransform.BuiltIn({
  name: 'structure-selections-label-3d',
  display: '3D Label',
  from: SO.Molecule.Structure.Selections,
  to: SO.Shape.Representation3D,
  params: () => ({
    ...LabelParams,
  }),
})({
  canAutoUpdate({ oldParams, newParams }) {
    return true;
  },
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Structure Label', async (ctx) => {
      const data = getLabelDataFromStructureSelections(a.data);
      const repr = LabelRepresentation(
        { webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes },
        () => LabelParams,
      );
      await repr.createOrUpdate(params, data).runInContext(ctx);

      // Support interactivity when needed
      const pickable = !!(params.snapshotKey?.trim() || params.tooltip?.trim());
      repr.setState({ pickable, markerActions: pickable ? MarkerActions.Highlighting : MarkerAction.None });

      return new SO.Shape.Representation3D({ repr, sourceData: data }, { label: `Label` });
    });
  },
  update({ a, b, oldParams, newParams }, plugin: PluginContext) {
    return Task.create('Structure Label', async (ctx) => {
      const props = { ...b.data.repr.props, ...newParams };
      const data = getLabelDataFromStructureSelections(a.data);
      await b.data.repr.createOrUpdate(props, data).runInContext(ctx);
      b.data.sourceData = data;

      // Update interactivity
      const pickable = !!(newParams.snapshotKey?.trim() || newParams.tooltip?.trim());
      b.data.repr.setState({ pickable, markerActions: pickable ? MarkerActions.Highlighting : MarkerAction.None });

      return StateTransformer.UpdateResult.Updated;
    });
  },
});

export { StructureSelectionsOrientation3D };
type StructureSelectionsOrientation3D = typeof StructureSelectionsOrientation3D;
const StructureSelectionsOrientation3D = PluginStateTransform.BuiltIn({
  name: 'structure-selections-orientation-3d',
  display: '3D Orientation',
  from: SO.Molecule.Structure.Selections,
  to: SO.Shape.Representation3D,
  params: () => ({
    ...OrientationParams,
  }),
})({
  canAutoUpdate({ oldParams, newParams }) {
    return true;
  },
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Structure Orientation', async (ctx) => {
      const data = getOrientationDataFromStructureSelections(a.data);
      const repr = OrientationRepresentation(
        { webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes },
        () => OrientationParams,
      );
      await repr.createOrUpdate(params, data).runInContext(ctx);
      return new SO.Shape.Representation3D({ repr, sourceData: data }, { label: `Orientation` });
    });
  },
  update({ a, b, oldParams, newParams }, plugin: PluginContext) {
    return Task.create('Structure Orientation', async (ctx) => {
      const props = { ...b.data.repr.props, ...newParams };
      const data = getOrientationDataFromStructureSelections(a.data);
      await b.data.repr.createOrUpdate(props, data).runInContext(ctx);
      b.data.sourceData = data;
      return StateTransformer.UpdateResult.Updated;
    });
  },
});

export { StructureSelectionsPlane3D };
type StructureSelectionsPlane3D = typeof StructureSelectionsPlane3D;
const StructureSelectionsPlane3D = PluginStateTransform.BuiltIn({
  name: 'structure-selections-plane-3d',
  display: '3D Plane',
  from: SO.Molecule.Structure.Selections,
  to: SO.Shape.Representation3D,
  params: () => ({
    ...PlaneParams,
  }),
})({
  canAutoUpdate({ oldParams, newParams }) {
    return true;
  },
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Structure Plane', async (ctx) => {
      const data = getPlaneDataFromStructureSelections(a.data);
      const repr = PlaneRepresentation(
        { webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes },
        () => PlaneParams,
      );
      await repr.createOrUpdate(params, data).runInContext(ctx);
      return new SO.Shape.Representation3D({ repr, sourceData: data }, { label: `Plane` });
    });
  },
  update({ a, b, oldParams, newParams }, plugin: PluginContext) {
    return Task.create('Structure Plane', async (ctx) => {
      const props = { ...b.data.repr.props, ...newParams };
      const data = getPlaneDataFromStructureSelections(a.data);
      await b.data.repr.createOrUpdate(props, data).runInContext(ctx);
      b.data.sourceData = data;
      return StateTransformer.UpdateResult.Updated;
    });
  },
});

export function getDistanceDataFromStructureSelections(
  s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>,
): DistanceData {
  const lociA = s[0].loci;
  const lociB = s[1].loci;
  return { pairs: [{ loci: [lociA, lociB] as const }] };
}

export function getAngleDataFromStructureSelections(
  s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>,
): AngleData {
  const lociA = s[0].loci;
  const lociB = s[1].loci;
  const lociC = s[2].loci;
  return { triples: [{ loci: [lociA, lociB, lociC] as const }] };
}

export function getDihedralDataFromStructureSelections(
  s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>,
): DihedralData {
  const lociA = s[0].loci;
  const lociB = s[1].loci;
  const lociC = s[2].loci;
  const lociD = s[3].loci;
  return { quads: [{ loci: [lociA, lociB, lociC, lociD] as const }] };
}

export function getLabelDataFromStructureSelections(
  s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>,
): LabelData {
  const loci = s[0].loci;
  return { infos: [{ loci }] };
}

export function getOrientationDataFromStructureSelections(
  s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>,
): OrientationData {
  return { locis: s.map((v) => v.loci) };
}

export function getPlaneDataFromStructureSelections(
  s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>,
): PlaneData {
  return { locis: s.map((v) => v.loci) };
}
