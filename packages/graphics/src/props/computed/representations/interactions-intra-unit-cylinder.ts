/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { chooseVisual } from '@molstar/graphics/repr/visual-factory';
import { Unit, Structure, StructureElement } from '@molstar/model/model/structure';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { type Loci, EmptyLoci } from '@molstar/model/model/loci';
import { Interval, OrderedSet, SortedArray } from '@molstar/core/data/int';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import type { VisualContext } from '@molstar/graphics/repr/visual';
import type { Theme } from '@molstar/graphics/theme/theme';
import { InteractionsProvider } from '@molstar/model/props/computed/interactions';
import {
  createLinkCylinderImpostors,
  createLinkCylinderMesh,
  type LinkBuilderProps,
  LinkCylinderParams,
  LinkStyle,
} from '@molstar/graphics/repr/structure/visual/util/link';
import {
  UnitsCylindersParams,
  UnitsCylindersVisual,
  UnitsMeshParams,
  type UnitsVisual,
  UnitsMeshVisual,
} from '@molstar/graphics/repr/structure/units-visual';
import { Cylinders } from '@molstar/graphics/geo/geometry/cylinders/cylinders';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { VisualUpdateState } from '@molstar/graphics/repr/util';
import { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import { Interactions } from '@molstar/model/props/computed/interactions/interactions';
import { FeatureType, InteractionFlag } from '@molstar/model/props/computed/interactions/common';
import { Sphere3D } from '@molstar/core/math/geometry';
import {
  type StructureGroup,
  checkCylinderImpostorSupport,
  isHydrogen,
} from '@molstar/graphics/repr/structure/visual/util/common';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import { impostorProps, InteractionsSharedParams, withDashedCylinderCount } from './shared.js';
import { InteractionType } from '@molstar/model/props/computed/interactions/common';
import { eachBondedAtom } from '@molstar/model/props/computed/chemistry/util';

/** The interactions of a unit as dashed links, for the impostor and the mesh builder; undefined if there are none */
function getIntraUnitInteractionsBuilderProps(
  unit: Unit,
  structure: Structure,
  theme: Theme,
  props: PD.Values<InteractionsIntraUnitParams>,
): LinkBuilderProps | undefined {
  if (!Unit.isAtomic(unit)) return;

  const { child } = structure;
  const childUnit = child?.unitMap.get(unit.id);
  if (child && !childUnit) return;

  const location = StructureElement.Location.create(structure, unit);

  const interactions = InteractionsProvider.get(structure).value!;
  const features = interactions.unitsFeatures.get(unit.id);
  const contacts = interactions.unitsContacts.get(unit.id);

  const { x, y, z, members, offsets, types } = features;
  const {
    edgeCount,
    a,
    b,
    edgeProps: { flag, type },
  } = contacts;
  const { ignoreHydrogens, ignoreHydrogensVariant, parentDisplay } = props;
  // as the impostors' uSizeFactor
  const sizeFactor = props.sizeFactor * props.sizeAspectRatio;

  if (!edgeCount) return;

  const { elements, conformation: c } = unit;
  const p = Vec3();
  /** A bonded hydrogen in this unit's invariant frame; with per-atom symmetry it can be in another unit */
  const hydrogenPosition = (u: Unit.Atomic, index: StructureElement.UnitIndex, out: Vec3) => {
    if (u === unit) return c.invariantPosition(u.elements[index], out);
    u.conformation.position(u.elements[index], out);
    return Vec3.transformMat4(out, out, c.operator.inverse);
  };
  const pA = Vec3();
  const pB = Vec3();

  return {
    linkCount: edgeCount * 2,
    position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
      const t = type[edgeIndex];
      if (
        (!ignoreHydrogens || ignoreHydrogensVariant !== 'all') &&
        (t === InteractionType.HydrogenBond ||
          (t === InteractionType.WeakHydrogenBond && ignoreHydrogensVariant !== 'non-polar'))
      ) {
        const idxA = members[offsets[a[edgeIndex]]];
        const idxB = members[offsets[b[edgeIndex]]];
        c.invariantPosition(elements[idxA], pA);
        c.invariantPosition(elements[idxB], pB);
        let minDistA = Vec3.distance(pA, pB);
        let minDistB = minDistA;
        Vec3.copy(posA, pA);
        Vec3.copy(posB, pB);
        const donorType =
          t === InteractionType.HydrogenBond ? FeatureType.HydrogenDonor : FeatureType.WeakHydrogenDonor;
        const isHydrogenDonorA = types[a[edgeIndex]] === donorType;

        if (isHydrogenDonorA) {
          eachBondedAtom(structure, unit, idxA, (u, idx) => {
            if (isHydrogen(structure, u, u.elements[idx], 'all')) {
              hydrogenPosition(u, idx, p);
              const dist = Vec3.distance(p, pB);
              if (dist < minDistA) {
                minDistA = dist;
                Vec3.copy(posA, p);
              }
            }
          });
        } else {
          eachBondedAtom(structure, unit, idxB, (u, idx) => {
            if (isHydrogen(structure, u, u.elements[idx], 'all')) {
              hydrogenPosition(u, idx, p);
              const dist = Vec3.distance(p, pA);
              if (dist < minDistB) {
                minDistB = dist;
                Vec3.copy(posB, p);
              }
            }
          });
        }
      } else {
        Vec3.set(posA, x[a[edgeIndex]], y[a[edgeIndex]], z[a[edgeIndex]]);
        Vec3.set(posB, x[b[edgeIndex]], y[b[edgeIndex]], z[b[edgeIndex]]);
      }
    },
    style: (edgeIndex: number) => LinkStyle.Dashed,
    radius: (edgeIndex: number) => {
      location.element = elements[members[offsets[a[edgeIndex]]]];
      const sizeA = theme.size.size(location);
      location.element = elements[members[offsets[b[edgeIndex]]]];
      const sizeB = theme.size.size(location);
      return Math.min(sizeA, sizeB) * sizeFactor;
    },
    ignore: (edgeIndex: number) => {
      if (flag[edgeIndex] === InteractionFlag.Filtered) return true;

      if (childUnit) {
        if (parentDisplay === 'stub') {
          const f = a[edgeIndex];
          for (let i = offsets[f], il = offsets[f + 1]; i < il; ++i) {
            const e = elements[members[i]];
            if (!SortedArray.has(childUnit.elements, e)) return true;
          }
        } else if (parentDisplay === 'full' || parentDisplay === 'between') {
          let flagA = false;
          let flagB = false;

          const fA = a[edgeIndex];
          for (let i = offsets[fA], il = offsets[fA + 1]; i < il; ++i) {
            const eA = elements[members[i]];
            if (!SortedArray.has(childUnit.elements, eA)) flagA = true;
          }

          const fB = b[edgeIndex];
          for (let i = offsets[fB], il = offsets[fB + 1]; i < il; ++i) {
            const eB = elements[members[i]];
            if (!SortedArray.has(childUnit.elements, eB)) flagB = true;
          }

          return parentDisplay === 'full' ? flagA && flagB : flagA === flagB;
        } else {
          assertUnreachable(parentDisplay);
        }
      }

      return false;
    },
  };
}

function setBoundingSphere(
  geometry: Mesh | Cylinders,
  boundingSphere: Sphere3D | undefined,
  unit: Unit,
  structure: Structure,
  sizeFactor: number,
) {
  if (boundingSphere) {
    geometry.setBoundingSphere(boundingSphere);
  } else if (geometry.kind === 'mesh' ? geometry.triangleCount > 0 : geometry.cylinderCount > 0) {
    const childUnit = structure.child?.unitMap.get(unit.id);
    const sphere = Sphere3D.expand(Sphere3D(), (childUnit ?? unit).boundary.sphere, 1 * sizeFactor);
    geometry.setBoundingSphere(sphere);
  }
}

async function createIntraUnitInteractionsCylinderImpostors(
  ctx: VisualContext,
  unit: Unit,
  structure: Structure,
  theme: Theme,
  props: PD.Values<InteractionsIntraUnitParams>,
  cylinders?: Cylinders,
) {
  const builderProps = getIntraUnitInteractionsBuilderProps(unit, structure, theme, props);
  if (!builderProps) return Cylinders.createEmpty(cylinders);

  const { cylinders: c, boundingSphere } = createLinkCylinderImpostors(
    ctx,
    withDashedCylinderCount(builderProps, props.dashCount),
    impostorProps(props),
    cylinders,
  );
  setBoundingSphere(c, boundingSphere, unit, structure, props.sizeFactor);
  return c;
}

async function createIntraUnitInteractionsCylinderMesh(
  ctx: VisualContext,
  unit: Unit,
  structure: Structure,
  theme: Theme,
  props: PD.Values<InteractionsIntraUnitParams>,
  mesh?: Mesh,
) {
  const builderProps = getIntraUnitInteractionsBuilderProps(unit, structure, theme, props);
  if (!builderProps) return Mesh.createEmpty(mesh);

  const { mesh: m, boundingSphere } = createLinkCylinderMesh(ctx, builderProps, props, mesh);
  setBoundingSphere(m, boundingSphere, unit, structure, props.sizeFactor);
  return m;
}

export const InteractionsIntraUnitParams = {
  ...UnitsMeshParams,
  ...UnitsCylindersParams,
  ...LinkCylinderParams,
  ...InteractionsSharedParams,
};
export type InteractionsIntraUnitParams = typeof InteractionsIntraUnitParams;

export const InteractionsIntraUnitVisual = chooseVisual(
  { impostor: InteractionsIntraUnitImpostorVisual, mesh: InteractionsIntraUnitMeshVisual },
  (props: PD.Values<InteractionsIntraUnitParams>, graphics?: GraphicsContext) =>
    props.tryUseImpostor && checkCylinderImpostorSupport(graphics) ? 'impostor' : 'mesh',
);

function setUpdateState(
  state: VisualUpdateState,
  newProps: PD.Values<InteractionsIntraUnitParams>,
  currentProps: PD.Values<InteractionsIntraUnitParams>,
  newTheme: Theme,
  currentTheme: Theme,
  newStructureGroup: StructureGroup,
  currentStructureGroup: StructureGroup,
) {
  state.createGeometry =
    newProps.dashCount !== currentProps.dashCount ||
    newProps.dashScale !== currentProps.dashScale ||
    newProps.dashCap !== currentProps.dashCap ||
    newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
    newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
    newProps.parentDisplay !== currentProps.parentDisplay;

  const interactionsHash = InteractionsProvider.get(newStructureGroup.structure).version;
  if ((state.info.interactionsHash as number) !== interactionsHash) {
    state.createGeometry = true;
    state.updateTransform = true;
    state.updateColor = true;
    // group sizes from a group-granularity size theme
    state.updateSize = true;
    state.info.interactionsHash = interactionsHash;
  }
}

/** The interactions version the new geometry is built from: the first update after a build does not rebuild it */
function initUpdateState(
  state: VisualUpdateState,
  newProps: PD.Values<InteractionsIntraUnitParams>,
  newTheme: Theme,
  newStructureGroup: StructureGroup,
) {
  state.info.interactionsHash = InteractionsProvider.get(newStructureGroup.structure).version;
}

/**
 * The mesh bakes the radius (`sizeFactor`, `sizeAspectRatio`) and the radial segments; impostors take them from
 * uniforms and the size theme, so a size slider rebuilds no impostor geometry (perf-audit-2 09/P16)
 */
function setMeshUpdateState(...args: Parameters<typeof setUpdateState>) {
  setUpdateState(...args);
  const [state, newProps, currentProps] = args;
  if (
    newProps.sizeFactor !== currentProps.sizeFactor ||
    newProps.radialSegments !== currentProps.radialSegments ||
    newProps.sizeAspectRatio !== currentProps.sizeAspectRatio
  ) {
    state.createGeometry = true;
  }
}

export function InteractionsIntraUnitImpostorVisual(materialId: number): UnitsVisual<InteractionsIntraUnitParams> {
  return UnitsCylindersVisual<InteractionsIntraUnitParams>(
    {
      defaultProps: PD.getDefaultValues(InteractionsIntraUnitParams),
      createGeometry: createIntraUnitInteractionsCylinderImpostors,
      createLocationIterator: createInteractionsIterator,
      getLoci: getInteractionLoci,
      eachLocation: eachInteraction,
      setUpdateState,
      initUpdateState,
      mustRecreate: (_: StructureGroup, props: PD.Values<InteractionsIntraUnitParams>, graphics?: GraphicsContext) =>
        !props.tryUseImpostor || !checkCylinderImpostorSupport(graphics),
    },
    materialId,
  );
}

export function InteractionsIntraUnitMeshVisual(materialId: number): UnitsVisual<InteractionsIntraUnitParams> {
  return UnitsMeshVisual<InteractionsIntraUnitParams>(
    {
      defaultProps: PD.getDefaultValues(InteractionsIntraUnitParams),
      createGeometry: createIntraUnitInteractionsCylinderMesh,
      createLocationIterator: createInteractionsIterator,
      getLoci: getInteractionLoci,
      eachLocation: eachInteraction,
      setUpdateState: setMeshUpdateState,
      initUpdateState,
      mustRecreate: (_: StructureGroup, props: PD.Values<InteractionsIntraUnitParams>, graphics?: GraphicsContext) =>
        props.tryUseImpostor && checkCylinderImpostorSupport(graphics),
    },
    materialId,
  );
}

function getInteractionLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
  const { objectId, instanceId, groupId } = pickingId;
  if (id === objectId) {
    const { structure, group } = structureGroup;
    const unit = structure.unitMap.get(group.units[instanceId].id);
    const interactions = InteractionsProvider.get(structure).value!;
    const { a, b } = interactions.unitsContacts.get(unit.id);
    return Interactions.Loci(structure, interactions, [
      { unitA: unit, indexA: a[groupId], unitB: unit, indexB: b[groupId] },
      { unitA: unit, indexA: b[groupId], unitB: unit, indexB: a[groupId] },
    ]);
  }
  return EmptyLoci;
}

const __contactIndicesSet = new Set<number>();

function eachInteraction(
  loci: Loci,
  structureGroup: StructureGroup,
  apply: (interval: Interval) => boolean,
  isMarking: boolean,
) {
  let changed = false;
  if (Interactions.isLoci(loci)) {
    const { structure, group } = structureGroup;
    if (!Structure.areEquivalent(loci.data.structure, structure)) return false;
    const interactions = InteractionsProvider.get(structure).value!;
    if (loci.data.interactions !== interactions) return false;
    const unit = group.units[0];
    const contacts = interactions.unitsContacts.get(unit.id);
    const groupCount = contacts.edgeCount * 2;
    for (const e of loci.elements) {
      if (e.unitA !== e.unitB) continue;
      const unitIdx = group.unitIndexMap.get(e.unitA.id);
      if (unitIdx !== undefined) {
        const idx = contacts.getDirectedEdgeIndex(e.indexA, e.indexB);
        if (idx !== -1) {
          if (apply(Interval.ofSingleton(unitIdx * groupCount + idx))) changed = true;
        }
      }
    }
  } else if (StructureElement.Loci.is(loci)) {
    const { structure, group } = structureGroup;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;

    const interactions = InteractionsProvider.get(structure).value;
    if (!interactions) return false;
    const unit = group.units[0];
    const contacts = interactions.unitsContacts.get(unit.id);
    const features = interactions.unitsFeatures.get(unit.id);
    const groupCount = contacts.edgeCount * 2;

    const { offset } = contacts;
    const { offsets: fOffsets, indices: fIndices } = features.elementsIndex;
    const { members, offsets } = features;

    for (const e of loci.elements) {
      const unitIdx = group.unitIndexMap.get(e.unit.id);
      if (unitIdx === undefined) continue;

      OrderedSet.forEach(e.indices, (v) => {
        for (let i = fOffsets[v], il = fOffsets[v + 1]; i < il; ++i) {
          const fI = fIndices[i];
          for (let j = offset[fI], jl = offset[fI + 1]; j < jl; ++j) {
            __contactIndicesSet.add(j);
          }
        }
      });

      __contactIndicesSet.forEach((i) => {
        if (isMarking) {
          const fA = contacts.a[i];
          for (let j = offsets[fA], jl = offsets[fA + 1]; j < jl; ++j) {
            if (!OrderedSet.has(e.indices, members[j])) return;
          }
          const fB = contacts.b[i];
          for (let j = offsets[fB], jl = offsets[fB + 1]; j < jl; ++j) {
            if (!OrderedSet.has(e.indices, members[j])) return;
          }
        }

        if (apply(Interval.ofSingleton(unitIdx * groupCount + i))) changed = true;
      });

      __contactIndicesSet.clear();
    }
  }
  return changed;
}

function createInteractionsIterator(structureGroup: StructureGroup): LocationIterator {
  const { structure, group } = structureGroup;
  const unit = group.units[0];
  const interactions = InteractionsProvider.get(structure).value!;
  const contacts = interactions.unitsContacts.get(unit.id);
  const groupCount = contacts.edgeCount * 2;
  const instanceCount = group.units.length;
  const location = Interactions.Location(interactions, structure);
  const { element } = location;
  const getLocation = (groupIndex: number, instanceIndex: number) => {
    const instanceUnit = group.units[instanceIndex];
    element.unitA = instanceUnit;
    element.indexA = contacts.a[groupIndex];
    element.unitB = instanceUnit;
    element.indexB = contacts.b[groupIndex];
    return location;
  };
  return LocationIterator(groupCount, instanceCount, 1, getLocation);
}
