/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { chooseVisual } from '@molstar/graphics/repr/visual-factory';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import type { VisualContext } from '@molstar/graphics/repr/visual';
import { Structure, StructureElement, type Unit } from '@molstar/model/model/structure';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import type { Vec3 } from '@molstar/core/math/linear-algebra';
import {
  createLinkCylinderImpostors,
  createLinkCylinderMesh,
  type LinkBuilderProps,
  LinkCylinderParams,
  LinkStyle,
} from '@molstar/graphics/repr/structure/visual/util/link';
import {
  ComplexCylindersParams,
  ComplexCylindersVisual,
  ComplexMeshParams,
  type ComplexVisual,
  ComplexMeshVisual,
} from '@molstar/graphics/repr/structure/complex-visual';
import { Cylinders } from '@molstar/graphics/geo/geometry/cylinders/cylinders';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { checkCylinderImpostorSupport } from '@molstar/graphics/repr/structure/visual/util/common';
import type { VisualUpdateState } from '@molstar/graphics/repr/util';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { EmptyLoci, type Loci } from '@molstar/model/model/loci';
import { NullLocation } from '@molstar/model/model/location';
import { type IntMap, Interval, OrderedSet, SortedArray } from '@molstar/core/data/int';
import { InteractionsProvider } from '@molstar/model/props/computed/interactions';
import { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import { type BridgeContacts, Bridges } from '@molstar/model/props/computed/interactions/interactions';
import { Sphere3D } from '@molstar/core/math/geometry';
import { impostorProps, InteractionsSharedParams, withDashedCylinderCount } from './shared.js';
import type { Features } from '@molstar/model/props/computed/interactions/features';

type CanonicalLegIndices = {
  endpointA: Int32Array;
  endpointB: Int32Array;
};

const CanonicalLegIndicesCache = new WeakMap<BridgeContacts, CanonicalLegIndices>();

function getCanonicalLegIndices(bridges: BridgeContacts): CanonicalLegIndices {
  const cached = CanonicalLegIndicesCache.get(bridges);
  if (cached) return cached;

  const n = bridges.length;
  const endpointA = new Int32Array(n);
  const endpointB = new Int32Array(n);

  const legA = new Map<string, number>();
  const legB = new Map<string, number>();

  for (let i = 0; i < n; i++) {
    const b = bridges[i];

    const kA = `${b.unitA}|${b.indexA}|${b.unitM}|${b.indexMA}`;
    const kB = `${b.unitM}|${b.indexMB}|${b.unitB}|${b.indexB}`;

    let ai = legA.get(kA);
    if (ai === undefined) {
      ai = i;
      legA.set(kA, i);
    }
    endpointA[i] = ai;

    let bi = legB.get(kB);
    if (bi === undefined) {
      bi = i;
      legB.set(kB, i);
    }
    endpointB[i] = bi;
  }

  const indices = { endpointA, endpointB };
  CanonicalLegIndicesCache.set(bridges, indices);
  return indices;
}

function getFeatureMember(features: Features, featureIndex: Features.FeatureIndex): StructureElement.UnitIndex {
  return features.members[features.offsets[featureIndex]] as StructureElement.UnitIndex;
}

function atomPosition(unit: Unit.Atomic, features: Features, featureIndex: Features.FeatureIndex, out: Vec3) {
  const atomLocalIdx = getFeatureMember(features, featureIndex);
  unit.conformation.position(unit.elements[atomLocalIdx], out);
}

function setFeatureLocation(
  structure: Structure,
  location: StructureElement.Location,
  unitId: number,
  features: Features,
  featureIndex: Features.FeatureIndex,
) {
  const unit = structure.unitMap.get(unitId) as Unit.Atomic;
  const atomLocalIdx = getFeatureMember(features, featureIndex);

  location.unit = unit;
  location.element = unit.elements[atomLocalIdx];
}

/** Whether every atom of the feature is in `child` */
function isFeatureInChild(structure: Structure, child: Structure, unitId: number, features: Features, index: number) {
  const childUnit = child.unitMap.get(unitId);
  if (!childUnit) return false;
  const { elements } = structure.unitMap.get(unitId);
  const { offsets, members } = features;
  for (let i = offsets[index], il = offsets[index + 1]; i < il; ++i) {
    if (!SortedArray.has(childUnit.elements, elements[members[i]])) return false;
  }
  return true;
}

/**
 * With `includeParent`, the drawn legs by `parentDisplay` and whether the end features A and B are in the child, as the
 * inter-unit visual does for contacts: 'stub' the legs of the ends in the child, 'full' both legs when an end is in the
 * child, 'between' both legs when exactly one is. A leg shared by bridges (drawn at its first bridge) is drawn when any
 * of them draws it (perf-audit-2 09/B6).
 */
function getChildLegs(
  structure: Structure,
  child: Structure,
  bridges: BridgeContacts,
  unitsFeatures: IntMap<Features>,
  canonical: CanonicalLegIndices,
  parentDisplay: PD.Values<BridgeParams>['parentDisplay'],
) {
  const n = bridges.length;
  const shownA = new Uint8Array(n);
  const shownB = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    const b = bridges[i];
    const a = isFeatureInChild(structure, child, b.unitA, unitsFeatures.get(b.unitA), b.indexA);
    const c = isFeatureInChild(structure, child, b.unitB, unitsFeatures.get(b.unitB), b.indexB);
    let legA: boolean, legB: boolean;
    if (parentDisplay === 'stub') {
      legA = a;
      legB = c;
    } else if (parentDisplay === 'full') {
      legA = legB = a || c;
    } else if (parentDisplay === 'between') {
      legA = legB = a !== c;
    } else {
      assertUnreachable(parentDisplay);
    }
    if (legA) shownA[canonical.endpointA[i]] = 1;
    if (legB) shownB[canonical.endpointB[i]] = 1;
  }
  return { shownA, shownB };
}

function applyLegA(
  bridgeIndex: number,
  bridgeCount: number,
  canonical: CanonicalLegIndices,
  apply: (interval: Interval) => boolean,
) {
  let changed = false;
  const i = canonical.endpointA[bridgeIndex];

  if (apply(Interval.ofSingleton(i))) changed = true;
  if (apply(Interval.ofSingleton(i + bridgeCount))) changed = true;

  return changed;
}

function applyLegB(
  bridgeIndex: number,
  bridgeCount: number,
  canonical: CanonicalLegIndices,
  apply: (interval: Interval) => boolean,
) {
  let changed = false;
  const i = canonical.endpointB[bridgeIndex];

  if (apply(Interval.ofSingleton(i + 2 * bridgeCount))) changed = true;
  if (apply(Interval.ofSingleton(i + 3 * bridgeCount))) changed = true;

  return changed;
}

/** The bridges as dashed links, for the impostor and the mesh builder; undefined if there are none */
function getBridgeBuilderProps(
  structure: Structure,
  theme: Theme,
  props: PD.Values<BridgeParams>,
): LinkBuilderProps | undefined {
  if (!structure.hasAtomic) return;

  const interactions = InteractionsProvider.get(structure).value;
  if (!interactions) return;

  const { bridges, unitsFeatures } = interactions;

  const n = bridges.length;
  if (!n) return;

  const l = StructureElement.Location.create(structure);
  // as the impostors' uSizeFactor
  const sizeFactor = props.sizeFactor * props.sizeAspectRatio;
  const canonical = getCanonicalLegIndices(bridges);
  const { child } = structure;
  const childLegs = child && getChildLegs(structure, child, bridges, unitsFeatures, canonical, props.parentDisplay);

  return {
    // Four half-cylinders per bridge; the link builders draw the A-side half per call:
    //   [0,   n): A→mediator, forward  (A side)
    //   [n,  2n): A→mediator, backward (mediator side)
    //   [2n, 3n): mediator→B, forward  (mediator side)
    //   [3n, 4n): mediator→B, backward (B side)
    //
    // When multiple bridges share the same physical leg, only the first
    // occurrence is drawn; later ones map back to the canonical edge index.
    linkCount: 4 * n,

    position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
      const b = bridges[edgeIndex % n];
      const uM = structure.unitMap.get(b.unitM) as Unit.Atomic;
      const fM = unitsFeatures.get(b.unitM);
      const leg = Math.floor(edgeIndex / n);

      if (leg === 0) {
        const uA = structure.unitMap.get(b.unitA) as Unit.Atomic;
        const fA = unitsFeatures.get(b.unitA);
        atomPosition(uA, fA, b.indexA, posA);
        atomPosition(uM, fM, b.indexMA, posB);
      } else if (leg === 1) {
        const uA = structure.unitMap.get(b.unitA) as Unit.Atomic;
        const fA = unitsFeatures.get(b.unitA);
        atomPosition(uM, fM, b.indexMA, posA);
        atomPosition(uA, fA, b.indexA, posB);
      } else if (leg === 2) {
        const uB = structure.unitMap.get(b.unitB) as Unit.Atomic;
        const fB = unitsFeatures.get(b.unitB);
        atomPosition(uM, fM, b.indexMB, posA);
        atomPosition(uB, fB, b.indexB, posB);
      } else {
        const uB = structure.unitMap.get(b.unitB) as Unit.Atomic;
        const fB = unitsFeatures.get(b.unitB);
        atomPosition(uB, fB, b.indexB, posA);
        atomPosition(uM, fM, b.indexMB, posB);
      }
    },

    ignore: (edgeIndex: number) => {
      const bi = edgeIndex % n;
      const leg = Math.floor(edgeIndex / n);

      if (leg <= 1) return canonical.endpointA[bi] !== bi || (!!childLegs && !childLegs.shownA[bi]);
      return canonical.endpointB[bi] !== bi || (!!childLegs && !childLegs.shownB[bi]);
    },

    style: (_edgeIndex: number) => LinkStyle.Dashed,

    radius: (edgeIndex: number) => {
      const b = bridges[edgeIndex % n];
      const leg = Math.floor(edgeIndex / n);
      const isLegA = leg <= 1;

      if (isLegA) {
        const fA = unitsFeatures.get(b.unitA);
        const fM = unitsFeatures.get(b.unitM);

        setFeatureLocation(structure, l, b.unitA, fA, b.indexA);
        const sizeA = theme.size.size(l);

        setFeatureLocation(structure, l, b.unitM, fM, b.indexMA);
        const sizeM = theme.size.size(l);

        return Math.min(sizeA, sizeM) * sizeFactor;
      } else {
        const fM = unitsFeatures.get(b.unitM);
        const fB = unitsFeatures.get(b.unitB);

        setFeatureLocation(structure, l, b.unitM, fM, b.indexMB);
        const sizeM = theme.size.size(l);

        setFeatureLocation(structure, l, b.unitB, fB, b.indexB);
        const sizeB = theme.size.size(l);

        return Math.min(sizeM, sizeB) * sizeFactor;
      }
    },
  };
}

function setBoundingSphere(
  geometry: Mesh | Cylinders,
  boundingSphere: Sphere3D | undefined,
  structure: Structure,
  sizeFactor: number,
) {
  if (boundingSphere) {
    geometry.setBoundingSphere(boundingSphere);
  } else if (geometry.kind === 'mesh' ? geometry.triangleCount > 0 : geometry.cylinderCount > 0) {
    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, sizeFactor);
    geometry.setBoundingSphere(sphere);
  }
}

function createBridgeCylinderImpostors(
  ctx: VisualContext,
  structure: Structure,
  theme: Theme,
  props: PD.Values<BridgeParams>,
  cylinders?: Cylinders,
) {
  const builderProps = getBridgeBuilderProps(structure, theme, props);
  if (!builderProps) return Cylinders.createEmpty(cylinders);

  const { cylinders: c, boundingSphere } = createLinkCylinderImpostors(
    ctx,
    withDashedCylinderCount(builderProps, props.dashCount),
    impostorProps(props),
    cylinders,
  );
  setBoundingSphere(c, boundingSphere, structure, props.sizeFactor);
  return c;
}

function createBridgeCylinderMesh(
  ctx: VisualContext,
  structure: Structure,
  theme: Theme,
  props: PD.Values<BridgeParams>,
  mesh?: Mesh,
) {
  const builderProps = getBridgeBuilderProps(structure, theme, props);
  if (!builderProps) return Mesh.createEmpty(mesh);

  const { mesh: m, boundingSphere } = createLinkCylinderMesh(ctx, builderProps, props, mesh);
  setBoundingSphere(m, boundingSphere, structure, props.sizeFactor);
  return m;
}

export const BridgeParams = {
  ...ComplexMeshParams,
  ...ComplexCylindersParams,
  ...LinkCylinderParams,
  ...InteractionsSharedParams,
};
export type BridgeParams = typeof BridgeParams;

export const BridgeVisual = chooseVisual(
  { impostor: BridgeImpostorVisual, mesh: BridgeMeshVisual },
  (props: PD.Values<BridgeParams>, graphics?: GraphicsContext) =>
    props.tryUseImpostor && checkCylinderImpostorSupport(graphics) ? 'impostor' : 'mesh',
);

function setUpdateState(
  state: VisualUpdateState,
  newProps: PD.Values<BridgeParams>,
  currentProps: PD.Values<BridgeParams>,
  newTheme: Theme,
  currentTheme: Theme,
  newStructure: Structure,
  _currentStructure: Structure,
) {
  state.createGeometry =
    newProps.dashCount !== currentProps.dashCount ||
    newProps.dashScale !== currentProps.dashScale ||
    newProps.dashCap !== currentProps.dashCap ||
    newProps.parentDisplay !== currentProps.parentDisplay;
  // a size change rebuilds the mesh (ComplexMeshVisual, by SizeTheme.areEqual) and updates the impostors' sizes

  const interactionsHash = InteractionsProvider.get(newStructure).version;
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
  newProps: PD.Values<BridgeParams>,
  newTheme: Theme,
  newStructure: Structure,
) {
  state.info.interactionsHash = InteractionsProvider.get(newStructure).version;
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

export function BridgeImpostorVisual(materialId: number): ComplexVisual<BridgeParams> {
  return ComplexCylindersVisual<BridgeParams>(
    {
      defaultProps: PD.getDefaultValues(BridgeParams),
      createGeometry: createBridgeCylinderImpostors,
      createLocationIterator: createBridgeIterator,
      getLoci: getBridgeLoci,
      eachLocation: eachBridgeInteraction,
      setUpdateState,
      initUpdateState,
      mustRecreate: (_: Structure, props: PD.Values<BridgeParams>, graphics?: GraphicsContext) =>
        !props.tryUseImpostor || !checkCylinderImpostorSupport(graphics),
    },
    materialId,
  );
}

export function BridgeMeshVisual(materialId: number): ComplexVisual<BridgeParams> {
  return ComplexMeshVisual<BridgeParams>(
    {
      defaultProps: PD.getDefaultValues(BridgeParams),
      createGeometry: createBridgeCylinderMesh,
      createLocationIterator: createBridgeIterator,
      getLoci: getBridgeLoci,
      eachLocation: eachBridgeInteraction,
      setUpdateState: setMeshUpdateState,
      initUpdateState,
      mustRecreate: (_: Structure, props: PD.Values<BridgeParams>, graphics?: GraphicsContext) =>
        props.tryUseImpostor && checkCylinderImpostorSupport(graphics),
    },
    materialId,
  );
}

function getBridgeLoci(pickingId: PickingId, structure: Structure, id: number) {
  const { objectId, groupId } = pickingId;
  if (id !== objectId) return EmptyLoci;

  const interactions = InteractionsProvider.get(structure).value;
  if (!interactions) return EmptyLoci;

  const { bridges, unitsFeatures } = interactions;
  const n = bridges.length;

  if (!n || groupId < 0 || groupId >= 4 * n) return EmptyLoci;

  const bridgeIndex = groupId % n;

  return Bridges.Loci({ structure, bridges, unitsFeatures }, [{ bridgeIndex }]);
}

const __unitMap = new Map<number, OrderedSet<StructureElement.UnitIndex>>();

function eachBridgeInteraction(
  loci: Loci,
  structure: Structure,
  apply: (interval: Interval) => boolean,
  _isMarking: boolean,
) {
  let changed = false;

  if (Bridges.isLoci(loci)) {
    if (!Structure.areEquivalent(loci.data.structure, structure)) return false;

    const interactions = InteractionsProvider.get(structure).value;
    if (!interactions) return false;

    const { bridges } = interactions;
    const n = bridges.length;
    if (!n) return false;

    const canonical = getCanonicalLegIndices(bridges);

    for (const e of loci.elements) {
      if (e.bridgeIndex < 0 || e.bridgeIndex >= n) continue;

      if (applyLegA(e.bridgeIndex, n, canonical, apply)) changed = true;
      if (applyLegB(e.bridgeIndex, n, canonical, apply)) changed = true;
    }
  } else if (StructureElement.Loci.is(loci)) {
    if (!Structure.areEquivalent(loci.structure, structure)) return false;

    const interactions = InteractionsProvider.get(structure).value;
    if (!interactions) return false;

    const { bridges, unitsFeatures } = interactions;
    const n = bridges.length;
    if (!n) return false;

    const canonical = getCanonicalLegIndices(bridges);

    __unitMap.clear();
    for (const e of loci.elements) {
      __unitMap.set(e.unit.id, e.indices);
    }

    for (let i = 0; i < n; i++) {
      const b = bridges[i];

      const indicesA = __unitMap.get(b.unitA);
      const indicesM = __unitMap.get(b.unitM);
      const indicesB = __unitMap.get(b.unitB);

      if (!indicesA && !indicesM && !indicesB) continue;

      let hitA = false;
      if (indicesA) {
        const fA = unitsFeatures.get(b.unitA);
        const mi = getFeatureMember(fA, b.indexA);
        hitA = OrderedSet.has(indicesA, mi);
      }

      let hitM = false;
      if (indicesM) {
        const fM = unitsFeatures.get(b.unitM);
        const miA = getFeatureMember(fM, b.indexMA);
        const miB = getFeatureMember(fM, b.indexMB);
        hitM = OrderedSet.has(indicesM, miA) || OrderedSet.has(indicesM, miB);
      }

      let hitB = false;
      if (indicesB) {
        const fB = unitsFeatures.get(b.unitB);
        const mi = getFeatureMember(fB, b.indexB);
        hitB = OrderedSet.has(indicesB, mi);
      }

      if (hitA || hitM) {
        if (applyLegA(i, n, canonical, apply)) changed = true;
      }

      if (hitB || hitM) {
        if (applyLegB(i, n, canonical, apply)) changed = true;
      }
    }

    __unitMap.clear();
  }

  return changed;
}

function createBridgeIterator(structure: Structure): LocationIterator {
  const interactions = InteractionsProvider.get(structure).value;
  if (!interactions) return LocationIterator(0, 1, 1, () => NullLocation, true);

  const { bridges, unitsFeatures } = interactions;

  const n = bridges.length;
  const groupCount = 4 * n;
  const instanceCount = 1;

  const data: Bridges.Data = { structure, bridges, unitsFeatures };
  const location = Bridges.Location(data);
  const { element } = location;

  const getLocation = (groupIndex: number) => {
    element.bridgeIndex = n === 0 ? 0 : groupIndex % n;
    return location;
  };

  return LocationIterator(groupCount, instanceCount, 1, getLocation, true);
}
