/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import {
  Representation,
  type RepresentationContext,
  type RepresentationParamsGetter,
} from '@molstar/graphics/repr/representation';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import type { Vec3 } from '@molstar/core/math/linear-algebra';
import { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { EmptyLoci, type Loci } from '@molstar/model/model/loci';
import { Interval } from '@molstar/core/data/int';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import type { VisualContext } from '@molstar/graphics/repr/visual';
import { createLinkCylinderMesh, LinkCylinderParams } from '@molstar/graphics/repr/structure/visual/util/link';
import {
  ComplexMeshParams,
  type ComplexVisual,
  ComplexMeshVisual,
} from '@molstar/graphics/repr/structure/complex-visual';
import type { VisualUpdateState } from '@molstar/graphics/repr/util';
import {
  ComplexRepresentation,
  type StructureRepresentation,
  StructureRepresentationStateBuilder,
  StructureRepresentationProvider,
} from '@molstar/graphics/repr/structure/representation';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import {
  CrossLinkRestraintProvider,
  CrossLinkRestraint,
} from '@molstar/model/props/integrative/cross-link-restraint/property';
import { Sphere3D } from '@molstar/core/math/geometry';

function createCrossLinkRestraintCylinderMesh(
  ctx: VisualContext,
  structure: Structure,
  theme: Theme,
  props: PD.Values<CrossLinkRestraintCylinderParams>,
  mesh?: Mesh,
) {
  const crossLinks = CrossLinkRestraintProvider.get(structure).value!;
  if (!crossLinks.count) return Mesh.createEmpty(mesh);
  const { sizeFactor } = props;

  const location = StructureElement.Location.create(structure);

  const builderProps = {
    linkCount: crossLinks.count,
    position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
      const b = crossLinks.pairs[edgeIndex];
      const uA = b.unitA,
        uB = b.unitB;
      uA.conformation.position(uA.elements[b.indexA], posA);
      uB.conformation.position(uB.elements[b.indexB], posB);
    },
    radius: (edgeIndex: number) => {
      const b = crossLinks.pairs[edgeIndex];
      location.unit = b.unitA;
      location.element = b.unitA.elements[b.indexA];
      return theme.size.size(location) * sizeFactor;
    },
  };

  const { mesh: m, boundingSphere } = createLinkCylinderMesh(ctx, builderProps, props, mesh);

  if (boundingSphere) {
    m.setBoundingSphere(boundingSphere);
  } else if (m.triangleCount > 0) {
    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, 1 * sizeFactor);
    m.setBoundingSphere(sphere);
  }

  return m;
}

export const CrossLinkRestraintCylinderParams = {
  ...ComplexMeshParams,
  ...LinkCylinderParams,
  sizeFactor: PD.Numeric(0.5, { min: 0, max: 10, step: 0.1 }),
};
export type CrossLinkRestraintCylinderParams = typeof CrossLinkRestraintCylinderParams;

export function CrossLinkRestraintVisual(materialId: number): ComplexVisual<CrossLinkRestraintCylinderParams> {
  return ComplexMeshVisual<CrossLinkRestraintCylinderParams>(
    {
      defaultProps: PD.getDefaultValues(CrossLinkRestraintCylinderParams),
      createGeometry: createCrossLinkRestraintCylinderMesh,
      createLocationIterator: createCrossLinkRestraintIterator,
      getLoci: getLinkLoci,
      eachLocation: eachCrossLink,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<CrossLinkRestraintCylinderParams>,
        currentProps: PD.Values<CrossLinkRestraintCylinderParams>,
      ) => {
        state.createGeometry =
          newProps.sizeFactor !== currentProps.sizeFactor ||
          newProps.radialSegments !== currentProps.radialSegments ||
          newProps.linkCap !== currentProps.linkCap;
      },
    },
    materialId,
  );
}

function createCrossLinkRestraintIterator(structure: Structure): LocationIterator {
  const crossLinkRestraints = CrossLinkRestraintProvider.get(structure).value!;
  const { pairs } = crossLinkRestraints;
  const groupCount = pairs.length;
  const instanceCount = 1;
  const location = CrossLinkRestraint.Location(crossLinkRestraints, structure);
  const getLocation = (groupIndex: number) => {
    location.element = groupIndex;
    return location;
  };
  return LocationIterator(groupCount, instanceCount, 1, getLocation, true);
}

function getLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
  const { objectId, groupId } = pickingId;
  if (id === objectId) {
    const crossLinkRestraints = CrossLinkRestraintProvider.get(structure).value!;
    const pair = crossLinkRestraints.pairs[groupId];
    if (pair) {
      return CrossLinkRestraint.Loci(structure, crossLinkRestraints, [groupId]);
    }
  }
  return EmptyLoci;
}

function eachCrossLink(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
  let changed = false;
  if (CrossLinkRestraint.isLoci(loci)) {
    if (!Structure.areEquivalent(loci.data.structure, structure)) return false;
    const crossLinkRestraints = CrossLinkRestraintProvider.get(structure).value!;
    if (loci.data.crossLinkRestraints !== crossLinkRestraints) return false;

    for (const e of loci.elements) {
      if (apply(Interval.ofSingleton(e))) changed = true;
    }
  }
  return changed;
}

//

const CrossLinkRestraintVisuals = {
  'cross-link-restraint': (
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<Structure, CrossLinkRestraintCylinderParams>,
  ) => ComplexRepresentation('Cross-link restraint', ctx, getParams, CrossLinkRestraintVisual),
};

export const CrossLinkRestraintParams = {
  ...CrossLinkRestraintCylinderParams,
};
export type CrossLinkRestraintParams = typeof CrossLinkRestraintParams;
export function getCrossLinkRestraintParams(ctx: ThemeRegistryContext, structure: Structure) {
  return PD.clone(CrossLinkRestraintParams);
}

export type CrossLinkRestraintRepresentation = StructureRepresentation<CrossLinkRestraintParams>;
export function CrossLinkRestraintRepresentation(
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<Structure, CrossLinkRestraintParams>,
): CrossLinkRestraintRepresentation {
  return Representation.createMulti(
    'CrossLinkRestraint',
    ctx,
    getParams,
    StructureRepresentationStateBuilder,
    CrossLinkRestraintVisuals as unknown as Representation.Def<Structure, CrossLinkRestraintParams>,
  );
}

export const CrossLinkRestraintRepresentationProvider = StructureRepresentationProvider({
  name: CrossLinkRestraint.Tag.CrossLinkRestraint,
  label: 'Cross Link Restraint',
  description: 'Displays cross-link restraints.',
  factory: CrossLinkRestraintRepresentation,
  getParams: getCrossLinkRestraintParams,
  defaultValues: PD.getDefaultValues(CrossLinkRestraintParams),
  defaultColorTheme: { name: 'cross-link' },
  defaultSizeTheme: { name: 'uniform' },
  isApplicable: (structure: Structure) => CrossLinkRestraint.isApplicable(structure),
  ensureCustomProperties: {
    attach: (ctx: CustomProperty.Context, structure: Structure) =>
      CrossLinkRestraintProvider.attach(ctx, structure, void 0, true),
    detach: (data) => CrossLinkRestraintProvider.ref(data, false),
  },
});
