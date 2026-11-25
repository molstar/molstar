/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { releaseGpuMeshVertices } from '@molstar/graphics/geo/geometry/mesh/gpu-vertices';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Structure, Unit, StructureElement, Bond } from '@molstar/model/model/structure';
import type { RepresentationProps } from '../representation.js';
import { Visual, type VisualContext } from '../visual.js';
import { Geometry, type GeometryUtils } from '@molstar/graphics/geo/geometry/geometry';
import type { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import { Theme } from '@molstar/graphics/theme/theme';
import { createUnitsTransform, includesUnitKind, type StructureGroup } from './visual/util/common.js';
import { getGroupLoci, isWholeStructureLoci } from './visual/util/loci-index.js';
import {
  createRenderObject,
  type GraphicsRenderObject,
  type RenderObjectValues,
} from '@molstar/graphics/gl/render-object';
import { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { type Loci, isEveryLoci, EmptyLoci } from '@molstar/model/model/loci';
import { Interval } from '@molstar/core/data/int';
import type { IntMap } from '@molstar/core/data/int/map';
import { type LocationCallback, VisualUpdateState } from '../util.js';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { createMarkers } from '@molstar/graphics/geo/geometry/marker-data';
import { type BaseGeometry, resolveInstanceGranularity } from '@molstar/graphics/geo/geometry/base';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ValueCell, deepEqual } from '@molstar/core/util';
import { createSizes } from '@molstar/graphics/geo/geometry/size-data';
import { createColors } from '@molstar/graphics/geo/geometry/color-data';
import type { Mat4 } from '@molstar/core/math/linear-algebra';
import type { Overpaint } from '@molstar/graphics/theme/overpaint';
import type { Transparency } from '@molstar/graphics/theme/transparency';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { Cylinders } from '@molstar/graphics/geo/geometry/cylinders/cylinders';
import { Points } from '@molstar/graphics/geo/geometry/points/points';
import { Lines } from '@molstar/graphics/geo/geometry/lines/lines';
import { Text } from '@molstar/graphics/geo/geometry/text/text';
import { DirectVolume } from '@molstar/graphics/geo/geometry/direct-volume/direct-volume';
import { TextureMesh } from '@molstar/graphics/geo/geometry/texture-mesh/texture-mesh';
import { Image } from '@molstar/graphics/geo/geometry/image/image';
import type { SizeValues } from '@molstar/graphics/gl/renderable/schema';
import {
  StructureParams,
  StructureMeshParams,
  StructureSpheresParams,
  StructurePointsParams,
  StructureLinesParams,
  StructureTextParams,
  StructureDirectVolumeParams,
  StructureTextureMeshParams,
  StructureCylindersParams,
  StructureImageParams,
} from './params.js';
import type { Clipping } from '@molstar/graphics/theme/clipping';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Substance } from '@molstar/graphics/theme/substance';
import type { Emissive } from '@molstar/graphics/theme/emissive';
import type { Wiggle } from '@molstar/graphics/theme/wiggle';
import { getGeometrySignatures, type SignatureContext, type SignatureDefines } from '@molstar/graphics/repr/signatures';

export interface UnitsVisual<P extends RepresentationProps = {}> extends Visual<StructureGroup, P> {
  /** the builder it is made of: a visual over many groups made of it alike (`UnitsMergedVisual`) */
  readonly builder?: UnitsVisualGeometryBuilder<any, any>;
}

function createUnitsRenderObject<G extends Geometry>(
  structureGroup: StructureGroup,
  geometry: G,
  locationIt: LocationIterator,
  theme: Theme,
  props: PD.Values<StructureParams & Geometry.Params<G>>,
  materialId: number,
) {
  const { createValues, createRenderableState } = Geometry.getUtils(geometry);
  const transform = createUnitsTransform(
    structureGroup,
    props.includeParent,
    geometry.boundingSphere,
    props.cellSize,
    props.batchSize,
  );
  const values = createValues(geometry, transform, locationIt, theme, props);
  const state = createRenderableState(props);
  return createRenderObject(geometry.kind, values, state, materialId);
}

interface UnitsVisualBuilder<P extends StructureParams, G extends Geometry> {
  /** define values the visual decides (e.g. varying groups), for its signatures, see `getGeometrySignatures` */
  signatureDefines?: (props: PD.Values<P>, graphics?: GraphicsContext) => SignatureDefines;
  defaultProps: PD.Values<P>;
  createGeometry(
    ctx: VisualContext,
    unit: Unit,
    structure: Structure,
    theme: Theme,
    props: PD.Values<P>,
    geometry?: G,
  ): Promise<G> | G;
  createLocationIterator(structureGroup: StructureGroup, props: PD.Values<P>): LocationIterator;
  getLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number): Loci;
  eachLocation(
    loci: Loci,
    structureGroup: StructureGroup,
    apply: (interval: Interval) => boolean,
    isMarking: boolean,
  ): boolean;
  setUpdateState(
    state: VisualUpdateState,
    newProps: PD.Values<P>,
    currentProps: PD.Values<P>,
    newTheme: Theme,
    currentTheme: Theme,
    newStructureGroup: StructureGroup,
    currentStructureGroup: StructureGroup,
  ): void;
  initUpdateState?: (
    state: VisualUpdateState,
    newProps: PD.Values<P>,
    newTheme: Theme,
    newStructureGroup: StructureGroup,
  ) => void;
  mustRecreate?: (structureGroup: StructureGroup, props: PD.Values<P>) => boolean;
  /**
   * The values after they are written (color smoothing); `graphics` none in a merged visual (`UnitsMergedVisual`):
   * smoothing on the CPU into data per vertex (a grid per object does not merge)
   */
  processValues?: (
    values: RenderObjectValues<G['kind']>,
    geometry: G,
    props: PD.Values<P>,
    theme: Theme,
    graphics?: GraphicsContext,
  ) => void;
  dispose?: (geometry: G) => void;
}

export interface UnitsVisualGeometryBuilder<P extends StructureParams, G extends Geometry>
  extends UnitsVisualBuilder<P, G> {
  geometryUtils: GeometryUtils<G>;
}

/**
 * The update state of a units visual's values `values` (of its geometry `geometry`) from the current props, theme and
 * group to the new: the builder's (`setUpdateState`) and what the theme, transforms, conformations and hierarchy
 * change. `foreign`: its values hold resources of another context than `graphics`
 */
export function setUnitsUpdateState<P extends StructureParams & Geometry.Params<G>, G extends Geometry>(
  updateState: VisualUpdateState,
  builder: UnitsVisualGeometryBuilder<P, G>,
  values: RenderObjectValues<G['kind']>,
  geometry: G | undefined,
  newProps: PD.Values<P>,
  currentProps: PD.Values<P>,
  newTheme: Theme,
  currentTheme: Theme,
  newStructureGroup: StructureGroup,
  currentStructureGroup: StructureGroup,
  graphics: GraphicsContext | undefined,
  foreign: boolean,
) {
  builder.setUpdateState(
    updateState,
    newProps,
    currentProps,
    newTheme,
    currentTheme,
    newStructureGroup,
    currentStructureGroup,
  );

  // after a backend switch or restore: textures and buffers of the previous context are rendered or computed again,
  // geometry relying on what the new context lacks (`GeometryUtils.isCompatible`) is created again, for the new context
  if (foreign) updateState.createGeometry = true;
  if (graphics && geometry && builder.geometryUtils.isCompatible?.(geometry, graphics) === false) {
    updateState.createGeometry = true;
  }

  if (!Structure.areHierarchiesEqual(currentStructureGroup.structure, newStructureGroup.structure)) {
    updateState.updateTransform = true;
    updateState.updateColor = true;
    updateState.updateSize = true;
  }

  if (!ColorTheme.areEqual(newTheme.color, currentTheme.color)) {
    updateState.updateColor = true;
  }

  if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) {
    updateState.updateSize = true;
  }

  if (currentStructureGroup.structure.child !== newStructureGroup.structure.child) {
    updateState.createGeometry = true;
    updateState.updateTransform = true;
  }

  // the representation keeps its props object while no prop changes (`shallowMerge2`): the same props, no prop changed
  const propsChanged = newProps !== currentProps;
  if (
    propsChanged &&
    (newProps.instanceGranularity !== currentProps.instanceGranularity ||
      newProps.cellSize !== currentProps.cellSize ||
      newProps.batchSize !== currentProps.batchSize)
  ) {
    updateState.updateTransform = true;
  }

  if (propsChanged && !deepEqual(newProps.unitKinds, currentProps.unitKinds)) {
    updateState.createGeometry = true;
  }

  if (newStructureGroup.group.transformHash !== currentStructureGroup.group.transformHash) {
    if (newStructureGroup.group.units.length !== currentStructureGroup.group.units.length || updateState.updateColor) {
      updateState.updateTransform = true;
    } else {
      updateState.updateMatrix = true;
    }
  }

  // check if the operator or conformation of unit has changed
  const newUnit = newStructureGroup.group.units[0];
  const currentUnit = currentStructureGroup.group.units[0];
  if (!Unit.areOperatorsEqual(newUnit, currentUnit)) {
    updateState.updateTransform = true;
  }
  if (!Unit.areConformationsEqual(newUnit, currentUnit)) {
    updateState.createGeometry = true;
  }

  if (updateState.updateTransform) {
    updateState.updateMatrix = true;
  }

  if (updateState.updateSize && !('uSize' in values)) {
    updateState.createGeometry = true;
  }

  if (updateState.createGeometry || updateState.updateTransform) {
    if (currentStructureGroup.structure.hashCode !== newStructureGroup.structure.hashCode) {
      updateState.updateColor = true;
      updateState.updateSize = true;
    }
    if (
      newTheme.color.granularity.startsWith('vertex') ||
      values.colorType.ref.value.startsWith('vertex') ||
      newTheme.color.granularity.startsWith('volume') ||
      values.colorType.ref.value.startsWith('volume')
    ) {
      updateState.updateColor = true;
    }
  }
}

/** A units visual's iterators: of its locations (theme data, marking) and of its geometry's positions */
export type UnitsIterators = { locationIt: LocationIterator; positionIt: LocationIterator };

/**
 * Updates the values `values` of a units visual by `updateState` (`setUnitsUpdateState`), not its renderable state:
 * markers, transforms, counts, bounds, sizes, colors and the values of the props; `its` made again where what they
 * iterate may differ. `newGeometry`: the geometry made again (`updateState.createGeometry`), else `geometry` stays
 */
export function updateUnitsValues<P extends StructureParams & Geometry.Params<G>, G extends Geometry>(
  builder: UnitsVisualGeometryBuilder<P, G>,
  updateState: VisualUpdateState,
  values: RenderObjectValues<G['kind']>,
  its: UnitsIterators,
  newGeometry: G | undefined,
  geometry: G,
  newProps: PD.Values<P>,
  currentProps: PD.Values<P>,
  newTheme: Theme,
  newStructureGroup: StructureGroup,
  currentStructureGroup: StructureGroup,
) {
  const { updateValues, updateBoundingSphere, createPositionIterator } = builder.geometryUtils;
  // a new iterator where what it iterates may differ: the transform, the group, the structure or the props (the
  // color and size writers reset it); a recolor keeps it
  if (
    updateState.updateTransform ||
    ((updateState.updateColor || updateState.updateSize) &&
      (newProps !== currentProps ||
        newStructureGroup.group !== currentStructureGroup.group ||
        newStructureGroup.structure !== currentStructureGroup.structure))
  ) {
    its.locationIt = builder.createLocationIterator(newStructureGroup, newProps);
  }
  const { locationIt } = its;

  if (updateState.updateTransform) {
    const { instanceCount, groupCount } = locationIt;
    if (resolveInstanceGranularity(newProps.instanceGranularity, groupCount, instanceCount)) {
      createMarkers(instanceCount, 'instance', values);
    } else {
      createMarkers(instanceCount * groupCount, 'groupInstance', values);
    }
  }

  if (updateState.updateMatrix) {
    createUnitsTransform(
      newStructureGroup,
      newProps.includeParent,
      newGeometry?.boundingSphere || values.invariantBoundingSphere.ref.value,
      newProps.cellSize,
      newProps.batchSize,
      values,
    );
    if ('lodLevels' in values) {
      // to trigger `uLod` update in `renderable.cull`
      ValueCell.update(values.lodLevels, values.lodLevels.ref.value);
    }
  }

  if (updateState.createGeometry) {
    if (newGeometry) {
      ValueCell.updateIfChanged(values.drawCount, Geometry.getDrawCount(newGeometry));
      ValueCell.updateIfChanged(values.uVertexCount, Geometry.getVertexCount(newGeometry));
      ValueCell.updateIfChanged(values.uGroupCount, locationIt.groupCount);
    } else {
      throw new Error('expected geometry to be given');
    }
  }

  if (updateState.updateTransform || updateState.createGeometry) {
    updateBoundingSphere(values, newGeometry || geometry);
    its.positionIt = createPositionIterator(newGeometry || geometry, values);
  }

  if (updateState.updateSize) {
    // not all geometries have size data, so check here
    if ('uSize' in values) {
      createSizes(locationIt, its.positionIt, newTheme.size, values as SizeValues);
    }
  }

  if (updateState.updateColor) {
    createColors(locationIt, its.positionIt, newTheme.color, values);
  }

  // the values from the props: again where the props or what the values derive from besides them changed (counts,
  // bounds, sizes); a recolor of many units keeps them
  if (
    newProps !== currentProps ||
    updateState.createGeometry ||
    updateState.updateTransform ||
    updateState.updateMatrix ||
    updateState.updateSize
  ) {
    updateValues(values, newProps);
    // the bounds' padding by the size factor (written by `updateValues`): again where it changed alone
    if (
      (newProps as { sizeFactor?: number }).sizeFactor !== (currentProps as { sizeFactor?: number }).sizeFactor &&
      !updateState.updateTransform &&
      !updateState.createGeometry
    ) {
      updateBoundingSphere(values, geometry);
    }
    return true;
  }
  return false;
}

/** Whether `loci` covers all of the structure of `structureGroup` (every location of a units visual) */
function isUnitsLociSuperset(loci: Loci, structureGroup: StructureGroup) {
  if (isEveryLoci(loci)) return true;
  if (Structure.isLoci(loci) && Structure.areRootsEquivalent(loci.structure, structureGroup.structure)) return true;
  if (StructureElement.Loci.is(loci) && Structure.areRootsEquivalent(loci.structure, structureGroup.structure)) {
    if (isWholeStructureLoci(loci)) return true;
  }
  return false;
}

function eachUnitsInstance(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
  let changed = false;
  if (Bond.isLoci(loci)) {
    const { structure, group } = structureGroup;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;
    for (const b of loci.bonds) {
      if (b.aUnit !== b.bUnit) continue;
      const unitIdx = group.unitIndexMap.get(b.aUnit.id);
      if (unitIdx !== undefined) {
        if (apply(Interval.ofSingleton(unitIdx))) changed = true;
      }
    }
  } else if (StructureElement.Loci.is(loci)) {
    const { structure, group } = structureGroup;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;
    for (const e of loci.elements) {
      const unitIdx = group.unitIndexMap.get(e.unit.id);
      if (unitIdx !== undefined) {
        if (apply(Interval.ofSingleton(unitIdx))) changed = true;
      }
    }
  }
  return changed;
}

/**
 * Applies `apply` to the intervals of a units visual's markers or layers `loci` reaches (by instance or by group and
 * instance, as its granularity resolves). `visualLoci`: `loci` on the group's units alone (`undefined`: none of them),
 * computed only where `loci` is not a superset
 */
export function applyUnitsLoci(
  eachLocation: UnitsVisualBuilder<any, any>['eachLocation'],
  instanceGranularity: PD.Values<BaseGeometry.Params>['instanceGranularity'],
  structureGroup: StructureGroup,
  groupCount: number,
  instanceCount: number,
  loci: Loci,
  getVisualLoci: (loci: Loci) => Loci | undefined,
  apply: (interval: Interval) => boolean,
  isMarking: boolean,
) {
  const byInstance = resolveInstanceGranularity(instanceGranularity, groupCount, instanceCount);
  if (isUnitsLociSuperset(loci, structureGroup)) {
    return apply(Interval.ofBounds(0, byInstance ? instanceCount : groupCount * instanceCount));
  }
  const visualLoci = getVisualLoci(loci);
  if (!visualLoci) return false;
  return byInstance
    ? eachUnitsInstance(visualLoci, structureGroup, apply)
    : eachLocation(visualLoci, structureGroup, apply, isMarking);
}

export function UnitsVisual<G extends Geometry, P extends StructureParams & Geometry.Params<G>>(
  builder: UnitsVisualGeometryBuilder<P, G>,
  materialId: number,
): UnitsVisual<P> {
  const {
    defaultProps,
    createGeometry,
    createLocationIterator,
    getLoci,
    eachLocation,
    initUpdateState,
    mustRecreate,
    processValues,
    dispose,
  } = builder;
  const { createEmpty: createEmptyGeometry, updateRenderableState, createPositionIterator } = builder.geometryUtils;
  const updateState = VisualUpdateState.create();
  const foreignResources = Visual.foreignResourcesCheck();
  let newGraphics: GraphicsContext | undefined;
  /** the context `geometry` was made with: a new render object over geometry of another context is not known clean */
  let geometryGraphics: GraphicsContext | undefined;
  const previousMark: Visual.PreviousMark = { loci: EmptyLoci, action: MarkerAction.None, status: -1 };

  let renderObject: GraphicsRenderObject<G['kind']> | undefined;

  let newProps: PD.Values<P> = Object.assign({}, defaultProps);
  let newTheme: Theme = Theme.createEmpty();
  let newStructureGroup: StructureGroup;

  let currentProps: PD.Values<P>;
  let currentTheme: Theme;
  let currentStructureGroup: StructureGroup;

  let geometry: G;
  let geometryVersion = -1;
  let locationIt: LocationIterator;
  let positionIt: LocationIterator;

  function prepareUpdate(
    theme: Theme,
    props: PD.Values<P>,
    structureGroup: StructureGroup,
    graphics?: GraphicsContext,
  ) {
    if (!structureGroup && !currentStructureGroup) {
      throw new Error('missing structureGroup');
    }

    newProps = props;
    newTheme = theme;
    newStructureGroup = structureGroup;

    VisualUpdateState.reset(updateState);
    newGraphics = graphics;

    if (!renderObject || !currentStructureGroup) {
      initUpdateState?.(updateState, newProps, newTheme, newStructureGroup);
      // console.log('create new');
      updateState.createNew = true;
      updateState.createGeometry = true;
      return;
    }

    setUnitsUpdateState(
      updateState,
      builder,
      renderObject.values,
      geometry,
      newProps,
      currentProps,
      newTheme,
      currentTheme,
      newStructureGroup,
      currentStructureGroup,
      graphics,
      !!graphics && foreignResources.hasForeign(renderObject, graphics),
    );
  }

  function update(newGeometry?: G) {
    if (updateState.createNew) {
      locationIt = createLocationIterator(newStructureGroup, newProps);
      if (newGeometry) {
        // the old object's mesh vertices on the GPU (the smoothing's splat), kept by its values
        if (renderObject) releaseGpuMeshVertices(renderObject.values);
        renderObject = createUnitsRenderObject(
          newStructureGroup,
          newGeometry,
          locationIt,
          newTheme,
          newProps,
          materialId,
        );
        if (geometry === undefined || geometryGraphics === newGraphics)
          foreignResources.created(renderObject, newGraphics);
        positionIt = createPositionIterator(newGeometry, renderObject.values);
      } else {
        throw new Error('expected geometry to be given');
      }
    } else {
      if (!renderObject) {
        throw new Error('expected renderObject to be available');
      }

      const its: UnitsIterators = { locationIt, positionIt };
      const propsUpdated = updateUnitsValues(
        builder,
        updateState,
        renderObject.values,
        its,
        newGeometry,
        geometry,
        newProps,
        currentProps,
        newTheme,
        newStructureGroup,
        currentStructureGroup,
      );
      locationIt = its.locationIt;
      positionIt = its.positionIt;
      if (propsUpdated) updateRenderableState(renderObject.state, newProps);
    }

    currentProps = newProps;
    currentTheme = newTheme;
    currentStructureGroup = newStructureGroup;
    if (newGeometry) {
      geometry = newGeometry;
      geometryGraphics = newGraphics;
      geometryVersion += 1;
    }
  }

  function _createGeometry(
    ctx: VisualContext,
    unit: Unit,
    structure: Structure,
    theme: Theme,
    props: PD.Values<P>,
    geometry?: G,
  ) {
    return includesUnitKind(props.unitKinds, unit)
      ? createGeometry(ctx, unit, structure, theme, props, geometry)
      : createEmptyGeometry(geometry);
  }

  /** the units of `currentStructureGroup.group` and its index, read once per group (see `getGroupLoci`) */
  const lociGroup: { group?: Unit.SymmetryGroup; units: ReadonlyArray<Unit>; unitIndexMap?: IntMap<number> } = {
    units: [],
  };

  /** `loci` on the units of this visual only, `undefined` when it has none of them (perf-audit 07 P3) */
  function getVisualLoci(loci: Loci): Loci | undefined {
    if (!StructureElement.Loci.is(loci)) return loci;
    const { group } = currentStructureGroup;
    if (lociGroup.group !== group) {
      lociGroup.group = group;
      lociGroup.units = group.units;
      lociGroup.unitIndexMap = group.unitIndexMap;
    }
    return getGroupLoci(loci, lociGroup.units, lociGroup.unitIndexMap!);
  }

  function lociApply(loci: Loci, apply: (interval: Interval) => boolean, isMarking: boolean) {
    return applyUnitsLoci(
      eachLocation,
      currentProps.instanceGranularity,
      currentStructureGroup,
      locationIt.groupCount,
      locationIt.instanceCount,
      loci,
      getVisualLoci,
      apply,
      isMarking,
    );
  }

  function finalize(ctx: VisualContext) {
    if (renderObject) {
      processValues?.(renderObject.values, geometry, currentProps, currentTheme, ctx.graphics);
    }
  }

  return {
    get groupCount() {
      return locationIt ? locationIt.count : 0;
    },
    get renderObject() {
      return locationIt && locationIt.count ? renderObject : undefined;
    },
    get geometryVersion() {
      return geometryVersion;
    },
    createOrUpdate(ctx: VisualContext, theme: Theme, props: PD.Values<P>, structureGroup?: StructureGroup) {
      prepareUpdate(theme, props, structureGroup || currentStructureGroup, ctx.graphics);
      if (updateState.createGeometry) {
        const newGeometry = _createGeometry(
          ctx,
          newStructureGroup.group.units[0],
          newStructureGroup.structure,
          newTheme,
          newProps,
          geometry,
        );
        return Visual.withPreparedGeometry(newGeometry, builder.geometryUtils, newTheme, ctx.graphics, (g) => {
          update(g);
          finalize(ctx);
        });
      }
      return Visual.withPreparedGeometry(geometry, builder.geometryUtils, newTheme, ctx.graphics, () => {
        update();
        finalize(ctx);
      });
    },
    getLoci(pickingId: PickingId) {
      if (!renderObject) return EmptyLoci;
      if (
        resolveInstanceGranularity(currentProps.instanceGranularity, locationIt.groupCount, locationIt.instanceCount)
      ) {
        pickingId = { ...pickingId, groupId: PickingId.Null };
      }
      return getLoci(pickingId, currentStructureGroup, renderObject.id);
    },
    eachLocation(cb: LocationCallback) {
      locationIt.reset();
      while (locationIt.hasNext) {
        const { location, isSecondary } = locationIt.move();
        cb(location, isSecondary);
      }
    },
    mark(loci: Loci, action: MarkerAction) {
      // on its own units only, so that `previousMark` compares (and keeps) just those; a loci of all of them stays as
      // it is, and a whole-structure one arrives as `EveryLoci` (`UnitsRepresentation.mark`)
      const visualLoci = getVisualLoci(loci);
      return visualLoci ? Visual.mark(renderObject, visualLoci, action, lociApply, previousMark) : false;
    },
    setVisibility(visible: boolean) {
      Visual.setVisibility(renderObject, visible);
    },
    setPickable(pickable: boolean) {
      Visual.setPickable(renderObject, pickable);
    },
    setColorOnly(colorOnly: boolean) {
      Visual.setColorOnly(renderObject, colorOnly);
    },
    setTransform(matrix?: Mat4, instanceMatrices?: Float32Array | null) {
      Visual.setTransform(renderObject, matrix, instanceMatrices);
    },
    setOverpaint(overpaint: Overpaint, graphics?: GraphicsContext, smooth = true) {
      const smoothing = smooth ? { geometry, props: currentProps, graphics } : undefined;
      Visual.setOverpaint(renderObject, overpaint, lociApply, true, smoothing);
    },
    setTransparency(transparency: Transparency, graphics?: GraphicsContext, smooth = true) {
      const smoothing = smooth ? { geometry, props: currentProps, graphics } : undefined;
      Visual.setTransparency(renderObject, transparency, lociApply, true, smoothing);
    },
    setEmissive(emissive: Emissive, graphics?: GraphicsContext, smooth = true) {
      const smoothing = smooth ? { geometry, props: currentProps, graphics } : undefined;
      Visual.setEmissive(renderObject, emissive, lociApply, true, smoothing);
    },
    setSubstance(substance: Substance, graphics?: GraphicsContext, smooth = true) {
      const smoothing = smooth ? { geometry, props: currentProps, graphics } : undefined;
      Visual.setSubstance(renderObject, substance, lociApply, true, smoothing);
    },
    setClipping(clipping: Clipping) {
      Visual.setClipping(renderObject, clipping, lociApply, true);
    },
    setWiggle(wiggle: Wiggle, graphics?: GraphicsContext) {
      Visual.setWiggle(renderObject, wiggle, lociApply, true);
    },
    setThemeStrength(strength: {
      overpaint: number;
      transparency: number;
      emissive: number;
      substance: number;
      wiggle: number;
    }) {
      Visual.setThemeStrength(renderObject, strength);
    },
    destroy() {
      dispose?.(geometry);
      if (renderObject) {
        releaseGpuMeshVertices(renderObject.values);
        renderObject.state.disposed = true;
        renderObject = undefined;
      }
    },
    getSignatures: (props: PD.Values<P>, ctx: SignatureContext) => {
      const p = { ...defaultProps, ...props };
      return getGeometrySignatures(builder.geometryUtils, p, ctx, builder.signatureDefines?.(p, ctx.graphics));
    },
    mustRecreate,
    builder,
  };
}

// mesh

export const UnitsMeshParams = { ...StructureMeshParams, ...StructureParams };
export type UnitsMeshParams = typeof UnitsMeshParams;
export interface UnitsMeshVisualBuilder<P extends UnitsMeshParams> extends UnitsVisualBuilder<P, Mesh> {}

export function UnitsMeshVisual<P extends UnitsMeshParams>(
  builder: UnitsMeshVisualBuilder<P>,
  materialId: number,
): UnitsVisual<P> {
  return UnitsVisual<Mesh, P>(
    {
      ...builder,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<P>,
        currentProps: PD.Values<P>,
        newTheme: Theme,
        currentTheme: Theme,
        newStructureGroup: StructureGroup,
        currentStructureGroup: StructureGroup,
      ) => {
        builder.setUpdateState(
          state,
          newProps,
          currentProps,
          newTheme,
          currentTheme,
          newStructureGroup,
          currentStructureGroup,
        );
        if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true;
      },
      geometryUtils: Mesh.Utils,
    },
    materialId,
  );
}

// spheres

export const UnitsSpheresParams = { ...StructureSpheresParams, ...StructureParams };
export type UnitsSpheresParams = typeof UnitsSpheresParams;
export interface UnitsSpheresVisualBuilder<P extends UnitsSpheresParams> extends UnitsVisualBuilder<P, Spheres> {}

export function UnitsSpheresVisual<P extends UnitsSpheresParams>(
  builder: UnitsSpheresVisualBuilder<P>,
  materialId: number,
): UnitsVisual<P> {
  return UnitsVisual<Spheres, P>(
    {
      ...builder,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<P>,
        currentProps: PD.Values<P>,
        newTheme: Theme,
        currentTheme: Theme,
        newStructureGroup: StructureGroup,
        currentStructureGroup: StructureGroup,
      ) => {
        builder.setUpdateState(
          state,
          newProps,
          currentProps,
          newTheme,
          currentTheme,
          newStructureGroup,
          currentStructureGroup,
        );
        if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
      },
      geometryUtils: Spheres.Utils,
    },
    materialId,
  );
}

// cylinders

export const UnitsCylindersParams = { ...StructureCylindersParams, ...StructureParams };
export type UnitsCylindersParams = typeof UnitsCylindersParams;
export interface UnitsCylindersVisualBuilder<P extends UnitsCylindersParams> extends UnitsVisualBuilder<P, Cylinders> {}

export function UnitsCylindersVisual<P extends UnitsCylindersParams>(
  builder: UnitsCylindersVisualBuilder<P>,
  materialId: number,
): UnitsVisual<P> {
  return UnitsVisual<Cylinders, P>(
    {
      ...builder,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<P>,
        currentProps: PD.Values<P>,
        newTheme: Theme,
        currentTheme: Theme,
        newStructureGroup: StructureGroup,
        currentStructureGroup: StructureGroup,
      ) => {
        builder.setUpdateState(
          state,
          newProps,
          currentProps,
          newTheme,
          currentTheme,
          newStructureGroup,
          currentStructureGroup,
        );
        if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
      },
      geometryUtils: Cylinders.Utils,
    },
    materialId,
  );
}

// points

export const UnitsPointsParams = { ...StructurePointsParams, ...StructureParams };
export type UnitsPointsParams = typeof UnitsPointsParams;
export interface UnitsPointVisualBuilder<P extends UnitsPointsParams> extends UnitsVisualBuilder<P, Points> {}

export function UnitsPointsVisual<P extends UnitsPointsParams>(
  builder: UnitsPointVisualBuilder<P>,
  materialId: number,
): UnitsVisual<P> {
  return UnitsVisual<Points, P>(
    {
      ...builder,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<P>,
        currentProps: PD.Values<P>,
        newTheme: Theme,
        currentTheme: Theme,
        newStructureGroup: StructureGroup,
        currentStructureGroup: StructureGroup,
      ) => {
        builder.setUpdateState(
          state,
          newProps,
          currentProps,
          newTheme,
          currentTheme,
          newStructureGroup,
          currentStructureGroup,
        );
        if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
      },
      geometryUtils: Points.Utils,
    },
    materialId,
  );
}

// lines

export const UnitsLinesParams = { ...StructureLinesParams, ...StructureParams };
export type UnitsLinesParams = typeof UnitsLinesParams;
export interface UnitsLinesVisualBuilder<P extends UnitsLinesParams> extends UnitsVisualBuilder<P, Lines> {}

export function UnitsLinesVisual<P extends UnitsLinesParams>(
  builder: UnitsLinesVisualBuilder<P>,
  materialId: number,
): UnitsVisual<P> {
  return UnitsVisual<Lines, P>(
    {
      ...builder,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<P>,
        currentProps: PD.Values<P>,
        newTheme: Theme,
        currentTheme: Theme,
        newStructureGroup: StructureGroup,
        currentStructureGroup: StructureGroup,
      ) => {
        builder.setUpdateState(
          state,
          newProps,
          currentProps,
          newTheme,
          currentTheme,
          newStructureGroup,
          currentStructureGroup,
        );
        if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
      },
      geometryUtils: Lines.Utils,
    },
    materialId,
  );
}

// text

export const UnitsTextParams = { ...StructureTextParams, ...StructureParams };
export type UnitsTextParams = typeof UnitsTextParams;
export interface UnitsTextVisualBuilder<P extends UnitsTextParams> extends UnitsVisualBuilder<P, Text> {}

export function UnitsTextVisual<P extends UnitsTextParams>(
  builder: UnitsTextVisualBuilder<P>,
  materialId: number,
): UnitsVisual<P> {
  return UnitsVisual<Text, P>(
    {
      ...builder,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<P>,
        currentProps: PD.Values<P>,
        newTheme: Theme,
        currentTheme: Theme,
        newStructureGroup: StructureGroup,
        currentStructureGroup: StructureGroup,
      ) => {
        builder.setUpdateState(
          state,
          newProps,
          currentProps,
          newTheme,
          currentTheme,
          newStructureGroup,
          currentStructureGroup,
        );
        if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
        if (newProps.background !== currentProps.background) state.createGeometry = true;
        if (newProps.backgroundMargin !== currentProps.backgroundMargin) state.createGeometry = true;
        if (newProps.tether !== currentProps.tether) state.createGeometry = true;
        if (newProps.tetherLength !== currentProps.tetherLength) state.createGeometry = true;
        if (newProps.tetherBaseWidth !== currentProps.tetherBaseWidth) state.createGeometry = true;
        if (newProps.attachment !== currentProps.attachment) state.createGeometry = true;

        if (newProps.fontFamily !== currentProps.fontFamily) state.createGeometry = true;
        if (newProps.fontQuality !== currentProps.fontQuality) state.createGeometry = true;
        if (newProps.fontStyle !== currentProps.fontStyle) state.createGeometry = true;
        if (newProps.fontVariant !== currentProps.fontVariant) state.createGeometry = true;
        if (newProps.fontWeight !== currentProps.fontWeight) state.createGeometry = true;
      },
      geometryUtils: Text.Utils,
    },
    materialId,
  );
}

// direct-volume

export const UnitsDirectVolumeParams = { ...StructureDirectVolumeParams, ...StructureParams };
export type UnitsDirectVolumeParams = typeof UnitsDirectVolumeParams;
export interface UnitsDirectVolumeVisualBuilder<P extends UnitsDirectVolumeParams>
  extends UnitsVisualBuilder<P, DirectVolume> {}

export function UnitsDirectVolumeVisual<P extends UnitsDirectVolumeParams>(
  builder: UnitsDirectVolumeVisualBuilder<P>,
  materialId: number,
): UnitsVisual<P> {
  return UnitsVisual<DirectVolume, P>(
    {
      ...builder,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<P>,
        currentProps: PD.Values<P>,
        newTheme: Theme,
        currentTheme: Theme,
        newStructureGroup: StructureGroup,
        currentStructureGroup: StructureGroup,
      ) => {
        builder.setUpdateState(
          state,
          newProps,
          currentProps,
          newTheme,
          currentTheme,
          newStructureGroup,
          currentStructureGroup,
        );
        if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true;
      },
      geometryUtils: DirectVolume.Utils,
    },
    materialId,
  );
}

// texture-mesh

export const UnitsTextureMeshParams = { ...StructureTextureMeshParams, ...StructureParams };
export type UnitsTextureMeshParams = typeof UnitsTextureMeshParams;
export interface UnitsTextureMeshVisualBuilder<P extends UnitsTextureMeshParams>
  extends UnitsVisualBuilder<P, TextureMesh> {}

export function UnitsTextureMeshVisual<P extends UnitsTextureMeshParams>(
  builder: UnitsTextureMeshVisualBuilder<P>,
  materialId: number,
): UnitsVisual<P> {
  return UnitsVisual<TextureMesh, P>(
    {
      ...builder,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<P>,
        currentProps: PD.Values<P>,
        newTheme: Theme,
        currentTheme: Theme,
        newStructureGroup: StructureGroup,
        currentStructureGroup: StructureGroup,
      ) => {
        builder.setUpdateState(
          state,
          newProps,
          currentProps,
          newTheme,
          currentTheme,
          newStructureGroup,
          currentStructureGroup,
        );
        if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true;
      },
      geometryUtils: TextureMesh.Utils,
    },
    materialId,
  );
}

// image

export const UnitsImageParams = { ...StructureImageParams, ...StructureParams };
export type UnitsImageParams = typeof UnitsImageParams;
export interface UnitsImageVisualBuilder<P extends UnitsImageParams> extends UnitsVisualBuilder<P, Image> {}

export function UnitsImageVisual<P extends UnitsImageParams>(
  builder: UnitsImageVisualBuilder<P>,
  materialId: number,
): UnitsVisual<P> {
  return UnitsVisual<Image, P>(
    {
      ...builder,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<P>,
        currentProps: PD.Values<P>,
        newTheme: Theme,
        currentTheme: Theme,
        newStructureGroup: StructureGroup,
        currentStructureGroup: StructureGroup,
      ) => {
        builder.setUpdateState(
          state,
          newProps,
          currentProps,
          newTheme,
          currentTheme,
          newStructureGroup,
          currentStructureGroup,
        );
        if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true;
      },
      geometryUtils: Image.Utils,
    },
    materialId,
  );
}
