/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import type { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { type SizeData, createGroupSize, createValueSize } from '@molstar/graphics/geo/geometry/size-data';
import { type Location, NullLocation } from '@molstar/model/model/location';
import { StructureElement, type Structure } from '@molstar/model/model/structure';
import type { Theme } from '@molstar/graphics/theme/theme';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import { PhysicalSizeTheme, getPhysicalRadius } from '@molstar/graphics/theme/size/physical';
import type { WebGLContext } from '@molstar/graphics/gl/webgl/context';
import type { VisualContext } from '../../../visual.js';
import { checkSphereImpostorSupport } from '../../../structure/visual/util/common.js';
import {
  ElementIterator,
  createStructureElementSphereImpostor,
  createStructureElementSphereMesh,
} from '../../../structure/visual/util/element.js';

/** The subset of representation props needed to build a structure target geometry. */
export interface StructureTargetProps {
  type: string;
  sizeFactor: number;
  detail: number;
  tryUseImpostor: boolean;
}

export function structureTargetUseImpostor(props: StructureTargetProps, webgl?: WebGLContext): boolean {
  return props.type !== 'blob-surface' && props.tryUseImpostor && checkSphereImpostorSupport(webgl);
}

/** The scale a size theme contributes to the target's own radii. */
export function getSizeScale(theme: Theme): number {
  const scale = (theme.size.props as { scale?: number }).scale;
  return typeof scale === 'number' ? scale : 1;
}

/**
 * The size theme applied to the elements of a target structure. The representation's size theme
 * assigns sizes to particles, not to the elements of their targets, so only its scale is used;
 * a uniform size theme however overrides the physical radii.
 */
export function getStructureTargetSizeTheme(theme: Theme): SizeTheme<any> {
  return theme.size.granularity === 'uniform' ? theme.size : PhysicalSizeTheme({}, { scale: getSizeScale(theme) });
}

/** Elements of a target structure are never filtered, the target is rendered as given. */
const ElementProps = {
  ignoreHydrogens: false,
  ignoreHydrogensVariant: 'all' as const,
  traceOnly: false,
};

/** Build the geometry for a structure target: impostor spheres or a sphere mesh. */
export function createStructureTargetGeometry(
  ctx: VisualContext,
  structure: Structure,
  theme: Theme,
  props: StructureTargetProps,
  webgl?: WebGLContext,
  existing?: Spheres | Mesh,
): Spheres | Mesh {
  const elementTheme = { color: theme.color, size: getStructureTargetSizeTheme(theme) };
  if (structureTargetUseImpostor(props, webgl)) {
    // Radii are supplied as size data (see `createStructureTargetSizes`), the geometry only holds centers.
    return createStructureElementSphereImpostor(
      ctx,
      structure,
      elementTheme,
      {
        ...ElementProps,
        sizeFactor: props.sizeFactor,
      },
      existing?.kind === 'spheres' ? existing : undefined,
    );
  } else {
    // Meshes carry no size data, so the effective radius has to be baked into the vertices.
    return createStructureElementSphereMesh(
      ctx,
      structure,
      elementTheme,
      {
        ...ElementProps,
        sizeFactor: props.sizeFactor,
        detail: props.detail,
      },
      existing?.kind === 'mesh' ? existing : undefined,
    );
  }
}

/**
 * Per-element sizes of a structure target come from the physical (vdW/coarse) radii of its
 * elements; a uniform size theme overrides them, any other theme only contributes its scale.
 */
export function createStructureTargetSizes(structure: Structure, theme: Theme, sizeData?: SizeData): SizeData {
  if (theme.size.granularity === 'uniform') {
    return createValueSize(theme.size.size(NullLocation), sizeData);
  }
  const scale = getSizeScale(theme);
  const size = (location: Location) => {
    if (!StructureElement.Location.is(location)) return scale;
    return getPhysicalRadius(location.unit, location.element) * scale;
  };
  return createGroupSize(ElementIterator.fromStructure(structure), size, sizeData);
}

/** Whether the sphere geometry must be rebuilt because relevant props changed. */
export function structureTargetGeometryPropsChanged(
  oldProps: StructureTargetProps,
  newProps: StructureTargetProps,
  webgl?: WebGLContext,
): boolean {
  if (structureTargetUseImpostor(newProps, webgl)) return false;
  return newProps.sizeFactor !== oldProps.sizeFactor || newProps.detail !== oldProps.detail;
}

/** Whether toggling impostor usage requires recreating the render object (geometry kind changes). */
export function structureTargetMustRecreate(oldProps: StructureTargetProps, newProps: StructureTargetProps): boolean {
  return oldProps.type !== newProps.type || oldProps.tryUseImpostor !== newProps.tryUseImpostor;
}

/** Whether the impostor sphere size factor changed (impostor spheres carry size as a uniform). */
export function structureTargetSizeFactorChanged(
  oldProps: StructureTargetProps,
  newProps: StructureTargetProps,
  webgl?: WebGLContext,
): boolean {
  return structureTargetUseImpostor(newProps, webgl) && newProps.sizeFactor !== oldProps.sizeFactor;
}
