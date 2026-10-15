/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Volume } from '@molstar/model/model/volume';
import type { Theme } from '@molstar/graphics/theme/theme';
import type { VisualContext } from '../../../visual.js';
import type { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import type { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { createVolumeIsosurfaceMesh, type VolumeIsosurfaceProps } from '../../../volume/isosurface.js';
import { createVolumeSphereImpostor, createVolumeSphereMesh, type VolumeSphereProps } from '../../../volume/dot.js';
import { checkSphereImpostorSupport } from '../../../structure/visual/util/common.js';

/** The subset of representation props needed to build a volume target geometry. */
export type VolumeTargetProps = VolumeIsosurfaceProps &
  Pick<VolumeSphereProps, 'perturbPositions' | 'sizeFactor' | 'detail' | 'lodLevels'> & {
    type: string;
    tryUseImpostor: boolean;
  };

export function volumeTargetUseImpostor(props: VolumeTargetProps, graphics?: GraphicsContext): boolean {
  return props.tryUseImpostor && checkSphereImpostorSupport(graphics);
}

/**
 * The geometry for a volume target is either an isosurface mesh of the volume, built on the CPU
 * (the GPU/texture-mesh variant cannot be instanced), or a dot per above-threshold cell. It is
 * instanced at each particle position.
 */
export function createVolumeTargetGeometry(
  ctx: VisualContext,
  volume: Volume,
  theme: Theme,
  props: VolumeTargetProps,
  graphics?: GraphicsContext,
  existing?: Spheres | Mesh,
): Spheres | Mesh | Promise<Mesh> {
  if (props.type === 'dot') {
    const dotProps = props as unknown as VolumeSphereProps;
    return volumeTargetUseImpostor(props, graphics)
      ? createVolumeSphereImpostor(ctx, volume, 0, theme, dotProps, existing?.kind === 'spheres' ? existing : undefined)
      : createVolumeSphereMesh(ctx, volume, 0, theme, dotProps, existing?.kind === 'mesh' ? existing : undefined);
  }
  return createVolumeIsosurfaceMesh(ctx, volume, 0, theme, props, existing?.kind === 'mesh' ? existing : undefined);
}

/** Whether the volume geometry must be rebuilt because relevant props changed. */
export function volumeTargetGeometryPropsChanged(
  volume: Volume,
  oldProps: VolumeTargetProps,
  newProps: VolumeTargetProps,
  graphics?: GraphicsContext,
): boolean {
  if (!Volume.IsoValue.areSame(oldProps.isoValue, newProps.isoValue, volume.grid.stats)) return true;
  if (newProps.type === 'dot') {
    return (
      newProps.perturbPositions !== oldProps.perturbPositions ||
      (!volumeTargetUseImpostor(newProps, graphics) &&
        (newProps.sizeFactor !== oldProps.sizeFactor || newProps.detail !== oldProps.detail))
    );
  }
  return newProps.wrap !== oldProps.wrap || newProps.floodfill !== oldProps.floodfill;
}

/** Whether switching the geometry kind requires recreating the render object. */
export function volumeTargetMustRecreate(oldProps: VolumeTargetProps, newProps: VolumeTargetProps): boolean {
  return oldProps.type !== newProps.type || oldProps.tryUseImpostor !== newProps.tryUseImpostor;
}

/** Whether the impostor dot size factor changed (impostor spheres carry size as a uniform). */
export function volumeTargetSizeFactorChanged(
  oldProps: VolumeTargetProps,
  newProps: VolumeTargetProps,
  graphics?: GraphicsContext,
): boolean {
  return (
    newProps.type === 'dot' &&
    volumeTargetUseImpostor(newProps, graphics) &&
    newProps.sizeFactor !== oldProps.sizeFactor
  );
}
