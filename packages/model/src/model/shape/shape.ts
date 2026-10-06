/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Color } from '@molstar/core/util/color';
import { UUID } from '@molstar/core/util';
import { OrderedSet } from '@molstar/core/data/int';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import type { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';

/** Data shared by a shape and its graphics integrations. `G` is deliberately
 * structural so the model package does not need to own graphics geometry. */
export interface ShapeGeometry {
  readonly kind: string;
  readonly boundingSphere: Sphere3D;
}

export interface Shape<G extends ShapeGeometry = ShapeGeometry> {
  readonly id: UUID;
  readonly name: string;
  readonly sourceData: unknown;
  readonly geometry: G;
  readonly transforms: Mat4[];
  readonly groupCount: number;
  getGroupBoundingSphere?: (groups: ShapeGroup.Loci['groups'], out?: Sphere3D) => Sphere3D;
  getColor(groupId: number, instanceId: number): Color;
  getSize(groupId: number, instanceId: number): number;
  getLabel(groupId: number, instanceId: number): string;
}

export namespace Shape {
  /** Create model-side shape data; rendering consumers should use the graphics Shape.create factory. */
  export function create<G extends ShapeGeometry>(
    name: string,
    sourceData: unknown,
    geometry: G,
    getColor: Shape['getColor'],
    getSize: Shape['getSize'],
    getLabel: Shape['getLabel'],
    groupCount: number,
    transforms?: Mat4[],
  ): Shape<G> {
    return {
      id: UUID.create22(),
      name,
      sourceData,
      geometry,
      transforms: transforms || [Mat4.identity()],
      groupCount,
      getColor,
      getSize,
      getLabel,
    };
  }

  export interface Loci {
    readonly kind: 'shape-loci';
    readonly shape: Shape;
  }
  export function Loci(shape: Shape): Loci {
    return { kind: 'shape-loci', shape };
  }
  export function isLoci(x: any): x is Loci {
    return !!x && x.kind === 'shape-loci';
  }
  export function areLociEqual(a: Loci, b: Loci) {
    return a.shape === b.shape;
  }
  export function isLociEmpty(loci: Loci) {
    return loci.shape.groupCount === 0;
  }
}

export namespace ShapeGroup {
  export interface Location {
    readonly kind: 'group-location';
    shape: Shape;
    group: number;
    instance: number;
  }

  export function Location(shape?: Shape, group = 0, instance = 0): Location {
    return { kind: 'group-location', shape: shape!, group, instance };
  }

  export function isLocation(x: any): x is Location {
    return !!x && x.kind === 'group-location';
  }

  export interface Loci {
    readonly kind: 'group-loci';
    readonly shape: Shape;
    readonly groups: ReadonlyArray<{
      readonly ids: OrderedSet<number>;
      readonly instance: number;
    }>;
  }

  export function Loci(shape: Shape, groups: Loci['groups']): Loci {
    return { kind: 'group-loci', shape, groups: groups as Loci['groups'] };
  }

  export function isLoci(x: any): x is Loci {
    return !!x && x.kind === 'group-loci';
  }

  export function areLociEqual(a: Loci, b: Loci) {
    if (a.shape !== b.shape) return false;
    if (a.groups.length !== b.groups.length) return false;
    for (let i = 0, il = a.groups.length; i < il; ++i) {
      const { ids: idsA, instance: instanceA } = a.groups[i];
      const { ids: idsB, instance: instanceB } = b.groups[i];
      if (instanceA !== instanceB) return false;
      if (!OrderedSet.areEqual(idsA, idsB)) return false;
    }
    return true;
  }

  export function isLociEmpty(loci: Loci) {
    return size(loci) === 0 ? true : false;
  }

  export function size(loci: Loci) {
    let size = 0;
    for (const group of loci.groups) size += OrderedSet.size(group.ids);
    return size;
  }
}
