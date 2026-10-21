/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StateObjectCell } from '@molstar/core/state';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Clip } from '@molstar/core/util/clip';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';

/**
 * The clip objects a representation cell holds in its params. Every geometry has `clip` (`BaseGeometry.Params`); a
 * representation transform keeps it at `type.params.clip` (structure, volume, particles: the type's params) or at
 * `clip` (a shape: the provider's own geometry params). The state tree is the record of what is clipped: these read
 * and write it there, keeping what they do not touch (other objects, the variant).
 */
export namespace RepresentationClip {
  export type Object = Clip.Props['objects'][number];

  /** The clip a representation's params carry, or undefined when its geometry has none */
  export function get(params: unknown): Clip.Props | undefined {
    const p = params as { type?: { params?: { clip?: Clip.Props } }; clip?: Clip.Props } | undefined;
    return p?.type?.params?.clip ?? p?.clip;
  }

  /** The params with `clip` in place of the one they carry (the same path); unchanged when they carry none */
  export function set<P>(params: P, clip: Clip.Props): P {
    const p = params as { type?: { name: string; params?: { clip?: Clip.Props } }; clip?: Clip.Props };
    if (p.type?.params?.clip) return { ...p, type: { ...p.type, params: { ...p.type.params, clip } } } as P;
    if (p.clip) return { ...p, clip } as P;
    return params;
  }

  /** A representation cell whose geometry can be clipped (its params carry `clip`) */
  export function isClippable(cell: StateObjectCell): boolean {
    return !!cell.obj && SO.isRepresentation3D(cell.obj) && !!get(cell.transform.params);
  }

  export function objectsEqual(a: Object, b: Object): boolean {
    return (
      a.type === b.type &&
      a.invert === b.invert &&
      Vec3.exactEquals(a.position, b.position) &&
      Vec3.exactEquals(a.scale, b.scale) &&
      Vec3.exactEquals(a.rotation.axis, b.rotation.axis) &&
      a.rotation.angle === b.rotation.angle &&
      Mat4.areEqual(a.transform, b.transform, 0) &&
      (a.tag ?? '') === (b.tag ?? '')
    );
  }

  /**
   * The same cut, maybe moved or resized: a tagged object (a clip manager's) by its tag, an untagged one by being
   * equal in every value
   */
  export function sameCut(a: Object, b: Object): boolean {
    const ta = a.tag ?? '';
    const tb = b.tag ?? '';
    return ta || tb ? ta === tb : objectsEqual(a, b);
  }

  export function areEqual(a: Clip.Props, b: Clip.Props): boolean {
    return PD.areEqual(Clip.Params, a, b);
  }

  /**
   * A cell's clip with the objects one owner put there (`previous`, the component manager's global objects) replaced by
   * its new ones (`next`) in their places: an object's index is what a clipping layer's group means (bit i, object i),
   * so the others keep theirs. The owner's first new object takes its first old one's place, and so on; more go at
   * the end, fewer leave their places out. Its variant is kept while it has objects of its own, else the owner's.
   */
  export function replaceOwned(clip: Clip.Props, previous: readonly Object[], next: Clip.Props): Clip.Props {
    const objects: Object[] = [];
    let k = 0;
    let own = 0;
    for (const o of clip.objects) {
      if (!previous.some((p) => objectsEqual(o, p))) {
        objects.push(o);
        own++;
      } else if (k < next.objects.length) objects.push(next.objects[k++]);
    }
    for (; k < next.objects.length; k++) objects.push(next.objects[k]);
    return { variant: own ? clip.variant : next.variant, objects };
  }
}
