/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
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
   * so the others keep theirs. The owner's first new object takes its first old one's place (or a place it left
   * holding), and so on; more go at the end. Fewer: a place before an object of the cell's own is held by a `none`
   * object (it clips nothing), so that one keeps its index; places at the end go. Its variant is kept while it has
   * objects of its own, else the owner's.
   */
  export function replaceOwned(clip: Clip.Props, previous: readonly Object[], next: Clip.Props): Clip.Props {
    const objects: (Object | undefined)[] = [];
    let k = 0;
    let own = 0;
    for (const o of clip.objects) {
      const owners = o.type === 'none' || previous.some((p) => objectsEqual(o, p));
      if (!owners) {
        objects.push(o);
        own++;
      } else objects.push(k < next.objects.length ? next.objects[k++] : undefined);
    }
    for (; k < next.objects.length; k++) objects.push(next.objects[k]);
    // places left at the end go; those before an object that stays are held
    while (objects.length && objects[objects.length - 1] === undefined) objects.pop();
    const held = objects.map((o) => o ?? placeholder());
    return { variant: own ? clip.variant : next.variant, objects: held };
  }

  /** A place held (`replaceOwned`): an object that clips nothing */
  function placeholder(): Object {
    return {
      type: 'none',
      invert: false,
      position: Vec3(),
      rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
      scale: Vec3.create(1, 1, 1),
      transform: Mat4.identity(),
    };
  }
}
