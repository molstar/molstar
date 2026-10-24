/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import { type Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { StateTransformer } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import type { Clip } from '@molstar/core/util/clip';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { stringToWords } from '@molstar/core/util/string';
import { type ClipObject, clipEdgesMesh, clipObjectFrame, clipWallsMesh } from '@molstar/graphics/geo/clip-volumes';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import * as ShapeUtils from '@molstar/graphics/geo/shape/shape';
import type { Representation } from '@molstar/graphics/repr/representation';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { registerShapeGeometry } from '@molstar/plugin/state/precompile';

// The volumes a representation's clip objects cut away, drawn in the scene as a shape (its walls translucent, its edges
// solid): a persistent view of the cuts, saved with the state, unlike a clip manager's gizmos. Read from the parent
// representation's live props, so a cut moved, turned or resized (an edit, an undo) is drawn where it is: the
// transform runs again with its parent. By tag (`Clip` objects' `tag`, a clip manager's cut), or all of the parent's.
// Its own representation is not clipped (a shape's `clip` is none by default).

/** Every cell of a clip volumes drawing (the shape and its representation) carries this tag */
export const ClipVolumesTag = 'clip-volumes';

/** The tag of the drawing of a tagged cut: one a cut, wherever it is (`clip-volume-<the cut's tag>`) */
export const clipVolumeTag = (cutTag: string) => `clip-volume-${cutTag}`;

export const ClipVolumesParams = {
  tags: PD.Value<readonly string[]>([], {
    isHidden: true,
    description: "The cuts drawn, by their tags (the clip objects' `tag`); empty: every clip object of the parent",
  }),
  walls: PD.Boolean(true, { description: 'The walls of what is cut away, translucent' }),
  wallColor: PD.Color(Color(0x8fb8ec), { hideIf: (p) => !p.walls }),
  wallOpacity: PD.Numeric(0.15, { min: 0, max: 1, step: 0.01 }, { hideIf: (p) => !p.walls }),
  edges: PD.Boolean(true, { description: 'Its edges, solid' }),
  edgeColor: PD.Color(Color(0x9db8e0), { hideIf: (p) => !p.edges }),
  edgeWidth: PD.Numeric(
    0,
    { min: 0, max: 5, step: 0.05 },
    { description: 'Å; 0: in proportion to what it cuts', hideIf: (p) => !p.edges },
  ),
  extent: PD.Numeric(
    1,
    { min: 0.25, max: 4, step: 0.05 },
    {
      description: "How far a plane's disk and a cone reach, relative to what it cuts (1: as far as it; less: shorter)",
    },
  ),
};
export type ClipVolumesParams = typeof ClipVolumesParams;
export type ClipVolumesProps = PD.Values<ClipVolumesParams>;

/** What a clip volumes shape is built from: the objects drawn, what they cut, how */
interface ClipVolumesData {
  readonly objects: readonly ClipObject[];
  readonly around: Sphere3D;
  readonly props: ClipVolumesProps;
}

type ClipVolumesShape3D = typeof ClipVolumesShape3D;
const ClipVolumesShape3D = PluginStateTransform.BuiltIn({
  name: 'clip-volumes-shape-3d',
  display: 'Clip Volumes',
  from: [
    SO.Molecule.Structure.Representation3D,
    SO.Volume.Representation3D,
    SO.Shape.Representation3D,
    SO.Particle.Representation3D,
  ],
  to: SO.Shape.Provider,
  params: ClipVolumesParams,
})({
  canAutoUpdate() {
    return true;
  },
  isApplicable(a) {
    return !!clipOf(a.data.repr);
  },
  apply({ a, params }) {
    return Task.create('Clip Volumes', async () => {
      const objects = (clipOf(a.data.repr)?.objects ?? []).filter(
        (o) => !params.tags.length || (!!o.tag && params.tags.includes(o.tag)),
      );
      // everything the mesh is built from is the shape's data: the representation keeps the `getShape` it was made
      // with and is given new data on an update (a closure over one mesh would draw only that one)
      const data: ClipVolumesData = { objects, around: boundsOf(a.data.repr), props: params };
      const label = objects.length === 1 ? `Cut: ${stringToWords(objects[0].type)}` : 'Cuts';
      return new SO.Shape.Provider(
        {
          label,
          data,
          params: Mesh.Params,
          // into the previous shape's mesh: the representation's buffers are that mesh's (an update gives new counts,
          // not new buffers), so a new mesh would be drawn with the old one's vertices
          getShape: (_, d: ClipVolumesData, __, previous) =>
            ShapeUtils.create(
              label,
              d,
              clipVolumesMesh(d.objects, d.around, d.props, previous?.geometry as Mesh | undefined),
              (g) => (g % 2 ? d.props.edgeColor : d.props.wallColor),
              () => 1,
              (g) => `${stringToWords(d.objects[g >> 1]?.type ?? '')} (${g % 2 ? 'edges' : 'walls'})`,
              undefined,
              Math.max(2 * d.objects.length, 1),
            ),
          geometryUtils: Mesh.Utils,
          // the walls translucent, the edges solid (applied by the shape representation, per group)
          transparency: (g) => (g % 2 ? 0 : 1 - params.wallOpacity),
        },
        { label, description: objects.length > 1 ? `${objects.length} cuts` : undefined },
      );
    });
  },
  update() {
    // the parent changed (a cut moved) or the params: drawn again
    return StateTransformer.UpdateResult.Recreate;
  },
});

export { ClipVolumesShape3D };

registerShapeGeometry(ClipVolumesShape3D, () => [Mesh.Utils]);

function clipOf(repr: Representation.Any | undefined): Clip.Props | undefined {
  return (repr?.props as { clip?: Clip.Props } | undefined)?.clip;
}

/** What a representation draws, its render objects' bounds together; a unit sphere when it draws nothing yet */
function boundsOf(repr: Representation.Any): Sphere3D {
  let out: Sphere3D | undefined;
  for (const ro of repr.renderObjects) {
    const s = ro.values.boundingSphere?.ref.value as Sphere3D | undefined;
    if (!s || s.radius <= 0) continue;
    out = out ? Sphere3D.expandBySphere(Sphere3D(), out, s) : Sphere3D.clone(s);
  }
  return out ?? Sphere3D.create(Vec3(), 1);
}

/** The objects' walls and edges in the world: object i's walls group 2i, its edges 2i + 1 */
export function clipVolumesMesh(
  objects: readonly ClipObject[],
  around: Sphere3D,
  params: ClipVolumesProps,
  into?: Mesh,
): Mesh {
  // relative to what it cuts: it stays in proportion when that changes
  const reach = Math.max(around.radius, 1) * (params.extent > 0 ? params.extent : 1);
  const width = params.edgeWidth > 0 ? params.edgeWidth : reach * 0.006;
  const state = MeshBuilder.createState(4096, 2048, into);
  objects.forEach((o, i) => {
    const frame: Mat4 = clipObjectFrame(o, around);
    if (params.walls) {
      state.currentGroup = 2 * i;
      MeshBuilder.addMesh(state, frame, clipWallsMesh(o, reach));
    }
    if (params.edges) {
      state.currentGroup = 2 * i + 1;
      MeshBuilder.addMesh(state, frame, clipEdgesMesh(o, reach, width));
    }
  });
  return MeshBuilder.getMesh(state);
}
