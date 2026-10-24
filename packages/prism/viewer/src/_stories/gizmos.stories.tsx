/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { twoFaced } from '@molstar/graphics/geo/clip-volumes';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { Box } from '@molstar/graphics/geo/primitive/box';
import { type GizmoBinding, GizmoController, type GizmoItem } from '@molstar/graphics/util/gizmo/controller';
import { type GizmoPose, poseGesture } from '@molstar/graphics/util/gizmo/gestures';
import { type TransformHandlesOptions, transformHandles } from '@molstar/graphics/util/gizmo/transform';
import { story } from '@molstar/prism-ui/testing/index';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { BehaviorSubject } from 'rxjs';
import type { PrismViewer } from '../viewer.js';
import { ExampleFiles, RealViewerStory } from '../testing/real-viewer.js';

// The gizmo library's standard transform handles (`@molstar/graphics/util/gizmo/transform`) on a box over 1CBS: drag an
// arrow to move along X (red), Y (green) or Z (blue), a square to move in its plane, a ring to turn about its axis, the
// center to move in the screen's plane; with scale knobs, past each arrow's tip. A drag runs live and its release keeps
// the box where it is (the story's own binding: nothing in the state; Esc cancels a drag).

const meta: Meta = { title: 'Viewer/Gizmos', parameters: { layout: 'fullscreen' } };
export default meta;

/** A box's walls of `scale`, in its frame (a pose's scale is the box's, the handles' size its own) */
function boxMesh(scale: Vec3) {
  const state = MeshBuilder.createState(64, 32);
  MeshBuilder.addPrimitive(state, Mat4.fromScaling(Mat4(), scale), Box());
  return twoFaced(MeshBuilder.getMesh(state));
}

/** A pose's frame: where it is, how it is turned (not its scale: the handles keep their size) */
function frameOf(pose: GizmoPose) {
  return Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), pose.position), Mat4.fromQuat(Mat4(), pose.rotation));
}

/** The box and its handles on the viewer's canvas, the release keeping the box where it was dragged */
function showTransformHandles(v: PrismViewer, options: TransformHandlesOptions) {
  const canvas = v.plugin.canvas3d;
  if (!canvas) return;
  // the structure's bounds (the canvas' may not be computed yet: nothing drawn)
  const structure = v.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
  const bounds = structure ? Sphere3D.clone(structure.boundary.sphere) : Sphere3D.create(Vec3(), 10);
  const r = Math.max(bounds.radius, 1);
  const items = new BehaviorSubject<readonly GizmoItem<GizmoPose>[]>([
    {
      key: 'box',
      label: 'Box',
      value: { position: bounds.center, rotation: Quat.identity(), scale: Vec3.create(r * 0.8, r * 0.5, r * 0.4) },
    },
  ]);
  const binding: GizmoBinding<GizmoPose> = {
    items,
    layout: (item) => ({
      body: boxMesh(item.value.scale),
      ...transformHandles(r * 0.6, options),
      frame: frameOf(item.value),
    }),
    sameGeometry: (a, b) => Vec3.exactEquals(a.scale, b.scale),
    equal: (a, b) =>
      Vec3.exactEquals(a.position, b.position) &&
      Quat.exactEquals(a.rotation, b.rotation) &&
      Vec3.exactEquals(a.scale, b.scale),
    apply: (from, handle, frame, dx, dy) => poseGesture(from, handle, frame, dx, dy) ?? from,
    // the gizmo itself is what moves: nothing else to preview, nothing to undo
    preview: async () => {},
    commit: async (key, _from, next) => {
      items.next(items.value.map((i) => (i.key === key ? { ...i, value: next } : i)));
      return true;
    },
    cancel: () => {},
  };
  const controller = new GizmoController(canvas, binding);
  // a click on a handle is the gizmo's (no focus, no select)
  v.plugin.managers.interactivity.clickOwners.claim((repr) => controller.owns(repr));
}

function transformStory(options: TransformHandlesOptions): StoryObj {
  return {
    render: () =>
      story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => showTransformHandles(v, options)} />),
  };
}

/** Arrows, plane squares, rings and the center */
export const TransformHandles: StoryObj = transformStory({});

/** With a scale knob past each arrow's tip */
export const WithScaleKnobs: StoryObj = transformStory({ scale: true });
