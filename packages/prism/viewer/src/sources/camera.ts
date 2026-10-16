/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { Camera } from '@molstar/graphics/canvas3d/camera';
import type { PluginContext } from '@molstar/plugin/context';
import { EMPTY, type Observable, Subject, type Subscription, switchMap } from 'rxjs';

// The camera as a view to follow and to set from outside (api.md §12): how a host keeps viewers' cameras together. Not
// an operation: a pointer-rate stream and a cut, neither of which belongs in the run path, the record or the history
// (`camera.get` and `camera.set` are for people and the agent: rounded, a result, a transition).

/** The view of a camera: where it is, what it looks at, which way is up, and its projection. Full precision, JSON. */
export interface CameraSnapshot {
  readonly mode: Camera.Mode;
  /** Field of view in radians (Mol*'s own). */
  readonly fov: number;
  readonly position: readonly [number, number, number];
  readonly target: readonly [number, number, number];
  readonly up: readonly [number, number, number];
  /** The radius around the target that is in focus: the clipping and fog follow it (a camera clamps it to its scene's). */
  readonly radius: number;
}

export interface SetCameraOptions {
  /**
   * The move is not announced on `changed`, so a host that mirrors the camera onto other viewers does not hear its own
   * setting back. Compared by value, not by timing (`changed`).
   */
  readonly silent?: boolean;
}

export interface ViewerCamera {
  /** The view now; undefined where there is no canvas (headless, or before it is mounted). */
  snapshot(): CameraSnapshot | undefined;
  /**
   * The view, after it moved: the person dragging, a call, a transition's frames, once per drawn frame. Not replaying. It
   * speaks only when the view differs from the one last announced or set, so a view set `silent` (and the frames that
   * redraw it) is never announced, and two viewers following each other settle at once instead of echoing.
   */
  readonly changed: Observable<CameraSnapshot>;
  /**
   * Cuts the camera to a view (no transition). False without a canvas. Never a history step, never recorded.
   */
  set(snapshot: CameraSnapshot, options?: SetCameraOptions): boolean;
}

const triple = (v: Vec3): [number, number, number] => [v[0], v[1], v[2]];

export function snapshotOf(
  state: Pick<Camera.Snapshot, 'mode' | 'fov' | 'position' | 'target' | 'up' | 'radius'>,
): CameraSnapshot {
  return {
    mode: state.mode,
    fov: state.fov,
    position: triple(state.position),
    target: triple(state.target),
    up: triple(state.up),
    radius: state.radius,
  };
}

const near = (a: readonly number[], b: readonly number[], eps: number) =>
  a[0] - b[0] <= eps &&
  b[0] - a[0] <= eps &&
  a[1] - b[1] <= eps &&
  b[1] - a[1] <= eps &&
  a[2] - b[2] <= eps &&
  b[2] - a[2] <= eps;

/** Whether two views are the same to what the eye sees: positions to a millionth of the distance, directions, fov and radius to 1e-6. */
export function sameView(a: CameraSnapshot, b: CameraSnapshot): boolean {
  if (
    a.mode !== b.mode ||
    Math.abs(a.fov - b.fov) > 1e-6 ||
    Math.abs(a.radius - b.radius) > 1e-6 * Math.max(1, a.radius)
  ) {
    return false;
  }
  const reach = Math.max(1, Vec3.distance(Vec3.create(...a.position), Vec3.create(...a.target)));
  return near(a.position, b.position, reach * 1e-6) && near(a.target, b.target, reach * 1e-6) && near(a.up, b.up, 1e-6);
}

/** The camera, and what releases it (with the viewer). */
export function createViewerCamera(plugin: PluginContext): { readonly camera: ViewerCamera; dispose(): void } {
  const out = new Subject<CameraSnapshot>();
  // The view last announced or set: a frame that shows it again is no change.
  let last: CameraSnapshot | undefined;
  const read = () => (plugin.canvas3d ? snapshotOf(plugin.canvas3d.camera.state) : undefined);

  // The canvas appears when the page mounts it; each frame it draws a different view is announced.
  const sub: Subscription = plugin.behaviors.canvas3d.initialized
    .pipe(switchMap((ready) => (ready ? (plugin.canvas3d?.camera?.changed ?? EMPTY) : EMPTY)))
    .subscribe(() => {
      const now = read();
      if (!now || (last && sameView(now, last))) return;
      last = now;
      out.next(now);
    });

  const camera: ViewerCamera = {
    snapshot: read,
    changed: out.asObservable(),
    set(snapshot, options) {
      const canvas = plugin.canvas3d;
      if (!canvas) return false;
      canvas.camera.setState(
        {
          mode: snapshot.mode,
          fov: snapshot.fov,
          position: Vec3.create(...snapshot.position),
          target: Vec3.create(...snapshot.target),
          up: Vec3.create(...snapshot.up),
          radius: snapshot.radius,
        },
        0,
      );
      if (options?.silent) last = snapshot;
      // The canvas draws when asked: the cut shows on the next frame.
      canvas.requestDraw();
      return true;
    },
  };
  return {
    camera,
    dispose() {
      sub.unsubscribe();
      out.complete();
    },
  };
}

/**
 * Keeps the cameras of viewers together: whichever moves, the others cut to its view (silently, so none echoes). Returns
 * what ends it. The first viewer's view is applied to the rest at once when `from` is given.
 */
export function syncCameras(
  viewers: readonly { readonly camera: ViewerCamera }[],
  from?: { readonly camera: ViewerCamera },
): () => void {
  const subs = viewers.map((source) =>
    source.camera.changed.subscribe((view) => {
      for (const other of viewers) if (other !== source) other.camera.set(view, { silent: true });
    }),
  );
  const first = from?.camera.snapshot();
  if (first) for (const other of viewers) if (other !== from) other.camera.set(first, { silent: true });
  return () => {
    for (const s of subs) s.unsubscribe();
  };
}
