/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { auditTime, distinctUntilChanged, EMPTY, map, merge, type Observable, of, startWith, switchMap } from 'rxjs';
import { type CameraValue, cameraValueOf } from '../../api/camera.js';

// The camera as the Export card's Viewport section reads it (v5's camera info, parity.md; state.md §4, a feature's
// source): cold, so it follows the canvas's redraws only while the section is shown and bound to it (prism AGENTS.md:
// hidden UI does no updates), at most a few times a second.

const EveryMs = 150;

const same3 = (a: readonly number[], b: readonly number[]) => a[0] === b[0] && a[1] === b[1] && a[2] === b[2];

const sameCamera = (a: CameraValue, b: CameraValue) =>
  a.mode === b.mode &&
  a.fov === b.fov &&
  same3(a.position, b.position) &&
  same3(a.target, b.target) &&
  same3(a.up, b.up) &&
  a.distance === b.distance &&
  a.radius === b.radius &&
  a.sceneRadius === b.sceneRadius &&
  a.showsAll === b.showsAll;

export function cameraUpdates(plugin: PluginContext): Observable<CameraValue | undefined> {
  return plugin.behaviors.canvas3d.initialized.pipe(
    switchMap((ready) => {
      const canvas = ready ? plugin.canvas3d : undefined;
      if (!canvas) return of(undefined);
      return canvas.didDraw.pipe(
        auditTime(EveryMs),
        startWith(0),
        map(() => cameraValueOf(canvas)),
        // A redraw with the camera still (accumulating samples, an animation, a hover) changes nothing to show.
        distinctUntilChanged(sameCamera),
      );
    }),
  );
}

/**
 * What changes the picture, for the Image section's preview (cold: followed while the section shows): the scene (the
 * state tree, a representation, its colors, a visibility) and the view (the camera moving, or the canvas resized, which
 * changes its projection). Not the canvas's redraws: an export's render starts the loop again with one, and the preview
 * would follow the exports (a preview renders on a camera of its own, so it moves nothing here).
 */
export function sceneChanges(plugin: PluginContext): Observable<unknown> {
  return merge(
    plugin.state.data.events.changed,
    plugin.state.data.events.cell.stateUpdated,
    plugin.behaviors.canvas3d.initialized.pipe(
      switchMap((ready) => (ready && plugin.canvas3d ? plugin.canvas3d.camera.changed : EMPTY)),
    ),
  );
}
