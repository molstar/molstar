/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import type { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import type { Viewport } from '@molstar/core/util/viewport';

/**
 * The camera the renderer draws with: what it reads of one, so that the renderer (`gl`) does not depend on the
 * canvas's cameras (`canvas3d/camera`'s `ICamera`, which has these and more)
 */
export interface RenderCamera {
  readonly viewport: Viewport;
  readonly view: Mat4;
  readonly projection: Mat4;
  readonly projectionView: Mat4;
  readonly near: number;
  readonly far: number;
  readonly fogNear: number;
  readonly fogFar: number;
  readonly headRotation: Mat4;
  readonly viewEye: Mat4;
  /** the scale of the model (`Camera.scale`), as the renderer's model matrix */
  readonly scale: number;
  readonly state: {
    readonly mode: 'perspective' | 'orthographic';
    /** fog strength, 0: none */
    readonly fog: number;
    readonly target: Vec3;
  };
  /** a part of a larger frame (`Camera.ViewOffset`) */
  readonly viewOffset: { readonly enabled: boolean; readonly offsetX: number; readonly offsetY: number };
  /**
   * The views drawn together (stereo's eyes, XR's views): culled as one, by a frustum enclosing all of them, and their
   * levels of detail by one camera plane; none: the camera's own
   */
  readonly cullView?: CullView;
}

/** The culling of views drawn together: a frustum enclosing all of them and the camera plane of their levels */
export interface CullView {
  readonly frustum: Frustum3D;
  readonly cameraPlane: Plane3D;
}
