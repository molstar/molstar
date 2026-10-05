/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';

export { Viewport };

type Viewport = {
  x: number;
  y: number;
  width: number;
  height: number;
};

function Viewport() {
  return Viewport.zero();
}

namespace Viewport {
  export function zero(): Viewport {
    return { x: 0, y: 0, width: 0, height: 0 };
  }
  export function create(x: number, y: number, width: number, height: number): Viewport {
    return { x, y, width, height };
  }
  export function clone(viewport: Viewport): Viewport {
    return { ...viewport };
  }
  export function copy(target: Viewport, source: Viewport): Viewport {
    return Object.assign(target, source);
  }
  export function set(viewport: Viewport, x: number, y: number, width: number, height: number): Viewport {
    viewport.x = x;
    viewport.y = y;
    viewport.width = width;
    viewport.height = height;
    return viewport;
  }

  export function toVec4(v4: Vec4, viewport: Viewport): Vec4 {
    v4[0] = viewport.x;
    v4[1] = viewport.y;
    v4[2] = viewport.width;
    v4[3] = viewport.height;
    return v4;
  }

  export function equals(a: Viewport, b: Viewport) {
    return a.x === b.x && a.y === b.y && a.width === b.width && a.height === b.height;
  }
}
