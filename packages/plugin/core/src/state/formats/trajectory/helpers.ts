/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import type { Trajectory } from '@molstar/model/model/structure';

export function trajectoryProps(trajectory: Trajectory) {
  const first = trajectory.representative;
  return {
    label: `${first.entry}`,
    description: `${trajectory.frameCount} model${trajectory.frameCount === 1 ? '' : 's'}`,
  };
}
