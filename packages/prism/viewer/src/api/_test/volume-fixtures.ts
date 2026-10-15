/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4, Tensor } from '@molstar/core/math/linear-algebra';
import type { Grid } from '@molstar/model/model/volume';
import { CCP4Writer } from '@molstar/model/formats/writer/ccp4/ccp4';

// Small data files for the volume and particle tests, made in memory (nothing is read from disk or the network).

/** A 12 × 12 × 12 map of one Gaussian blob at 1.5 Å per voxel, as a CCP4 file named `name`. */
export function mapFile(name = 'test.ccp4'): File {
  const n = 12;
  const space = Tensor.Space([n, n, n], [2, 1, 0], Float32Array);
  const data = new Float32Array(n * n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      for (let k = 0; k < n; k++) {
        const d = (i - 5.5) ** 2 + (j - 5.5) ** 2 + (k - 5.5) ** 2;
        data[space.dataOffset(i, j, k)] = Math.exp(-d / 8);
      }
    }
  }
  const grid: Grid = {
    transform: { kind: 'matrix', matrix: Mat4.fromScaling(Mat4(), [1.5, 1.5, 1.5] as never) },
    cells: Tensor.create(space, Tensor.Data1(data)),
    stats: { min: 0, max: 1, mean: 0, sigma: 1 },
  };
  return new File([CCP4Writer.writeMrc(grid, data)], name);
}

/** Five oriented particles in a CryoET Data Portal NDJSON file. */
export function particlesFile(name = 'particles.ndjson'): File {
  const lines = [0, 1, 2, 3, 4].map((i) =>
    JSON.stringify({
      type: 'orientedPoint',
      location: { x: i * 10, y: i * 5, z: 1 },
      xyz_rotation_matrix: [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
      ],
    }),
  );
  return new File([lines.join('\n')], name);
}
