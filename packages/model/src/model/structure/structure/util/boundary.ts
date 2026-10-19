/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import type { Structure } from '../structure.js';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import type { Boundary } from '@molstar/core/math/geometry/boundary';

const tmpSphere = Sphere3D();

// 14 directions (refined to 98 from their extrema) for every unit count: the radius pass visits every unit sphere's
// extrema, so the sphere encloses the atoms as with 98 directions; unit spheres keep 98 (05-model-structure P2)
const boundaryHelper = new BoundaryHelper('14');

export function computeStructureBoundary(s: Structure): Boundary {
  const { units } = s;

  boundaryHelper.reset();

  for (let i = 0, _i = units.length; i < _i; i++) {
    const u = units[i];
    const invariantBoundary = u.boundary;
    const o = u.conformation.operator;

    if (o.isIdentity) {
      boundaryHelper.includeSphere(invariantBoundary.sphere);
    } else {
      Sphere3D.transform(tmpSphere, invariantBoundary.sphere, o.matrix);
      boundaryHelper.includeSphere(tmpSphere);
    }
  }

  boundaryHelper.finishedIncludeStep();

  for (let i = 0, _i = units.length; i < _i; i++) {
    const u = units[i];
    const invariantBoundary = u.boundary;
    const o = u.conformation.operator;

    if (o.isIdentity) {
      boundaryHelper.radiusSphere(invariantBoundary.sphere);
    } else {
      Sphere3D.transform(tmpSphere, invariantBoundary.sphere, o.matrix);
      boundaryHelper.radiusSphere(tmpSphere);
    }
  }

  return { box: boundaryHelper.getBox(), sphere: boundaryHelper.getSphere() };
}
