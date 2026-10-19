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

// 98 directions up to 10k points (the unit spheres' extrema), as for unit spheres (`getBoundary`): the sphere keeps
// 98 extrema for spheres merged from it, e.g. the camera's. Above, 14 directions refined to 98 from their extrema: the
// radius pass still visits every unit sphere's extrema, so the sphere encloses the atoms as before, but it keeps only
// 14 distinct extrema (05-model-structure P2).
const boundaryHelperCoarse = new BoundaryHelper('14');
const boundaryHelperFine = new BoundaryHelper('98');
const FinePointCount = 10_000;

export function computeStructureBoundary(s: Structure): Boundary {
  const { units } = s;

  let pointCount = 0;
  for (let i = 0, _i = units.length; i < _i; i++) {
    const { extrema } = units[i].boundary.sphere;
    pointCount += extrema && extrema.length > 1 ? extrema.length : 1;
  }
  const boundaryHelper = pointCount > FinePointCount ? boundaryHelperCoarse : boundaryHelperFine;
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
