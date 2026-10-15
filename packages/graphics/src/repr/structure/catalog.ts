/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { namedCatalog } from '@molstar/graphics/util/named-catalog';
import { BallAndStickRepresentationProvider } from './representation/ball-and-stick.js';
import { BlobSurfaceRepresentationProvider } from './representation/blob-surface.js';
import { CarbohydrateRepresentationProvider } from './representation/carbohydrate.js';
import { CartoonRepresentationProvider } from './representation/cartoon.js';
import { EllipsoidRepresentationProvider } from './representation/ellipsoid.js';
import { GaussianSurfaceRepresentationProvider } from './representation/gaussian-surface.js';
import { LabelRepresentationProvider } from './representation/label.js';
import { MolecularSurfaceRepresentationProvider } from './representation/molecular-surface.js';
import { OrientationRepresentationProvider } from './representation/orientation.js';
import { PointRepresentationProvider } from './representation/point.js';
import { PuttyRepresentationProvider } from './representation/putty.js';
import { SpacefillRepresentationProvider } from './representation/spacefill.js';
import { LineRepresentationProvider } from './representation/line.js';
import { GaussianVolumeRepresentationProvider } from './representation/gaussian-volume.js';
import { BackboneRepresentationProvider } from './representation/backbone.js';
import { PolyhedronRepresentationProvider } from './representation/polyhedron.js';
import { PlaneRepresentationProvider } from './representation/plane.js';

export const BuiltInStructureRepresentations = namedCatalog({
  cartoon: CartoonRepresentationProvider,
  backbone: BackboneRepresentationProvider,
  'ball-and-stick': BallAndStickRepresentationProvider,
  'blob-surface': BlobSurfaceRepresentationProvider,
  carbohydrate: CarbohydrateRepresentationProvider,
  ellipsoid: EllipsoidRepresentationProvider,
  'gaussian-surface': GaussianSurfaceRepresentationProvider,
  'gaussian-volume': GaussianVolumeRepresentationProvider,
  label: LabelRepresentationProvider,
  line: LineRepresentationProvider,
  'molecular-surface': MolecularSurfaceRepresentationProvider,
  orientation: OrientationRepresentationProvider,
  plane: PlaneRepresentationProvider,
  point: PointRepresentationProvider,
  putty: PuttyRepresentationProvider,
  spacefill: SpacefillRepresentationProvider,
  polyhedron: PolyhedronRepresentationProvider,
});
