/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Registers every built-in transformer id (for snapshot restore and tools such as state-docs).
// Importing this module has side effects only; it exports nothing.

import '@molstar/plugin/state/formats/cif';
import '@molstar/plugin/state/formats/coordinates/dcd';
import '@molstar/plugin/state/formats/coordinates/lammps';
import '@molstar/plugin/state/formats/coordinates/nctraj';
import '@molstar/plugin/state/formats/coordinates/trr';
import '@molstar/plugin/state/formats/coordinates/xtc';
import '@molstar/plugin/state/formats/particles/em';
import '@molstar/plugin/state/formats/particles/mmcif-assembly';
import '@molstar/plugin/state/formats/particles/ndjson';
import '@molstar/plugin/state/formats/particles/simularium';
import '@molstar/plugin/state/formats/particles/star';
import '@molstar/plugin/state/formats/particles/tbl';
import '@molstar/plugin/state/formats/shape/obj';
import '@molstar/plugin/state/formats/shape/ply';
import '@molstar/plugin/state/formats/shape/vtp';
import '@molstar/plugin/state/formats/topology/prmtop';
import '@molstar/plugin/state/formats/topology/psf';
import '@molstar/plugin/state/formats/topology/top';
import '@molstar/plugin/state/formats/trajectory/cif-core';
import '@molstar/plugin/state/formats/trajectory/gro';
import '@molstar/plugin/state/formats/trajectory/lammps';
import '@molstar/plugin/state/formats/trajectory/mmcif';
import '@molstar/plugin/state/formats/trajectory/mol';
import '@molstar/plugin/state/formats/trajectory/mol2';
import '@molstar/plugin/state/formats/trajectory/pdb';
import '@molstar/plugin/state/formats/trajectory/sdf';
import '@molstar/plugin/state/formats/trajectory/xyz';
import '@molstar/plugin/state/formats/volume/ccp4';
import '@molstar/plugin/state/formats/volume/cube';
import '@molstar/plugin/state/formats/volume/density-server';
import '@molstar/plugin/state/formats/volume/dsn6';
import '@molstar/plugin/state/formats/volume/dx';
import '@molstar/plugin/state/formats/volume/mtz';
import '@molstar/plugin/state/formats/volume/segmentation';
import '@molstar/plugin/state/formats/volume/structure-factors';
import '@molstar/plugin/state/transforms/data/fetch';
import '@molstar/plugin/state/transforms/data/json';
import '@molstar/plugin/state/transforms/misc/group';
import '@molstar/plugin/state/transforms/particles/ops';
import '@molstar/plugin/state/transforms/particles/representation';
import '@molstar/plugin/state/transforms/particles/unitcell';
import '@molstar/plugin/state/transforms/shape/box';
import '@molstar/plugin/state/transforms/shape/representation';
import '@molstar/plugin/state/transforms/structure/animation';
import '@molstar/plugin/state/transforms/structure/bounding-box';
import '@molstar/plugin/state/transforms/structure/effects/clipping';
import '@molstar/plugin/state/transforms/structure/effects/emissive';
import '@molstar/plugin/state/transforms/structure/effects/overpaint';
import '@molstar/plugin/state/transforms/structure/effects/substance';
import '@molstar/plugin/state/transforms/structure/effects/theme-strength';
import '@molstar/plugin/state/transforms/structure/effects/transparency';
import '@molstar/plugin/state/transforms/structure/effects/wiggle';
import '@molstar/plugin/state/transforms/structure/hierarchy';
import '@molstar/plugin/state/transforms/structure/measurement';
import '@molstar/plugin/state/transforms/structure/representation';
import '@molstar/plugin/state/transforms/structure/selection';
import '@molstar/plugin/state/transforms/structure/unitcell';
import '@molstar/plugin/state/transforms/volume/ops';
import '@molstar/plugin/state/transforms/volume/representation';
