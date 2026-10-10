/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { all, current } from './basic.js';
import { disulfideBridges, nosBridges } from './bond.js';
import {
  complement,
  covalentlyBonded,
  covalentlyBondedComponent,
  covalentlyOrMetallicBonded,
  surroundingAtoms,
  surroundingLigands,
  surroundings,
  wholeResidues,
} from './manipulate.js';
import { aromaticRing, nonStandardPolymer, ring } from './residue.js';
import { backbone, beta, helix, sidechain, sidechainWithTrace, trace } from './structure.js';
import {
  branched,
  branchedConnectedOnly,
  branchedPlusConnected,
  coarse,
  connectedOnly,
  ion,
  ligand,
  ligandConnectedOnly,
  ligandPlusConnected,
  lipid,
  nucleic,
  polymer,
  protein,
  water,
} from './type.js';

export const StructureSelectionQueries = {
  all,
  current,
  polymer,
  trace,
  backbone,
  sidechain,
  sidechainWithTrace,
  protein,
  nucleic,
  helix,
  beta,
  water,
  ion,
  lipid,
  branched,
  branchedPlusConnected,
  branchedConnectedOnly,
  ligand,
  ligandPlusConnected,
  ligandConnectedOnly,
  connectedOnly,
  disulfideBridges,
  nosBridges,
  nonStandardPolymer,
  coarse,
  ring,
  aromaticRing,
  surroundings,
  surroundingLigands,
  surroundingAtoms,
  complement,
  covalentlyBonded,
  covalentlyOrMetallicBonded,
  covalentlyBondedComponent,
  wholeResidues,
};
