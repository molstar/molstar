/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Fred Ludlow <Fred.Ludlow@astx.com>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 *
 * based in part on NGL (https://github.com/arose/ngl)
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Structure, Unit, StructureElement } from '@molstar/model/model/structure';
import { calcAngles } from '../chemistry/geometry.js';
import { type FeaturesBuilder, Features } from './features.js';
import type { ElementSymbol } from '@molstar/model/model/structure/model/types';
import { typeSymbol, eachBondedAtom } from '../chemistry/util.js';
import { Elements } from '@molstar/model/model/structure/model/properties/atomic/types';
import { degToRad } from '@molstar/core/math/misc';
import { FeatureType, FeatureGroup, InteractionType } from './common.js';
import type { ContactProvider } from './contacts.js';
import { HalogenBondsParams } from './params.js';

type HalogenBondsProps = PD.Values<HalogenBondsParams>;

const halBondElements = [Elements.CL, Elements.BR, Elements.I, Elements.AT] as ElementSymbol[];

/**
 * Halogen bond donors (X-C, with X one of Cl, Br, I or At) not F!
 */
function addUnitHalogenDonors(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
  const { elements } = unit;
  const { x, y, z } = unit.model.atomicConformation;

  for (let i = 0 as StructureElement.UnitIndex, il = elements.length; i < il; ++i) {
    const element = typeSymbol(unit, i);
    if (halBondElements.includes(element)) {
      builder.add(FeatureType.HalogenDonor, FeatureGroup.None, x[elements[i]], y[elements[i]], z[elements[i]], i);
    }
  }
}

const X = [Elements.N, Elements.O, Elements.S] as ElementSymbol[];
const Y = [Elements.C, Elements.N, Elements.P, Elements.S] as ElementSymbol[];

/**
 * Halogen bond acceptors (Y-{O|N|S}, with Y=C,P,N,S)
 */
function addUnitHalogenAcceptors(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
  const { elements } = unit;
  const { x, y, z } = unit.model.atomicConformation;

  for (let i = 0 as StructureElement.UnitIndex, il = elements.length; i < il; ++i) {
    const element = typeSymbol(unit, i);
    if (X.includes(element)) {
      let flag = false;
      eachBondedAtom(structure, unit, i, (unitB, indexB) => {
        if (Y.includes(typeSymbol(unitB, indexB))) {
          flag = true;
        }
      });
      if (flag) {
        builder.add(FeatureType.HalogenAcceptor, FeatureGroup.None, x[elements[i]], y[elements[i]], z[elements[i]], i);
      }
    }
  }
}

function isHalogenBond(ti: FeatureType, tj: FeatureType) {
  return (
    (ti === FeatureType.HalogenAcceptor && tj === FeatureType.HalogenDonor) ||
    (ti === FeatureType.HalogenDonor && tj === FeatureType.HalogenAcceptor)
  );
}

// http://www.pnas.org/content/101/48/16789.full
const OptimalHalogenAngle = degToRad(180); // adjusted from 165 to account for spherical statistics
const OptimalAcceptorAngle = degToRad(120);

function getOptions(props: HalogenBondsProps) {
  return {
    angleMax: degToRad(props.angleMax),
  };
}
type Options = ReturnType<typeof getOptions>;

function testHalogenBond(
  structure: Structure,
  infoA: Features.Info,
  infoB: Features.Info,
  opts: Options,
): InteractionType | undefined {
  const typeA = infoA.types[infoA.feature];
  const typeB = infoB.types[infoB.feature];

  if (!isHalogenBond(typeA, typeB)) return;

  const [don, acc] = typeA === FeatureType.HalogenDonor ? [infoA, infoB] : [infoB, infoA];

  const donIndex = don.members[don.offsets[don.feature]];
  const accIndex = acc.members[acc.offsets[acc.feature]];

  const [halogenAngles] = calcAngles(structure, don.unit, donIndex, acc.unit, accIndex);
  // Singly bonded halogen only (not bromide ion for example)
  if (halogenAngles.length !== 1) return;
  if (OptimalHalogenAngle - halogenAngles[0] > opts.angleMax) return;

  const [acceptorAngles] = calcAngles(structure, acc.unit, accIndex, don.unit, donIndex);
  // Angle must be defined. Excludes water as acceptor. Debatable
  if (acceptorAngles.length === 0) return;
  if (acceptorAngles.some((acceptorAngle) => OptimalAcceptorAngle - acceptorAngle > opts.angleMax)) return;

  return InteractionType.HalogenBond;
}

//

export const HalogenDonorProvider = Features.Provider([FeatureType.HalogenDonor], addUnitHalogenDonors);
export const HalogenAcceptorProvider = Features.Provider([FeatureType.HalogenAcceptor], addUnitHalogenAcceptors);

export const HalogenBondsProvider: ContactProvider<HalogenBondsParams> = {
  name: 'halogen-bonds',
  params: HalogenBondsParams,
  createTester: (props: HalogenBondsProps) => {
    const opts = getOptions(props);
    return {
      maxDistance: props.distanceMax,
      requiredFeatures: new Set([FeatureType.HalogenDonor, FeatureType.HalogenAcceptor]),
      getType: (structure, infoA, infoB) => testHalogenBond(structure, infoA, infoB, opts),
    };
  },
};
