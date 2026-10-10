/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Fred Ludlow <Fred.Ludlow@astx.com>
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * Parameter definitions of the interactions computation. This module holds no computation code, so UI and
 * state code can describe the interaction params without loading the interaction providers.
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

export const ContactsParams = {
  lineOfSightDistFactor: PD.Numeric(1.0, { min: 0, max: 3, step: 0.1 }),
};
export type ContactsParams = typeof ContactsParams;
export type ContactsProps = PD.Values<ContactsParams>;

export const IonicParams = {
  distanceMax: PD.Numeric(5.0, { min: 0, max: 8, step: 0.1 }),
};
export type IonicParams = typeof IonicParams;

export const PiStackingParams = {
  distanceMax: PD.Numeric(5.5, { min: 1, max: 8, step: 0.1 }),
  offsetMax: PD.Numeric(2.0, { min: 0, max: 4, step: 0.1 }),
  angleDevMax: PD.Numeric(30, { min: 0, max: 180, step: 1 }),
};
export type PiStackingParams = typeof PiStackingParams;

export const CationPiParams = {
  distanceMax: PD.Numeric(6.0, { min: 1, max: 8, step: 0.1 }),
  offsetMax: PD.Numeric(2.0, { min: 0, max: 4, step: 0.1 }),
};
export type CationPiParams = typeof CationPiParams;

export const HalogenBondsParams = {
  distanceMax: PD.Numeric(4.0, { min: 1, max: 5, step: 0.1 }),
  angleMax: PD.Numeric(30, { min: 0, max: 60, step: 1 }),
};
export type HalogenBondsParams = typeof HalogenBondsParams;

export const GeometryParams = {
  distanceMax: PD.Numeric(3.5, { min: 1, max: 5, step: 0.1 }),
  backbone: PD.Boolean(true, { description: 'Include backbone-to-backbone hydrogen bonds' }),
  accAngleDevMax: PD.Numeric(
    45,
    { min: 0, max: 180, step: 1 },
    { description: 'Max deviation from ideal acceptor angle' },
  ),
  ignoreHydrogens: PD.Boolean(false, { description: 'Ignore explicit hydrogens in geometric constraints' }),
  donAngleDevMax: PD.Numeric(
    45,
    { min: 0, max: 180, step: 1 },
    { description: 'Max deviation from ideal donor angle' },
  ),
  accOutOfPlaneAngleMax: PD.Numeric(90, { min: 0, max: 180, step: 1 }),
  donOutOfPlaneAngleMax: PD.Numeric(45, { min: 0, max: 180, step: 1 }),
};
export type GeometryParams = typeof GeometryParams;

export const HydrogenBondsParams = {
  ...GeometryParams,
  water: PD.Boolean(false, { description: 'Include water-to-water hydrogen bonds' }),
  sulfurDistanceMax: PD.Numeric(4.1, { min: 1, max: 5, step: 0.1 }),
};
export type HydrogenBondsParams = typeof HydrogenBondsParams;

export const WeakHydrogenBondsParams = {
  ...GeometryParams,
};
export type WeakHydrogenBondsParams = typeof WeakHydrogenBondsParams;

export const HydrophobicParams = {
  distanceMax: PD.Numeric(4.0, { min: 1, max: 5, step: 0.1 }),
};
export type HydrophobicParams = typeof HydrophobicParams;

export const MetalCoordinationParams = {
  distanceMax: PD.Numeric(3.0, { min: 1, max: 5, step: 0.1 }),
};
export type MetalCoordinationParams = typeof MetalCoordinationParams;
export type MetalCoordinationProps = PD.Values<MetalCoordinationParams>;

export const WaterBridgesParams = {
  backbone: PD.Boolean(true, { description: 'Include backbone hydrogen bonds' }),
  ignoreHydrogens: PD.Boolean(true, { description: 'Ignore explicit hydrogens in geometric constraints' }),
  legDistMin: PD.Numeric(2.5, { min: 1, max: 4, step: 0.1 }, { description: 'Minimum leg distance (Å)' }),
  legDistMax: PD.Numeric(4.1, { min: 1, max: 6, step: 0.1 }, { description: 'Maximum leg distance (Å)' }),
  donAngleDevMax: PD.Numeric(
    80,
    { min: 0, max: 180, step: 1 },
    { description: 'Max deviation from ideal donor angle' },
  ),
  accAngleDevMax: PD.Numeric(
    50,
    { min: 0, max: 180, step: 1 },
    { description: 'Max deviation from ideal acceptor angle' },
  ),
  donOutOfPlaneAngleMax: PD.Numeric(45, { min: 0, max: 180, step: 1 }),
  accOutOfPlaneAngleMax: PD.Numeric(90, { min: 0, max: 180, step: 1 }),
  omegaMin: PD.Numeric(71, { min: 0, max: 180, step: 1 }, { description: 'Minimum A–W–B angle (°)' }),
  omegaMax: PD.Numeric(140, { min: 0, max: 180, step: 1 }, { description: 'Maximum A–W–B angle (°)' }),
};
export type WaterBridgesParams = typeof WaterBridgesParams;
export type WaterBridgesProps = PD.Values<WaterBridgesParams>;

const ContactParamDefinitions = {
  ionic: IonicParams,
  'pi-stacking': PiStackingParams,
  'cation-pi': CationPiParams,
  'halogen-bonds': HalogenBondsParams,
  'hydrogen-bonds': HydrogenBondsParams,
  'weak-hydrogen-bonds': WeakHydrogenBondsParams,
  hydrophobic: HydrophobicParams,
  'metal-coordination': MetalCoordinationParams,
};
type ContactParamDefinitions = typeof ContactParamDefinitions;

function getProvidersParams(defaultOn: string[] = []) {
  const params: {
    [k in keyof ContactParamDefinitions]: PD.Mapped<
      PD.NamedParamUnion<{
        on: PD.Group<ContactParamDefinitions[k]>;
        off: PD.Group<{}>;
      }>
    >;
  } = Object.create(null);

  Object.keys(ContactParamDefinitions).forEach((k) => {
    (params as any)[k] = PD.MappedStatic(
      defaultOn.includes(k) ? 'on' : 'off',
      {
        on: PD.Group(ContactParamDefinitions[k as keyof ContactParamDefinitions]),
        off: PD.Group({}),
      },
      { cycle: true },
    );
  });
  return params;
}
export const ContactProviderParams = getProvidersParams([
  // 'ionic',
  'cation-pi',
  'pi-stacking',
  'hydrogen-bonds',
  'halogen-bonds',
  // 'hydrophobic',
  'metal-coordination',
  // 'weak-hydrogen-bonds',
]);

const BridgeParamDefinitions = {
  'water-bridges': WaterBridgesParams,
};
type BridgeParamDefinitions = typeof BridgeParamDefinitions;

function getBridgeProviderParams(defaultOn: string[] = []) {
  const params: {
    [k in keyof BridgeParamDefinitions]: PD.Mapped<
      PD.NamedParamUnion<{
        on: PD.Group<BridgeParamDefinitions[k]>;
        off: PD.Group<{}>;
      }>
    >;
  } = Object.create(null);

  Object.keys(BridgeParamDefinitions).forEach((k) => {
    (params as any)[k] = PD.MappedStatic(
      defaultOn.includes(k) ? 'on' : 'off',
      {
        on: PD.Group(BridgeParamDefinitions[k as keyof BridgeParamDefinitions]),
        off: PD.Group({}),
      },
      { cycle: true },
    );
  });
  return params;
}
export const BridgeProviderParams = getBridgeProviderParams([]);

export const InteractionsParams = {
  providers: PD.Group(ContactProviderParams, { isFlat: true }),
  bridges: PD.Group(BridgeProviderParams, { isFlat: true }),
  contacts: PD.Group(ContactsParams, { label: 'Advanced Options' }),
};
export type InteractionsParams = typeof InteractionsParams;
export type InteractionsProps = PD.Values<InteractionsParams>;
