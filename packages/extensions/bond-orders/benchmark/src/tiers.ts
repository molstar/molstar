/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

/** The tiers of plan §5. */
export interface Tier {
  name: string;
  coordinates: 'ideal' | 'model';
  /** H kept in the input: always in tier 0; in the perception tiers, by `--hydrogens kept` (`withHydrogens`) */
  hydrogens: boolean;
  /** run perception (`forceCompute`); tier 0 keeps the CCD's bonds */
  perceive: boolean;
  /** tier C: the linked copy in the representative entry (plan §3.2); `coordinates` then unused */
  entry?: boolean;
}

export const Tiers: Record<string, Tier> = {
  '0': { name: '0', coordinates: 'ideal', hydrogens: true, perceive: false },
  A: { name: 'A', coordinates: 'ideal', hydrogens: false, perceive: true },
  B: { name: 'B', coordinates: 'model', hydrogens: false, perceive: true },
  C: { name: 'C', coordinates: 'model', hydrogens: false, perceive: true, entry: true },
};

/**
 * How a run treats H (plan §3.3): `stripped`, perception on H-suppressed input (the default); `kept`, the CCD's H in the
 * input, to measure the double bonds placed with H present (NAD+ against NADH). Protonation and charge are judged in
 * neither: the compared layers are the same.
 */
export type Hydrogens = 'stripped' | 'kept';

export function withHydrogens(tier: Tier, hydrogens: Hydrogens): Tier {
  return tier.perceive ? { ...tier, hydrogens: hydrogens === 'kept' } : tier;
}
