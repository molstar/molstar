/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type CapabilityKey } from './capabilities.js';
import { type VisualTestCase } from './test-case.js';

export interface Variant<T> {
  readonly name: string;
  readonly value: T;
  readonly requires?: readonly CapabilityKey[];
  readonly tags?: readonly string[];
}

export function variants<T>(
  values: Record<string, T>,
  extra?: Record<string, Pick<Variant<T>, 'requires' | 'tags'>>,
): Variant<T>[] {
  return Object.keys(values).map((name) => ({ name, value: values[name], ...extra?.[name] }));
}

/** Cartesian product of two variant lists; `make` may return undefined to prune a combination */
export function cartesian<A, B>(
  a: readonly Variant<A>[],
  b: readonly Variant<B>[],
  make: (a: Variant<A>, b: Variant<B>) => VisualTestCase | undefined,
): VisualTestCase[] {
  const cases: VisualTestCase[] = [];
  for (const va of a) {
    for (const vb of b) {
      const c = make(va, vb);
      if (c) cases.push(c);
    }
  }
  return cases;
}

export function mergeRequires(...lists: (readonly CapabilityKey[] | undefined)[]): CapabilityKey[] {
  const set = new Set<CapabilityKey>();
  for (const l of lists) l?.forEach((k) => set.add(k));
  return Array.from(set);
}
