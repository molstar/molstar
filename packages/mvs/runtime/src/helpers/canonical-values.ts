/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Objects made once per content, e.g. the color maps of an animation (`animation.md §3`): two of them made for the same
 * place are equal only when they are the same object, so comparisons can stop at identity instead of walking them.
 * Two from different places may have equal content; taking them as different costs only reuse (a new transform with
 * the same version), never a wrong result. Never changed after they are made.
 */
export const CanonicalValues = new WeakSet<object>();

/** Different canonical objects: different content, without walking them. */
export function areDistinctCanonical(a: unknown, b: unknown): boolean {
  return (
    a !== b &&
    typeof a === 'object' &&
    a !== null &&
    typeof b === 'object' &&
    b !== null &&
    CanonicalValues.has(a) &&
    CanonicalValues.has(b)
  );
}
