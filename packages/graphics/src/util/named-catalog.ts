/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Compile-time check that every catalog key equals the `name` of the provider stored under it. */
export function namedCatalog<T extends { [K in keyof T]: { readonly name: K } }>(t: T): T {
  return t;
}
