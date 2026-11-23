/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The extent of Keys (_spec/README.md "Queries", Keys): an annotation's items by their ranges; a structure group's
// (`groups`' keys, _spec/README.md "Queries", Language) by its owner: a system's position (a residue, an atom), a
// complex's unit by its id, or its entity's units.

import { type Ctx, fail, stage, tick } from './form-kit.js';
import type { Keys } from './keys.js';
import { isPosition, type Location, location } from './location.js';
import type { QueryModel } from './model.js';

/**
 * A structure group's ranges by its owner, or undefined when the owner is no system or complex: on a system, each key
 * a position; on a complex, each key a unit id (its whole system) or an entity id (its units', whole).
 */
export function structureKeyRanges(
  model: QueryModel,
  owner: string,
  keys: ReadonlySet<string | number>,
): { system: string; ranges: number[] }[] | undefined {
  const units = model.units?.(owner);
  if (units) {
    const out: { system: string; ranges: number[] }[] = [];
    for (const u of units) {
      if (!keys.has(u.id) && !(u.entity !== undefined && keys.has(u.entity))) continue;
      const s = model.system(u.address);
      if (s) out.push({ system: u.address, ranges: [0, s.length] });
    }
    return out;
  }
  const s = model.system(owner);
  if (!s) return undefined;
  const ranges: number[] = [];
  for (const k of [...keys].filter((k): k is number => typeof k === 'number').sort((x, y) => x - y))
    if (isPosition(k) && k < s.length) ranges.push(k, k + 1);
  return [{ system: owner, ranges }];
}

/**
 * The extent of keys: an annotation's items by their ranges on their systems (a pair's two ends on two), a structure
 * group's by what it names (`structureKeyRanges`). A key whose owner is neither is `seq.query.form.prefix`.
 */
export function keysExtent(ctx: Ctx, e: Keys): Location {
  const byOwner = new Map<string, Set<string | number>>();
  for (const k of e.keys) {
    let s = byOwner.get(k.owner);
    if (!s) byOwner.set(k.owner, (s = new Set()));
    s.add(k.key);
  }
  const parts: { system: string; ranges: number[] }[] = [];
  for (const [owner, keys] of byOwner) {
    const a = ctx.model.annotation(owner);
    if (!a) {
      const s = structureKeyRanges(ctx.model, owner, keys);
      if (!s)
        return fail('seq.query.form.prefix', 'A key’s owner is no annotation, system or complex', { address: owner });
      parts.push(...s);
      continue;
    }
    const read = a.items();
    stage(ctx, read.pending);
    // by system: a pair's two ends may be on two
    const bySystem = new Map<string, number[]>([[a.target, []]]);
    for (const it of read.items) {
      tick(ctx);
      if (!keys.has(it.key)) continue;
      if (!isPosition(it.start) || !isPosition(it.end) || it.start > it.end)
        fail('seq.query.position.invalid', 'An item’s extent is not 32-bit positions with start ≤ end', {
          address: owner,
        });
      const system = it.system ?? a.target;
      let ranges = bySystem.get(system);
      if (!ranges) bySystem.set(system, (ranges = []));
      ranges.push(it.start, it.end);
    }
    for (const [system, ranges] of bySystem) parts.push({ system, ranges });
  }
  return location({ parts });
}
