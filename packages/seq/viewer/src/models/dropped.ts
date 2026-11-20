/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { elements, type Location, location } from '@molstar/seqx-core/query';
import type { Emphasized } from './selection.js';

// What a host's update removed, taken out of session state (_spec/README.md "StoreModel", Host updates): the parts of
// a selection on a node gone (a unit, a transcript) or under it, as when a patch removes that node.

/** Whether an address is one of `gone` or under one. */
export function isUnder(address: string, gone: readonly string[]): boolean {
  for (const g of gone) if (address === g || address.startsWith(`${g}/`)) return true;
  return false;
}

function without(l: Location, gone: readonly string[]): Location {
  const parts = l.parts.filter((p) => !isUnder(p.system, gone));
  const rects = l.rects.filter((r) => !isUnder(r.a.system, gone) && !isUnder(r.b.system, gone));
  return parts.length === l.parts.length && rects.length === l.rects.length ? l : location({ parts, rects });
}

/** The value less what is on or under `gone`; the same object when nothing is (its address dropped with a part). */
export function withoutGone(value: Emphasized, gone: readonly string[]): Emphasized {
  const keys = value.elements.keys.filter((k) => !isUnder(k.owner, gone));
  const loc = without(value.location, gone);
  const shown = without(value.shown, gone);
  if (keys.length === value.elements.keys.length && loc === value.location && shown === value.shown) return value;
  return {
    elements: elements(keys),
    location: loc,
    shown,
    ...(value.style !== undefined ? { style: value.style } : {}),
  };
}
