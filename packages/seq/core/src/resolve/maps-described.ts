/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Node, parseAddress, printAddress, type Value } from '@molstar/core/repr/tree';

// A spec map with an end a source describes (_spec/README.md "Overlay", Described ends; Structures, `toUnits`): a map
// onto a unit of a complex a source's description gives (SIFTS per chain) names it by its dynamic address,
// `<source>/(node :key <key>)/<unit>`. Resolve keeps it aside, never a missing ref; the overlay resolves it into its
// steps once the description has arrived, and until then it is pending (EC-169), never a map that pairs nothing.

/** A spec map with an end a source describes: its node and its two ends' addresses (canonical). */
export interface DescribedMap {
  readonly address: string;
  readonly node: Node;
  readonly from: string;
  readonly to: string;
}

/**
 * The canonical address `v` names when it is a node a spec source describes (`<source>/(node :key <key>)…`, `isSource`
 * saying the source is the spec's); undefined otherwise.
 */
export function describedAddress(v: Value | undefined, isSource: (address: string) => boolean): string | undefined {
  if (typeof v !== 'string') return undefined;
  const p = parseAddress(v);
  if (p.issues.length || !p.pattern || p.pattern.subtree) return undefined;
  const segments = p.pattern.segments;
  const i = segments.findIndex((s) => typeof s !== 'string');
  const form = segments[i];
  if (i < 1 || typeof form === 'string' || form?.kind !== 'node') return undefined;
  return isSource(printAddress(segments.slice(0, i))) ? printAddress(segments) : undefined;
}

/**
 * The described maps whose ends are not both held yet: pending (EC-169) until the description arrives (`held`: the
 * static model's systems and the overlay's).
 */
export function pendingMaps(
  maps: Iterable<DescribedMap> | undefined,
  held: (system: string) => boolean,
): readonly string[] {
  const out: string[] = [];
  for (const d of maps ?? []) if (!held(d.from) || !held(d.to)) out.push(d.address);
  return out;
}
