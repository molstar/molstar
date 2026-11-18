/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// User-defined refs (_spec/README.md, Language, Refs): a `ref` under the `refs` container names what an address, or a
// set of addresses, names; `/refs/<id>` stands wherever an address does, a prefix of one too (`/refs/p53/(range :from
// 1 :to 9)` when the ref names `/data/p53`). Following a ref is a rewrite of the address's leading ids, bounded by a
// chain depth; a cycle, a chain past the bound and a set followed by more segments are refusals, never a guess.

import { type Node, parseAddress, printAddress, resolveAddress, type Segment } from '@molstar/core/repr/tree';

/** The reserved id of the refs container: a ref's address is `/refs/<id>`. */
export const RefsId = 'refs';

/** How many refs a chain may pass through before it is refused (`/refs/a` naming `/refs/b` is two). */
export const DefaultRefDepth = 8;

/** The addresses a ref node names, as written: one, or each of a set. */
export function refAddresses(node: Node): readonly string[] {
  const to = node.params?.to;
  if (typeof to === 'string') return [to];
  if (Array.isArray(to)) return to.filter((v): v is string => typeof v === 'string');
  return [];
}

/** Why a ref could not be followed. */
export type RefRefusal =
  /** The chain comes back to a ref it passed. */
  | 'cycle'
  /** The chain is longer than the bound. */
  | 'depth'
  /** A set (or a subtree) is followed by more segments: they would apply to each, which an address cannot say. */
  | 'set'
  /** The ref names no address (its `to` is not one). */
  | 'empty';

export type FollowedRef =
  | {
      /** The addresses the chain ends at: the address itself when it names no ref; one, or each of a set. */
      readonly addresses: readonly string[];
      /** The refs passed, in order (`/refs/dbd`), none when the address names no ref. */
      readonly via: readonly string[];
    }
  | { readonly refused: RefRefusal; readonly via: readonly string[] };

/** The leading ids of a pattern's segments. */
function idRun(segments: readonly Segment[]): string[] {
  const ids: string[] = [];
  for (const s of segments) {
    if (typeof s !== 'string') break;
    ids.push(s);
  }
  return ids;
}

/**
 * The ref an address starts with: `/refs/<id>`, a node of kind `ref` in the spec, and the segments after it; undefined
 * when the address does not start with one (or does not parse).
 */
export function refPrefix(
  spec: Node,
  address: string,
):
  | { readonly ref: Node; readonly address: string; readonly rest: readonly Segment[]; readonly subtree: boolean }
  | undefined {
  const pattern = parseAddress(address).pattern;
  if (!pattern || pattern.segments[0] !== RefsId || typeof pattern.segments[1] !== 'string') return undefined;
  const ids = idRun(pattern.segments).slice(0, 2);
  const ref = resolveAddress(spec, ids);
  if (ref?.kind !== 'ref') return undefined;
  return { ref, address: printAddress(ids), rest: pattern.segments.slice(2), subtree: pattern.subtree };
}

/**
 * What an address names once its refs are followed: an address that starts with a ref has that ref's target in place
 * of `/refs/<id>`, again until no ref leads it (at most `maxDepth` refs). A set of addresses is followed member by
 * member; one followed by more segments is refused (`set`).
 */
export function followRef(spec: Node, address: string, maxDepth = DefaultRefDepth): FollowedRef {
  /** The chain being followed (a cycle is judged along it alone), and every ref passed, for the caller. */
  const chain: string[] = [];
  const passed: string[] = [];
  const out: string[] = [];
  const go = (a: string): RefRefusal | undefined => {
    const r = refPrefix(spec, a);
    if (!r) {
      out.push(a);
      return undefined;
    }
    if (!passed.includes(r.address)) passed.push(r.address);
    if (chain.includes(r.address)) return 'cycle';
    if (chain.length >= maxDepth) return 'depth';
    const targets = refAddresses(r.ref);
    if (!targets.length) return 'empty';
    const more = r.rest.length > 0 || r.subtree;
    if (targets.length > 1 && more) return 'set';
    chain.push(r.address);
    for (const t of targets) {
      const p = parseAddress(t).pattern;
      if (!p) return 'empty';
      // a subtree names many nodes: a segment after it would apply to each, as after a set
      if (p.subtree && r.rest.length) return 'set';
      const why = go(more ? printAddress([...p.segments, ...r.rest], r.subtree || p.subtree) : t);
      if (why) return why;
    }
    chain.pop();
    return undefined;
  };
  const why = go(address);
  return why ? { refused: why, via: passed } : { addresses: out, via: passed };
}
