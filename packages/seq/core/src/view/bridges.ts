/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A view's bridges (_spec/README.md "Resolved view", Bridges): a band of pairs between two neighbouring sections in the
// stack, its sides sections in the tracks layout (blocks in a row), its layers `pairs` (another type drawn as its chain's
// first `pairs`, said), each pair's ends placed at layout by the systems its sides show (no path: a pair is never
// projected). The stack is the view's sections and bridges in order.

import { exactAddress, type Issue, type Node, type Value } from '@molstar/core/repr/tree';
import type { ViewIssues } from './issues.js';
import type { ResolvedLayer, ResolvedSection } from './view.js';

/**
 * A bridge (_spec/README.md "Resolved view", Bridges): its two sides (section addresses), its params with their defaults
 * where the figure has one, its layers, and whether it is drawn (its sides neighbours in the stack and rows; false with
 * an issue). `above` names which of its sides is the section above it.
 */
export interface ResolvedBridge {
  readonly address: string;
  /** Its sides as written: the sections' addresses (`a` the data's `a` end's side when both show a system). */
  readonly a: string;
  readonly b: string;
  /** Its `label`, else "a to b" by the sides' ids. */
  readonly label: string;
  /** Its height in CSS px as written; the render config's `bridges.height` when absent. */
  readonly height?: number;
  readonly lock: 'none' | 'px';
  readonly hidden: boolean;
  /** Its `pairs` layers, resolved on their data's systems (no section: a bridge has no axis). */
  readonly layers: readonly ResolvedLayer[];
  /** Which side is the section above it in the stack (`a` or `b`); `a` while not placed. */
  readonly above: 'a' | 'b';
  /** Whether it is drawn: its sides neighbours in the stack, each in the tracks layout with blocks in a row. */
  readonly placed: boolean;
  /** What resolving it said (in the view's issues too). */
  readonly issues: readonly Issue[];
}

/** A place in the view's stack: a section or a bridge, by its index in `sections` or `bridges`. */
export interface StackEntry {
  readonly kind: 'section' | 'bridge';
  readonly index: number;
}

/** What resolving a bridge reads of the view being built. */
export interface BridgeInput {
  readonly view: Node;
  readonly address: string;
  readonly sections: readonly ResolvedSection[];
  /** Resolves one layer of a bridge (its type chain and data, `bridgeLayer` in view.ts). */
  layer(
    l: Node,
    index: number,
    at: string,
    say: (code: BridgeCode, message: string, at: string) => void,
  ): ResolvedLayer;
}

/** The issue codes a bridge says (ViewIssues). */
export type BridgeCode = keyof typeof ViewIssues.codes;

const param = (n: Node, key: string): Value | undefined => n.params?.[key];

const idOf = (address: string) => address.slice(address.lastIndexOf('/') + 1);

/**
 * The view's stack and bridges (_spec/README.md "Resolved view", Bridges): sections and bridges in the order written,
 * each bridge's sides checked against its neighbours, their layouts, and its layers resolved.
 */
export function resolveBridges(
  input: BridgeInput,
  report: (code: BridgeCode, message: string, at: string) => Issue,
): { bridges: ResolvedBridge[]; stack: StackEntry[] } {
  const { view, address, sections } = input;
  const bridges: ResolvedBridge[] = [];
  const stack: StackEntry[] = [];
  const sectionIndex = new Map(sections.map((s, i) => [s.address, i]));
  const children = (view.children ?? []).filter(
    (c) => (c.kind === 'section' || c.kind === 'bridge') && c.id !== undefined,
  );
  for (let ci = 0; ci < children.length; ci++) {
    const c = children[ci];
    const at = `${address}/${c.id}`;
    if (c.kind === 'section') {
      const i = sectionIndex.get(at);
      if (i !== undefined) stack.push({ kind: 'section', index: i });
      continue;
    }
    const issues: Issue[] = [];
    const say = (code: BridgeCode, message: string, where: string) => issues.push(report(code, message, where));
    const a = exactAddress(param(c, 'a')) ?? '';
    const b = exactAddress(param(c, 'b')) ?? '';
    // its neighbours in the stack: the section just above and the one just below (a bridge between: none)
    const above = ci > 0 && children[ci - 1].kind === 'section' ? `${address}/${children[ci - 1].id}` : undefined;
    const below =
      ci + 1 < children.length && children[ci + 1].kind === 'section' ? `${address}/${children[ci + 1].id}` : undefined;
    let placed = a !== '' && b !== '' && a !== b;
    let top: 'a' | 'b' = 'a';
    if (placed && !((above === a && below === b) || (above === b && below === a))) {
      say(
        'seq.view.bridge.neighbours',
        `The bridge joins ${idOf(a)} and ${idOf(b)}, which are not the sections just above and below it`,
        at,
      );
      placed = false;
    } else if (placed) top = above === a ? 'a' : 'b';
    // each side a row: the tracks layout, its blocks (if any) in a row
    for (const side of placed ? [a, b] : []) {
      const s = sections[sectionIndex.get(side)!];
      const row = s.layout === 'tracks' && !(s.blocks && s.params.arrange === 'stack');
      if (!row) {
        say(
          'seq.view.bridge.layout',
          `The side ${idOf(side)} is not a row of tracks (a ${s.blocks && s.params.arrange === 'stack' ? 'stack of blocks' : s.layout})`,
          at,
        );
        placed = false;
      }
    }
    const lock = param(c, 'lock') === 'px' ? 'px' : 'none';
    if (placed && lock === 'px') {
      const sa = sections[sectionIndex.get(a)!];
      const sb = sections[sectionIndex.get(b)!];
      if (sa.follow?.leader === b || sb.follow?.leader === a)
        say('seq.view.bridge.lock', 'The sides already move together by a follow: the lock moves nothing more', at);
    }
    const layers: ResolvedLayer[] = [];
    (c.children ?? []).forEach((l, index) => {
      if (l.kind === 'layer') layers.push(input.layer(l, index, at, say));
    });
    const label = param(c, 'label');
    const height = param(c, 'height');
    stack.push({ kind: 'bridge', index: bridges.length });
    bridges.push({
      address: at,
      a,
      b,
      label: typeof label === 'string' ? label : `${idOf(a)} to ${idOf(b)}`,
      ...(typeof height === 'number' && Number.isFinite(height) && height >= 0 ? { height } : {}),
      lock,
      hidden: param(c, 'hidden') === true,
      layers,
      above: top,
      placed,
      issues,
    });
  }
  return { bridges, stack };
}
