/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MVSData } from '@molstar/mvs-builder';
import type { ComponentExpressionT, HexColorT } from '@molstar/mvs-builder/tree/mvs/param-types';

// The fake structure viewer's state is MolViewSpec: a structure of an entry with one representation and the `color`
// nodes written over it, validated by `@molstar/mvs-builder` before the fake takes it. The host builds it (the colors
// of a track), the fake reads its color nodes back.

/** What the representation is drawn with before a `color` node says otherwise. */
export const BaseColor: HexColorT = '#c4c8cc';

/** A `color` node: the cells a selector names (`all` when absent) take the color. */
export interface ColorNode {
  readonly selector?: readonly ComponentExpressionT[];
  readonly color: HexColorT;
}

/** The MolViewSpec state of an entry colored by `colors`, in order (a later node over an earlier). */
export function structureState(entry: string, colors: readonly ColorNode[], title?: string): MVSData {
  const b = MVSData.createBuilder();
  const rep = b
    .download({ url: `https://files.rcsb.org/download/${entry}.cif` })
    .parse({ format: 'mmcif' })
    .modelStructure()
    .component({ selector: 'all' })
    .representation({ type: 'cartoon' });
  rep.color({ color: BaseColor });
  for (const c of colors) if (c.selector?.length) rep.color({ selector: [...c.selector], color: c.color });
  return b.getState(title ? { title } : {});
}

/** The issues of a state, or undefined when it is MolViewSpec as the builder's validation reads it. */
export function stateIssues(mvs: MVSData): string[] | undefined {
  return MVSData.validationIssues(mvs);
}

/** A state's `color` nodes in order. */
export function colorNodesOf(mvs: MVSData): ColorNode[] {
  const out: ColorNode[] = [];
  const walk = (n: { kind: string; params?: unknown; children?: unknown[] }) => {
    if (n.kind === 'color') {
      const p = n.params as { selector?: unknown; color: HexColorT };
      const sel = p.selector;
      out.push({
        color: p.color,
        ...(sel === undefined || sel === 'all'
          ? {}
          : { selector: (Array.isArray(sel) ? sel : [sel]) as ComponentExpressionT[] }),
      });
    }
    for (const c of n.children ?? []) walk(c as never);
  };
  walk(mvs.kind === 'multiple' ? mvs.snapshots[0].root : mvs.root);
  return out;
}
