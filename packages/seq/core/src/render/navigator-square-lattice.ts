/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node, Value } from '@molstar/core/repr/tree';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import {
  blockSpan,
  type ResolvedNavigator,
  type ResolvedSection,
  type ResolvedView,
  resolveViewNode,
} from '@molstar/seqx-core/view';

// The square's lattice (_spec/README.md "Render", Navigators, The square): the assembly's blocks as one section, what
// a bound section shows of it, its washed runs and its names' shared prefix. Pure helpers of `SquareNavigator`.

/**
 * The washed cells (a cell washed when a bound section does not show its x block or its y block) as rectangles of
 * block indices, `[i0, i1, j0, j1]` each: a column run of x blocks not shown (every row), then for each row run of y
 * blocks not shown, a rectangle over each column run of x blocks shown (no cell twice: a wash is translucent).
 */
export function washRuns(shownX: readonly boolean[], shownY: readonly boolean[]): Int32Array {
  const runs = (shown: readonly boolean[], want: boolean): number[] => {
    const out: number[] = [];
    for (let k = 0; k < shown.length; ) {
      if (shown[k] !== want) {
        k++;
        continue;
      }
      const a = k;
      while (k < shown.length && shown[k] === want) k++;
      out.push(a, k);
    }
    return out;
  };
  const n = shownY.length;
  const hiddenX = runs(shownX, false);
  const shownXRuns = runs(shownX, true);
  const hiddenY = runs(shownY, false);
  const out: number[] = [];
  for (let r = 0; r < hiddenX.length; r += 2) out.push(hiddenX[r], hiddenX[r + 1], 0, n);
  for (let s = 0; s < hiddenY.length; s += 2)
    for (let r = 0; r < shownXRuns.length; r += 2)
      out.push(shownXRuns[r], shownXRuns[r + 1], hiddenY[s], hiddenY[s + 1]);
  return Int32Array.from(out);
}

/** The letters every name starts with, never a whole name (`chr` of `chr1`…`chrX`; none of `chr1` alone). */
export function sharedPrefix(names: readonly string[]): string {
  if (names.length < 2) return '';
  let p = names[0];
  for (const n of names) while (p && !n.startsWith(p)) p = p.slice(0, -1);
  while (p && names.some((n) => n.length === p.length)) p = p.slice(0, -1);
  // letters only (never a digit of the name: `chr1` and `chr12` share `chr1`)
  return p.replace(/[0-9]+$/, '');
}

/** A section's units end to end: its last block's end. */
export function sectionTotal(s: ResolvedSection): number {
  let t = 0;
  for (const b of s.blocks ?? []) t = Math.max(t, (b.offset ?? 0) + blockSpan(b));
  return t;
}

/** Whether a bound section shows a system (one of its blocks is on it). */
export function shows(s: ResolvedSection | undefined, system: string): boolean {
  return !!s?.blocks?.some((b) => b.axis.system === system);
}

/** The lattice: a section of the navigator's assembly (its `include` alone), else of the first bound section's `of`. */
export function latticeOf(model: StaticModel, view: ResolvedView, nav: ResolvedNavigator): ResolvedSection | undefined {
  const of = latticeAssembly(model, view, nav);
  if (of === undefined) return undefined;
  const node: Node = {
    kind: 'view',
    id: nav.id,
    params: { layout: 'tracks' },
    children: [{ kind: 'section', id: 'square', params: { of } }],
  };
  return resolveViewNode(model, node, nav.address).sections[0];
}

/** What the lattice shows: the navigator's assembly, else the first bound section's `of`. */
export function latticeAssembly(model: StaticModel, view: ResolvedView, nav: ResolvedNavigator): Value | undefined {
  const first = view.sections.find((s) => s.address === nav.of[0]);
  let of: Value | undefined = nav.assemblies
    ? nav.assemblies.length === 1
      ? nav.assemblies[0]
      : (nav.assemblies as Value)
    : undefined;
  if (of === undefined && first) {
    let n: Node | undefined = model.spec;
    for (const id of first.address.split('/').filter(Boolean)) n = n?.children?.find((c) => c.id === id);
    of = n?.params?.of as Value | undefined;
  }
  return of;
}
