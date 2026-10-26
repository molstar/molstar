/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { EveryLoci, type Loci } from '@molstar/model/model/loci';
import type { Representation } from '../../repr/representation.js';
import { MarkIndex } from '../mark-index.js';

/** a representation that marks the structure loci of roots with `hashes` (`undefined`: of any root) */
const repr = (hashes?: number[]) =>
  ({ getStructureRootHashes: hashes === undefined ? undefined : () => hashes }) as unknown as Representation.Any;

/** the loci kinds `MarkIndex` reads: the root's `hashCode` */
const elementLoci = (hash: number) =>
  ({ kind: 'element-loci', structure: { root: { hashCode: hash } } }) as unknown as Loci;
const structureLoci = (hash: number) =>
  ({ kind: 'structure-loci', structure: { root: { hashCode: hash } } }) as unknown as Loci;

const visited = (index: MarkIndex, loci: Loci, action = MarkerAction.Highlight) => {
  const out: Representation.Any[] = [];
  index.forEach(loci, action, (r) => out.push(r));
  return out;
};

describe('MarkIndex', () => {
  it('visits the representations of a structure loci root and those that may mark any', () => {
    const index = new MarkIndex();
    const a = repr([1]);
    const b = repr([2]);
    const ab = repr([1, 2]);
    const any = repr();
    const none = repr([]);
    for (const r of [a, b, ab, any, none]) index.set(r);

    expect(new Set(visited(index, elementLoci(1)))).toEqual(new Set([a, ab, any]));
    expect(new Set(visited(index, structureLoci(2)))).toEqual(new Set([b, ab, any]));
    expect(visited(index, elementLoci(3))).toEqual([any]);
    // other loci: all, each once
    expect(new Set(visited(index, EveryLoci))).toEqual(new Set([a, b, ab, any, none]));
    expect(visited(index, EveryLoci).length).toBe(5);
  });

  it('takes a representation again when its root changed, and forgets it', () => {
    const index = new MarkIndex();
    let hashes = [1];
    const r = { getStructureRootHashes: () => hashes } as unknown as Representation.Any;
    index.set(r);
    expect(visited(index, elementLoci(1))).toEqual([r]);

    hashes = [2];
    index.set(r);
    expect(visited(index, elementLoci(1))).toEqual([]);
    expect(visited(index, elementLoci(2))).toEqual([r]);

    index.delete(r);
    expect(visited(index, elementLoci(2))).toEqual([]);
    expect(visited(index, EveryLoci)).toEqual([]);
  });

  it('deselects all on the representations markers were added to, also after a remove and add', () => {
    const index = new MarkIndex();
    const [a, b, c] = [repr([1]), repr([2]), repr()];
    for (const r of [a, b, c]) index.set(r);
    // nothing marked: nothing to remove
    expect(visited(index, EveryLoci, MarkerAction.Deselect)).toEqual([]);

    index.note(a, MarkerAction.Select);
    index.note(b, MarkerAction.RemoveHighlight); // removing: not noted
    index.note(c, MarkerAction.Highlight);
    expect(new Set(visited(index, EveryLoci, MarkerAction.Deselect))).toEqual(new Set([a, c]));
    // kept: a deselect leaves highlights, a remove of highlights selections
    expect(new Set(visited(index, EveryLoci, MarkerAction.RemoveHighlight))).toEqual(new Set([a, c]));
    // adding actions on `EveryLoci`: all
    expect(new Set(visited(index, EveryLoci, MarkerAction.Select))).toEqual(new Set([a, b, c]));

    // removed and added again: still holds its markers
    index.delete(a);
    expect(visited(index, EveryLoci, MarkerAction.Deselect)).toEqual([c]);
    index.set(a);
    expect(new Set(visited(index, EveryLoci, MarkerAction.Deselect))).toEqual(new Set([a, c]));

    // noted before it was added
    const d = repr([3]);
    index.note(d, MarkerAction.Toggle);
    index.set(d);
    expect(visited(index, EveryLoci, MarkerAction.Deselect)).toContain(d);

    // a clear removes all markers: none left to remove
    expect(new Set(visited(index, EveryLoci, MarkerAction.Clear))).toEqual(new Set([a, c, d]));
    expect(visited(index, EveryLoci, MarkerAction.Deselect)).toEqual([]);
  });
});
