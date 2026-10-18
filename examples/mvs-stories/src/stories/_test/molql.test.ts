/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { buildStory } from '../molql.js';

describe('MolQL story', () => {
  it('builds a valid story from MolScriptBuilder and a PyMOL-transpiled expression', () => {
    const story = buildStory();
    expect(MVSData.validationIssues(story)).toEqual(undefined);
    expect(story.snapshots).toHaveLength(4);
    expect(findNodes(story, 'camera')).toHaveLength(0);

    const componentSelectors = findNodes(story, 'component').map((node) => node.params.selector);
    expect(componentSelectors.some((selector) => selector?.molql)).toEqual(true);

    const colorSelectors = findNodes(story, 'color').map((node) => node.params.selector);
    expect(colorSelectors.filter((selector) => selector?.molql)).toHaveLength(2);

    const primitivePositions = findNodes(story, 'primitive').flatMap((node) => [
      node.params.position,
      node.params.start,
      node.params.end,
    ]);
    expect(primitivePositions.filter((position) => position?.molql)).toHaveLength(2);

    const serialized = MVSData.toMVSJ(story);
    expect(serialized).toContain('"kind":"component"');
    expect(serialized).toContain('"kind":"color"');
    expect(serialized).toContain('"kind":"primitive"');
    expect(serialized).toContain('"selector":{"molql"');
  });
});

function findNodes(tree: any, kind: string): any[] {
  const nodes: any[] = [];
  const visit = (node: any) => {
    if (node.kind === kind) nodes.push(node);
    for (const child of node.children ?? []) visit(child);
  };

  for (const snapshot of tree.snapshots ?? []) visit(snapshot.root);
  return nodes;
}
