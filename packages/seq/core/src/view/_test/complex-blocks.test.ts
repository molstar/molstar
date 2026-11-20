/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section of a complex (_spec/README.md "Resolved view", Blocks): a block a unit, polymers in the complex's order
// then each class of `classes` (water left out by default, EC-101), copies merged into one block `×N` (EC-65, EC-67),
// each block's form (`whole` for a molecule, `residues` for a polymer, `atoms` on the atoms system) and its label (a
// molecule's name, else its one component, else its id).

import { resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { moleculesSpec } from '../../render/_test/molecules-fixtures.js';
import { resolveView } from '../index.js';

const blocksOf = (spec: Parameters<typeof resolveSpec>[0]) => {
  const model = resolveSpec(spec);
  const view = resolveView(model, '/views/v')!;
  return { view, blocks: view.sections[0].blocks!, model };
};

describe('a section of a complex', () => {
  it('a block a unit: polymers in order, then branched, ligand, ion; water left out (EC-101); copies merged (EC-65)', () => {
    const { blocks, view } = blocksOf(moleculesSpec());
    expect(blocks.map((b) => [b.label, b.form, b.members.length])).toEqual([
      ['A ×2', 'residues', 2],
      ['Glycan', 'whole', 1],
      ['ATP', 'whole', 1],
      ['NA', 'whole', 1],
    ]);
    expect(blocks[0].systems).toEqual(['/data/c/A', '/data/c/B']);
    expect(view.issues).toEqual([]);
  });

  it('`classes` orders and filters; water shown when listed; `copies: each` a block a unit', () => {
    const { blocks } = blocksOf(moleculesSpec({ classes: ['water', 'polymer'], copies: 'each' }));
    expect(blocks.map((b) => b.label)).toEqual(['HOH', 'A', 'B']);
  });

  it('200 copies are one block (EC-67)', () => {
    const { blocks } = blocksOf(moleculesSpec({ chains: 200 }));
    expect(blocks[0].label).toBe('A ×200');
    expect(blocks[0].members.length).toBe(200);
  });

  it('`form: atoms` shows the atoms system; explicit blocks win', () => {
    const { blocks } = blocksOf(moleculesSpec({ atpForm: 'atoms' }));
    const l = blocks.find((b) => b.address.endsWith('/l'))!;
    expect([l.form, l.axis.system, l.axis.ranges]).toEqual(['atoms', '/data/c/l/(atoms)', [0, 4]]);
  });
});
