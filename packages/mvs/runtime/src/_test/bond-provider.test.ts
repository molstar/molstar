/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { structureProps } from '../load-helpers.js';

// A Mol*-specific bond provider on a structure node, through `custom.molstar_bond_provider` (bond-providers plan §4.2).

const node = (params: object, custom?: object) => ({ kind: 'structure', params, custom }) as any;

describe('custom.molstar_bond_provider', () => {
  it('without it, the structure props are as before: no bonds param', () => {
    expect(structureProps(node({ type: 'model' }))).toEqual({ type: { name: 'model', params: {} } });
    expect(structureProps(node({ type: 'assembly', assembly_id: '1' }))).toEqual({
      type: { name: 'assembly', params: { id: '1' } },
    });
  });

  it('becomes the bonds param of every structure type, with params completed', () => {
    const none = { molstar_bond_provider: { name: 'none' } };
    expect(structureProps(node({ type: 'model' }, none)).type.params).toEqual({ bonds: { name: 'none', params: {} } });
    expect(structureProps(node({ type: 'assembly', assembly_id: '1' }, none)).type.params).toEqual({
      id: '1',
      bonds: { name: 'none', params: {} },
    });
    const withParams = { molstar_bond_provider: { name: 'my-bonds', params: { cutoff: 2 } } };
    expect(structureProps(node({ type: 'symmetry_mates', radius: 5 }, withParams)).type.params).toEqual({
      radius: 5,
      bonds: { name: 'my-bonds', params: { cutoff: 2 } },
    });
  });

  it('ignores a value without a name', () => {
    expect(structureProps(node({ type: 'model' }, { molstar_bond_provider: { params: {} } }))).toEqual({
      type: { name: 'model', params: {} },
    });
  });
});
