/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Table } from '@molstar/core/data/db';
import { BasicSchema, createBasic } from '../schema.js';
import { getModelGroupName, getModelGroupNames } from '../util.js';

describe('model group names', () => {
  it("are the first link's first group's name", () => {
    const data = createBasic({
      ihm_model_group: Table.ofArrays(BasicSchema.ihm_model_group, {
        id: [1, 2, 2, 3],
        name: ['cluster 1', 'cluster 2', 'cluster 2 again', 'cluster 3'],
      }),
      ihm_model_group_link: Table.ofArrays(BasicSchema.ihm_model_group_link, {
        model_id: [10, 11, 12, 11, 13, 14],
        group_id: [1, 2, 2, 1, 4, 3],
      }),
    });
    const expected = ['cluster 1', 'cluster 2', 'cluster 2', '', 'cluster 3', ''];
    expect([10, 11, 12, 13, 14, 15].map(getModelGroupNames(data))).toEqual(expected);
    expect([10, 11, 12, 13, 14, 15].map((id) => getModelGroupName(id, data))).toEqual(expected);
  });

  it('are empty without groups', () => {
    const data = createBasic({});
    expect(getModelGroupName(1, data)).toBe('');
    expect(getModelGroupNames(data)(1)).toBe('');
  });
});
