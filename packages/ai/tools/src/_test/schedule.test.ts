/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { toolSchedule } from '../index.js';

const call = (callId: string, tool: string, order: number) => ({ callId, tool: { id: tool }, order });
const reads = (id: string) => id.startsWith('get_');

describe('tool schedule', () => {
  it("edits one at a time in the model's order; reads between them together", () => {
    const calls = [
      call('a', 'get_scene', 0),
      call('b', 'get_elements', 1),
      call('c', 'add_scene', 2),
      call('d', 'set_look', 3),
      call('e', 'get_scene', 4),
      call('f', 'set_camera', 5),
    ];
    expect(toolSchedule(calls, reads)).toEqual([['a', 'b'], ['c'], ['d'], ['e'], ['f']]);
  });

  it("follows the model's order, not the list's", () => {
    expect(toolSchedule([call('x', 'set_look', 1), call('y', 'add_scene', 0)], reads)).toEqual([['y'], ['x']]);
  });

  it('calls without an order keep the list order', () => {
    expect(
      toolSchedule(
        [
          { callId: 'a', tool: { id: 'get_a' } },
          { callId: 'b', tool: { id: 'get_b' } },
          { callId: 'c', tool: { id: 'edit' } },
        ],
        reads,
      ),
    ).toEqual([['a', 'b'], ['c']]);
    expect(toolSchedule([], reads)).toEqual([]);
  });
});
