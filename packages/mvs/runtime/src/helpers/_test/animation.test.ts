/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { generateStateTransition } from '../animation.js';

function snapshot(transitions: object[]) {
  return {
    root: {
      kind: 'root',
      children: [
        {
          kind: 'primitives',
          children: [
            {
              kind: 'primitive',
              ref: 'lines',
              params: { kind: 'lines', vertices: [0, 0, 0, 1, 1, 1], indices: [0, 1], group_colors: { 1: '#000000' } },
            },
          ],
        },
      ],
    },
    animation: {
      kind: 'animation',
      params: { frame_time_ms: 100 },
      children: transitions.map((params) => ({ kind: 'interpolate', params })),
    },
    metadata: {},
  } as any;
}

const groupColors = (tree: any) => tree.children[0].children[0].params.group_colors;

describe('generateStateTransition, colors', () => {
  it('leaves a color without a start value as it is (a key the color map does not have, no start)', async () => {
    const result = await generateStateTransition(
      SyncRuntimeContext,
      snapshot([
        { kind: 'color', target_ref: 'lines', property: ['group_colors', 0], start_ms: 0, duration_ms: 200 },
        {
          kind: 'color',
          target_ref: 'lines',
          property: ['group_colors', 1],
          start_ms: 0,
          duration_ms: 200,
          end: '#ffffff',
        },
      ]),
      0,
      1,
    );
    const frames = result!.frames;
    expect(frames.length).toBe(3);
    expect(groupColors(frames[0][0])).toEqual({ 1: '#000000' });
    expect(groupColors(frames[2][0])).toEqual({ 1: '#ffffff' });
    for (const [tree] of frames) expect('0' in groupColors(tree)).toBe(false);
  });

  it('does not change a color without a start value when an end is given', async () => {
    const result = await generateStateTransition(
      SyncRuntimeContext,
      snapshot([
        {
          kind: 'color',
          target_ref: 'lines',
          property: ['group_colors', 0],
          start_ms: 0,
          duration_ms: 200,
          end: '#ff0000',
        },
      ]),
      0,
      1,
    );
    expect(result!.frames.length).toBe(1);
    expect(groupColors(result!.frames[0][0])).toEqual({ 1: '#000000' });
  });
});
