/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Task } from '@molstar/core/task';
import { Representation } from '../representation.js';

describe('Representation.createUpdateQueue', () => {
  it('runs the updates of a representation one after another, also when started together', async () => {
    const serialize = Representation.createUpdateQueue();
    const log: string[] = [];
    const update = (name: string, ms: number) =>
      serialize(
        Task.create(name, async () => {
          log.push(`${name} start`);
          await new Promise((resolve) => setTimeout(resolve, ms));
          log.push(`${name} end`);
        }),
      );
    // e.g. a state update in progress and Canvas3D's update after a context change
    await Promise.all([update('a', 20).run(), update('b', 0).run()]);
    expect(log).toEqual(['a start', 'a end', 'b start', 'b end']);
  });

  it('runs the next update after one failed', async () => {
    const serialize = Representation.createUpdateQueue();
    const failing = serialize(
      Task.create('a', async () => {
        throw new Error('failed');
      }),
    );
    const next = serialize(Task.create('b', async () => 'b'));
    await expect(failing.run()).rejects.toThrow('failed');
    await expect(next.run()).resolves.toBe('b');
  });
});
