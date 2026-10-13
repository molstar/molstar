/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { DefaultDragAndDrop } from '../../../default-registry.js';
import { DragAndDropManager, type PluginDragAndDropHandler } from '../drag-and-drop.js';

function create() {
  const runTask = jest.fn();
  const applyAction = jest.fn(() => 'apply-action');
  const dispatch = jest.fn(() => Promise.resolve());
  const plugin = { runTask, commands: { dispatch }, state: { data: { applyAction } } } as unknown as PluginContext;
  return { manager: new DragAndDropManager(plugin), runTask, applyAction, dispatch };
}

function file(name: string) {
  return { name } as File;
}

function handler(log: string[], name: string, result = false): PluginDragAndDropHandler {
  return jest.fn((_files: File[]) => {
    log.push(name);
    return result;
  });
}

describe('DragAndDropManager registration', () => {
  it('counts the same name and handle', async () => {
    const { manager } = create();
    const log: string[] = [];
    const h = handler(log, 'a', true);
    manager.addHandler('a', h);
    manager.addHandler('a', h);
    manager.addEntry({ name: 'a', handle: h });

    manager.removeHandler('a');
    manager.removeHandler('a');
    await manager.handle([file('x.txt')]);
    expect(log).toEqual(['a']);

    manager.removeHandler('a');
    log.length = 0;
    await manager.handle([file('x.txt')]);
    expect(log).toEqual([]);
  });

  it('throws for a different function under an existing name and keeps the original', async () => {
    const { manager } = create();
    const log: string[] = [];
    manager.addHandler('a', handler(log, 'first', true));
    expect(() => manager.addHandler('a', handler(log, 'second', true))).toThrow(/DragAndDropManager.*'a'/);

    await manager.handle([file('x.txt')]);
    expect(log).toEqual(['first']);
  });

  it('lists entries in registration order', () => {
    const { manager } = create();
    const h = jest.fn(() => true);
    manager.addHandler('a', h);
    manager.addHandler('b', h, { fallback: true });
    manager.addHandler('a', h);
    expect(manager.list()).toEqual([
      { name: 'a', handle: h },
      { name: 'b', handle: h, fallback: true },
    ]);
  });

  it('findConflict reports without changing anything', async () => {
    const { manager } = create();
    const log: string[] = [];
    const h = handler(log, 'a', true);
    const other = handler(log, 'b', true);
    expect(manager.findConflict({ name: 'a', handle: h })).toBeUndefined();
    manager.addEntry({ name: 'a', handle: h });
    expect(manager.findConflict({ name: 'a', handle: h })).toBeUndefined();
    expect(manager.findConflict({ name: 'a', handle: other })).toMatch(/'a'/);
    expect(manager.findConflict({ name: 'b', handle: other })).toBeUndefined();

    // the lookup did not count anything
    manager.removeHandler('a');
    await manager.handle([file('x.txt')]);
    expect(log).toEqual([]);
  });

  it('treats entries with the same name and handle as the same provider', () => {
    const { manager } = create();
    const h = jest.fn(() => true);
    manager.addEntry({ name: 'a', handle: h });
    expect(() => manager.addEntry({ name: 'a', handle: h, fallback: true })).not.toThrow();
  });

  it('ignores unknown removals', async () => {
    const { manager } = create();
    const log: string[] = [];
    const h = handler(log, 'a', true);
    manager.removeHandler('missing');
    manager.removeEntry({ name: 'missing', handle: h });
    manager.addHandler('a', h);
    // same name, different function
    manager.removeEntry({ name: 'a', handle: handler(log, 'b') });
    await manager.handle([file('x.txt')]);
    expect(log).toEqual(['a']);

    manager.removeEntry({ name: 'a', handle: h });
    log.length = 0;
    await manager.handle([file('x.txt')]);
    expect(log).toEqual([]);
  });
});

describe('DragAndDropManager.handle', () => {
  it('tries non-fallback handlers in reverse registration order and stops at the first that handles', async () => {
    const { manager } = create();
    const log: string[] = [];
    manager.addHandler('a', handler(log, 'a', true));
    manager.addHandler('b', handler(log, 'b', true));
    manager.addHandler('c', handler(log, 'c', false));
    await manager.handle([file('x.txt')]);
    expect(log).toEqual(['c', 'b']);
  });

  it('runs session handling after non-fallback and before fallback handlers', async () => {
    const { manager, runTask, dispatch } = create();
    const log: string[] = [];
    manager.addHandler('fallback', handler(log, 'fallback', true), { fallback: true });
    manager.addHandler('normal', handler(log, 'normal', false));

    await manager.handle([file('x.molx')]);
    expect(log).toEqual(['normal']);
    expect(dispatch).toHaveBeenCalledTimes(1);
    expect(runTask).not.toHaveBeenCalled();

    log.length = 0;
    await manager.handle([file('x.MOLJ')]);
    expect(dispatch).toHaveBeenCalledTimes(2);
  });

  it('a non-fallback handler takes precedence over session handling', async () => {
    const { manager, dispatch } = create();
    const log: string[] = [];
    manager.addHandler('normal', handler(log, 'normal', true));
    await manager.handle([file('x.molx')]);
    expect(log).toEqual(['normal']);
    expect(dispatch).not.toHaveBeenCalled();
  });

  it('runs fallback handlers last, in reverse registration order, whatever the registration order', async () => {
    const { manager, runTask } = create();
    const log: string[] = [];
    manager.addHandler('f1', handler(log, 'f1', false), { fallback: true });
    manager.addHandler('n1', handler(log, 'n1', false));
    manager.addHandler('f2', handler(log, 'f2', false), { fallback: true });
    manager.addHandler('n2', handler(log, 'n2', false));
    await manager.handle([file('x.txt')]);
    expect(log).toEqual(['n2', 'n1', 'f2', 'f1']);
    // nothing handled and no default entry: the files are not opened
    expect(runTask).not.toHaveBeenCalled();
  });

  it('stops at a fallback handler that handles the files', async () => {
    const { manager, runTask } = create();
    const log: string[] = [];
    manager.addHandler('f1', handler(log, 'f1', true), { fallback: true });
    manager.addHandler('f2', handler(log, 'f2', true), { fallback: true });
    await manager.handle([file('x.txt')]);
    expect(log).toEqual(['f2']);
    expect(runTask).not.toHaveBeenCalled();
  });

  it('does not open unrecognized drops without the default entry', async () => {
    const { manager, runTask, applyAction } = create();
    await manager.handle([file('x.txt')]);
    expect(applyAction).not.toHaveBeenCalled();
    expect(runTask).not.toHaveBeenCalled();
  });

  describe('with DefaultDragAndDrop', () => {
    function createDefault() {
      const ctx = create();
      for (const entry of DefaultDragAndDrop.dragAndDrop!) ctx.manager.addEntry(entry);
      return ctx;
    }

    it('opens unrecognized drops with the open-files fallback', async () => {
      const { manager, runTask, applyAction } = createDefault();
      await manager.handle([file('x.txt')]);
      expect(applyAction).toHaveBeenCalledTimes(1);
      expect(runTask).toHaveBeenCalledWith('apply-action');
    });

    it('keeps session handling ahead of the open-files fallback', async () => {
      const { manager, runTask, dispatch } = createDefault();
      await manager.handle([file('x.molx')]);
      expect(dispatch).toHaveBeenCalledTimes(1);
      expect(runTask).not.toHaveBeenCalled();
    });

    it('runs the open-files fallback after other handlers that do not take the files', async () => {
      const { manager, runTask } = createDefault();
      const log: string[] = [];
      manager.addHandler('other', handler(log, 'other', false));
      await manager.handle([file('x.txt')]);
      expect(log).toEqual(['other']);
      expect(runTask).toHaveBeenCalledTimes(1);
    });

    it('does not open files another handler took', async () => {
      const { manager, runTask } = createDefault();
      const log: string[] = [];
      manager.addHandler('other', handler(log, 'other', true));
      await manager.handle([file('x.txt')]);
      expect(runTask).not.toHaveBeenCalled();
    });
  });

  it('dispose drops all handlers', async () => {
    const { manager } = create();
    const log: string[] = [];
    manager.addHandler('a', handler(log, 'a', true));
    manager.dispose();
    await manager.handle([file('x.txt')]);
    expect(log).toEqual([]);
  });
});
