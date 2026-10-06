/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { MarkdownExtensionManager, type MarkdownExtension } from '../markdown-extensions.js';
import { BuiltInMarkdownExtension } from '../../markdown/catalog.js';

function create() {
  return new MarkdownExtensionManager({} as unknown as PluginContext);
}

function ext(name: string, result?: string): MarkdownExtension & { execute: jest.Mock } {
  return { name, execute: jest.fn(), reactRenderFn: () => result };
}

describe('MarkdownExtensionManager extensions', () => {
  it('counts the same extension object', () => {
    const manager = create();
    const e = ext('test-ext', 'rendered');
    manager.registerExtension(e);
    manager.registerExtension(e);

    manager.removeExtension(e);
    expect(manager.tryRender({}, [])).toBe('rendered');
    manager.removeExtension(e);
    expect(manager.tryRender({}, [])).toBe(null);
  });

  it('throws for a different object under an existing name and keeps the original', () => {
    const manager = create();
    manager.registerExtension(ext('test-ext', 'first'));
    expect(() => manager.registerExtension(ext('test-ext', 'second'))).toThrow(/MarkdownExtensionManager.*'test-ext'/);
    expect(manager.tryRender({}, [])).toBe('first');
  });

  it('findConflict reports without changing anything', () => {
    const manager = create();
    const e = ext('test-ext', 'first');
    expect(manager.findConflict(e)).toBeUndefined();
    manager.registerExtension(e);
    expect(manager.findConflict(e)).toBeUndefined();
    expect(manager.findConflict(ext('test-ext'))).toMatch(/'test-ext'/);

    // the lookups did not count anything
    manager.removeExtension(e);
    expect(manager.tryRender({}, [])).toBe(null);
  });

  it('removes by name, decrementing, and ignores unknown removals', () => {
    const manager = create();
    const e = ext('test-ext', 'rendered');
    manager.removeExtension('missing');
    manager.removeExtension(e);

    manager.registerExtension(e);
    manager.registerExtension(e);
    // a different object with the same name is not the registered one
    manager.removeExtension(ext('test-ext'));
    manager.removeExtension('test-ext');
    expect(manager.tryRender({}, [])).toBe('rendered');
    manager.removeExtension('test-ext');
    expect(manager.tryRender({}, [])).toBe(null);
  });

  it('executes registered extensions', () => {
    const manager = create();
    const e = ext('test-ext');
    manager.registerExtension(e);
    manager.tryExecute('click', { a: 'b' });
    expect(e.execute).toHaveBeenCalledTimes(1);
    manager.removeExtension(e);
    manager.tryExecute('click', { a: 'b' });
    expect(e.execute).toHaveBeenCalledTimes(1);
  });

  it('starts empty', () => {
    const manager = create();
    expect((manager as any).extension as MarkdownExtension[]).toEqual([]);
    expect(manager.tryRender({}, [])).toBe(null);
  });

  it('registers the same built-in extensions again by counting', () => {
    const manager = create();
    const builtIn = (manager as any).extension as MarkdownExtension[];
    for (const e of BuiltInMarkdownExtension) manager.registerExtension(e);
    expect(builtIn.length).toBeGreaterThan(0);
    const before = [...builtIn];
    for (const e of before) manager.registerExtension(e);
    expect([...builtIn]).toEqual(before);
    for (const e of before) manager.removeExtension(e);
    expect([...builtIn]).toEqual(before);
    for (const e of before) manager.removeExtension(e);
    expect(builtIn.length).toBe(0);
  });
});
