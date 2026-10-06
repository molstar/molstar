/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { BehaviorSubject } from 'rxjs';
import type { PluginStateAnimation } from '../../animation/model.js';
import { PluginAnimationManager } from '../animation.js';

function create() {
  const context = { behaviors: { state: { isAnimating: new BehaviorSubject(false) } } } as unknown as PluginContext;
  return new PluginAnimationManager(context);
}

function animation(name: string, extra?: Partial<PluginStateAnimation>): PluginStateAnimation {
  return {
    name,
    display: { name: `Display ${name}` },
    params: () => ({}),
    initialState: () => ({}),
    apply: async () => ({ kind: 'finished' }),
    ...extra,
  } as PluginStateAnimation;
}

function names(manager: PluginAnimationManager) {
  return manager.animations.map((a) => a.name);
}

describe('PluginAnimationManager registration', () => {
  it('counts the same animation object', () => {
    const manager = create();
    const a = animation('a');
    manager.register(a);
    manager.register(a);
    expect(names(manager)).toEqual(['a']);

    manager.unregister(a);
    expect(names(manager)).toEqual(['a']);
    manager.unregister(a);
    expect(names(manager)).toEqual([]);
    expect(manager.isEmpty).toBe(true);
  });

  it('throws for a different object under an existing name and keeps the original', () => {
    const manager = create();
    const a = animation('a');
    manager.register(a);
    expect(() => manager.register(animation('a'))).toThrow(/PluginAnimationManager.*'a'/);
    expect(manager.animations).toEqual([a]);
  });

  it('findConflict reports without changing anything', () => {
    const manager = create();
    const a = animation('a');
    expect(manager.findConflict(a)).toBeUndefined();
    manager.register(a);
    expect(manager.findConflict(a)).toBeUndefined();
    expect(manager.findConflict(animation('a'))).toMatch(/'a'/);

    // the lookups did not count anything
    manager.unregister(a);
    expect(manager.isEmpty).toBe(true);
  });

  it('ignores unknown removals', () => {
    const manager = create();
    const a = animation('a');
    manager.unregister(a);
    manager.register(a);
    manager.unregister(animation('a'));
    manager.unregister(animation('other'));
    expect(manager.animations).toEqual([a]);
  });

  it('selects the first registered animation as the current one', () => {
    const manager = create();
    manager.register(animation('a'));
    manager.register(animation('b'));
    expect(manager.current.anim.name).toBe('a');
    expect(manager.state.params.current).toBe('a');
  });

  it('removing a non-current animation keeps the current one and resets cached params', () => {
    const manager = create();
    const a = animation('a');
    const b = animation('b');
    manager.register(a);
    manager.register(b);
    expect(manager.getParams().current).toBeDefined();

    manager.unregister(b);
    expect(manager.current.anim).toBe(a);
    expect(JSON.stringify(manager.getParams())).not.toContain('Display b');
  });

  it('removing the current animation selects the first remaining one', () => {
    const manager = create();
    const a = animation('a');
    const b = animation('b');
    const c = animation('c');
    manager.register(a);
    manager.register(b);
    manager.register(c);

    manager.unregister(a);
    expect(manager.current.anim).toBe(b);
    expect(manager.state.params.current).toBe('b');
    expect(JSON.stringify(manager.getParams())).not.toContain('Display a');
  });

  it('removing the last animation leaves no current one', () => {
    const manager = create();
    const a = animation('a');
    manager.register(a);
    manager.unregister(a);
    expect(manager.current).toBeUndefined();
    expect(manager.state.params.current).toBe('');
    expect(manager.isAnimatingStateTransition).toBe(false);

    // registering again selects it
    manager.register(a);
    expect(manager.current.anim).toBe(a);
  });

  it('stops and tears down a playing current animation when it is removed', async () => {
    const manager = create();
    const teardown = jest.fn();
    const a = animation('a', { teardown });
    const b = animation('b');
    manager.register(a);
    manager.register(b);
    await manager.start();
    expect(manager.isAnimating).toBe(true);

    manager.unregister(a);
    expect(manager.isAnimating).toBe(false);
    expect(teardown).toHaveBeenCalledTimes(1);
    expect(manager.current.anim).toBe(b);
  });

  it('adopts an animation added on demand by play without removing it on unregister', async () => {
    const manager = create();
    const a = animation('a');
    await manager.play(a, {});
    expect(manager.animations).toEqual([a]);
    await manager.stop();

    // not counted, so unregister is a no-op
    manager.unregister(a);
    expect(manager.animations).toEqual([a]);

    // a counted registration adopts it
    manager.register(a);
    manager.unregister(a);
    expect(manager.animations).toEqual([]);
  });
});
