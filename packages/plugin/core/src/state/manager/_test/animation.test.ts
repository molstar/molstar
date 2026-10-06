/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { BehaviorSubject } from 'rxjs';
import type { PluginStateAnimation } from '../../animation/model.js';
import { PluginAnimationManager } from '../animation.js';

const warnings: string[] = [];

function create() {
  const context = {
    behaviors: { state: { isAnimating: new BehaviorSubject(false) } },
    log: { warn: (message: string) => warnings.push(message) },
  } as unknown as PluginContext;
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
    expect(manager.current?.anim).toBe(a);
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
    expect(manager.current?.anim).toBe(b);
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
    expect(manager.current?.anim).toBe(a);
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
    expect(manager.current?.anim).toBe(b);
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

describe('PluginAnimationManager current and snapshots', () => {
  beforeEach(() => {
    warnings.length = 0;
  });

  it('has no current animation until one is registered', () => {
    const manager = create();
    expect(manager.current).toBeUndefined();
    expect(manager.isAnimatingStateTransition).toBe(false);
    expect(manager.getSnapshot().current).toBeUndefined();
    manager.updateCurrentParams({ x: 1 });
    expect(manager.current).toBeUndefined();

    const a = animation('a');
    manager.register(a);
    expect(manager.current?.anim).toBe(a);
    manager.unregister(a);
    expect(manager.current).toBeUndefined();
  });

  it('restores the state of a registered animation', () => {
    const manager = create();
    manager.register(animation('a'));
    manager.register(animation('b'));
    manager.setSnapshot({
      state: { params: { current: 'b' }, animationState: 'stopped' },
      current: { paramValues: { x: 1 }, state: { y: 2 } },
    });
    expect(manager.current?.anim.name).toBe('b');
    expect(manager.current?.paramValues).toEqual({ x: 1 });
    expect(manager.current?.state).toEqual({ y: 2 });
    expect(warnings).toEqual([]);
  });

  it('skips the state of an unregistered animation and warns', () => {
    const manager = create();
    manager.register(animation('a'));
    const before = manager.current!.paramValues;
    manager.setSnapshot({
      state: { params: { current: 'missing' }, animationState: 'stopped' },
      current: { paramValues: { x: 1 }, state: { y: 2 } },
    });
    expect(manager.current?.anim.name).toBe('a');
    expect(manager.current?.paramValues).toBe(before);
    expect(warnings).toHaveLength(1);
    expect(warnings[0]).toContain("'missing'");
  });
});
