/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ColorTheme } from '../color.js';
import { SizeTheme } from '../size.js';
import { ThemeRegistry } from '../theme.js';

type Registry = ColorTheme.Registry;

function provider(name: string, label = name, category = 'cat'): ColorTheme.Provider<{}> {
  return { ...ColorTheme.EmptyProvider, name, label, category };
}

function create(): Registry {
  return new ThemeRegistry(ColorTheme.EmptyProvider);
}

const names = (r: Registry) => r.list.map((e) => e.name);

describe('ThemeRegistry', () => {
  it('starts empty and does not preload built-ins', () => {
    const r = create();
    expect(r.list.length).toBe(0);
    expect(r.default).toBeUndefined();
    expect(r.get('uniform')).toBe(ColorTheme.EmptyProvider);
  });

  it('counts the same object registered twice', () => {
    const r = create();
    const a = provider('a');
    r.add(a);
    r.add(a);
    expect(names(r)).toEqual(['a']);
    r.remove(a);
    expect(r.has(a)).toBe(true);
    r.remove(a);
    expect(r.has(a)).toBe(false);
    expect(names(r)).toEqual([]);
  });

  it('throws for a different object under the same name', () => {
    const r = create();
    const a = provider('a');
    r.add(a);
    expect(() => r.add(provider('a'))).toThrow(/'a'/);
    expect(r.get('a')).toBe(a);
    expect(names(r)).toEqual(['a']);
  });

  it('findConflict reports the add error and changes nothing', () => {
    const r = create();
    const a = provider('a');
    expect(r.findConflict(a)).toBeUndefined();
    expect(r.has(a)).toBe(false);
    r.add(a);
    expect(r.findConflict(a)).toBeUndefined();
    const other = provider('a');
    const message = r.findConflict(other);
    expect(typeof message).toBe('string');
    expect(() => r.add(other)).toThrow(message);
    r.remove(a);
    expect(r.has('a')).toBe(false);
  });

  it('remove of an unknown provider is a no-op and keeps other entries', () => {
    const r = create();
    const a = provider('a');
    const b = provider('b');
    r.add(a);
    r.add(b);
    r.remove(provider('c'));
    r.remove(provider('b')); // same name, not the registered object
    expect(names(r)).toEqual(['a', 'b']);
    expect(r.has(a)).toBe(true);
    expect(r.has(b)).toBe(true);
  });

  it('remove on an empty registry is a no-op', () => {
    const r = create();
    expect(() => r.remove(provider('a'))).not.toThrow();
    expect(names(r)).toEqual([]);
  });

  it('removes only the given provider', () => {
    const r = create();
    const a = provider('a');
    const b = provider('b');
    const c = provider('c');
    r.add(a);
    r.add(b);
    r.add(c);
    r.remove(b);
    expect(names(r)).toEqual(['a', 'c']);
    expect(r.has('b')).toBe(false);
  });

  it('clear drops providers and counts', () => {
    const r = create();
    const a = provider('a');
    r.add(a);
    r.add(a);
    r.add(provider('b'));
    r.clear();
    expect(names(r)).toEqual([]);
    expect(r.default).toBeUndefined();
    expect(r.has(a)).toBe(false);
    r.add(a);
    r.clear();
    r.remove(a);
    r.add(a);
    r.remove(a);
    expect(r.has(a)).toBe(false);
    expect(() => r.getName(a)).toThrow();
  });

  it('has works by name and by provider', () => {
    const r = create();
    const a = provider('a');
    expect(r.has('a')).toBe(false);
    expect(r.has(a)).toBe(false);
    r.add(a);
    expect(r.has('a')).toBe(true);
    expect(r.has(a)).toBe(true);
    expect(r.has('b')).toBe(false);
    expect(r.has(provider('a'))).toBe(false);
  });

  it('keeps entries sorted by category and label', () => {
    const r = create();
    const x = provider('x', 'Zed', 'B');
    r.add(x);
    r.add(provider('y', 'Alpha', 'B'));
    r.add(provider('z', 'Mid', 'A'));
    expect(names(r)).toEqual(['z', 'y', 'x']);
    expect(r.default?.name).toBe('z');
    // counting a provider again does not reorder or duplicate
    r.add(x);
    expect(names(r)).toEqual(['z', 'y', 'x']);
  });

  it('createRegistry preloads the built-in catalogs', () => {
    const color = ColorTheme.createRegistry();
    expect(color.has('uniform')).toBe(true);
    expect(color.list.length).toBeGreaterThan(10);
    const size = SizeTheme.createRegistry();
    expect(size.has('uniform')).toBe(true);
    // adding a preloaded built-in again counts instead of throwing
    const uniform = color.get('uniform');
    expect(() => color.add(uniform)).not.toThrow();
    color.remove(uniform);
    expect(color.has('uniform')).toBe(true);
  });
});
