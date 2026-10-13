/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RepresentationRegistry, type RepresentationProvider } from '../representation.js';
import type { Representation } from '../representation.js';

type Registry = RepresentationRegistry<any, Representation.State>;

function provider(name: string): RepresentationProvider<any, {}, Representation.State> {
  return {
    name,
    label: name,
    description: name,
    factory: () => undefined as any,
    getParams: () => ({}),
    defaultValues: {},
    defaultColorTheme: { name: '' },
    defaultSizeTheme: { name: '' },
    isApplicable: () => true,
  };
}

const names = (r: Registry) => r.list.map((e) => e.name);

describe('RepresentationRegistry', () => {
  it('default is undefined when empty', () => {
    const r: Registry = new RepresentationRegistry();
    expect(r.default).toBeUndefined();
    r.add(provider('a'));
    expect(r.default?.name).toBe('a');
  });

  it('default is the first registered entry', () => {
    const r: Registry = new RepresentationRegistry();
    r.add(provider('b'));
    r.add(provider('a'));
    expect(r.default?.name).toBe('b');
  });

  it('counts the same object registered twice', () => {
    const r: Registry = new RepresentationRegistry();
    const a = provider('a');
    r.add(a);
    r.add(a);
    expect(names(r)).toEqual(['a']);
    r.remove(a);
    expect(r.has(a)).toBe(true);
    expect(names(r)).toEqual(['a']);
    r.remove(a);
    expect(r.has(a)).toBe(false);
    expect(names(r)).toEqual([]);
  });

  it('throws for a different object under the same name', () => {
    const r: Registry = new RepresentationRegistry();
    const a = provider('a');
    r.add(a);
    expect(() => r.add(provider('a'))).toThrow(/'a'/);
    expect(r.get('a')).toBe(a);
    expect(names(r)).toEqual(['a']);
  });

  it('findConflict reports the add error and changes nothing', () => {
    const r: Registry = new RepresentationRegistry();
    const a = provider('a');
    expect(r.findConflict(a)).toBeUndefined();
    expect(r.has(a)).toBe(false);
    r.add(a);
    expect(r.findConflict(a)).toBeUndefined();
    const other = provider('a');
    const message = r.findConflict(other);
    expect(typeof message).toBe('string');
    expect(() => r.add(other)).toThrow(message);
    // the count of `a` is unchanged: a single remove drops it
    r.remove(a);
    expect(r.has('a')).toBe(false);
  });

  it('remove of an unknown provider is a no-op and keeps other entries', () => {
    const r: Registry = new RepresentationRegistry();
    const a = provider('a');
    const b = provider('b');
    r.add(a);
    r.add(b);
    r.remove(provider('c'));
    r.remove(provider('a')); // same name, not the registered object
    expect(names(r)).toEqual(['a', 'b']);
    expect(r.has(a)).toBe(true);
    expect(r.has(b)).toBe(true);
  });

  it('remove on an empty registry is a no-op', () => {
    const r: Registry = new RepresentationRegistry();
    expect(() => r.remove(provider('a'))).not.toThrow();
    expect(names(r)).toEqual([]);
  });

  it('removes only the given provider', () => {
    const r: Registry = new RepresentationRegistry();
    const a = provider('a');
    const b = provider('b');
    const c = provider('c');
    r.add(a);
    r.add(b);
    r.add(c);
    r.remove(b);
    expect(names(r)).toEqual(['a', 'c']);
    expect(r.has('b')).toBe(false);
    expect(r.get('b').name).toBe('');
  });

  it('can register another object under a name after it is removed', () => {
    const r: Registry = new RepresentationRegistry();
    r.add(provider('a'));
    r.remove(r.get('a'));
    const a2 = provider('a');
    r.add(a2);
    expect(r.get('a')).toBe(a2);
    expect(r.getName(a2)).toBe('a');
  });

  it('clear drops providers and counts', () => {
    const r: Registry = new RepresentationRegistry();
    const a = provider('a');
    r.add(a);
    r.add(a);
    r.add(provider('b'));
    r.clear();
    expect(names(r)).toEqual([]);
    expect(r.default).toBeUndefined();
    expect(r.has(a)).toBe(false);
    // a later remove of a cleared provider is a no-op
    r.add(a);
    r.clear();
    r.remove(a);
    r.add(a);
    r.remove(a);
    expect(r.has(a)).toBe(false);
    expect(() => r.getName(a)).toThrow();
  });

  it('has works by name and by provider', () => {
    const r: Registry = new RepresentationRegistry();
    const a = provider('a');
    expect(r.has('a')).toBe(false);
    expect(r.has(a)).toBe(false);
    r.add(a);
    expect(r.has('a')).toBe(true);
    expect(r.has(a)).toBe(true);
    expect(r.has('b')).toBe(false);
    expect(r.has(provider('a'))).toBe(false);
  });
});
