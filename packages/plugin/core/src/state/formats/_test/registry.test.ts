/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DataFormatProvider } from '../provider.js';
import { DataFormatRegistry } from '../registry.js';

function format(name: string, extra: Partial<DataFormatProvider> = {}) {
  return DataFormatProvider({
    name,
    label: name.toUpperCase(),
    description: name,
    parse: async () => undefined,
    ...extra,
  });
}

function empty() {
  return new DataFormatRegistry();
}

const names = (registry: DataFormatRegistry) => registry.list.map((e) => e.name);

describe('DataFormatProvider.withName', () => {
  it('returns the provider when the name matches', () => {
    const a = format('a');
    expect(DataFormatProvider.withName(a, 'a')).toBe(a);
  });

  it('returns a memoized copy otherwise', () => {
    const a = format('a', { stringExtensions: ['abc'] });
    const copy = DataFormatProvider.withName(a, 'b');
    expect(copy).not.toBe(a);
    expect(copy.name).toBe('b');
    expect(copy.stringExtensions).toEqual(['abc']);
    expect(a.name).toBe('a');
    expect(DataFormatProvider.withName(a, 'b')).toBe(copy);
    expect(DataFormatProvider.withName(a, 'c')).not.toBe(copy);
  });

  it('names an unnamed provider', () => {
    const { name: _, ...unnamed } = format('a');
    const named = DataFormatProvider.withName(unnamed, 'x');
    expect(named.name).toBe('x');
    expect(DataFormatProvider.withName(unnamed, 'x')).toBe(named);
  });
});

describe('DataFormatRegistry', () => {
  it('starts empty', () => {
    const registry = new DataFormatRegistry();
    expect(registry.list).toHaveLength(0);
    expect(registry.has('mmcif')).toBe(false);
    expect(registry.extensions.size).toBe(0);
    expect(() => registry.get('mmcif')).toThrow(/not registered/);
  });

  it('clear drops everything', () => {
    const registry = new DataFormatRegistry();
    registry.add(format('a', { stringExtensions: ['abc'] }));
    expect(registry.extensions.size).toBe(1);
    registry.clear();
    expect(registry.list).toHaveLength(0);
    expect(registry.has('a')).toBe(false);
    expect(registry.extensions.size).toBe(0);
  });

  it('counts repeated adds of the same provider', () => {
    const registry = empty();
    const a = format('a');
    registry.add(a);
    registry.add(a);
    expect(names(registry)).toEqual(['a']);
    registry.remove(a);
    expect(registry.has('a')).toBe(true);
    registry.remove(a);
    expect(registry.has('a')).toBe(false);
    expect(registry.list).toHaveLength(0);
  });

  it('throws for a different provider under an existing name', () => {
    const registry = empty();
    const a = format('a');
    registry.add(a);
    const other = format('a');
    expect(() => registry.add(other)).toThrow(/data format registry.*'a'/i);
    expect(registry.get('a')).toBe(a);
    // the failed add did not change the count
    registry.remove(a);
    expect(registry.has('a')).toBe(false);
  });

  it('findConflict reports without changing anything', () => {
    const registry = empty();
    const a = format('a');
    registry.add(a);
    expect(registry.findConflict(a)).toBeUndefined();
    expect(registry.findConflict(format('b'))).toBeUndefined();
    const other = format('a');
    const message = registry.findConflict(other);
    expect(message).toMatch(/'a'/);
    expect(() => registry.add(other)).toThrow(message);
    expect(names(registry)).toEqual(['a']);
    registry.remove(a);
    expect(registry.has('a')).toBe(false);
  });

  it('legacy add with a matching name registers the provider itself', () => {
    const registry = empty();
    const a = format('a');
    registry.add('a', a);
    expect(registry.get('a')).toBe(a);
    registry.add(a);
    registry.remove(a);
    expect(registry.has('a')).toBe(true);
  });

  it('legacy add with a different name registers a stable named copy and warns once', () => {
    const warn = jest.spyOn(console, 'warn').mockImplementation(() => {});
    try {
      const registry = empty();
      const a = format('a');
      registry.add('b', a);
      const copy = registry.get('b')!;
      expect(copy).not.toBe(a);
      expect(copy.name).toBe('b');
      expect(copy).toBe(DataFormatProvider.withName(a, 'b'));
      expect(warn).toHaveBeenCalledTimes(1);

      registry.add('b', a);
      expect(registry.get('b')).toBe(copy);
      expect(warn).toHaveBeenCalledTimes(1);

      // an original and a renamed copy under the same name conflict
      const original = format('b');
      expect(() => registry.add(original)).toThrow(/'b'/);
      expect(registry.get('b')).toBe(copy);
    } finally {
      warn.mockRestore();
    }
  });

  it('removes by provider and by name', () => {
    const registry = empty();
    const a = format('a');
    const b = format('b');
    registry.add(a);
    registry.add(b);
    registry.add(b);

    registry.remove(a);
    expect(names(registry)).toEqual(['b']);

    registry.remove('b');
    expect(registry.has('b')).toBe(true);
    registry.remove('b');
    expect(registry.has('b')).toBe(false);
  });

  it('remove of an unknown provider or name is a no-op', () => {
    const registry = empty();
    const a = format('a');
    registry.add(a);
    registry.remove(format('b'));
    registry.remove('b');
    // a different object with a registered name does not decrement either
    registry.remove(format('a'));
    expect(names(registry)).toEqual(['a']);
    registry.remove(a);
    expect(registry.list).toHaveLength(0);
  });

  it('has and get', () => {
    const registry = empty();
    const a = format('a');
    registry.add(a);
    expect(registry.has('a')).toBe(true);
    expect(registry.has('b')).toBe(false);
    expect(registry.get('a')).toBe(a);
    expect(() => registry.get('b')).toThrow("Data format 'b' is not registered in this plugin.");
  });

  it('points an empty registry at DefaultRegistry', () => {
    expect(() => empty().get('mmcif')).toThrow(/No data formats are registered.*DefaultRegistry/);
    const registry = empty();
    registry.add(format('a'));
    expect(() => registry.get('mmcif')).not.toThrow(/DefaultRegistry/);
  });

  it('keeps priority-then-registration order for auto', () => {
    const registry = empty();
    const info = { ext: 'x' } as any;
    const data = { data: '' } as any;
    const first = format('first', { stringExtensions: ['x'] });
    const second = format('second', { stringExtensions: ['x'] });
    const urgent = format('urgent', { stringExtensions: ['x'], priority: 1 });
    registry.add(first);
    registry.add(second);
    expect(registry.auto(info, data)).toBe(first);
    registry.add(urgent);
    expect(registry.auto(info, data)).toBe(urgent);
    registry.remove(urgent);
    expect(registry.auto(info, data)).toBe(first);
    registry.remove(first);
    expect(registry.auto(info, data)).toBe(second);
  });

  it('keeps list, extensions and options correct after removals', () => {
    const registry = empty();
    const a = format('a', { stringExtensions: ['aa'], category: 'c1' });
    const b = format('b', { binaryExtensions: ['bb'], category: 'c2' });
    const c = format('c', { stringExtensions: ['cc'] });
    registry.add(a);
    registry.add(b);
    registry.add(c);
    expect([...registry.extensions]).toEqual(['aa', 'bb', 'cc']);
    expect(registry.options.map((o) => o[0])).toEqual(['a', 'b', 'c']);

    registry.remove(b);
    expect(registry.list).toEqual([
      { name: 'a', provider: a },
      { name: 'c', provider: c },
    ]);
    expect([...registry.extensions]).toEqual(['aa', 'cc']);
    expect([...registry.binaryExtensions]).toEqual([]);
    expect(registry.options).toEqual([
      ['a', 'A', 'c1'],
      ['c', 'C', ''],
    ]);
    expect(registry.types).toEqual([
      ['a', 'A'],
      ['c', 'C'],
    ]);
  });
});
