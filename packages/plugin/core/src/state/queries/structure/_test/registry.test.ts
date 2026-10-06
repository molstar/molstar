/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptBuilder as MS } from '@molstar/model/script/language/builder';
import { StructureSelectionQuery } from '../query.js';
import { StructureSelectionQueryRegistry } from '../registry.js';
import { AminoAcidSelectionQueries, NucleicBaseSelectionQueries, ResidueSelectionQueries } from '../residue.js';

function query(label: string) {
  return StructureSelectionQuery(label, MS.struct.generator.all(), { category: 'Test' });
}

describe('StructureSelectionQueryRegistry', () => {
  it('preloads the built-in queries without duplicates', () => {
    const registry = new StructureSelectionQueryRegistry();
    expect(registry.list.length).toBeGreaterThan(0);
    expect(new Set(registry.list).size).toBe(registry.list.length);
    expect(registry.options.length).toBe(registry.list.length);
    expect(registry.options.map((o) => o[0])).toEqual(registry.list);
  });

  it('preloads the module-scope residue queries, so registering them again only counts', () => {
    const registry = new StructureSelectionQueryRegistry();
    const size = registry.list.length;
    const version = registry.version;
    for (const q of ResidueSelectionQueries) {
      expect(registry.list).toContain(q);
      registry.add(q);
    }
    expect(registry.list.length).toBe(size);
    expect(registry.version).toBe(version);
    expect(registry.list.filter((q) => AminoAcidSelectionQueries.includes(q)).length).toBe(
      AminoAcidSelectionQueries.length,
    );
    expect(registry.list.filter((q) => NucleicBaseSelectionQueries.includes(q)).length).toBe(
      NucleicBaseSelectionQueries.length,
    );

    // the built-ins are still there after one matching removal each
    for (const q of ResidueSelectionQueries) registry.remove(q);
    expect(registry.list.length).toBe(size);
  });

  it('counts the same query', () => {
    const registry = new StructureSelectionQueryRegistry();
    const size = registry.list.length;
    const q = query('a');
    registry.add(q);
    registry.add(q);
    expect(registry.list.length).toBe(size + 1);
    expect(registry.options.filter((o) => o[0] === q).length).toBe(1);

    registry.remove(q);
    expect(registry.list).toContain(q);
    registry.remove(q);
    expect(registry.list).not.toContain(q);
    expect(registry.list.length).toBe(size);
    expect(registry.options.length).toBe(size);
  });

  it('ignores unknown removals', () => {
    const registry = new StructureSelectionQueryRegistry();
    const size = registry.list.length;
    const version = registry.version;
    registry.remove(query('missing'));
    expect(registry.list.length).toBe(size);
    expect(registry.version).toBe(version);
  });

  it('has no conflicts for queries', () => {
    const registry = new StructureSelectionQueryRegistry();
    const q = query('a');
    expect(registry.findConflict(q)).toBeUndefined();
    registry.add(q);
    expect(registry.findConflict(q)).toBeUndefined();
    expect(registry.findConflict(query('a'))).toBeUndefined();
  });

  it('bumps the version only when the list changes', () => {
    const registry = new StructureSelectionQueryRegistry();
    const q = query('a');
    const v0 = registry.version;
    registry.add(q);
    const v1 = registry.version;
    expect(v1).toBeGreaterThan(v0);
    registry.add(q);
    expect(registry.version).toBe(v1);
    registry.remove(q);
    expect(registry.version).toBe(v1);
    registry.remove(q);
    expect(registry.version).toBeGreaterThan(v1);
  });
});
