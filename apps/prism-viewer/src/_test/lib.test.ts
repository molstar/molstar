/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { lib } from '../lib.js';

// `molstar.lib` (app.md §1): v5's keys where the app has them, and the Prism additions.

describe('molstar.lib', () => {
  it("keeps v5's keys, without the classic UI's", () => {
    expect(Object.keys(lib).sort()).toEqual([
      'extensions',
      'loci',
      'math',
      'plugin',
      'prism',
      'shape',
      'structure',
      'volume',
    ]);
    expect(lib.structure.Structure).toBeDefined();
    expect(lib.structure.StructureElement).toBeDefined();
    expect(typeof lib.loci.Loci.isEmpty).toBe('function');
    expect(typeof lib.math.LinearAlgebra.Mat4.identity).toBe('function');
    expect(Object.keys(lib.plugin)).not.toContain('PluginUIContext');
    expect(Object.keys(lib.plugin)).not.toContain('DefaultPluginUISpec');
    expect(typeof lib.plugin.DefaultPluginSpec).toBe('function');
    expect(Object.keys(lib.extensions).sort()).toEqual(['modelArchive', 'mvs', 'plugin', 'wwPDBStructConn']);
    expect(typeof lib.extensions.mvs.createBuilder).toBe('function');
  });

  it("has the Prism viewer's entry points, the node selector builder and MolScript", () => {
    const { prism } = lib;
    expect(typeof prism.createPrismViewer).toBe('function');
    expect(prism.DefaultPrismViewerSpec().plugin).toBeDefined();
    expect(prism.FullPrismViewerSpec().extensions?.length).toBeGreaterThan(
      prism.DefaultPrismViewerSpec().extensions?.length ?? 0,
    );
    expect(typeof prism.node).toBe('function');
    expect(prism.node.ref('x').text).toBe('#x');
    expect(prism.MS.struct.generator.all()).toBeDefined();
  });
});
