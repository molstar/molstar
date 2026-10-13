/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginState } from '@molstar/plugin/state';
import { StructureRepresentation3D } from '../transforms/structure/representation.js';
import { VolumeRepresentation3D } from '../transforms/volume/representation.js';
import { ParticlesRepresentation3D } from '../transforms/particles/representation.js';
import { createFakePlugin } from './fake-plugin.js';

interface Names {
  type?: string;
  color?: string;
  size?: string;
}

function node(transformer: string, ref: string, names: Names): any {
  return {
    parent: 'ref-root',
    transformer,
    ref,
    version: '1',
    params: {
      type: { name: names.type ?? 'cartoon', params: {} },
      colorTheme: { name: names.color ?? 'uniform', params: {} },
      sizeTheme: { name: names.size ?? 'uniform', params: {} },
    },
  };
}

function tree(...transforms: any[]): any {
  return { tree: { transforms } };
}

function report(snapshot: Partial<PluginState.Snapshot>, empty: Parameters<typeof createFakePlugin>[0] = []) {
  const { plugin, warnings } = createFakePlugin(empty);
  PluginState.prototype.reportUnregisteredNames.call({ plugin } as any, { id: 'id' as any, ...snapshot });
  return warnings;
}

const Structure = StructureRepresentation3D.id;
const Volume = VolumeRepresentation3D.id;
const Particles = ParticlesRepresentation3D.id;

describe('PluginState.reportUnregisteredNames', () => {
  it('is silent for registered names and empty snapshots', () => {
    expect(report({})).toEqual([]);
    expect(report({ data: tree() })).toEqual([]);
    const data = tree(
      node(Structure, 'a', { type: 'cartoon', color: 'chain-id', size: 'uniform' }),
      node(Volume, 'b', { type: 'isosurface', color: 'uniform', size: 'uniform' }),
      node(Particles, 'c', { type: 'spacefill', color: 'uniform', size: 'uniform' }),
    );
    expect(report({ data })).toEqual([]);
  });

  it('warns for missing names in the data tree, naming the kind and scope', () => {
    const warnings = report({
      data: tree(
        node(Structure, 'a', { type: 'missing-repr', color: 'missing-color', size: 'missing-size' }),
        node(Volume, 'b', { type: 'isosurface', color: 'missing-volume-color' }),
        node(Particles, 'c', { type: 'spacefill', size: 'missing-particles-size' }),
      ),
    });
    expect(warnings).toEqual([
      "Structure representation 'missing-repr' is not registered in this plugin; the registry default is used",
      "Structure color theme 'missing-color' is not registered in this plugin; the registry default is used",
      "Structure size theme 'missing-size' is not registered in this plugin; the registry default is used",
      "Volume color theme 'missing-volume-color' is not registered in this plugin; the registry default is used",
      "Particles size theme 'missing-particles-size' is not registered in this plugin; the registry default is used",
    ]);
  });

  it('warns for missing names in transition frames', () => {
    const warnings = report({
      data: tree(node(Structure, 'a', {})),
      transition: {
        frames: [
          { durationInMs: 1, data: tree(node(Structure, 'a', {})) },
          { durationInMs: 1, data: tree(node(Volume, 'b', { type: 'isosurface', color: 'frame-color' })) },
        ],
      },
    });
    expect(warnings).toEqual([
      "Volume color theme 'frame-color' is not registered in this plugin; the registry default is used",
    ]);
  });

  it('warns once per scope, kind, and name', () => {
    const missing = { color: 'missing-color' };
    const warnings = report({
      data: tree(
        node(Structure, 'a', missing),
        node(Structure, 'b', missing),
        node(Volume, 'c', { ...missing, type: 'isosurface' }),
      ),
      transition: { frames: [{ durationInMs: 1, data: tree(node(Structure, 'a', missing)) }] },
    });
    expect(warnings.length).toBe(2);
  });

  it('checks the name in the matching scope only', () => {
    // 'isosurface' is a volume representation, not a structure one
    expect(report({ data: tree(node(Structure, 'a', { type: 'isosurface' })) }).length).toBe(1);
    expect(report({ data: tree(node(Volume, 'a', { type: 'isosurface' })) })).toEqual([]);
  });

  it('ignores other transformers and missing params', () => {
    const data = tree(
      { parent: 'ref-root', transformer: 'test.other', ref: 'a', version: '1', params: { type: { name: 'x' } } },
      { parent: 'ref-root', transformer: Structure, ref: 'b', version: '1' },
    );
    expect(report({ data })).toEqual([]);
  });

  it('warns for every name when the scope registries are empty', () => {
    const warnings = report({ data: tree(node(Structure, 'a', {})) }, [
      'representations',
      'color-themes',
      'size-themes',
    ]);
    expect(warnings.length).toBe(3);
  });
});
