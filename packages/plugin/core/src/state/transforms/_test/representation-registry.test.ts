/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Structure } from '@molstar/model/model/structure';
import { Volume } from '@molstar/model/model/volume';
import { ParticleList } from '@molstar/model/model/particles/particle-list';
import { createFakePlugin, type RegistryKind } from '../../_test/fake-plugin.js';
import { StructureRepresentation3D } from '../structure/representation.js';
import { VolumeRepresentation3D } from '../volume/representation.js';
import { ParticlesRepresentation3D } from '../particles/representation.js';

const kinds: RegistryKind[][] = [
  [],
  ['representations'],
  ['color-themes'],
  ['size-themes'],
  ['representations', 'color-themes', 'size-themes'],
];

const scopes = [
  { name: 'structure', transformer: StructureRepresentation3D, data: Structure.Empty },
  { name: 'volume', transformer: VolumeRepresentation3D, data: Volume.One },
  { name: 'particles', transformer: ParticlesRepresentation3D, data: {} as ParticleList },
] as const;

const emptyParams = {
  type: { name: '', params: {} },
  colorTheme: { name: '', params: {} },
  sizeTheme: { name: '', params: {} },
};

describe('representation transformer param definitions', () => {
  for (const { name, transformer, data } of scopes) {
    for (const empty of kinds) {
      it(`${name}: do not throw with empty [${empty.join(', ')}]`, () => {
        const { plugin } = createFakePlugin(empty);
        const params = transformer.definition.params as any;
        expect(() => PDValues(params(undefined, plugin))).not.toThrow();
        expect(() => PDValues(params({ data }, plugin))).not.toThrow();
      });
    }

    it(`${name}: fall back to empty mapped params when nothing is registered`, () => {
      const { plugin } = createFakePlugin(['representations']);
      const params = (transformer.definition.params as any)(undefined, plugin);
      expect(PDValues(params)).toEqual(emptyParams);
      expect(params.type.select.options.map((o: string[]) => o[0])).toEqual(['']);
    });

    it(`${name}: keep the registry defaults for populated registries`, () => {
      const { plugin } = createFakePlugin();
      const values = PDValues((transformer.definition.params as any)(undefined, plugin));
      expect(values.type.name).toBe(plugin.representation[name].registry.default!.name);
      expect(values.colorTheme.name).not.toBe('');
      expect(values.sizeTheme.name).not.toBe('');
    });
  }
});

describe('representation transformers on empty registries', () => {
  for (const { name, transformer, data } of scopes) {
    const a = { data } as any;
    const params = emptyParams as any;

    it(`${name}: apply throws without representations`, () => {
      const { plugin } = createFakePlugin(['representations']);
      expect(() => (transformer.definition.apply as any)({ a, params }, plugin)).toThrow(
        `No ${name} representations are registered in this plugin`,
      );
    });

    it(`${name}: update throws without representations`, () => {
      const { plugin } = createFakePlugin(['representations']);
      expect(() =>
        (transformer.definition.update as any)({ a, b: {}, oldParams: params, newParams: params }, plugin),
      ).toThrow(`No ${name} representations are registered in this plugin`);
    });

    it(`${name}: apply and update throw without color themes`, () => {
      const { plugin } = createFakePlugin(['color-themes']);
      const message = `No ${name} color themes are registered in this plugin`;
      expect(() => (transformer.definition.apply as any)({ a, params }, plugin)).toThrow(message);
      expect(() =>
        (transformer.definition.update as any)({ a, b: {}, oldParams: params, newParams: params }, plugin),
      ).toThrow(message);
    });

    it(`${name}: apply and update throw without size themes`, () => {
      const { plugin } = createFakePlugin(['size-themes']);
      const message = `No ${name} size themes are registered in this plugin`;
      expect(() => (transformer.definition.apply as any)({ a, params }, plugin)).toThrow(message);
      expect(() =>
        (transformer.definition.update as any)({ a, b: {}, oldParams: params, newParams: params }, plugin),
      ).toThrow(message);
    });
  }
});

function PDValues(params: PD.Params) {
  return PD.getDefaultValues(params);
}
