/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomProperties } from '@molstar/model/model/custom-property';
import {
  Particle,
  type ParticleList,
  type ParticleTarget,
  ParticleValues,
} from '@molstar/model/model/particles/particle-list';
import { ParticleSlots } from '@molstar/model/model/particles/particle-slots';
import { ParticleTileView, type ParticleTileViewSource } from '@molstar/model/model/particles/particle-tile-view';
import { type ParticleTileSet, ParticleTiles } from '@molstar/model/model/particles/particle-tiles';
import type { ColorTheme } from '../../color.js';
import { ParticleAttributeColorThemeProvider } from '../particle-attribute.js';
import { ParticleCompartmentColorThemeProvider } from '../particle-compartment.js';
import { ParticleEntityColorThemeProvider } from '../particle-entity.js';
import { ParticleHierarchyColorThemeProvider } from '../particle-hierarchy.js';
import { ParticleIndexColorThemeProvider } from '../particle-index.js';
import { fillParticleColors } from '../particle-pages.js';

const Count = 120;

/** particles over 0-100 Å, of 3 targets, 4 entities and 2 compartments, an attribute */
function createList(): ParticleList {
  const coordinates = new Float32Array(Count * 3);
  for (let i = 0; i < Count; ++i) coordinates.set([(i * 37) % 101, (i * 53) % 97, (i * 11) % 89], i * 3);
  const sites = (n: number): ParticleTarget => ({
    kind: 'sites',
    sites: { count: n, positions: new Float32Array(n * 3), radii: new Float32Array(n).fill(1) },
  });
  return {
    label: 'test',
    count: Count,
    targets: Int32Array.from({ length: Count }, (_, i) => 1 + (Math.floor(i / 10) % 3)),
    targetInfo: new Map([
      [1, { name: 'one' }],
      [2, { name: 'two' }],
      [3, { name: 'three' }],
    ]),
    targetMapping: new Map([
      [1, sites(1)],
      [2, sites(1)],
      [3, sites(1)],
    ]),
    entities: Int32Array.from({ length: Count }, (_, i) => i % 4),
    entityInfo: new Map([0, 1, 2, 3].map((e) => [e, { name: `e${e}` }])),
    compartments: Int32Array.from({ length: Count }, (_, i) => (i < 60 ? 0 : 1)),
    compartmentInfo: new Map([
      [0, { name: 'c0' }],
      [1, { name: 'c1' }],
    ]),
    coordinates,
    attributes: new Map([
      [
        'score',
        {
          column: Column.ofFloatArray(Float32Array.from({ length: Count }, (_, i) => i / 10)),
          label: 'Score',
          min: 0,
          max: 12,
        },
      ],
    ]),
    getParticleLabel: (i) => `particle ${i}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

const sourceOf = (set: ParticleTileSet): ParticleTileViewSource => ({
  getTile: (id, frame) => set.provider.getTile(id, frame) as ParticleList,
  getTarget: (id) => set.provider.getTarget(id) as ParticleTarget,
});

/** `provider`'s colors of `list`, by index (a slot list's: by slot) */
function colorsOf(provider: ColorTheme.Provider<any, any>, list: ParticleList) {
  const ctx = { particles: list };
  const theme = provider.factory(ctx, PD.getDefaultValues(provider.getParams(ctx)));
  return (i: number) => theme.color(Particle.Location(list, i), false);
}

describe('particle color themes of a slot list', () => {
  const list = createList();
  const set = ParticleTiles.fromList(list, [2, 2, 2]);
  const source = sourceOf(set);
  const ids = [...set.tiles.keys()].sort((a, b) => a - b);
  const flat = ParticleTileView.create(set, source, ids, 0).list;
  const slots = ParticleTileView.createSlots(set, source, ids, 0).list;

  it('applies the themes to a slot list as to its flat equal', () => {
    for (const p of [
      ParticleEntityColorThemeProvider,
      ParticleCompartmentColorThemeProvider,
      ParticleHierarchyColorThemeProvider,
      ParticleAttributeColorThemeProvider,
      ParticleIndexColorThemeProvider,
    ]) {
      expect(p.isApplicable!({ particles: slots }), p.name).toBe(p.isApplicable!({ particles: flat }));
    }
  });

  it('colors a slot as the flat view its particle: by entity, compartment, hierarchy, attribute and index', () => {
    for (const p of [
      ParticleEntityColorThemeProvider,
      ParticleCompartmentColorThemeProvider,
      ParticleHierarchyColorThemeProvider,
      ParticleAttributeColorThemeProvider,
      ParticleIndexColorThemeProvider,
    ]) {
      const f = colorsOf(p, flat);
      const s = colorsOf(p, slots);
      for (let i = 0; i < flat.count; ++i) expect(s(flat.keys![i]), `${p.name} ${i}`).toBe(f(i));
    }
  });

  it('colors a slot by its compartment and hierarchy alike in views of other tiles', () => {
    const fewer = ParticleTileView.createSlots(set, source, ids.slice(0, 4), 0).list;
    for (const p of [ParticleCompartmentColorThemeProvider, ParticleHierarchyColorThemeProvider]) {
      const all = colorsOf(p, slots);
      const some = colorsOf(p, fewer);
      // a particle's color by what it is: the same (compartment, entity) the same color
      const byKind = new Map<string, number>();
      for (const id of ids) {
        const tile = set.tiles.get(id)!;
        for (let k = 0; k < tile.count; ++k) {
          const slot = tile.base + k;
          const local = set.provider.getTile(id, 0) as ParticleList;
          const kind = `${ParticleValues.get(local.compartments!, k)}/${ParticleValues.get(local.entities!, k)}`;
          const c = all(slot);
          if (byKind.has(kind)) expect(c).toBe(byKind.get(kind));
          else byKind.set(kind, c);
        }
      }
      // the compartments by their info: the same in a view of fewer tiles
      if (p === ParticleCompartmentColorThemeProvider) {
        for (const id of ids.slice(0, 4)) {
          const tile = set.tiles.get(id)!;
          for (let k = 0; k < tile.count; ++k) expect(some(tile.base + k)).toBe(all(tile.base + k));
        }
      }
    }
  });

  it('colors by compartment and hierarchy in the catalog order: a view of some tiles, flat or of slots, as the whole', () => {
    const some = ids.filter((_, k) => k % 3 !== 1);
    const flatSome = ParticleTileView.create(set, source, some, 0).list;
    const slotsSome = ParticleTileView.createSlots(set, source, some, 0).list;
    // the catalog's compartments hold their entities
    expect(set.catalog.compartmentInfo!.get(0)!.entities).toEqual([0, 1, 2, 3]);
    for (const p of [ParticleCompartmentColorThemeProvider, ParticleHierarchyColorThemeProvider]) {
      const whole = colorsOf(p, slots);
      const f = colorsOf(p, flatSome);
      const s = colorsOf(p, slotsSome);
      for (let i = 0; i < flatSome.count; ++i) {
        const key = flatSome.keys![i];
        expect(f(i), `${p.name} ${i}`).toBe(whole(key));
        expect(s(key), `${p.name} ${i}`).toBe(whole(key));
      }
    }
  });

  it('colors a gap as no particle', () => {
    const some = ParticleTileView.createSlots(set, source, ids.slice(1), 0).list;
    const gap = set.tiles.get(ids[0])!.base;
    const entity = colorsOf(ParticleEntityColorThemeProvider, some);
    const defaultColor = colorsOf(ParticleEntityColorThemeProvider, some)(-1);
    expect(entity(gap)).toBe(defaultColor);
  });
});

describe('particle color themes: colors all at once', () => {
  const list = createList();
  const set = ParticleTiles.fromList(list, [2, 2, 2]);
  const source = sourceOf(set);
  const ids = [...set.tiles.keys()].sort((a, b) => a - b);
  const slots = ParticleTileView.createSlots(set, source, ids, 0).list;
  const providers = [
    ParticleEntityColorThemeProvider,
    ParticleCompartmentColorThemeProvider,
    ParticleHierarchyColorThemeProvider,
    ParticleAttributeColorThemeProvider,
    ParticleIndexColorThemeProvider,
  ];
  const themeOf = (p: ColorTheme.Provider<any, any>, l: ParticleList) =>
    p.factory({ particles: l }, PD.getDefaultValues(p.getParams({ particles: l })));
  const rgbaOf = (color: number) => [(color >> 16) & 255, (color >> 8) & 255, color & 255, 255];

  it("fills a slot list's page as its slots' colors", () => {
    for (const p of providers) {
      const theme = themeOf(p, slots);
      expect(theme.fillParticleColors, p.name).toBeDefined();
      for (const page of slots.pages!.pages) {
        const out = new Uint8Array(page.list.count * 4 + 8);
        fillParticleColors(theme, page.list, out, 8);
        for (let i = 0; i < page.list.count; ++i) {
          const c = theme.color(Particle.Location(slots, page.start + i), false);
          expect(Array.from(out.subarray(8 + i * 4, 12 + i * 4)), `${p.name} ${page.start + i}`).toEqual(rgbaOf(c));
        }
      }
    }
  });

  it("fills a flat list as its particles' colors, another list a location each", () => {
    for (const p of providers) {
      const theme = themeOf(p, list);
      for (const l of [list, set.provider.getTile(ids[0], 0) as ParticleList]) {
        const out = new Uint8Array(l.count * 4);
        fillParticleColors(theme, l, out, 0);
        for (let i = 0; i < l.count; ++i) {
          const c = theme.color(Particle.Location(l, i), false);
          expect(Array.from(out.subarray(i * 4, i * 4 + 4)), `${p.name} ${i}`).toEqual(rgbaOf(c));
        }
      }
    }
  });
});

describe('particle index colors of many slots', () => {
  it('colors a slot list of more slots than a table holds by the table, its fill as its locations', () => {
    // 200,000 slots, a page of 1,000 at 150,000 (its keys its slots)
    const base = createList();
    const n = 1000;
    const page = Particle.subset(
      { ...base, count: Count },
      Int32Array.from({ length: n }, (_, i) => i % Count),
      { keys: Int32Array.from({ length: n }, (_, i) => 150_000 + i) },
    );
    const slots = ParticleSlots.create(base, 200_000, [{ start: 150_000, list: page }]);
    const ctx = { particles: slots };
    const theme = ParticleIndexColorThemeProvider.factory(
      ctx,
      PD.getDefaultValues(ParticleIndexColorThemeProvider.getParams(ctx)),
    );
    const out = new Uint8Array(n * 4);
    fillParticleColors(theme, page, out, 0);
    const words = new Uint32Array(out.buffer);
    const rgba = new Uint8Array(4);
    const word = new Uint32Array(rgba.buffer);
    for (let i = 0; i < n; ++i) {
      const c = theme.color(Particle.Location(slots, 150_000 + i), false);
      rgba.set([(c >> 16) & 255, (c >> 8) & 255, c & 255, 255]);
      expect(words[i]).toBe(word[0]);
    }
    // neighbors in the whole about alike, the ends apart
    const first = theme.color(Particle.Location(slots, 150_000), false);
    expect(theme.color(Particle.Location(slots, 150_001), false)).not.toBe(-1);
    expect(first).not.toBe(theme.color(Particle.Location(slots, 0), false));
  });

  it("gives a target's hierarchy color where its info holds its compartment and entity: its particles' (colors by member)", () => {
    // entities and compartments a target's (a pack's): its info holds them
    const byTarget = new Map([
      [1, { entity: 0, compartment: 0 }],
      [2, { entity: 1, compartment: 1 }],
      [3, { entity: undefined, compartment: 1 }],
    ]);
    const base = createList();
    const targets = base.targets as Int32Array;
    const own: ParticleList = {
      ...base,
      entities: Int32Array.from(targets, (t) => byTarget.get(t)!.entity ?? -1),
      compartments: Int32Array.from(targets, (t) => byTarget.get(t)!.compartment),
      targetInfo: new Map([...base.targetInfo].map(([t, info]) => [t, { ...info, ...byTarget.get(t)! }])),
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    const ownSet = ParticleTiles.fromList(own, [2, 2, 2]);
    const ownIds = [...ownSet.tiles.keys()].sort((a, b) => a - b);
    const ownSlots = ParticleTileView.createSlots(ownSet, sourceOf(ownSet), ownIds, 0).list;
    for (const by of ['entity', 'function'] as const) {
      const ctx = { particles: ownSlots };
      const theme = ParticleHierarchyColorThemeProvider.factory(ctx, {
        ...PD.getDefaultValues(ParticleHierarchyColorThemeProvider.getParams(ctx)),
        by,
      });
      expect(theme.targetColor).toBeDefined();
      const flat = ParticleTileView.create(ownSet, sourceOf(ownSet), ownIds, 0).list;
      for (let i = 0; i < flat.count; ++i) {
        const slot = flat.keys![i];
        const t = ParticleValues.get(flat.targets, i);
        expect(theme.targetColor!(t), `${by} ${i}`).toBe(theme.color(Particle.Location(ownSlots, slot), false));
      }
    }
    // without a target's compartment in its info: none (by particle)
    const plainSet = ParticleTiles.fromList(base, [2, 2, 2]);
    const plainIds = [...plainSet.tiles.keys()].sort((a, b) => a - b);
    const ctx = { particles: ParticleTileView.createSlots(plainSet, sourceOf(plainSet), plainIds, 0).list };
    const plain = ParticleHierarchyColorThemeProvider.factory(
      ctx,
      PD.getDefaultValues(ParticleHierarchyColorThemeProvider.getParams(ctx)),
    );
    expect(plain.targetColor?.(1)).toBeUndefined();
  });
});
