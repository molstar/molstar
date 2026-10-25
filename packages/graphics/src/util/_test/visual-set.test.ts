/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it, vi } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { RuntimeContext, Task } from '@molstar/core/task';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { Color } from '@molstar/core/util/color';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Shape } from '@molstar/model/model/shape';
import { Structure } from '@molstar/model/model/structure';
import { addSphere } from '@molstar/graphics/geo/geometry/mesh/builder/sphere';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Representation } from '@molstar/graphics/repr/representation';
import { CartoonRepresentationProvider } from '@molstar/graphics/repr/structure/representation/cartoon';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { VisualSet, type VisualSetCanvas } from '../visual-set.js';

const graphics = { capabilities: { flatInterpolation: true } } as unknown as GraphicsContext;

function fakeCanvas() {
  const added = new Set<Representation.Any>();
  const canvas = {
    graphics,
    add: vi.fn((r: Representation.Any) => void added.add(r)),
    remove: vi.fn((r: Representation.Any) => void added.delete(r)),
    update: vi.fn(),
    requestDraw: vi.fn(),
  } satisfies VisualSetCanvas;
  return { canvas, added };
}

function sphereShape(_: unknown, n: number) {
  const state = MeshBuilder.createState(256, 256);
  for (let g = 0; g < n; g++) {
    state.currentGroup = g;
    addSphere(state, Vec3.create(g * 3, 0, 0), 1, 1);
  }
  return Shape.create(
    'gizmo',
    {},
    MeshBuilder.getMesh(state),
    () => Color(0x3f7fbf),
    () => 1,
    () => 'handle',
    n,
  );
}

/** A representation that only records what is done to it */
function stubRepresentation(build?: () => Promise<void>) {
  const state = Representation.createState();
  return {
    label: 'stub',
    updated: { next() {} } as any,
    groupCount: 0,
    renderObjects: [],
    geometryVersion: 0,
    props: {},
    params: {},
    get state() {
      return state;
    },
    theme: {} as any,
    createOrUpdate: () => Task.create('stub', async () => build?.()),
    setState: (s: Partial<Representation.State>) => Object.assign(state, s),
    setTheme: () => {},
    getLoci: () => ({ kind: 'empty-loci' }) as any,
    getAllLoci: () => [],
    eachLocation: () => {},
    mark: () => false,
    destroy: vi.fn(),
  } satisfies Representation.Any;
}

describe('VisualSet', () => {
  it('adds a shape built with the canvas graphics, auto-synced, and tells its own loci', async () => {
    const { canvas, added } = fakeCanvas();
    const set = new VisualSet(canvas);
    const repr = await set.addShape(sphereShape, Mesh.Utils, 3, {}, { pickable: false });

    expect(added.has(repr)).toBe(true);
    expect(repr.renderObjects.length).toBeGreaterThan(0);
    expect(repr.state.syncManually).toBe(false);
    expect(repr.state.pickable).toBe(false);
    expect(set.owns({ loci: repr.getAllLoci()[0], repr })).toBe(true);
    expect(set.owns({ loci: repr.getAllLoci()[0], repr: stubRepresentation() })).toBe(false);
    expect(set.owns({ loci: repr.getAllLoci()[0] })).toBe(false);
    set.dispose();
  });

  it('moves a representation through the canvas, its bounds with it, without rebuilding it', async () => {
    const { canvas } = fakeCanvas();
    const set = new VisualSet(canvas);
    const repr = await set.addShape(sphereShape, Mesh.Utils, 1);
    const ids = repr.renderObjects.map((ro) => ro.id);
    const m = Mat4.fromTranslation(Mat4(), Vec3.create(5, 0, 0));

    set.setTransform(repr, m);

    expect(Mat4.areEqual(repr.state.transform, m, 0)).toBe(true);
    // the scene's bounds recomputed (`keepBoundingSphere` not set): a moved gizmo is not culled or clipped by the camera
    expect(canvas.update).toHaveBeenCalledWith(repr);
    const sphere = repr.renderObjects[0].values.boundingSphere.ref.value;
    expect(sphere.center[0]).toBeCloseTo(5, 3);
    expect(repr.renderObjects.map((ro) => ro.id)).toEqual(ids);
  });

  it('removes and destroys one, or all at once; disposed, it takes no more', async () => {
    const { canvas, added } = fakeCanvas();
    const set = new VisualSet(canvas);
    const a = set.add(stubRepresentation());
    const b = set.add(stubRepresentation());
    const c = set.add(stubRepresentation());

    set.remove(a);
    expect(added.has(a)).toBe(false);
    expect(a.destroy).toHaveBeenCalledOnce();
    expect(set.size).toBe(2);

    set.clear();
    expect(added.size).toBe(0);
    expect(b.destroy).toHaveBeenCalledOnce();
    expect(c.destroy).toHaveBeenCalledOnce();

    set.add(stubRepresentation());
    expect(set.size).toBe(1);
    set.dispose();
    expect(set.size).toBe(0);
    expect(() => set.add(stubRepresentation())).toThrow();
    // a representation of someone else is not removed
    const other = stubRepresentation();
    new VisualSet(canvas).remove(other);
    expect(canvas.remove).not.toHaveBeenCalledWith(other);
  });

  it('does not add what a build finishes after a clear, and destroys it', async () => {
    const { canvas, added } = fakeCanvas();
    const set = new VisualSet(canvas);
    let finish!: () => void;
    const repr = stubRepresentation(() => new Promise<void>((resolve) => (finish = resolve)));
    const pending = set.create(() => repr);

    set.clear();
    finish();
    await pending;

    expect(added.size).toBe(0);
    expect(set.size).toBe(0);
    expect(repr.destroy).toHaveBeenCalledOnce();
  });

  it('destroys what fails to build and passes the error on', async () => {
    const { canvas, added } = fakeCanvas();
    const set = new VisualSet(canvas);
    const repr = stubRepresentation(async () => {
      throw new Error('no geometry');
    });

    await expect(set.create(() => repr)).rejects.toThrow('no geometry');
    expect(added.size).toBe(0);
    expect(repr.destroy).toHaveBeenCalledOnce();
  });

  it('draws a molecule outside the state tree (a shadow molecule) and releases its properties', async () => {
    const data = fs.readFileSync(path.resolve(import.meta.dirname, '../../../../../data/examples/1crn.cif'), 'utf8');
    const parsed = await CIF.parse(data).run();
    if (parsed.isError) throw new Error(parsed.message);
    const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
    const structure = Structure.ofModel(trajectory.representative);

    const colorThemeRegistry = ColorTheme.createRegistry();
    for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
    const sizeThemeRegistry = SizeTheme.createRegistry();
    for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);

    const { canvas, added } = fakeCanvas();
    const set = new VisualSet(canvas);
    const detach = vi.spyOn(CartoonRepresentationProvider.ensureCustomProperties!, 'detach');
    const repr = await set.fromProvider(CartoonRepresentationProvider, structure, {
      themes: { colorThemeRegistry, sizeThemeRegistry },
      themeData: { structure },
      properties: { runtime: RuntimeContext.Synchronous, assetManager: new DefaultAssetManager() },
      colorTheme: { name: 'uniform', params: { value: Color(0xcccccc) } },
      state: { pickable: false },
    });

    expect(added.has(repr)).toBe(true);
    expect(repr.renderObjects.length).toBeGreaterThan(0);
    expect(repr.state.pickable).toBe(false);
    expect(repr.theme.color.factory).toBe(colorThemeRegistry.get('uniform').factory);

    set.remove(repr);
    expect(detach).toHaveBeenCalledWith(structure);
  });
});
