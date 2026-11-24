/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it, vi } from 'vitest';
import { setDebugMode } from '@molstar/core/util/debug';
import { Interval } from '@molstar/core/data/int';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Clip } from '@molstar/core/util/clip';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ValueCell } from '@molstar/core/util/value-cell';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { isEmptyLoci } from '@molstar/model/model/loci';
import { Structure, StructureElement, StructureSymmetry, type Unit } from '@molstar/model/model/structure';
import {
  type GraphicsRenderObject,
  getMergedMembers,
  isMergedRenderObject,
  type MergedGraphicsRenderObject,
} from '@molstar/graphics/gl/render-object';
import type { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { Color } from '@molstar/core/util/color';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { Theme } from '@molstar/graphics/theme/theme';
import { Overpaint } from '@molstar/graphics/theme/overpaint';
import { UnitsRepresentation } from '../units-representation.js';
import { UnitsVisual, type UnitsVisualGeometryBuilder } from '../units-visual.js';
import { ElementSphereImpostorVisual, ElementSphereParams } from '../visual/element-sphere.js';
import { MolecularSurfaceMeshParams, MolecularSurfaceMeshVisual } from '../visual/molecular-surface-mesh.js';

// The merged units visual's parts: groups of another merge key (here a define their geometry decides, as a group's
// counts decide the instance granularity) merged into objects of their own, picked and marked as per-group visuals

const registries = (() => {
  const colorThemeRegistry = ColorTheme.createRegistry();
  for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
  const sizeThemeRegistry = SizeTheme.createRegistry();
  for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
  return { colorThemeRegistry, sizeThemeRegistry };
})();

async function load1cbs() {
  const data = fs.readFileSync(
    path.resolve(import.meta.dirname, '../../../../../../data/examples/1cbs_updated.cif'),
    'utf8',
  );
  const parsed = await CIF.parse(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

/** the element spheres' impostors, a define by the group's size: groups of more than 500 atoms (the protein) lit apart */
function splitBuilder() {
  const base = ElementSphereImpostorVisual(-1).builder as UnitsVisualGeometryBuilder<ElementSphereParams, Spheres>;
  return {
    ...base,
    processValues: (values: any, geometry: Spheres) => {
      ValueCell.updateIfChanged(values.dIgnoreLight, geometry.sphereCount > 500);
    },
  };
}

async function createRepr(
  structure: Structure,
  renderMerged: boolean,
  props: Partial<PD.Values<ElementSphereParams>> = {},
  builder: UnitsVisualGeometryBuilder<ElementSphereParams, Spheres> = splitBuilder(),
) {
  // the capabilities the impostors' choice reads (no device)
  const graphics = { capabilities: { fragDepth: true } } as any;
  const repr = UnitsRepresentation(
    'split spheres',
    { graphics, ...registries } as any,
    () => ElementSphereParams,
    (id) => UnitsVisual(builder, id),
  );
  repr.setTheme(
    Theme.create(
      registries as any,
      { structure },
      { colorTheme: { name: 'uniform' }, sizeTheme: { name: 'physical' } },
    ),
  );
  await repr.createOrUpdate({ ...PD.getDefaultValues(ElementSphereParams), renderMerged, ...props }, structure).run();
  return repr;
}

/** each values' color type before `setUnmergeable` first changed it */
const colorTypes = new WeakMap<object, string>();
/** `values` made not to merge (`canMergeValues`: a color type the merged object lacks, the grid's), or back */
function setUnmergeable(values: any, on: boolean) {
  if (!colorTypes.has(values)) colorTypes.set(values, values.colorType.ref.value);
  ValueCell.updateIfChanged(values.colorType, on ? 'volume' : colorTypes.get(values)!);
}

/** the markers of each group's unit by unit id: its record's (merged) or its render object's (per group) */
function markersByUnit(repr: Awaited<ReturnType<typeof createRepr>>) {
  const out: Record<number, string> = {};
  for (const o of repr.renderObjects) {
    const sources = isMergedRenderObject(o) ? (o as MergedGraphicsRenderObject).merged.records : [o.values];
    for (const v of sources) {
      const marker = v.tMarker!.ref.value as { array: Uint8Array; count: number };
      out[v.instanceCount.ref.value * 1e6 + v.uGroupCount.ref.value] =
        `${v.markerStatus.ref.value}/${Array.from(marker.array.subarray(0, marker.count)).join('')}`;
    }
  }
  return out;
}

describe('the merged units visual', async () => {
  const structure = await load1cbs();
  const protein = structure.units.find((u) => u.elements.length > 500)!;
  const unitLoci = (u: Unit) =>
    StructureElement.Loci(structure, [{ unit: u, indices: Interval.ofLength(u.elements.length) }]);

  it('merges the groups of another merge key into a part of their own', async () => {
    const repr = await createRepr(structure, true);
    const merged = repr.renderObjects.filter((o) => isMergedRenderObject(o)) as MergedGraphicsRenderObject[];
    expect(repr.renderObjects.length).toBe(2);
    expect(merged.length).toBe(2);
    const lit = merged.map((o) => o.values.dIgnoreLight.ref.value).sort();
    expect(lit).toEqual([false, true]);
    expect(merged.reduce((n, o) => n + o.merged.records.length, 0)).toBe(structure.unitSymmetryGroups.length);
    repr.destroy();
  });

  it('gives its records out as members, read as the per-group visuals’ objects (geometry export)', async () => {
    const merged = await createRepr(structure, true);
    const single = await createRepr(structure, false);
    const read = (o: GraphicsRenderObject) => {
      const v = o.values as any;
      return {
        type: o.type,
        drawCount: v.drawCount.ref.value,
        instanceCount: v.instanceCount.ref.value,
        groupCount: v.uGroupCount.ref.value,
        transform: Array.from(v.aTransform.ref.value as Float32Array),
        positions: Array.from((v.tPositionGroup.ref.value as { array: Float32Array }).array),
        ignoreLight: v.dIgnoreLight.ref.value,
      };
    };
    const byCount = (a: ReturnType<typeof read>, b: ReturnType<typeof read>) =>
      a.groupCount - b.groupCount || a.instanceCount - b.instanceCount || a.transform[12] - b.transform[12];
    const members = (merged.renderObjects as MergedGraphicsRenderObject[]).flatMap((o) => [...getMergedMembers(o)]);
    expect(members.length).toBe(structure.unitSymmetryGroups.length);
    // the same views while the records are
    const first = merged.renderObjects[0] as MergedGraphicsRenderObject;
    expect(getMergedMembers(first)).toBe(getMergedMembers(first));
    expect(members.map(read).sort(byCount)).toEqual(single.renderObjects.map(read).sort(byCount));
    merged.destroy();
    single.destroy();
  });

  it('picks the groups of each part', async () => {
    const repr = await createRepr(structure, true);
    for (const o of repr.renderObjects) {
      const loci = repr.getLoci({ objectId: o.id, instanceId: 0, groupId: 0 });
      expect(isEmptyLoci(loci)).toBe(false);
    }
    repr.destroy();
  });

  it('makes its objects again for a new context, drawn: the old ones disposed alone', async () => {
    const repr = await createRepr(structure, true);
    const before = repr.renderObjects.slice();
    // a backend switch or a restore: another context, the objects made again for it
    repr.setGraphics({ capabilities: { fragDepth: true } } as any);
    await repr.createOrUpdate().run();
    expect(repr.renderObjects.length).toBe(before.length);
    for (const o of repr.renderObjects) {
      expect(before).not.toContain(o);
      expect(o.state.disposed).toBe(false);
    }
    for (const o of before) expect(o.state.disposed).toBe(true);
    repr.destroy();
  });

  it('moves its clip as the per-group visuals do, its objects kept', async () => {
    const merged = await createRepr(structure, true);
    const single = await createRepr(structure, false);
    const objects = merged.renderObjects.slice();
    const clipOf = (repr: typeof merged) =>
      repr.renderObjects.map((o) => {
        const v = o.values as any;
        return [v.dClipObjects.ref.value, v.uClipObjectCount.ref.value, Array.from(v.uClipObjectFrame.ref.value)];
      });
    for (const x of [0, 5, -3]) {
      const plane: Clip.Props['objects'][number] = {
        type: 'plane',
        invert: false,
        position: Vec3.create(x, 0, 0),
        rotation: { axis: Vec3.create(0, 0, 1), angle: 30 },
        scale: Vec3.create(1, 1, 1),
        transform: Mat4.identity(),
      };
      const clip: Clip.Props = { variant: 'pixel', objects: [plane] };
      await merged.createOrUpdate({ clip } as any).run();
      await single.createOrUpdate({ clip } as any).run();
      const [first] = clipOf(merged);
      expect(first[0]).toBe(true);
      // every object as the per-group visuals' (all alike)
      for (const c of [...clipOf(merged), ...clipOf(single)]) expect(c).toEqual(first);
      expect(merged.renderObjects).toEqual(objects);
      // no record written: a sync reads the members of changed observed cells alone (`Merged.recordsUnchanged`)
      for (const o of objects as MergedGraphicsRenderObject[]) {
        o.merged.sync();
        expect(o.merged.recordsUnchanged!()).toBe(true);
      }
    }
    // a theme: the records written
    merged.setTheme(
      Theme.create(
        registries as any,
        { structure },
        { colorTheme: { name: 'element-symbol' }, sizeTheme: { name: 'physical' } },
      ),
    );
    await merged.createOrUpdate().run();
    for (const o of objects as MergedGraphicsRenderObject[]) expect(o.merged.recordsUnchanged!()).toBe(false);
    merged.destroy();
    single.destroy();
  });

  it('updates its props as a merged visual made with them', async () => {
    const repr = await createRepr(structure, true);
    const objects = repr.renderObjects.slice() as MergedGraphicsRenderObject[];
    /** the objects' values of the props, and every record's cells (but `omit`) */
    const valuesOf = (r: typeof repr, omit: string[] = []) =>
      (r.renderObjects as MergedGraphicsRenderObject[]).map((o) => ({
        emissive: o.values.uEmissive.ref.value,
        sizeFactor: (o.values as any).uSizeFactor.ref.value,
        records: o.merged.records.map((record) => {
          const cells = record as unknown as { [k: string]: { ref: { value: unknown } } };
          const keys = Object.keys(cells).filter((k) => !omit.includes(k));
          return Object.fromEntries(keys.map((k) => [k, cells[k].ref.value]));
        }),
      }));
    for (const o of objects) o.merged.sync();

    // the shared values alone: the records not written
    await repr.createOrUpdate({ emissive: 0.4 }).run();
    for (const o of objects) expect(o.merged.recordsUnchanged!()).toBe(true);
    const emissive = await createRepr(structure, true, { emissive: 0.4 });
    expect(valuesOf(repr)).toEqual(valuesOf(emissive));

    // the records' levels again (by their spheres' sizes): every segment's
    await repr.createOrUpdate({ sizeFactor: 2 }).run();
    for (const o of objects) expect(o.merged.recordsUnchanged!()).toBe(false);
    const sized = await createRepr(structure, true, { emissive: 0.4, sizeFactor: 2 });
    // the bounds of the geometry kept (its grid's too), padded by the new size factor: one made with it is bounded otherwise
    const bounds = ['boundingSphere', 'invariantBoundingSphere', 'uInvariantBoundingSphere', 'instanceGrid'];
    expect(valuesOf(repr, bounds)).toEqual(valuesOf(sized, bounds));
    expect(repr.renderObjects).toEqual(objects);
    for (const r of [repr, emissive, sized]) r.destroy();
  });

  it("keeps each group's update info: a version its geometry was built of, every group's made again when it moves", async () => {
    // as the cartoon's secondary structure version or the interactions': seeded by `initUpdateState`, read and set by
    // `setUpdateState`
    let version = 0;
    const made: Unit.SymmetryGroup[] = [];
    const base = splitBuilder();
    const builder: UnitsVisualGeometryBuilder<ElementSphereParams, Spheres> = {
      ...base,
      createGeometry: (ctx, unit, structure, theme, props, geometry) => {
        made.push(unit as any);
        return base.createGeometry(ctx, unit, structure, theme, props, geometry);
      },
      setUpdateState: (state, ...rest) => {
        base.setUpdateState(state, ...rest);
        if (state.info.version !== version) {
          state.createGeometry = true;
          state.info.version = version;
        }
      },
      initUpdateState: (state) => {
        state.info.version = version;
      },
    };
    const repr = await createRepr(structure, true, {}, builder);
    const groups = structure.unitSymmetryGroups.length;
    expect(made.length).toBe(groups);

    // the version as built: nothing made again
    made.length = 0;
    await repr.createOrUpdate({ emissive: 0.2 }).run();
    expect(made.length).toBe(0);

    // moved: every group's geometry, not the first's alone
    version = 1;
    await repr.createOrUpdate({ emissive: 0.3 }).run();
    expect(made.length).toBe(groups);

    // moved with nothing else changed (the props, theme and structure the same objects: interactions computed again
    // on the same structure): read by each group's update state too
    made.length = 0;
    version = 2;
    await repr.createOrUpdate({}, structure).run();
    expect(made.length).toBe(groups);
    repr.destroy();
  });

  it('does not merge groups whose values do not, its per-group visuals updated and the merged one ended early', async () => {
    let made = 0;
    const base = splitBuilder();
    const builder: UnitsVisualGeometryBuilder<ElementSphereParams, Spheres> = {
      ...base,
      createGeometry: (...args) => {
        made += 1;
        return base.createGeometry(...args);
      },
      // a color type the merged object lacks does not merge (`canMergeValues`)
      processValues: (values: any, geometry: Spheres) => {
        base.processValues!(values, geometry, undefined as any, undefined as any);
        setUnmergeable(values, true);
      },
    };
    const repr = await createRepr(structure, true, {}, builder);
    const groups = structure.unitSymmetryGroups.length;
    expect(repr.renderObjects.some((o) => isMergedRenderObject(o))).toBe(false);
    // the merged visual's first group, then a visual each
    expect(made).toBe(1 + groups);
    const objects = repr.renderObjects.slice();

    // the props again: the per-group visuals updated, not made again, their values deciding (no geometry made to try)
    made = 0;
    await repr.createOrUpdate({ emissive: 0.2 }).run();
    expect(made).toBe(0);
    expect(repr.renderObjects).toEqual(objects);
    expect((repr.renderObjects[0].values as any).uEmissive.ref.value).toBe(0.2);
    // a new structure alike
    await repr.createOrUpdate({}, Structure.ofModel(structure.model)).run();
    expect(made).toBe(0);
    repr.destroy();
  });

  it('merges again where the per-group visuals’ values do after an update', async () => {
    let unmergeable = true;
    const base = splitBuilder();
    const builder: UnitsVisualGeometryBuilder<ElementSphereParams, Spheres> = {
      ...base,
      processValues: (values: any, geometry: Spheres) => {
        base.processValues!(values, geometry, undefined as any, undefined as any);
        setUnmergeable(values, unmergeable);
      },
    };
    const repr = await createRepr(structure, true, {}, builder);
    expect(repr.renderObjects.some((o) => isMergedRenderObject(o))).toBe(false);
    unmergeable = false;
    await repr.createOrUpdate({ emissive: 0.2 }).run();
    expect(repr.renderObjects.length).toBeGreaterThan(0);
    expect(repr.renderObjects.every((o) => isMergedRenderObject(o))).toBe(true);
    repr.destroy();
  });

  it('merges again after a theme change alone (the props the same)', async () => {
    const base = splitBuilder();
    const builder: UnitsVisualGeometryBuilder<ElementSphereParams, Spheres> = {
      ...base,
      // a color type the merged object lacks under one theme
      processValues: (values: any, geometry: Spheres, props, theme) => {
        base.processValues!(values, geometry, props, theme);
        setUnmergeable(values, theme.color.factory === BuiltInColorThemes['element-symbol'].factory);
      },
    };
    const theme = (name: string) =>
      Theme.create(registries as any, { structure }, { colorTheme: { name }, sizeTheme: { name: 'physical' } });
    const repr = await createRepr(structure, true, {}, builder);
    repr.setTheme(theme('element-symbol'));
    await repr.createOrUpdate().run();
    expect(repr.renderObjects.some((o) => isMergedRenderObject(o))).toBe(false);
    repr.setTheme(theme('uniform'));
    await repr.createOrUpdate().run();
    expect(repr.renderObjects.length).toBeGreaterThan(0);
    expect(repr.renderObjects.every((o) => isMergedRenderObject(o))).toBe(true);
    repr.destroy();
  });

  /** `builder` counting the geometries it makes and disposes (by identity: one updated in place counted once) */
  function countingBuilder(builder: UnitsVisualGeometryBuilder<ElementSphereParams, Spheres>) {
    const made = new Set<Spheres>();
    const disposed = new Set<Spheres>();
    const counted: UnitsVisualGeometryBuilder<ElementSphereParams, Spheres> = {
      ...builder,
      createGeometry: (...args) => {
        const geometry = builder.createGeometry(...args) as Spheres;
        made.add(geometry);
        return geometry;
      },
      dispose: (geometry) => {
        disposed.add(geometry);
      },
    };
    const undisposed = () => [...made].filter((g) => !disposed.has(g)).length;
    return { builder: counted, undisposed };
  }

  it('disposes the group it made where its values do not merge, each update', async () => {
    const base = splitBuilder();
    const { builder, undisposed } = countingBuilder({
      ...base,
      processValues: (values: any, geometry: Spheres) => {
        base.processValues!(values, geometry, undefined as any, undefined as any);
        setUnmergeable(values, true);
      },
    });
    const repr = await createRepr(structure, true, {}, builder);
    await repr.createOrUpdate({ emissive: 0.2 }).run();
    await repr.createOrUpdate({ emissive: 0.3 }).run();
    // the per-group visuals' alone not disposed while drawn; the merged visual's tries each disposed
    expect(undisposed()).toBe(structure.unitSymmetryGroups.length);
    repr.destroy();
    expect(undisposed()).toBe(0);
  });

  it('keeps its groups through an update that fails, disposing what that update made', async () => {
    let fail = false;
    let updates = 0;
    const base = splitBuilder();
    const { builder, undisposed } = countingBuilder({
      ...base,
      setUpdateState: (state, ...rest) => {
        base.setUpdateState(state, ...rest);
        // the second group's update of a failing update throws (an abort mid-way)
        if (fail && ++updates === 2) throw new Error('aborted');
      },
    });
    const repr = await createRepr(structure, true, {}, builder);
    // the same groups of another structure (found by group, not the same list): the first group's segment reused
    const again = Structure.ofModel(structure.model);
    fail = true;
    await expect(repr.createOrUpdate({}, again).run()).rejects.toThrow('aborted');
    fail = false;
    // the next update finds every group as it was: none made again and lost
    await repr.createOrUpdate({}, again).run();
    expect(repr.renderObjects.some((o) => isMergedRenderObject(o))).toBe(true);
    repr.destroy();
    expect(undisposed()).toBe(0);
  });

  it('brings each group from what it was written with after an update that failed: its props back drawn as made', async () => {
    let fail = false;
    let updates = 0;
    const base = splitBuilder();
    const builder: UnitsVisualGeometryBuilder<ElementSphereParams, Spheres> = {
      ...base,
      setUpdateState: (state, ...rest) => {
        base.setUpdateState(state, ...rest);
        // the third group's update throws: two written with the failing update's props
        if (fail && ++updates === 3) throw new Error('aborted');
      },
    };
    const drawCounts = (r: Awaited<ReturnType<typeof createRepr>>) =>
      (r.renderObjects as MergedGraphicsRenderObject[]).map((o) => o.merged.records.map((v) => v.drawCount.ref.value));
    const repr = await createRepr(structure, true, {}, builder);
    const made = drawCounts(repr);
    fail = true;
    await expect(repr.createOrUpdate({ stride: 2 }).run()).rejects.toThrow('aborted');
    fail = false;
    // the props it was made with again: the groups written with the other ones made again too
    await repr.createOrUpdate({ stride: 1 }).run();
    expect(drawCounts(repr)).toEqual(made);
    repr.destroy();
  });

  it('keeps what it draws through a rebuild that fails (another instance granularity): its objects not disposed', async () => {
    let fail = false;
    const base = splitBuilder();
    const builder: UnitsVisualGeometryBuilder<ElementSphereParams, Spheres> = {
      ...base,
      createGeometry: (...args) => {
        if (fail) throw new Error('aborted');
        return base.createGeometry(...args);
      },
    };
    const repr = await createRepr(structure, true, {}, builder);
    const objects = repr.renderObjects.slice();
    fail = true;
    await expect(repr.createOrUpdate({ instanceGranularity: 'groupInstance' } as any).run()).rejects.toThrow('aborted');
    fail = false;
    expect(repr.renderObjects).toEqual(objects);
    expect(objects.every((o) => !o.state.disposed)).toBe(true);
    // and made again by the next update
    await repr.createOrUpdate({ instanceGranularity: 'groupInstance' } as any).run();
    expect(repr.renderObjects.every((o) => isMergedRenderObject(o) && !o.state.disposed)).toBe(true);
    expect(objects.every((o) => o.state.disposed)).toBe(true);
    repr.destroy();
  });

  it('checks each group against a new context alone (a backend switch, nothing else changed)', async () => {
    let made = 0;
    const base = splitBuilder();
    const builder: UnitsVisualGeometryBuilder<ElementSphereParams, Spheres> = {
      ...base,
      createGeometry: (...args) => {
        made += 1;
        return base.createGeometry(...args);
      },
      // its geometry not drawn alike by another context (as a mesh's group ids on WebGL1)
      geometryUtils: { ...base.geometryUtils, isCompatible: (_g, graphics) => !(graphics as any).other },
    };
    const repr = await createRepr(structure, true, {}, builder);
    made = 0;
    repr.setGraphics({ capabilities: { fragDepth: true }, other: true } as any);
    await repr.createOrUpdate().run();
    expect(made).toBe(structure.unitSymmetryGroups.length);
    repr.destroy();
  });

  it("adjusts every group's levels to a new size factor, not its part's first group's alone", async () => {
    const levels = (r: Awaited<ReturnType<typeof createRepr>>) =>
      (r.renderObjects as MergedGraphicsRenderObject[])
        .flatMap((o) => o.merged.records.map((rec) => rec.lodLevels!.ref.value as number[][]))
        .map((l) => l.map((level) => level[4]));
    const repr = await createRepr(structure, true, { lodPixelSize: 3 });
    await repr.createOrUpdate({ sizeFactor: 0.25 }).run();
    const sized = await createRepr(structure, true, { lodPixelSize: 3, sizeFactor: 0.25 });
    // the strides adjusted for the size factor (coarser where the spheres are smaller)
    expect(levels(repr)).toEqual(levels(sized));
    expect(levels(repr)).not.toEqual(levels(await createRepr(structure, true, { lodPixelSize: 3 })));
    repr.destroy();
    sized.destroy();
  });

  it('merges the groups it makes again under layers into the parts of their merge keys', async () => {
    const repr = await createRepr(structure, true);
    const loci = unitLoci(protein);
    repr.setState({ overpaint: Overpaint(loci.kind, [{ loci, color: Color(0xff0000), clear: false }]) });
    // its segments all made again (another instance granularity), the overpaint kept
    await repr.createOrUpdate({ instanceGranularity: 'groupInstance' } as any).run();
    const merged = repr.renderObjects as MergedGraphicsRenderObject[];
    expect(merged.length).toBe(2);
    expect(merged.every((o) => isMergedRenderObject(o) && o.values.dOverpaint.ref.value)).toBe(true);
    expect(merged.reduce((n, o) => n + o.merged.records.length, 0)).toBe(structure.unitSymmetryGroups.length);
    repr.destroy();
  });

  const sizedTheme = (s: Structure, size: string) =>
    Theme.create(registries as any, { structure: s }, { colorTheme: { name: 'uniform' }, sizeTheme: { name: size } });
  /** each merged object's sizes as made: their type and the texels of each */
  const sizesOf = (r: Awaited<ReturnType<typeof createRepr>>) =>
    (r.renderObjects as MergedGraphicsRenderObject[])
      .map((o) => `${(o.values as any).dSizeType.ref.value}:${(o.values as any).tSize.ref.value.count}`)
      .sort();

  it('follows a live segment where its part outlives its first: its sizes by a new size theme as made with it', async () => {
    // the second group the first of its part (the protein a part of its own, `splitBuilder`), then gone
    const subset = Structure.create(structure.units.filter((u) => u !== structure.unitSymmetryGroups[1].units[0]));
    const repr = await createRepr(structure, true);
    repr.setTheme(sizedTheme(structure, 'uniform'));
    await repr.createOrUpdate({}, structure).run();
    repr.setTheme(sizedTheme(subset, 'uniform'));
    await repr.createOrUpdate({}, subset).run();
    repr.setTheme(sizedTheme(subset, 'physical'));
    await repr.createOrUpdate({}, subset).run();
    const fresh = await createRepr(subset, true);
    fresh.setTheme(sizedTheme(subset, 'physical'));
    await fresh.createOrUpdate({}, subset).run();
    expect(sizesOf(repr)).toEqual(sizesOf(fresh));
    repr.destroy();
    fresh.destroy();
  });

  it('moves a segment whose levels another size theme moves to the part of its key, as made with it', async () => {
    // levels by each group's largest size (`lodRadius` 0) on screen: one part under a uniform size
    const plain = ElementSphereImpostorVisual(-1).builder as UnitsVisualGeometryBuilder<ElementSphereParams, Spheres>;
    const props = { lodPixelSize: 2, lodRadius: 0 } as Partial<PD.Values<ElementSphereParams>>;
    const repr = await createRepr(structure, true, props, plain);
    repr.setTheme(sizedTheme(structure, 'uniform'));
    await repr.createOrUpdate({}, structure).run();
    expect(repr.renderObjects.length).toBe(1);
    repr.setTheme(sizedTheme(structure, 'physical'));
    await repr.createOrUpdate({}, structure).run();
    const fresh = await createRepr(structure, true, props, plain);
    expect(fresh.renderObjects.length).toBeGreaterThan(1);
    const recordsOf = (r: typeof repr) =>
      (r.renderObjects as MergedGraphicsRenderObject[]).map((o) => o.merged.records.length).sort();
    expect(recordsOf(repr)).toEqual(recordsOf(fresh));
    repr.destroy();
    fresh.destroy();
  });

  it('gives views of its members each with its own instance count (the geometry export)', async () => {
    // groups of different instance counts: the first group's every copy, one of the others'
    const range = await StructureSymmetry.buildSymmetryRange(
      structure,
      Vec3.create(-1, 0, 0),
      Vec3.create(1, 0, 0),
    ).run();
    const firstGroup = range.unitSymmetryGroups[0];
    const s = Structure.create(
      range.unitSymmetryGroups.flatMap((g) => (g === firstGroup ? [...g.units] : [g.units[0]])),
    );
    // one part (not split by size, `splitBuilder`): groups of other instance counts in it
    const plain = ElementSphereImpostorVisual(-1).builder as UnitsVisualGeometryBuilder<ElementSphereParams, Spheres>;
    const repr = await createRepr(s, true, {}, plain);
    // every segment's transform made again (another cell size, `createUnitsTransform`): written through its record
    await repr.createOrUpdate({ cellSize: 13 } as any).run();
    const counts = new Set<number>();
    for (const o of repr.renderObjects as MergedGraphicsRenderObject[]) {
      for (const m of getMergedMembers(o)) {
        expect(m.values.uInstanceCount.ref.value).toBe(m.values.instanceCount.ref.value);
        counts.add(m.values.instanceCount.ref.value);
      }
    }
    expect(counts.size).toBeGreaterThan(1);
    repr.destroy();
  });

  it('tells in debug mode a cell its segments write apart through their shared values', async () => {
    const plain = ElementSphereImpostorVisual(-1).builder as UnitsVisualGeometryBuilder<ElementSphereParams, Spheres>;
    // a cell of the shared values (not a segment's own) written by each segment as its own
    const builder: UnitsVisualGeometryBuilder<ElementSphereParams, Spheres> = {
      ...plain,
      processValues: (values: any, geometry: Spheres) => {
        ValueCell.updateIfChanged(values.uEmissive, geometry.sphereCount);
      },
    };
    const error = vi.spyOn(console, 'error').mockImplementation(() => {});
    setDebugMode(true);
    try {
      const repr = await createRepr(structure, true, {}, builder);
      await repr.createOrUpdate({ sizeFactor: 2 }).run();
      expect(error.mock.calls.some(([m]) => String(m).includes("'uEmissive'"))).toBe(true);
      repr.destroy();
    } finally {
      setDebugMode(false);
      error.mockRestore();
    }
  });

  it('makes its parts again where an update gives segments of one part other keys (by its builder, per segment)', async () => {
    // lit apart by size under one theme alone: a part's segments of two keys after a theme change
    const plain = ElementSphereImpostorVisual(-1).builder as UnitsVisualGeometryBuilder<ElementSphereParams, Spheres>;
    const builder: UnitsVisualGeometryBuilder<ElementSphereParams, Spheres> = {
      ...plain,
      processValues: (values: any, geometry: Spheres, _props, theme) => {
        const split = theme.color.factory === BuiltInColorThemes['element-symbol'].factory;
        ValueCell.updateIfChanged(values.dIgnoreLight, split && geometry.sphereCount > 500);
      },
    };
    const theme = (name: string) =>
      Theme.create(registries as any, { structure }, { colorTheme: { name }, sizeTheme: { name: 'physical' } });
    const lit = (r: Awaited<ReturnType<typeof createRepr>>) =>
      (r.renderObjects as MergedGraphicsRenderObject[]).map((o) => o.values.dIgnoreLight.ref.value).sort();
    const repr = await createRepr(structure, true, {}, builder);
    expect(lit(repr)).toEqual([false]);
    repr.setTheme(theme('element-symbol'));
    await repr.createOrUpdate({}, structure).run();
    const fresh = await createRepr(structure, true, {}, builder);
    fresh.setTheme(theme('element-symbol'));
    await fresh.createOrUpdate({}, structure).run();
    expect(lit(fresh)).toEqual([false, true]);
    expect(lit(repr)).toEqual(lit(fresh));
    repr.destroy();
    fresh.destroy();
  });

  it('merges a molecular surface whose colors are smoothed: per vertex (no grid an object)', async () => {
    const graphics = { capabilities: { fragDepth: true } } as any;
    const repr = UnitsRepresentation(
      'surface',
      { graphics, ...registries } as any,
      () => MolecularSurfaceMeshParams,
      (id) => MolecularSurfaceMeshVisual(id),
    );
    repr.setTheme(
      Theme.create(
        registries as any,
        { structure },
        { colorTheme: { name: 'chain-id' }, sizeTheme: { name: 'physical' } },
      ),
    );
    const props = {
      ...PD.getDefaultValues(MolecularSurfaceMeshParams),
      renderMerged: true,
      resolution: 1,
      smoothColors: { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } },
    } as PD.Values<MolecularSurfaceMeshParams>;
    await repr.createOrUpdate(props, structure).run();
    const objects = repr.renderObjects;
    expect(objects.length).toBeGreaterThan(0);
    expect(objects.every((o) => isMergedRenderObject(o))).toBe(true);
    expect(objects.some((o) => o.values.colorType.ref.value === 'vertex')).toBe(true);
    repr.destroy();
  });

  it('marks as the per-group visuals do', async () => {
    const merged = await createRepr(structure, true);
    const single = await createRepr(structure, false);
    expect(single.renderObjects.some((o) => isMergedRenderObject(o))).toBe(false);
    for (const [loci, action] of [
      [unitLoci(protein), MarkerAction.Highlight],
      [unitLoci(protein), MarkerAction.RemoveHighlight],
      [unitLoci(structure.units[structure.units.length - 1]), MarkerAction.Select],
    ] as const) {
      expect(merged.mark(loci, action)).toBe(single.mark(loci, action));
      expect(markersByUnit(merged)).toEqual(markersByUnit(single));
    }
    merged.destroy();
    single.destroy();
  });
});
