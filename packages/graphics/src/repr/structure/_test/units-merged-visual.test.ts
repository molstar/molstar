/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { Interval } from '@molstar/core/data/int';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Clip } from '@molstar/core/util/clip';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ValueCell } from '@molstar/core/util/value-cell';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { isEmptyLoci } from '@molstar/model/model/loci';
import { Structure, StructureElement, type Unit } from '@molstar/model/model/structure';
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
    expect(merged.every((o) => o.members.length === 0)).toBe(true);
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
    const bounds = ['boundingSphere', 'invariantBoundingSphere', 'instanceGrid'];
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
      // a palette's colors do not merge (`canMergeValues`)
      processValues: (values: any, geometry: Spheres) => {
        base.processValues!(values, geometry, undefined as any, undefined as any);
        ValueCell.updateIfChanged(values.dUsePalette, true);
      },
    };
    const repr = await createRepr(structure, true, {}, builder);
    const groups = structure.unitSymmetryGroups.length;
    expect(repr.renderObjects.some((o) => isMergedRenderObject(o))).toBe(false);
    // the merged visual's first group, then a visual each
    expect(made).toBe(1 + groups);
    const objects = repr.renderObjects.slice();

    // the props again: the merged visual tried (its first group), the per-group visuals updated, not made again
    made = 0;
    await repr.createOrUpdate({ emissive: 0.2 }).run();
    expect(made).toBe(1);
    expect(repr.renderObjects).toEqual(objects);
    expect((repr.renderObjects[0].values as any).uEmissive.ref.value).toBe(0.2);
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
