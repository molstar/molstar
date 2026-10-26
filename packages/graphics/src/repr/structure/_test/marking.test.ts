/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { Interval, SortedArray } from '@molstar/core/data/int';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { RuntimeContext } from '@molstar/core/task';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { Color } from '@molstar/core/util/color';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { SegmentedArray } from '@molstar/graphics/gl/segmented-array';
import { EmptyLoci, EveryLoci, type Loci } from '@molstar/model/model/loci';
import { Bond, Structure, StructureElement, StructureProperties, Unit } from '@molstar/model/model/structure';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import {
  canMergeRenderObjects,
  isMergedRenderObject,
  type GraphicsRenderObject,
} from '@molstar/graphics/gl/render-object';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { Clipping } from '@molstar/graphics/theme/clipping';
import { Overpaint } from '@molstar/graphics/theme/overpaint';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { Theme } from '@molstar/graphics/theme/theme';
import { UnitsRepresentation } from '../units-representation.js';
import { GaussianSurfaceMeshParams, GaussianSurfaceMeshVisual } from '../visual/gaussian-surface-mesh.js';
import { Transparency } from '@molstar/graphics/theme/transparency';
import { Wiggle } from '@molstar/graphics/theme/wiggle';
import type { Representation } from '../../representation.js';
import { MarkIndex } from '@molstar/graphics/canvas3d/mark-index';
import { BallAndStickRepresentationProvider } from '../representation/ball-and-stick.js';
import { CartoonRepresentationProvider } from '../representation/cartoon.js';
import { GaussianSurfaceRepresentationProvider } from '../representation/gaussian-surface.js';

// Marking and layers over the visuals of a representation: the marker, overpaint and transparency arrays each visual
// ends up with, for loci of the shapes the interactive paths see (perf-audit 07 P1, P3-P5). Pinned from the
// implementation before the loci index and the single-pass mark; re-pinned where a removal on an all-clear visual
// left a mixed status (the marker arrays unchanged, the statuses exact).

const registries = (() => {
  const colorThemeRegistry = ColorTheme.createRegistry();
  for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
  const sizeThemeRegistry = SizeTheme.createRegistry();
  for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
  return { colorThemeRegistry, sizeThemeRegistry };
})();
const propCtx: CustomProperty.Context = {
  runtime: RuntimeContext.Synchronous,
  assetManager: new DefaultAssetManager(),
};

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

async function createRepr(provider: any, structure: Structure, props: object = {}) {
  await provider.ensureCustomProperties?.attach(propCtx, structure);
  const themeProps = {
    colorTheme: { name: provider.defaultColorTheme.name, params: provider.defaultColorTheme.props ?? {} },
    sizeTheme: { name: provider.defaultSizeTheme.name, params: provider.defaultSizeTheme.props ?? {} },
  };
  await Theme.ensureDependencies(propCtx, registries, { structure }, themeProps);
  const repr = provider.factory({ ...registries }, provider.getParams) as Representation<Structure>;
  repr.setTheme(Theme.create(registries, { structure }, themeProps));
  const params = PD.getDefaultValues(provider.getParams(registries, structure));
  await repr.createOrUpdate({ ...params, ...props }, structure).run();
  return repr;
}

/** loci over the atoms `pick` accepts, a full unit as an interval */
function lociOf(structure: Structure, pick: (l: StructureElement.Location) => boolean) {
  const elements: { unit: Unit; indices: StructureElement.UnitIndex[] | Interval }[] = [];
  const l = StructureElement.Location.create(structure);
  for (const unit of structure.units) {
    l.unit = unit;
    const indices: number[] = [];
    for (let i = 0; i < unit.elements.length; ++i) {
      l.element = unit.elements[i];
      if (pick(l)) indices.push(i);
    }
    if (!indices.length) continue;
    elements.push({
      unit,
      indices: indices.length === unit.elements.length ? Interval.ofLength(indices.length) : (indices as any),
    });
  }
  return StructureElement.Loci(
    structure,
    elements.map(({ unit, indices }) => ({
      unit,
      indices: (Interval.is(indices) ? indices : SortedArray.ofSortedArray(indices)) as any,
    })),
  );
}

function makeLoci(structure: Structure): Record<string, Loci> {
  const seqId = StructureProperties.residue.label_seq_id;
  const symbol = StructureProperties.atom.type_symbol;
  const protein = structure.units.find((u) => Unit.isAtomic(u) && u.polymerElements.length > 0)!;
  const atomic = (l: StructureElement.Location) => Unit.isAtomic(l.unit);
  const bondLocations: Bond.Location[] = [];
  const { bonds } = protein as Unit.Atomic;
  for (let i = 0; i < bonds.edgeCount && bondLocations.length < 20; i += 7) {
    const a = bonds.a[i] as StructureElement.UnitIndex;
    const b = bonds.b[i] as StructureElement.UnitIndex;
    bondLocations.push(Bond.Location(structure, protein, a, structure, protein, b));
    bondLocations.push(Bond.Location(structure, protein, b, structure, protein, a));
  }
  return {
    'one chain': StructureElement.Loci(structure, [
      { unit: protein, indices: Interval.ofLength(protein.elements.length) as any },
    ]),
    'one residue': lociOf(structure, (l) => l.unit === protein && seqId(l) === 42),
    'scattered residues': lociOf(structure, (l) => atomic(l) && seqId(l) % 5 === 0),
    'residue range': lociOf(structure, (l) => l.unit.id % 2 === 0 && seqId(l) >= 20 && seqId(l) < 60),
    carbons: lociOf(structure, (l) => atomic(l) && symbol(l) === 'C'),
    'all but one atom': lociOf(structure, (l) => !(l.unit === protein && l.element === protein.elements[0])),
    'whole structure': lociOf(structure, () => true),
    'structure loci': Structure.Loci(structure),
    bonds: Bond.Loci(structure, bondLocations),
    'empty element loci': StructureElement.Loci(structure, []),
    'empty loci': EmptyLoci,
  };
}

function hash(array: ArrayLike<number> | ArrayBufferView) {
  const bytes = ArrayBuffer.isView(array)
    ? new Uint8Array(array.buffer, array.byteOffset, array.byteLength)
    : Uint8Array.from(array);
  return hashFnv32a(bytes as unknown as number[]).toString(16);
}

function members(repr: Representation<Structure>): GraphicsRenderObject[] {
  const out: GraphicsRenderObject[] = [];
  for (const o of repr.renderObjects) {
    if (isMergedRenderObject(o)) {
      o.merged.sync();
      out.push(...o.members);
    } else {
      out.push(o);
    }
  }
  return out;
}

/** markers of every visual: status, average and elements; merged objects also their own (synced) copy */
function markers(repr: Representation<Structure>) {
  const parts: string[] = [];
  for (const o of [...members(repr), ...repr.renderObjects.filter((o) => isMergedRenderObject(o))]) {
    const v = o.values;
    parts.push(
      `${v.markerStatus.ref.value}/${v.markerAverage.ref.value}/${hash(layerData(v.tMarker.ref.value.array).subarray(0, v.tMarker.ref.value.count))}`,
    );
  }
  return (
    hash(
      parts
        .join(' ')
        .split('')
        .map((c) => c.charCodeAt(0)),
    ) + ` ${parts.length}`
  );
}

/** A layer's data: merged objects' layers are their members' as parts (`SegmentedArray`), composed here */
function layerData<T extends Uint8Array | Float32Array>(array: T): T {
  return SegmentedArray.is(array) ? (SegmentedArray.flatten(array) as T) : array;
}

function isZero(array: Uint8Array) {
  for (let i = 0; i < array.length; ++i) if (array[i]) return false;
  return true;
}

/** the texels a layer of `type` reads: textures are reused, what lies beyond is never read */
function usedTexels(values: GraphicsRenderObject['values'], type: string) {
  const instances = values.instanceCount.ref.value;
  if (type === 'instance') return instances;
  if (type === 'vertexInstance') return values.uVertexCount.ref.value * instances;
  return values.uGroupCount.ref.value * instances;
}

/** a layer as rendered: off, all clear (whatever the layer type) or the type and its data */
function layer(values: GraphicsRenderObject['values'], name: 'Overpaint' | 'Transparency', itemSize: number) {
  if (!values[`d${name}`].ref.value) return 'off';
  const type = values[`${name.toLowerCase()}Type`].ref.value;
  const array = layerData(values[`t${name}`].ref.value.array).subarray(0, usedTexels(values, type) * itemSize);
  return isZero(array) ? 'clear' : `${type}:${hash(array)}`;
}

function overpaints(repr: Representation<Structure>) {
  return members(repr).map((o) => layer(o.values, 'Overpaint', 4));
}

function transparencies(repr: Representation<Structure>) {
  return members(repr).map((o) => {
    const v = o.values;
    const t = layer(v, 'Transparency', 1);
    return `${t} ${v.transparencyAverage.ref.value.toFixed(4)}/${v.transparencyMin.ref.value.toFixed(4)}`;
  });
}

/** the marker state after a hover, a select under a hover, and the reverse actions */
function markSequence(repr: Representation<Structure>, loci: Loci, other: Loci) {
  const steps: string[] = [];
  const step = (name: string, changed: boolean) => steps.push(`${name} ${changed ? '+' : '-'} ${markers(repr)}`);
  step('highlight', repr.mark(loci, MarkerAction.Highlight));
  step('remove highlight', repr.mark(loci, MarkerAction.RemoveHighlight));
  step('select', repr.mark(loci, MarkerAction.Select));
  step('highlight other', repr.mark(other, MarkerAction.Highlight));
  step('toggle', repr.mark(loci, MarkerAction.Toggle));
  step('remove other', repr.mark(other, MarkerAction.RemoveHighlight));
  step('select again', repr.mark(loci, MarkerAction.Select));
  step('deselect', repr.mark(loci, MarkerAction.Deselect));
  step('clear', repr.mark(loci, MarkerAction.Clear));
  return steps;
}

function layerSequence(repr: Representation<Structure>, loci: Loci, other: Loci) {
  const red = Color(0xff0000);
  const blue = Color(0x0000ff);
  const out: Record<string, string[]> = {};
  repr.setState({ overpaint: Overpaint(loci.kind, [{ loci, color: red, clear: false }]) });
  out.overpaint = overpaints(repr);
  repr.setState({
    overpaint: Overpaint(loci.kind, [
      { loci, color: red, clear: false },
      { loci: other, color: blue, clear: false },
    ]),
  });
  out['overpaint + other'] = overpaints(repr);
  repr.setState({ overpaint: Overpaint.Empty });
  out['overpaint cleared'] = overpaints(repr);
  repr.setState({ transparency: Transparency(loci.kind, [{ loci, value: 0.5 }]) });
  out.transparency = transparencies(repr);
  repr.setState({ transparency: Transparency.Empty });
  out['transparency cleared'] = transparencies(repr);
  return out;
}

describe('representation marking and layers', async () => {
  const model = await load1cbs();
  // two unit cells along a: 8 copies of each unit, instanced per symmetry group
  const cells = await StructureSymmetry.buildSymmetryRange(model, Vec3.create(0, 0, 0), Vec3.create(1, 0, 0)).run();

  const reprs = [
    ['ball-and-stick', BallAndStickRepresentationProvider, {}],
    ['ball-and-stick merged', BallAndStickRepresentationProvider, { renderMerged: true }],
    ['cartoon', CartoonRepresentationProvider, {}],
    [
      'gaussian surface',
      GaussianSurfaceRepresentationProvider,
      { smoothColors: { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } } },
    ],
  ] as const;

  for (const [structureName, structure] of [
    ['model', model],
    ['cells', cells],
  ] as const) {
    const loci = makeLoci(structure);
    const other = loci['one residue'];

    for (const [reprName, provider, props] of reprs) {
      // the gaussian surface's layers are smoothed on the CPU, seconds on the cells: its marks on the model are the
      // quick lane's representative, the rest is medium (AGENTS.md, Tests)
      const surface = reprName === 'gaussian surface';
      const marksLane = surface && structureName === 'cells' ? { tags: ['medium'] } : {};
      const layersLane = surface ? { tags: ['medium'] } : {};
      it(`marks ${reprName} on ${structureName}`, marksLane, async () => {
        const repr = await createRepr(provider, structure, props);
        const result: Record<string, string[]> = {};
        for (const [name, l] of Object.entries(loci)) result[name] = markSequence(repr, l, other);
        expect(result).toMatchSnapshot();
        repr.destroy();
      });

      it(`applies layers on ${reprName} on ${structureName}`, layersLane, async () => {
        const repr = await createRepr(provider, structure, props);
        const result: Record<string, Record<string, string[]>> = {};
        for (const name of ['one chain', 'scattered residues', 'carbons', 'whole structure', 'bonds', 'empty loci']) {
          result[name] = layerSequence(repr, loci[name], other);
        }
        expect(result).toMatchSnapshot();
        repr.destroy();
      });
    }
  }
});

describe('layers on surfaces', async () => {
  const model = await load1cbs();
  const props = { smoothColors: { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } } };
  const protein = model.units.find((u) => Unit.isAtomic(u) && u.polymerElements.length > 0)!;
  const chain = StructureElement.Loci(model, [{ unit: protein, indices: Interval.ofLength(protein.elements.length) }]);
  const overpaint = Overpaint(chain.kind, [{ loci: chain, color: Color(0xff0000), clear: false }]);

  /** the layer type of each member and whether it has the color */
  function layerTypes(repr: Representation<Structure>) {
    return members(repr).map((o) => {
      const v = o.values;
      return `${v.overpaintType.ref.value} ${isZero(layerData(v.tOverpaint.ref.value.array)) ? 'clear' : 'set'}`;
    });
  }

  it('smooths the layer only on the visuals it reaches', async () => {
    const repr = await createRepr(GaussianSurfaceRepresentationProvider, model, { ...props, renderMerged: false });
    expect(repr.renderObjects.length).toBeGreaterThan(1);
    repr.setState({ overpaint });
    const types = layerTypes(repr);
    expect(types.filter((t) => t === 'vertexInstance set')).toHaveLength(1);
    expect(types.filter((t) => t !== 'vertexInstance set').every((t) => t === 'groupInstance clear')).toBe(true);
    repr.destroy();
  });

  it('does not smooth members of a merged object, which keep one layer type', async () => {
    const repr = await createRepr(GaussianSurfaceRepresentationProvider, model, { ...props, renderMerged: true });
    const merged = repr.renderObjects.find((o) => isMergedRenderObject(o));
    if (!merged || !isMergedRenderObject(merged)) throw new Error('expected a merged render object');
    repr.setState({ overpaint });
    const types = layerTypes(repr);
    expect(types.filter((t) => t === 'groupInstance set')).toHaveLength(1);
    expect(types.filter((t) => t !== 'groupInstance set').every((t) => t === 'groupInstance clear')).toBe(true);
    expect(canMergeRenderObjects(merged.members)).toBe(true);
    // the merged copy shows the layer
    merged.merged.sync();
    expect(merged.values.overpaintType.ref.value).toBe('groupInstance');
    expect(isZero(layerData(merged.values.tOverpaint.ref.value.array))).toBe(false);
    repr.destroy();
  });

  it('merges visuals recreated with a layer set (perf-audit-2 07/P16)', async () => {
    // gaussian surface visuals that can be made to recreate (as after a GPU <-> CPU switch)
    let recreate = false;
    const visualCtor = (materialId: number) => {
      const visual = GaussianSurfaceMeshVisual(materialId);
      visual.mustRecreate = () => recreate;
      return visual;
    };
    const params = () => GaussianSurfaceMeshParams as any;
    const repr = UnitsRepresentation('test', { ...registries }, params, visualCtor as any);
    repr.setTheme(
      Theme.create(
        registries,
        { structure: model },
        { colorTheme: { name: 'uniform' }, sizeTheme: { name: 'uniform' } },
      ),
    );
    const values = { ...PD.getDefaultValues(GaussianSurfaceMeshParams), ...props, renderMerged: true };
    await repr.createOrUpdate(values, model).run();
    repr.setState({ overpaint });
    expect(repr.renderObjects.some((o) => isMergedRenderObject(o))).toBe(true);

    recreate = true;
    await repr.createOrUpdate({}, model).run();
    expect(repr.renderObjects.some((o) => isMergedRenderObject(o))).toBe(true);
    const types = layerTypes(repr);
    expect(types.filter((t) => t === 'groupInstance set')).toHaveLength(1);
    expect(types.filter((t) => t !== 'groupInstance set').every((t) => t === 'groupInstance clear')).toBe(true);
    repr.destroy();
  });

  it('smooths the layers again when the surface is no longer merged', async () => {
    const repr = await createRepr(GaussianSurfaceRepresentationProvider, model, { ...props, renderMerged: true });
    repr.setState({ overpaint });
    expect(layerTypes(repr).filter((t) => t === 'groupInstance set')).toHaveLength(1);
    await repr.createOrUpdate({ renderMerged: false }).run();
    expect(repr.renderObjects.some((o) => isMergedRenderObject(o))).toBe(false);
    const types = layerTypes(repr);
    expect(types.filter((t) => t === 'vertexInstance set')).toHaveLength(1);
    expect(types.filter((t) => t !== 'vertexInstance set').every((t) => t === 'groupInstance clear')).toBe(true);
    repr.destroy();
  });
});

describe('layer updates and picking', async () => {
  const model = await load1cbs();
  const units = model.units.filter((u) => Unit.isAtomic(u));
  const unitLoci = (u: Unit) =>
    StructureElement.Loci(model, [{ unit: u, indices: Interval.ofLength(u.elements.length) }]);

  it('leaves the layer textures of visuals a layer change does not reach as they are', async () => {
    const repr = await createRepr(BallAndStickRepresentationProvider, model);
    const versions = () =>
      members(repr).map((o) => [o.values.tOverpaint.ref.version, o.values.tTransparency.ref.version].join('/'));
    const [a, b] = units;
    const red = Color(0xff0000);
    repr.setState({ overpaint: Overpaint('element-loci', [{ loci: unitLoci(a), color: red, clear: false }]) });
    repr.setState({ transparency: Transparency('element-loci', [{ loci: unitLoci(a), value: 0.5 }]) });
    const before = versions();
    const op = overpaints(repr);
    // the same layers again: only the visuals of `a` are set again
    repr.setState({ overpaint: Overpaint('element-loci', [{ loci: unitLoci(a), color: red, clear: false }]) });
    repr.setState({ transparency: Transparency('element-loci', [{ loci: unitLoci(a), value: 0.5 }]) });
    const after = versions();
    const reached = members(repr).map((o) => !isZero(layerData(o.values.tOverpaint.ref.value.array)));
    expect(reached.some((r) => r)).toBe(true);
    expect(reached.some((r) => !r)).toBe(true);
    for (let i = 0; i < before.length; ++i) {
      if (reached[i]) expect(after[i]).not.toBe(before[i]);
      else expect(after[i]).toBe(before[i]);
    }
    expect(overpaints(repr)).toEqual(op);
    // a layer on `b`: the visuals of `a` are cleared, those of `b` set
    repr.setState({ overpaint: Overpaint('element-loci', [{ loci: unitLoci(b), color: red, clear: false }]) });
    const moved = members(repr).map((o) => !isZero(layerData(o.values.tOverpaint.ref.value.array)));
    expect(moved).not.toEqual(reached);
    expect(moved.some((r) => r)).toBe(true);
    repr.destroy();
  });

  it('leaves the clipping and wiggle textures of visuals a change does not reach as they are (perf-audit-2 07/P14)', async () => {
    const repr = await createRepr(BallAndStickRepresentationProvider, model);
    const versions = () =>
      members(repr).map((o) => [o.values.tClipping.ref.version, o.values.tWiggle.ref.version].join('/'));
    const [a, b] = units;
    const groups = Clipping.Groups.create(Clipping.Groups.Flag.One);
    const set = (u: Unit) => {
      repr.setState({ clipping: Clipping('element-loci', [{ loci: unitLoci(u), groups }]) });
      repr.setState({ wiggle: Wiggle('element-loci', [{ loci: unitLoci(u), value: 0.5 }]) });
    };
    set(a);
    const before = versions();
    const reached = members(repr).map((o) => !isZero(o.values.tClipping.ref.value.array));
    expect(reached.some((r) => r)).toBe(true);
    expect(reached.some((r) => !r)).toBe(true);
    // the same layers again: only the visuals of `a` are set again
    set(a);
    const after = versions();
    for (let i = 0; i < before.length; ++i) {
      if (reached[i]) expect(after[i]).not.toBe(before[i]);
      else expect(after[i]).toBe(before[i]);
    }
    expect(members(repr).map((o) => !isZero(o.values.tClipping.ref.value.array))).toEqual(reached);
    // a layer on `b`: the visuals of `a` are cleared, those of `b` set
    set(b);
    const moved = members(repr).map((o) => !isZero(o.values.tClipping.ref.value.array));
    const wiggled = members(repr).map((o) => !isZero(new Uint8Array(o.values.tWiggle.ref.value.array.buffer)));
    expect(moved).not.toEqual(reached);
    expect(moved.some((r) => r)).toBe(true);
    expect(wiggled).toEqual(moved);
    repr.destroy();
  });

  it('resolves a pick by its render object', async () => {
    const repr = await createRepr(BallAndStickRepresentationProvider, model);
    for (const o of repr.renderObjects) {
      const loci = repr.getLoci({ objectId: o.id, instanceId: 0, groupId: 0 });
      expect(loci.kind).not.toBe('empty-loci');
    }
    const unknown = Math.max(...repr.renderObjects.map((o) => o.id)) + 1000;
    expect(repr.getLoci({ objectId: unknown, instanceId: 0, groupId: 0 }).kind).toBe('empty-loci');
    repr.destroy();
  });
});

describe('marks through the index of Canvas3D', async () => {
  const root = await load1cbs();
  const seqId = StructureProperties.residue.label_seq_id;
  /** a subset of `parent` with `parent` as its root, as the plugin's components */
  const subset = (parent: Structure, pick: (l: StructureElement.Location) => boolean) =>
    Structure.create(StructureElement.Loci.toStructure(lociOf(parent, pick)).units, { parent });

  it('visits the representations a mark changes: includeParent (`asParent`), a child, not another root', async () => {
    const child = subset(root, (l) => seqId(l) < 50);
    // the child's root is `root`; the proxy `asParent` makes is a root of its own, equivalent to `root` (equal hashes)
    expect(child.root).toBe(root);
    expect(child.asParent().root).not.toBe(root);
    // ball-and-stick's `getData` with `includeParent`: the parent, `child` as its child
    const withParent = await createRepr(BallAndStickRepresentationProvider, child.asParent(), { includeParent: true });
    const ofChild = await createRepr(CartoonRepresentationProvider, child);
    // another root: a structure of its own (no parent), not equivalent to `root`
    const other = await createRepr(CartoonRepresentationProvider, Structure.create(child.units));
    const reprs = [withParent, ofChild, other];
    const index = new MarkIndex();
    for (const r of reprs) index.set(r);

    for (const loci of [lociOf(child, (l) => seqId(l) === 10), lociOf(root, (l) => seqId(l) === 20)]) {
      const visited = new Set<Representation.Any>();
      index.forEach(loci, MarkerAction.Highlight, (r) => visited.add(r));
      expect(visited).toEqual(new Set([withParent, ofChild]));
      // exact: those visited change, the one left out would not
      for (const r of reprs) {
        expect(r.mark(loci, MarkerAction.Highlight)).toBe(visited.has(r));
        r.mark(loci, MarkerAction.RemoveHighlight);
      }
    }
    for (const r of reprs) r.destroy();
  });

  it('deselects all on the representations a select reached, the others hold no selection', async () => {
    const a = await createRepr(CartoonRepresentationProvider, root);
    // a root of its own, not equivalent to `root`
    const b = await createRepr(
      CartoonRepresentationProvider,
      Structure.create(subset(root, (l) => seqId(l) < 50).units),
    );
    const index = new MarkIndex();
    for (const r of [a, b]) index.set(r);
    // as `Canvas3D.applyMark`
    const mark = (loci: Loci, action: MarkerAction) => {
      const visited: Representation.Any[] = [];
      index.forEach(loci, action, (r) => {
        visited.push(r);
        if (r.mark(loci, action)) index.note(r, action);
      });
      return visited;
    };
    expect(
      mark(
        lociOf(root, (l) => seqId(l) === 10),
        MarkerAction.Select,
      ),
    ).toEqual([a]);
    expect(mark(EveryLoci, MarkerAction.Deselect)).toEqual([a]);
    expect(a.renderObjects.every((o) => o.values.markerStatus.ref.value === 0)).toBe(true);
    // the one left out: nothing to deselect
    expect(b.mark(EveryLoci, MarkerAction.Deselect)).toBe(false);
    for (const r of [a, b]) r.destroy();
  });
});

describe('marks on all-clear visuals', async () => {
  const root = await load1cbs();
  const seqId = StructureProperties.residue.label_seq_id;
  const statusOf = (r: Representation<Structure>) => r.renderObjects.map((o) => o.values.markerStatus.ref.value);

  it('reports no change for a removal on a clear visual, and keeps its status clear', async () => {
    const repr = await createRepr(CartoonRepresentationProvider, root);
    const loci = lociOf(root, (l) => seqId(l) >= 25 && seqId(l) <= 35);
    expect(repr.mark(loci, MarkerAction.Highlight)).toBe(true);
    // the highlight removed by a mark of all (as one aimed at the representation), the loci's removal comes after
    expect(repr.mark(EveryLoci, MarkerAction.RemoveHighlight)).toBe(true);
    expect(statusOf(repr).every((s) => s === 0)).toBe(true);
    const versions = repr.renderObjects.map((o) => o.values.tMarker.ref.version);
    expect(repr.mark(loci, MarkerAction.RemoveHighlight)).toBe(false);
    expect(repr.mark(loci, MarkerAction.Deselect)).toBe(false);
    expect(statusOf(repr).every((s) => s === 0)).toBe(true);
    expect(repr.renderObjects.map((o) => o.values.tMarker.ref.version)).toEqual(versions);
    // still marks after
    expect(repr.mark(loci, MarkerAction.Select)).toBe(true);
    expect(statusOf(repr).some((s) => s !== 0)).toBe(true);
    repr.destroy();
  });
});
