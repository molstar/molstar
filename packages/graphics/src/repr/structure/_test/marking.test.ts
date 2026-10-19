/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { Interval, SortedArray } from '@molstar/core/data/int';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { RuntimeContext } from '@molstar/core/task';
import { AssetManager } from '@molstar/core/util/assets';
import { Color } from '@molstar/core/util/color';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { EmptyLoci, type Loci } from '@molstar/model/model/loci';
import { Bond, Structure, StructureElement, StructureProperties, Unit } from '@molstar/model/model/structure';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import { isMergedRenderObject, type GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { Overpaint } from '@molstar/graphics/theme/overpaint';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { Theme } from '@molstar/graphics/theme/theme';
import { Transparency } from '@molstar/graphics/theme/transparency';
import type { Representation } from '../../representation.js';
import { BallAndStickRepresentationProvider } from '../representation/ball-and-stick.js';
import { CartoonRepresentationProvider } from '../representation/cartoon.js';
import { GaussianSurfaceRepresentationProvider } from '../representation/gaussian-surface.js';

// Marking and layers over the visuals of a representation: the marker, overpaint and transparency arrays each visual
// ends up with, for loci of the shapes the interactive paths see (perf-audit 07 P1, P3-P5). Pinned from the
// implementation before the loci index and the single-pass mark.

const registries = (() => {
  const colorThemeRegistry = ColorTheme.createRegistry();
  for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
  const sizeThemeRegistry = SizeTheme.createRegistry();
  for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
  return { colorThemeRegistry, sizeThemeRegistry };
})();
const propCtx: CustomProperty.Context = { runtime: RuntimeContext.Synchronous, assetManager: new AssetManager() };

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

/** markers of every visual: status, average and array; merged objects also their own (synced) copy */
function markers(repr: Representation<Structure>) {
  const parts: string[] = [];
  for (const o of [...members(repr), ...repr.renderObjects.filter((o) => isMergedRenderObject(o))]) {
    const v = o.values;
    parts.push(`${v.markerStatus.ref.value}/${v.markerAverage.ref.value}/${hash(v.tMarker.ref.value.array)}`);
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
  const array = values[`t${name}`].ref.value.array.subarray(0, usedTexels(values, type) * itemSize);
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
      it(`marks ${reprName} on ${structureName}`, async () => {
        const repr = await createRepr(provider, structure, props);
        const result: Record<string, string[]> = {};
        for (const [name, l] of Object.entries(loci)) result[name] = markSequence(repr, l, other);
        expect(result).toMatchSnapshot();
        repr.destroy();
      });

      it(`applies layers on ${reprName} on ${structureName}`, async () => {
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
