/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { RuntimeContext } from '@molstar/core/task';
import { CIF, type CifFile } from '@molstar/io/reader/cif';
import { CifWriter } from '@molstar/io/writer/cif';
import {
  createParticleListFromMmcifAssembly,
  getAssemblyIdsFromMmcif,
  type MmcifVariant,
} from '@molstar/model/formats/particles/mmcif';
import { createParticleListFromPack, readParticlePack, writeParticlePack } from '@molstar/model/formats/particles/pack';
import {
  type ParticleList,
  ParticleRotations,
  ParticleSites,
  ParticleValues,
} from '@molstar/model/model/particles/particle-list';
import { AtomicNumbers, VdwRadius } from '@molstar/model/model/structure/model/properties/atomic/measures';
import { type ElementSymbol, getElementFromAtomicNumber } from '@molstar/model/model/structure/model/types';
import { ParticleTiles } from '@molstar/model/model/particles/particle-tiles';

export interface MmcifToPackOptions {
  /** the assembly the particles are (the file's first by default) */
  readonly assemblyId?: string;
  readonly variant?: MmcifVariant;
  /** tile-ordered for a grid of this many tiles an axis (a power of two: `WriteParticlePackOptions.tiles`) */
  readonly tiles?: number;
  /** with `tiles`, the coarse levels (`WriteParticlePackOptions.levels`) */
  readonly levels?: boolean;
}

/** The mmCIF (text or binary) in `data` parsed */
export async function parseMmcif(data: Uint8Array, binary: boolean): Promise<CifFile> {
  const parsed = binary ? await CIF.parseBinary(data).run() : await CIF.parseText(new TextDecoder().decode(data)).run();
  if (parsed.isError) throw new Error(parsed.message);
  return parsed.result;
}

/** The particle list of an mmCIF assembly (cellPACK, PetWorld) as a pack */
export async function mmcifToPack(cif: CifFile, options: MmcifToPackOptions = {}): Promise<Uint8Array> {
  const assemblyId = options.assemblyId ?? getAssemblyIdsFromMmcif(cif)[0];
  if (assemblyId === undefined) throw new Error('the mmCIF has no assembly');
  const list = await createParticleListFromMmcifAssembly(RuntimeContext.Synchronous, cif, {
    assemblyId,
    variant: options.variant,
  });
  return writeParticlePack(list, { tiles: options.tiles, levels: options.levels });
}

/** The pack in `data` written again (tile-ordered with `tiles`, else by target), in the current version */
export function packToPack(
  data: Uint8Array,
  options: { readonly tiles?: number; readonly levels?: boolean } = {},
): Uint8Array {
  return writeParticlePack(createParticleListFromPack(readParticlePack(data)), {
    tiles: options.tiles,
    levels: options.levels,
  });
}

/**
 * The pack in `data` as standard mmCIF, what the mmCIF particle reader (`createParticleListFromMmcifAssembly`) reads
 * back as the same list: a target's sites as its atoms (the site's element where the pack stores elements, else an
 * element whose van der Waals radius is the site's), a

 * particle an operator (its pose placing the target's center: `R`, `p - R c`), a target a gen row (a fiber one each).
 * cellPACK's layout where every target is one chain (its entity's description its compartment and name, the fibers
 * kept), PetWorld's otherwise (a model a target, `PDB_model_num`). `unmatchedRadii`: sites of a radius no element has,
 * written as carbon
 */
export function packToMmcif(data: Uint8Array, binary: boolean): { data: Uint8Array | string; unmatchedRadii: number } {
  let list = createParticleListFromPack(readParticlePack(data));
  // tile-ordered: its particles by tile, a target's not consecutive; written again sorted by target
  if (ParticleTiles.getDirectory(list)) list = createParticleListFromPack(readParticlePack(writeParticlePack(list)));
  return writeParticleListMmcif(list, binary);
}

/** A target as the writer lays it out */
interface TargetRows {
  readonly id: number;
  readonly sites?: ParticleSites;
  readonly center: Vec3;
  /** its chains' labels, unique in the file */
  readonly chains: string[];
  readonly entity: number;
  readonly first: number;
  readonly count: number;
}

export function writeParticleListMmcif(
  list: ParticleList,
  binary: boolean,
): { data: Uint8Array | string; unmatchedRadii: number } {
  const n = list.count;
  // the targets as ranges of the list (a pack's are): their particles consecutive
  const targets: TargetRows[] = [];
  const usedChains = new Set<string>();
  const byId = new Map<number, { first: number; count: number }>();
  // by the targets' runs (a pack's: a run a target)
  ParticleValues.forEachRun(list.targets, n, (start, end, t) => {
    const r = byId.get(t);
    if (!r) byId.set(t, { first: start, count: end - start });
    else if (r.first + r.count === start) r.count += end - start;
    else throw new Error('a target whose particles are not consecutive (write the mmCIF of a pack)');
  });
  for (const [id, { first, count }] of byId) {
    const target = list.targetMapping?.get(id);
    const sites = target?.kind === 'sites' ? target.sites : undefined;
    const center = target?.kind === 'sites' && target.center ? target.center : Vec3();
    const labels = sites?.chains ? sites.chains.labels : [list.targetInfo.get(id)?.name ?? String(id)];
    const chains = labels.map((l) => uniqueChain(l, usedChains));
    targets.push({
      id,
      sites,
      center,
      chains,
      entity: list.entities ? ParticleValues.get(list.entities, first) : -1,
      first,
      count,
    });
  }
  const cellpack = targets.every((t) => t.chains.length === 1);

  // the fibers: a particle's fiber, a gen row each
  const fiberOf = new Int32Array(n).fill(-1);
  if (list.fibers) {
    const { offsets, indices } = list.fibers;
    for (let f = 0; f < list.fibers.count; ++f) {
      for (let k = offsets[f]; k < offsets[f + 1]; ++k) fiberOf[indices[k]] = f;
    }
  }

  const entityName = (e: number) => list.entityInfo?.get(e)?.name;
  const compartmentName = (i: number) => {
    const c = list.compartments ? ParticleValues.get(list.compartments, i) : -1;
    return c >= 0 ? list.compartmentInfo?.get(c)?.name : undefined;
  };

  const encoder = CifWriter.createEncoder({ binary, encoderName: 'particle-pack' });
  encoder.startDataBlock(cellpack ? `${list.entryId ?? 'particles'}_cellpack` : (list.entryId ?? 'particles'));

  // entities: cellPACK's description `compartment.name` (a fiber's compartment ends in `fibers`), a target each
  const entityRows = targets.map((t, k) => {
    const name = entityName(t.entity) ?? list.targetInfo.get(t.id)?.name ?? `target ${t.id}`;
    const compartment = t.count > 0 ? compartmentName(t.first) : undefined;
    return {
      id: String(k + 1),
      description: cellpack && compartment ? `${compartment}.${name}` : name,
      details: list.entityInfo?.get(t.entity)?.function ?? '',
    };
  });
  encoder.writeCategory({
    name: 'entity',
    instance: () => ({
      fields: CifWriter.fields<number, typeof entityRows>()
        .str('id', (i, d) => d[i].id)
        .str('type', () => 'polymer')
        .str('pdbx_description', (i, d) => d[i].description)
        .str('details', (i, d) => d[i].details)
        .getFields(),
      source: [{ data: entityRows, rowCount: entityRows.length }],
    }),
  });

  if (!cellpack) {
    // PetWorld: a model a target, its name the entity's
    encoder.writeCategory({
      name: 'pdbx_model',
      instance: () => ({
        fields: CifWriter.fields<number, typeof entityRows>()
          .int('id', (i) => i + 1)
          .str('name', (i, d) => entityName(targets[i].entity) ?? d[i].description)
          .getFields(),
        source: [{ data: entityRows, rowCount: entityRows.length }],
      }),
    });
  }

  // the sites as atoms, a target's chains in order
  const radiusSymbol = getRadiusSymbols();
  let unmatchedRadii = 0;
  const atoms: { target: number; chain: number; site: number }[] = [];
  for (let k = 0; k < targets.length; ++k) {
    const t = targets[k];
    if (!t.sites) continue;
    const offsets = t.sites.chains?.offsets;
    for (let c = 0; c < t.chains.length; ++c) {
      const start = offsets ? offsets[c] : 0;
      const end = offsets ? offsets[c + 1] : t.sites.count;
      for (let s = start; s < end; ++s) atoms.push({ target: k, chain: c, site: s });
    }
  }
  // a site's element as stored, by its radius where the elements are not known (the first element with it)
  const symbolOf = (a: { target: number; site: number }) => {
    const sites = targets[a.target].sites!;
    const element = ParticleSites.getElement(sites, a.site);
    if (element > 0) return getElementFromAtomicNumber(element) as string;
    const r = Math.fround(ParticleSites.getRadius(sites, a.site));
    const symbol = radiusSymbol.get(r);
    if (symbol === undefined) unmatchedRadii++;
    return symbol ?? 'C';
  };
  const symbols = atoms.map(symbolOf);
  const sitePosition = Vec3();
  const position = (a: { target: number; site: number }, k: number) =>
    ParticleSites.getPosition(targets[a.target].sites!, a.site, sitePosition)[k];
  const float32 = { typedArray: Float32Array, encoder: CifWriter.Encoding.by(CifWriter.Encoding.byteArray) };
  encoder.writeCategory({
    name: 'atom_site',
    instance: () => ({
      fields: CifWriter.fields<number, typeof atoms>()
        .str('group_PDB', () => 'ATOM')
        .index('id')
        .str('type_symbol', (i) => symbols[i])
        .str('label_atom_id', (i) => symbols[i])
        .str('label_comp_id', () => 'UNK')
        .str('label_asym_id', (i, d) => targets[d[i].target].chains[d[i].chain])
        .str('label_entity_id', (i, d) => entityRows[d[i].target].id)
        .int('label_seq_id', () => 1)
        .float('Cartn_x', (i, d) => position(d[i], 0), float32)
        .float('Cartn_y', (i, d) => position(d[i], 1), float32)
        .float('Cartn_z', (i, d) => position(d[i], 2), float32)
        .int('pdbx_PDB_model_num', (i, d) => (cellpack ? 1 : d[i].target + 1))
        .getFields(),
      source: [{ data: atoms, rowCount: atoms.length }],
    }),
  });

  // the operators: a particle's pose placing its target's center, in the list's order (their ids 1 to n)
  const op = Mat4();
  const q = Quat();
  const rowsOf = new Float32Array(n * 12);
  for (let i = 0; i < n; ++i) {
    const t = targets[targetIndexOf(targets, i)];
    ParticleRotations.get(list.rotations, i, q);
    Mat4.fromQuat(op, q);
    const p = GridPositions.get(list.coordinates, i, Vec3());
    const rc = Vec3.transformMat4(Vec3(), t.center, op);
    const o = i * 12;
    for (let r = 0; r < 3; ++r) {
      rowsOf[o + r * 4] = op[r];
      rowsOf[o + r * 4 + 1] = op[4 + r];
      rowsOf[o + r * 4 + 2] = op[8 + r];
      rowsOf[o + r * 4 + 3] = p[r] - rc[r];
    }
  }
  const opFields = CifWriter.fields<number, Float32Array>()
    .str('id', (i) => String(i + 1))
    .str('type', () => 'general operation');
  for (let r = 0; r < 3; ++r) {
    for (let c = 0; c < 3; ++c) opFields.float(`matrix[${r + 1}][${c + 1}]`, (i, d) => d[i * 12 + r * 4 + c], float32);
    opFields.float(`vector[${r + 1}]`, (i, d) => d[i * 12 + r * 4 + 3], float32);
  }
  encoder.writeCategory({
    name: 'pdbx_struct_oper_list',
    instance: () => ({ fields: opFields.getFields(), source: [{ data: rowsOf, rowCount: n }] }),
  });

  // the gen rows: a target's particles runs of operator ids, a fiber a row of its own where it starts (its particles
  // consecutive in a pack), in the list's order: the reader's particles in the same order
  const gen: { expression: string; chains: string; model: number }[] = [];
  for (let k = 0; k < targets.length; ++k) {
    const t = targets[k];
    const chains = t.chains.join(',');
    let runStart = -1;
    const flush = (end: number) => {
      if (runStart < 0) return;
      gen.push({ expression: `(${runStart + 1}-${end})`, chains, model: k + 1 });
      runStart = -1;
    };
    for (let i = t.first; i < t.first + t.count; ++i) {
      const f = fiberOf[i];
      if (f < 0) {
        if (runStart < 0) runStart = i;
        continue;
      }
      flush(i);
      const { offsets, indices } = list.fibers!;
      if (indices[offsets[f]] !== i) continue;
      const ids: number[] = [];
      for (let j = offsets[f]; j < offsets[f + 1]; ++j) ids.push(indices[j] + 1);
      gen.push({ expression: `(${ids.join(',')})`, chains, model: k + 1 });
    }
    flush(t.first + t.count);
  }
  const genFields = CifWriter.fields<number, typeof gen>()
    .str('assembly_id', () => '1')
    .str('oper_expression', (i, d) => d[i].expression)
    .str('asym_id_list', (i, d) => d[i].chains);
  if (!cellpack) genFields.int('PDB_model_num', (i, d) => d[i].model);
  encoder.writeCategory({
    name: 'pdbx_struct_assembly_gen',
    instance: () => ({ fields: genFields.getFields(), source: [{ data: gen, rowCount: gen.length }] }),
  });
  encoder.writeCategory({
    name: 'pdbx_struct_assembly',
    instance: () => ({
      fields: CifWriter.fields<number, undefined>()
        .str('id', () => '1')
        .str('details', () => list.label ?? 'particles')
        .getFields(),
      source: [{ rowCount: 1 }],
    }),
  });

  encoder.encode();
  return { data: encoder.getData(), unmatchedRadii };
}

/** The target `targets` (ranges, in order) holding particle `i` */
function targetIndexOf(targets: readonly TargetRows[], i: number) {
  let lo = 0;
  let hi = targets.length - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (targets[mid].first <= i) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

/** `label`, or one made unique by a suffix where taken */
function uniqueChain(label: string, used: Set<string>) {
  let l = label;
  for (let k = 2; used.has(l); ++k) l = `${label}_${k}`;
  used.add(l);
  return l;
}

/** An element symbol by its van der Waals radius (float32): the first in the table with it */
function getRadiusSymbols() {
  const out = new Map<number, string>();
  for (const symbol of Object.keys(AtomicNumbers)) {
    const r = Math.fround(VdwRadius(symbol as unknown as ElementSymbol));
    if (!out.has(r)) out.set(r, symbol);
  }
  return out;
}

export { writeParticlePack };
