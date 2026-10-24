/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Structure files read headless with Mol*'s parsers (seqspec:libraries#the-mol-adapter): one parse per file, kept as a
// `StructureSource` that the `molstar:structure` complex provider, the structure providers and a host's 3D view share.
// No canvas, no plugin.

import { CIF, type CifFrame } from '@molstar/io/reader/cif';
import { parseMol2 } from '@molstar/io/reader/mol2/parser';
import { parsePDB } from '@molstar/io/reader/pdb/parser';
import { parseSdf } from '@molstar/io/reader/sdf/parser';
import { type Model, Structure, type Trajectory } from '@molstar/model/model/structure';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import { ArrayTrajectory } from '@molstar/model/model/structure/trajectory';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { trajectoryFromMol2 } from '@molstar/model/formats/structure/mol2';
import { trajectoryFromPDB } from '@molstar/model/formats/structure/pdb';
import { trajectoryFromSdf } from '@molstar/model/formats/structure/sdf';
import { Task } from '@molstar/core/task';

export type StructureFormat = 'mmcif' | 'bcif' | 'pdb' | 'sdf' | 'mol2';
export const STRUCTURE_FORMATS: readonly StructureFormat[] = ['mmcif', 'bcif', 'pdb', 'sdf', 'mol2'];

/** A parsed structure file: its trajectory, its CIF block (for categories) and its structures, each made once. */
export interface StructureSource {
  readonly kind: 'molstar:structure-source';
  readonly format: StructureFormat;
  readonly trajectory: Trajectory;
  /** The CIF data block, null for PDB, SDF and MOL2. */
  readonly frame: CifFrame | null;
  /** Structures by model and assembly ("0", "0@1"), made on first use. */
  readonly structures: Map<string, Promise<Structure>>;
  /**
   * The Structure the host's 3D view draws over this trajectory, when the host says: what its computed properties are
   * read from (the secondary structure its cartoon shows, `secondary.ts`). Set by the host as it changes.
   */
  drawn?: Structure;
}

const isSource = (v: unknown): v is StructureSource =>
  typeof v === 'object' && v !== null && (v as StructureSource).kind === 'molstar:structure-source';
export const StructureSource = {
  is: isSource,
  /** A source over a trajectory a host parsed itself (its 3D view's), so `Loaders.parsed` can hand it on. */
  of(trajectory: Trajectory, format: StructureFormat, frame: CifFrame | null = null): StructureSource {
    return { kind: 'molstar:structure-source', format, trajectory, frame, structures: new Map() };
  },
};

async function frameOf(bytes: Uint8Array, binary: boolean): Promise<CifFrame> {
  const parsed = await (binary ? CIF.parseBinary(bytes) : CIF.parseText(new TextDecoder().decode(bytes))).run();
  if (parsed.isError) throw new Error(`not a CIF file: ${parsed.message}`);
  const block = parsed.result.blocks[0];
  if (!block) throw new Error('the CIF file has no data block');
  return block;
}

const parses = new WeakMap<Uint8Array, Map<StructureFormat, Promise<StructureSource>>>();

/** Reads a structure file once per bytes and format (a format entry is asked once per request kind). */
export function readStructure(bytes: Uint8Array, format: StructureFormat): Promise<StructureSource> {
  const byFormat = parses.get(bytes) ?? new Map<StructureFormat, Promise<StructureSource>>();
  parses.set(bytes, byFormat);
  const hit = byFormat.get(format);
  if (hit) return hit;
  const p = parse(bytes, format);
  byFormat.set(format, p);
  return p;
}

async function parse(bytes: Uint8Array, format: StructureFormat): Promise<StructureSource> {
  const text = () => new TextDecoder().decode(bytes);
  if (format === 'mmcif' || format === 'bcif') {
    const frame = await frameOf(bytes, format === 'bcif');
    return StructureSource.of(await trajectoryFromMmCIF(frame).run(), format, frame);
  }
  if (format === 'pdb') {
    const parsed = await parsePDB(text(), 'structure').run();
    if (parsed.isError) throw new Error(`not a PDB file: ${parsed.message}`);
    return StructureSource.of(await trajectoryFromPDB(parsed.result).run(), format);
  }
  if (format === 'sdf') {
    const parsed = await parseSdf(text()).run();
    if (parsed.isError) throw new Error(`not an SDF file: ${parsed.message}`);
    const models: Model[] = [];
    for (const compound of parsed.result.compounds) {
      const t = await trajectoryFromSdf(compound).run();
      for (let i = 0; i < t.frameCount; i++) models.push(await modelAt(t, i));
    }
    if (!models.length) throw new Error('the SDF file has no compound');
    return StructureSource.of(new ArrayTrajectory(models), format);
  }
  const parsed = await parseMol2(text(), 'structure').run();
  if (parsed.isError) throw new Error(`not a MOL2 file: ${parsed.message}`);
  return StructureSource.of(await trajectoryFromMol2(parsed.result).run(), format);
}

async function modelAt(t: Trajectory, i: number): Promise<Model> {
  const f = t.getFrameAtIndex(i);
  return Task.is(f) ? await f.run() : f;
}

/** Model `model` (0-based) as a structure, or its assembly `assembly`; made once per source. */
export function structureOf(source: StructureSource, model: number, assembly: string | null): Promise<Structure> {
  if (model < 0 || model >= source.trajectory.frameCount) {
    return Promise.reject(new Error(`no model ${model}: the file has ${source.trajectory.frameCount}`));
  }
  const key = assembly === null ? `${model}` : `${model}@${assembly}`;
  const hit = source.structures.get(key);
  if (hit) return hit;
  const p = (async () => {
    const s = Structure.ofModel(await modelAt(source.trajectory, model));
    return assembly === null ? s : await StructureSymmetry.buildAssembly(s, assembly).run();
  })();
  source.structures.set(key, p);
  return p;
}

/**
 * A CIF category as a table (MVS's annotation files): each column's cells as strings, `.` and `?` as null, so the
 * kernel types them by their column's role. Null: no such category.
 */
export function categoryTable(
  frame: CifFrame,
  name: string,
): { columns: Record<string, (string | null)[]>; length: number } | null {
  const cat = frame.categories[name];
  if (!cat) return null;
  const columns: Record<string, (string | null)[]> = {};
  for (const f of cat.fieldNames) {
    const field = cat.getField(f);
    if (!field) continue;
    const col: (string | null)[] = new Array(cat.rowCount);
    for (let i = 0; i < cat.rowCount; i++) col[i] = field.valueKind(i) === 0 ? field.str(i) : null;
    columns[f] = col;
  }
  return { columns, length: cat.rowCount };
}
