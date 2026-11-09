/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import type {
  ParticleCompartmentInfo,
  ParticleEntityInfo,
  ParticleList,
  ParticleSites,
  ParticleTarget,
  ParticleTargetInfo,
} from '@molstar/model/model/particles/particle-list';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { ElementSymbol } from '@molstar/model/model/structure/model/types';
import { AtomNumber, VdwRadiusOfNumber } from '@molstar/model/model/structure/model/properties/atomic/measures';
import { type CifBlock, type CifField, type CifFile, readCifField } from '@molstar/io/reader/cif/data-model';
import { toDatabase } from '@molstar/io/reader/cif/schema';
import { mmCIF_Schema } from '@molstar/io/reader/cif/schema/mmcif';
import {
  type OperatorCombinations,
  composeOperatorCombination,
  parseOperatorCombinations,
  formatOperatorCombination,
  readParticleOperators,
} from './operators.js';
import type { ModelFormat } from '../format.js';
import type { RuntimeContext } from '@molstar/core/task';

export type MmcifVariant = 'auto' | 'cellpack' | 'petworld' | 'standard';

// The content sniff lives apart (mmcif-sniff.ts): a host can ask it without loading the particle code.
export { looksLikeMmcifParticles } from './mmcif-sniff.js';

export interface MmcifAssemblyParticleListOptions {
  /** Assembly identifier from `_pdbx_struct_assembly.id`. */
  readonly assemblyId: string;
  /**
   * Individual chain IDs (`label_asym_id`) to include.
   * An empty or omitted array means "include all chains for this assembly".
   */
  readonly asymIds?: ReadonlyArray<string>;
  readonly label?: string;
  /**
   * Which mmCIF variant to use for building the particle list.
   * - `'auto'` (default): detect automatically from the CIF block header and categories.
   * - `'cellpack'`: interpret `_entity.pdbx_description` as a dot-separated `compartment.entityName` path;
   *   entities below a `.fibers` compartment use assembly-gen rows as ordered fibers.
   * - `'standard'`: like cellpack but treat `_entity.pdbx_description` as a plain entity name only.
   * - `'petworld'`: read the non-standard `_pdbx_struct_assembly_gen.PDB_model_num` field;
   *   one particle per operator combo per gen row.
   */
  readonly variant?: MmcifVariant;
  /**
   * Build the per-target reference geometry from the same mmCIF data and attach it
   * as `ParticleList.targetMapping`. Defaults to `true`.
   */
  readonly resolveTargets?: boolean;
}

/**
 * Return all unique assembly IDs found in `_pdbx_struct_assembly_gen` of the first CIF block.
 * Results are sorted lexicographically.
 */
export function getAssemblyIdsFromMmcif(cifFile: CifFile): string[] {
  const block = cifFile.blocks[0];
  if (!block) return [];
  const db = toDatabase(mmCIF_Schema, block);
  const { assembly_id } = db.pdbx_struct_assembly_gen;
  const ids = new Set<string>();
  for (let i = 0, il = db.pdbx_struct_assembly_gen._rowCount; i < il; i++) {
    ids.add(assembly_id.value(i));
  }
  return Array.from(ids).sort();
}

/**
 * Return the sorted unique individual chain IDs for the given assembly in the first CIF block.
 */
export function getAsymIdsFromMmcif(cifFile: CifFile, assemblyId: string): string[] {
  const block = cifFile.blocks[0];
  if (!block) return [];
  const db = toDatabase(mmCIF_Schema, block);
  const { assembly_id, asym_id_list } = db.pdbx_struct_assembly_gen;
  const ids = new Set<string>();
  for (let i = 0, il = db.pdbx_struct_assembly_gen._rowCount; i < il; i++) {
    if (assembly_id.value(i) !== assemblyId) continue;
    for (const chain of asym_id_list.value(i)) ids.add(chain);
  }
  return Array.from(ids).sort();
}

export async function createParticleListFromMmcifAssembly(
  ctx: RuntimeContext,
  cifFile: CifFile,
  options: MmcifAssemblyParticleListOptions,
): Promise<ParticleList> {
  const block = cifFile.blocks[0];
  if (!block) throw new Error('CIF file contains no data blocks.');
  const variant = resolveVariant(block, options.variant);

  // the targets as their sites (`ParticleTarget` kind `sites`), read from `_atom_site`: not a model and a structure
  // each (millions of atoms in a mesoscale model's targets)
  const built =
    variant === 'petworld'
      ? await buildPetworldParticleList(ctx, cifFile, block, options)
      : await buildCellpackStandardParticleList(ctx, cifFile, block, options, variant);

  if (built.targetMapping.size === 0) return built.list;
  return { ...built.list, targetMapping: built.targetMapping };
}

interface BuiltMmcifParticleList {
  readonly list: ParticleList;
  /** Target ID → the reference structure instanced at each particle with that ID. */
  readonly targetMapping: ReadonlyMap<number, ParticleTarget>;
}

/**
 * A group of `_atom_site` rows (a target's): their centroid, their maximum distance from it, their sites (positions and
 * van der Waals radii, `ParticleSites`)
 */
interface AtomSiteGroup {
  readonly center: Vec3;
  readonly radius: number;
  readonly sites: ParticleSites;
}

/** A row's chain as a model's (`findHierarchyOffsets`, `sortAtomSite`): its entity and its `label_asym_id` */
interface AtomSiteChain {
  readonly entity?: CifField;
  readonly asym: CifField;
}

/** Whether rows `a` and `b` are of the same chain: the fields' values compared, not read as strings (a string a row) */
function isSameChain(chain: AtomSiteChain, a: number, b: number) {
  return chain.asym.areValuesEqual(a, b) && (!chain.entity || chain.entity.areValuesEqual(a, b));
}

/**
 * The rows of `_atom_site` grouped by `key` (`onFirst` called with the first row of each), each group's sites by chain
 * (`chain`; the group one chain, labelled by its key, without), the chains in a model's order (by entity, then chain,
 * each in the order it first appears: `sortAtomSite`): in row order where they are already, else reordered (a stable
 * counting sort). Applying the operator of a particle to the centroid of its reference atoms gives a meaningful position even
 * when the operator has zero translation (pure rotation). The columns read once (`readCifField`): not kept decoded
 * with the file (no model shares them)
 */
async function getAtomSiteGroups<K>(
  ctx: RuntimeContext,
  block: CifBlock,
  key: (row: number) => K,
  onFirst?: (row: number, key: K) => void,
  chain?: AtomSiteChain,
): Promise<Map<K, AtomSiteGroup>> {
  const category = block.categories['atom_site'];
  const rowCount = category?.rowCount ?? 0;
  const groups = new Map<K, AtomSiteGroup>();
  if (!category || rowCount === 0) return groups;
  // as Float32, as a model's conformation decodes them
  const coordinate = (name: string) =>
    readCifField(category, name)?.toFloatArray({ array: Float32Array }) ?? new Float32Array(rowCount);
  const x = coordinate('Cartn_x');
  const y = coordinate('Cartn_y');
  const z = coordinate('Cartn_z');
  const typeSymbol = readCifField(category, 'type_symbol');
  const updateChunk = 10000;

  interface Accum {
    key: K;
    x: number;
    y: number;
    z: number;
    n: number;
    radiusSq: number;
    positions: Float32Array;
    radii: Float32Array;
    /** the chains in the order they first appear (by entity and `label_asym_id`): labels, entities' ranks, counts */
    chainLabels: string[];
    chainIndex: Map<string, number>;
    chainEntity: number[];
    chainCounts: number[];
    entityRank: Map<string, number>;
    /** the chain of the group's last row */
    lastChain: number;
    /** a chain's rows not consecutive, or the chains not in a model's order: the sites reordered */
    reorder: boolean;
    /** the next site to fill: the group's in row order, each chain's where reordered */
    next: number;
    chainNext: number[];
    /** the chains in a model's order */
    order: number[];
  }
  const accums: Accum[] = [];
  const byKey = new Map<K, number>();
  // each row's group, for the second pass: the keys read once (a string a row, millions of atoms); a row's key looked
  // up only where it differs from the previous row's (a chain's or a model's rows are consecutive)
  const rowGroup = new Int32Array(rowCount);

  let prevKey: K | undefined;
  let g = -1;
  let c = 0;
  for (let i = 0; i < rowCount; i++) {
    const k = key(i);
    let groupChanged = false;
    if (g < 0 || k !== prevKey) {
      prevKey = k;
      groupChanged = true;
      const known = byKey.get(k);
      if (known === undefined) {
        g = accums.length;
        byKey.set(k, g);
        accums.push({
          key: k,
          x: 0,
          y: 0,
          z: 0,
          n: 0,
          radiusSq: 0,
          positions: EmptyF32,
          radii: EmptyF32,
          chainLabels: [],
          chainIndex: new Map(),
          chainEntity: [],
          chainCounts: [],
          entityRank: new Map(),
          lastChain: -1,
          reorder: false,
          next: 0,
          chainNext: [],
          order: [],
        });
        onFirst?.(i, k);
      } else {
        g = known;
      }
    }
    rowGroup[i] = g;
    const acc = accums[g];
    // a row's chain read and looked up only where it or the group differs from the previous row's
    if (groupChanged || (chain && !isSameChain(chain, i - 1, i))) {
      const asym = chain ? chain.asym.str(i) : '';
      const entity = chain?.entity ? chain.entity.str(i) : '';
      const chainKey = `${entity}\u0000${asym}`;
      const known = acc.chainIndex.get(chainKey);
      if (known === undefined) {
        c = acc.chainLabels.length;
        acc.chainIndex.set(chainKey, c);
        acc.chainLabels.push(chain ? asym : String(k));
        let rank = acc.entityRank.get(entity);
        if (rank === undefined) {
          rank = acc.entityRank.size;
          acc.entityRank.set(entity, rank);
        }
        acc.chainEntity.push(rank);
        acc.chainCounts.push(0);
      } else {
        c = known;
        if (c !== acc.lastChain) acc.reorder = true;
      }
    }
    acc.lastChain = c;
    acc.chainCounts[c]++;
    acc.x += x[i];
    acc.y += y[i];
    acc.z += z[i];
    acc.n++;
    if (i % updateChunk === 0 && ctx.shouldUpdate) {
      await ctx.update({ message: 'Reading particle targets', current: i, max: rowCount * 2 });
    }
  }

  for (const acc of accums) {
    acc.x /= acc.n;
    acc.y /= acc.n;
    acc.z /= acc.n;
    acc.positions = new Float32Array(acc.n * 3);
    acc.radii = new Float32Array(acc.n);
    // by entity, then by first appearance (a stable sort): out of row order where an entity's chains are not consecutive
    const chainCount = acc.chainLabels.length;
    acc.order = Array.from({ length: chainCount }, (_, i) => i).sort((a, b) => acc.chainEntity[a] - acc.chainEntity[b]);
    for (let i = 0; i < chainCount; i++) if (acc.order[i] !== i) acc.reorder = true;
    if (!acc.reorder) continue;
    // each chain's first site, the chains consecutive in their order
    acc.chainNext = new Array(chainCount);
    let offset = 0;
    for (const ci of acc.order) {
      acc.chainNext[ci] = offset;
      offset += acc.chainCounts[ci];
    }
  }

  // the van der Waals radius by element symbol as a model's (`ElementSymbol`, `AtomNumber`), once a symbol
  const radiusOf = new Map<string, number>();
  const radius = (row: number) => {
    const symbol = typeSymbol ? typeSymbol.str(row) : '';
    let r = radiusOf.get(symbol);
    if (r === undefined) {
      r = VdwRadiusOfNumber(AtomNumber(ElementSymbol(symbol)));
      radiusOf.set(symbol, r);
    }
    return r;
  };

  let prevGroup = -1;
  for (let i = 0; i < rowCount; i++) {
    const gi = rowGroup[i];
    const acc = accums[gi];
    const dx = x[i] - acc.x;
    const dy = y[i] - acc.y;
    const dz = z[i] - acc.z;
    const d2 = dx * dx + dy * dy + dz * dz;
    if (d2 > acc.radiusSq) acc.radiusSq = d2;
    let s: number;
    if (acc.reorder) {
      // a row's chain read again only in a group reordered, where it differs from the previous row's
      if (gi !== prevGroup || !isSameChain(chain!, i - 1, i)) {
        const entity = chain!.entity ? chain!.entity.str(i) : '';
        c = acc.chainIndex.get(`${entity}\u0000${chain!.asym.str(i)}`)!;
      }
      s = acc.chainNext[c]++;
    } else {
      s = acc.next++;
    }
    prevGroup = gi;
    acc.positions[s * 3] = x[i];
    acc.positions[s * 3 + 1] = y[i];
    acc.positions[s * 3 + 2] = z[i];
    acc.radii[s] = radius(i);
    if (i % updateChunk === 0 && ctx.shouldUpdate) {
      await ctx.update({ message: 'Reading particle targets', current: rowCount + i, max: rowCount * 2 });
    }
  }

  for (const acc of accums) {
    const { order } = acc;
    const offsets = new Int32Array(order.length + 1);
    for (let k = 0; k < order.length; k++) offsets[k + 1] = offsets[k] + acc.chainCounts[order[k]];
    groups.set(acc.key, {
      center: Vec3.create(acc.x, acc.y, acc.z),
      radius: Math.sqrt(acc.radiusSq),
      sites: {
        count: acc.n,
        positions: acc.positions,
        radii: acc.radii,
        label: String(acc.key),
        chains: { offsets, labels: order.map((ci) => acc.chainLabels[ci]) },
      },
    });
  }
  return groups;
}

const EmptyF32 = new Float32Array(0);

function detectMmcifVariant(block: CifBlock): 'cellpack' | 'petworld' | 'standard' {
  if (block.header.toLowerCase().includes('cellpack')) return 'cellpack';
  if (block.categories['pdbx_struct_assembly_gen']?.getField('PDB_model_num')?.isDefined) return 'petworld';
  return 'standard';
}

function resolveVariant(block: CifBlock, requested?: MmcifVariant): 'cellpack' | 'petworld' | 'standard' {
  if (!requested || requested === 'auto') return detectMmcifVariant(block);
  return requested;
}

async function buildCellpackStandardParticleList(
  ctx: RuntimeContext,
  cifFile: CifFile,
  block: CifBlock,
  options: MmcifAssemblyParticleListOptions,
  variant: 'cellpack' | 'standard',
): Promise<BuiltMmcifParticleList> {
  const db = toDatabase(mmCIF_Schema, block);
  const { pdbx_struct_assembly_gen, pdbx_struct_oper_list } = db;

  if (pdbx_struct_oper_list._rowCount === 0) {
    throw new Error('CIF file contains no _pdbx_struct_oper_list entries.');
  }

  const operators = readParticleOperators(block);

  const { assembly_id, oper_expression, asym_id_list } = pdbx_struct_assembly_gen;
  const asymFilter = options.asymIds?.length ? new Set<string>(options.asymIds) : undefined;

  // Collect matching gen rows and their expanded operator combinations.
  interface GenEntry {
    genIndex: number;
    combinations: OperatorCombinations;
    /** the row's chains (`asym_id_list`, filtered), read once */
    chains: string[];
  }

  const entries: GenEntry[] = [];
  let totalCount = 0;

  for (let i = 0, il = pdbx_struct_assembly_gen._rowCount; i < il; i++) {
    if (assembly_id.value(i) !== options.assemblyId) continue;

    const combinations = parseOperatorCombinations(oper_expression.value(i), operators);
    const chains = asymFilter ? asym_id_list.value(i).filter((c) => asymFilter.has(c)) : asym_id_list.value(i);
    entries.push({ genIndex: i, combinations, chains });
    totalCount += chains.length * combinations.count;
  }

  if (totalCount === 0) {
    if (entries.length === 0) {
      throw new Error(`No _pdbx_struct_assembly_gen rows found for assembly '${options.assemblyId}'.`);
    }
    throw new Error(
      asymFilter !== undefined
        ? `No chains matched assembly '${options.assemblyId}' and asym IDs '${options.asymIds!.join(', ')}'.`
        : `Assembly '${options.assemblyId}' has no expanded operators.`,
    );
  }

  // In the 'standard' variant, the full description is used as the entity name with no compartment.
  // _entity.details holds the functional annotation of the entity, if any.
  // In the 'cellpack' variant, descriptions are dot-separated paths:
  //   compartment = all segments except the last  (e.g. "root.mge.surface.proteins")
  //   entity name = last segment                  (e.g. "MG_191_192_NAP")
  // CellPack entities whose penultimate segment is "fibers" use each assembly-gen row as
  // a separate ordered fiber. In the 'standard' variant, the full description is used as
  // the entity name with no compartment.
  const entityToCompartment = new Map<string, string>();
  const entityToName = new Map<string, string>();
  const fiberEntityIds = new Set<string>();
  const entityToFunction = new Map<string, string>();
  for (let i = 0, il = db.entity._rowCount; i < il; i++) {
    const entityId = db.entity.id.value(i);
    const desc = db.entity.pdbx_description.value(i).join(',');
    if (variant === 'cellpack') {
      const segments = desc.split('.');
      if (segments.length > 1) {
        entityToCompartment.set(entityId, segments.slice(0, -1).join('.'));
        entityToName.set(entityId, segments[segments.length - 1]);
        if (segments[segments.length - 2] === 'fibers') fiberEntityIds.add(entityId);
      } else if (desc) {
        entityToName.set(entityId, desc);
      }
    } else {
      if (desc) entityToName.set(entityId, desc);
    }
    const details = db.entity.details.value(i);
    if (details) entityToFunction.set(entityId, details);
  }

  // Build unique compartment name → index map (populated lazily below).
  const compartmentNameToIdx = new Map<string, number>();
  // Build unique entity name → index map (populated lazily below).
  const entityNameToIdx = new Map<string, number>();
  const entityNameToFunction = new Map<string, string>();

  // read once, not kept decoded with the file (`readCifField`)
  const atomSite = block.categories['atom_site'];
  const siteAsymId = atomSite && readCifField(atomSite, 'label_asym_id');
  const siteEntityId = atomSite && readCifField(atomSite, 'label_entity_id');
  // a chain's entity from its first row, in the same pass
  const chainToEntityId = new Map<string, string>();
  const chainBounds = await getAtomSiteGroups(
    ctx,
    block,
    (i) => siteAsymId?.str(i) ?? '',
    (row, chain) => chainToEntityId.set(chain, siteEntityId?.str(row) ?? ''),
  );
  const updateChunk = 10000;

  let fiberCount = 0;
  let fiberPointCount = 0;
  for (const entry of entries) {
    if (entry.combinations.count < 2) continue;
    for (const chain of entry.chains) {
      const entityId = chainToEntityId.get(chain);
      if (entityId !== undefined && fiberEntityIds.has(entityId)) {
        fiberCount++;
        fiberPointCount += entry.combinations.count;
      }
    }
  }

  const keys = new Int32Array(totalCount);
  const targets = new Int32Array(totalCount);
  const compartments = new Int32Array(totalCount).fill(-1);
  const entities = new Int32Array(totalCount).fill(-1);
  const coordinates = new Float32Array(totalCount * 3);
  const rotations = new Float32Array(totalCount * 4);
  const radii = new Float32Array(totalCount);
  const fiberOffsets = new Int32Array(fiberCount + 1);
  const fiberIndices = new Int32Array(fiberPointCount);

  // Per-particle metadata for labels, kept as indices so no string is allocated per particle.
  const labelEntry = new Int32Array(totalCount);
  const labelCombo = new Int32Array(totalCount);

  const combined = Mat4();
  const quaternion = Quat();
  const centroid = Vec3();
  const position = Vec3();

  // Assign a unique target index per distinct chain ID.
  const chainToTargetIdx = new Map<string, number>();
  let nextTargetIdx = 0;
  for (let ei = 0, eil = entries.length; ei < eil; ei++) {
    for (const chain of entries[ei].chains) {
      if (!chainToTargetIdx.has(chain)) chainToTargetIdx.set(chain, nextTargetIdx++);
    }
  }

  let count = 0;
  let fiberIndex = 0;
  let fiberPosition = 0;
  for (let ei = 0, eil = entries.length; ei < eil; ei++) {
    const entry = entries[ei];

    for (const chain of entry.chains) {
      const bounds = chainBounds.get(chain);
      Vec3.copy(centroid, bounds ? bounds.center : Vec3.origin);
      const radius = bounds ? bounds.radius : 0;
      const entityId = chainToEntityId.get(chain);
      const isFiber = entry.combinations.count >= 2 && entityId !== undefined && fiberEntityIds.has(entityId);

      // the chain's target, compartment and entity, once a chain and not a particle (millions of particles); their
      // indices assigned in the order of the chains' first particles, as before
      const target = chainToTargetIdx.get(chain)!;
      let compartment = -1;
      let entity = -1;
      if (entry.combinations.count > 0) {
        const compartmentName = entityId !== undefined ? entityToCompartment.get(entityId) : undefined;
        if (compartmentName !== undefined) {
          if (!compartmentNameToIdx.has(compartmentName)) {
            compartmentNameToIdx.set(compartmentName, compartmentNameToIdx.size);
          }
          compartment = compartmentNameToIdx.get(compartmentName)!;
        }
        const entityName = entityId !== undefined ? entityToName.get(entityId) : undefined;
        if (entityName !== undefined) {
          if (!entityNameToIdx.has(entityName)) {
            entityNameToIdx.set(entityName, entityNameToIdx.size);
            const fn = entityId !== undefined ? entityToFunction.get(entityId) : undefined;
            if (fn !== undefined) entityNameToFunction.set(entityName, fn);
          }
          entity = entityNameToIdx.get(entityName)!;
        }
      }

      for (let ci = 0, cil = entry.combinations.count; ci < cil; ci++) {
        composeOperatorCombination(operators, entry.combinations, ci, combined);

        Vec3.transformMat4(position, centroid, combined);
        Quat.normalize(quaternion, Quat.fromMat4(quaternion, combined));

        const cOffset = count * 3;
        coordinates[cOffset + 0] = position[0];
        coordinates[cOffset + 1] = position[1];
        coordinates[cOffset + 2] = position[2];

        const qOffset = count * 4;
        rotations[qOffset + 0] = quaternion[0];
        rotations[qOffset + 1] = quaternion[1];
        rotations[qOffset + 2] = quaternion[2];
        rotations[qOffset + 3] = quaternion[3];

        radii[count] = radius * Mat4.getMaxScaleOnAxis(combined);

        keys[count] = count;
        targets[count] = target;
        compartments[count] = compartment;
        entities[count] = entity;

        labelEntry[count] = ei;
        labelCombo[count] = ci;
        if (isFiber) fiberIndices[fiberPosition++] = count;
        count++;
        if (count % updateChunk === 0 && ctx.shouldUpdate) {
          await ctx.update({ message: 'Building particles', current: count, max: totalCount });
        }
      }
      if (isFiber) fiberOffsets[++fiberIndex] = fiberPosition;
    }
  }
  if (ctx.shouldUpdate) {
    await ctx.update({ message: 'Building particles', current: count, max: totalCount });
  }

  // Build the chain mapping: target index → chain ID, and the sites of each target.
  const chainMapping = new Map<number, string>();
  const targetMapping = new Map<number, ParticleTarget>();
  for (const [chain, idx] of chainToTargetIdx) {
    chainMapping.set(idx, chain);
    const group = options.resolveTargets === false ? undefined : chainBounds.get(chain);
    // the centroid the particles were placed by: the atoms end up where the operators put them
    if (group) targetMapping.set(idx, { kind: 'sites', sites: group.sites, center: group.center });
  }

  // Build compartmentInfo: compartment index → compartment info.
  const compartmentInfo = new Map<number, ParticleCompartmentInfo>();
  for (const [name, idx] of compartmentNameToIdx) {
    compartmentInfo.set(idx, { name });
  }

  // Build entityInfo: entity index → entity info.
  const entityInfo = new Map<number, ParticleEntityInfo>();
  for (const [name, idx] of entityNameToIdx) {
    entityInfo.set(idx, { name, function: entityNameToFunction.get(name) });
  }

  // Build targetInfo: target index → chain name and its (single, by construction) entity index.
  const targetInfo = new Map<number, ParticleTargetInfo>();
  for (const [chain, idx] of chainToTargetIdx) {
    const entityId = chainToEntityId.get(chain);
    const entityName = entityId !== undefined ? entityToName.get(entityId) : undefined;
    const entity = entityName !== undefined ? entityNameToIdx.get(entityName) : undefined;
    targetInfo.set(idx, { name: chain, entity });
  }

  const assemblyId = options.assemblyId;
  // Only the ids are captured below, so the label closure does not retain the matrix array.
  const operatorIds = operators.ids;
  const combinations = entries.map((e) => e.combinations);

  return {
    list: {
      label: buildMmcifLabel(options.label, assemblyId),
      count,
      keys,
      targets,
      targetInfo,
      compartments: compartmentInfo.size > 0 ? compartments : undefined,
      compartmentInfo: compartmentInfo.size > 0 ? compartmentInfo : undefined,
      entities: entityInfo.size > 0 ? entities : undefined,
      entityInfo: entityInfo.size > 0 ? entityInfo : undefined,
      coordinates,
      rotations,
      radii,
      fibers: fiberCount > 0 ? { count: fiberCount, offsets: fiberOffsets, indices: fiberIndices } : undefined,
      getParticleLabel: (index: number) => {
        const entity = entityInfo.get(entities[index])?.name;
        const chain = chainMapping.get(targets[index]);
        const opCombo = formatOperatorCombination(operatorIds, combinations[labelEntry[index]], labelCombo[index]);
        return `Particle ${index + 1} · ${entity} · chain ${chain} · ops ${opCombo}`;
      },
      sourceData: MmcifParticleFormat.create(cifFile),
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    },
    targetMapping,
  };
}

async function buildPetworldParticleList(
  ctx: RuntimeContext,
  cifFile: CifFile,
  block: CifBlock,
  options: MmcifAssemblyParticleListOptions,
): Promise<BuiltMmcifParticleList> {
  const db = toDatabase(mmCIF_Schema, block);
  const { pdbx_struct_assembly_gen, pdbx_struct_oper_list } = db;

  if (pdbx_struct_oper_list._rowCount === 0) {
    throw new Error('CIF file contains no _pdbx_struct_oper_list entries.');
  }

  const operators = readParticleOperators(block);
  const { assembly_id, oper_expression } = pdbx_struct_assembly_gen;

  // PDB_model_num is a non-standard field absent from the typed mmCIF_Schema; read it via raw categories.
  const pdbModelNumField = block.categories['pdbx_struct_assembly_gen'].getField('PDB_model_num')!;

  // Read _pdbx_model.name for entity names (non-standard category; rows are in model order, 1-indexed).
  const pdbxModelName = block.categories['pdbx_model']?.getField('name');
  const modelNumToName = new Map<number, string>();
  if (pdbxModelName) {
    for (let i = 0, il = pdbxModelName.rowCount; i < il; i++) {
      modelNumToName.set(i + 1, pdbxModelName.str(i));
    }
  }

  // Collect matching gen rows.
  // One particle per operator combo per gen row = one particle per complex instance.
  interface GenEntry {
    genIndex: number;
    modelNum: number;
    combinations: OperatorCombinations;
  }

  const entries: GenEntry[] = [];
  let totalCount = 0;

  for (let i = 0, il = pdbx_struct_assembly_gen._rowCount; i < il; i++) {
    if (assembly_id.value(i) !== options.assemblyId) continue;
    const combinations = parseOperatorCombinations(oper_expression.value(i), operators);
    const modelNum = pdbModelNumField.int(i);
    entries.push({ genIndex: i, modelNum, combinations });
    totalCount += combinations.count;
  }

  if (totalCount === 0) {
    if (entries.length === 0) {
      throw new Error(`No _pdbx_struct_assembly_gen rows found for assembly '${options.assemblyId}'.`);
    }
    throw new Error(`Assembly '${options.assemblyId}' has no expanded operators.`);
  }

  // Map distinct PDB_model_num values (ascending) → 0-based trajectory model index.
  // mol* `TrajectoryFromMmCif` creates one model per distinct pdbx_PDB_model_num in ascending order,
  // so the rank of a model number equals its trajectory frame index.
  const distinctModelNums = Array.from(new Set(entries.map((e) => e.modelNum))).sort((a, b) => a - b);
  const modelNumToIndex = new Map<number, number>();
  distinctModelNums.forEach((mn, idx) => modelNumToIndex.set(mn, idx));

  const atomSite = block.categories['atom_site'];
  // read once, not kept decoded with the file (`readCifField`)
  const siteModelNum = atomSite && readCifField(atomSite, 'pdbx_PDB_model_num');
  // a model's chains (its units): by entity and `label_asym_id`, else `auth_asym_id`
  const siteLabelAsymId = atomSite && readCifField(atomSite, 'label_asym_id');
  const siteAsym = siteLabelAsymId?.isDefined ? siteLabelAsymId : atomSite && readCifField(atomSite, 'auth_asym_id');
  const siteEntity = atomSite && readCifField(atomSite, 'label_entity_id');
  const modelBounds = await getAtomSiteGroups(
    ctx,
    block,
    (i) => siteModelNum?.int(i) ?? 0,
    undefined,
    siteAsym ? { asym: siteAsym, entity: siteEntity } : undefined,
  );

  // targetMapping: target ID (= the model number's rank) → the model's sites, centered by the centroid its particles
  // were placed by: the atoms end up where the operators put them
  const targetMapping = new Map<number, ParticleTarget>();
  if (options.resolveTargets !== false) {
    for (const [modelNum, idx] of modelNumToIndex) {
      const group = modelBounds.get(modelNum);
      if (group) targetMapping.set(idx, { kind: 'sites', sites: group.sites, center: group.center });
    }
  }

  const keys = new Int32Array(totalCount);
  const targets = new Int32Array(totalCount);
  const entities = new Int32Array(totalCount).fill(-1);
  const coordinates = new Float32Array(totalCount * 3);
  const rotations = new Float32Array(totalCount * 4);
  const radii = new Float32Array(totalCount);

  const labelEntry = new Int32Array(totalCount);
  const labelCombo = new Int32Array(totalCount);

  const entityNameToIdx = new Map<string, number>();

  const combined = Mat4();
  const quaternion = Quat();
  const centroid = Vec3();
  const position = Vec3();
  const updateChunk = 10000;

  let count = 0;
  for (let ei = 0, eil = entries.length; ei < eil; ei++) {
    const entry = entries[ei];
    const targetIdx = modelNumToIndex.get(entry.modelNum)!; // 0-based trajectory model index
    const modelName = modelNumToName.get(entry.modelNum) || `Model ${entry.modelNum}`;

    const bounds = modelBounds.get(entry.modelNum);
    Vec3.copy(centroid, bounds ? bounds.center : Vec3.origin);

    for (let ci = 0, cil = entry.combinations.count; ci < cil; ci++) {
      composeOperatorCombination(operators, entry.combinations, ci, combined);

      Vec3.transformMat4(position, centroid, combined);
      Quat.normalize(quaternion, Quat.fromMat4(quaternion, combined));

      const cOffset = count * 3;
      coordinates[cOffset + 0] = position[0];
      coordinates[cOffset + 1] = position[1];
      coordinates[cOffset + 2] = position[2];

      const qOffset = count * 4;
      rotations[qOffset + 0] = quaternion[0];
      rotations[qOffset + 1] = quaternion[1];
      rotations[qOffset + 2] = quaternion[2];
      rotations[qOffset + 3] = quaternion[3];

      radii[count] = (bounds ? bounds.radius : 0) * Mat4.getMaxScaleOnAxis(combined);

      keys[count] = count;
      targets[count] = targetIdx;

      if (!entityNameToIdx.has(modelName)) entityNameToIdx.set(modelName, entityNameToIdx.size);
      entities[count] = entityNameToIdx.get(modelName)!;

      labelEntry[count] = ei;
      labelCombo[count] = ci;
      count++;
      if (count % updateChunk === 0 && ctx.shouldUpdate) {
        await ctx.update({ message: 'Building particles', current: count, max: totalCount });
      }
    }
  }
  if (ctx.shouldUpdate) {
    await ctx.update({ message: 'Building particles', current: count, max: totalCount });
  }

  // Build entityInfo: entity index → model name.
  const entityInfo = new Map<number, ParticleEntityInfo>();
  for (const [name, idx] of entityNameToIdx) {
    entityInfo.set(idx, { name });
  }

  // Build targetInfo: target index → model name and its (single, by construction) entity index.
  const targetInfo = new Map<number, ParticleTargetInfo>();
  for (const [modelNum, idx] of modelNumToIndex) {
    const modelName = modelNumToName.get(modelNum) || `Model ${modelNum}`;
    targetInfo.set(idx, { name: modelName, entity: entityNameToIdx.get(modelName) });
  }

  const assemblyId = options.assemblyId;
  // Only the ids are captured below, so the label closure does not retain the matrix array.
  const operatorIds = operators.ids;
  const combinations = entries.map((e) => e.combinations);

  return {
    list: {
      label: buildMmcifLabel(options.label, assemblyId),
      count,
      keys,
      targets,
      targetInfo,
      entities: entityInfo.size > 0 ? entities : undefined,
      entityInfo: entityInfo.size > 0 ? entityInfo : undefined,
      coordinates,
      rotations,
      radii,
      getParticleLabel: (index: number) => {
        const entity = entityInfo.get(entities[index])!.name;
        const opCombo = formatOperatorCombination(operatorIds, combinations[labelEntry[index]], labelCombo[index]);
        return `Particle ${index + 1} · ${entity} · ops ${opCombo}`;
      },
      sourceData: MmcifParticleFormat.create(cifFile),
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    },
    targetMapping,
  };
}

function buildMmcifLabel(label?: string, assemblyId?: string): string {
  const base = label || 'Particles';
  return assemblyId ? `${base} (assembly ${assemblyId})` : base;
}

//

export { MmcifParticleFormat };

type MmcifParticleFormat = ModelFormat<CifFile>;

namespace MmcifParticleFormat {
  export function is(x?: ModelFormat): x is MmcifParticleFormat {
    return x?.kind === 'mmcif-particle';
  }

  export function create(cifFile: CifFile): MmcifParticleFormat {
    return { kind: 'mmcif-particle', name: 'mmcif-particle', data: cifFile };
  }
}
