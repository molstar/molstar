/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Coordinate systems by name over a set (seqspec:libraries#the-system-registry): a `Translator` holding a host's systems
// and the maps registered with them, which also names a `SeqSet`'s own systems by key, `seqstar:<ref>` (a molecule, a
// unit, a complex, an alignment's columns, a row `<alignment>/<k>`). No viewer: a host gives it the set (`update`) and
// hears of every registration through a plain callback (`onChange`); the viewer's `ctx.systems` is one with a reactive
// version. Own descriptors are built on first use, once per system and lineage, and cached by key. Edge rate but the
// reads (`system`, `descriptor`, `Coordinates.toExtents` over it).

import {
  type CoordinateMap,
  type CoordinateSystem,
  CoordinateSystemKind as K,
  Coords,
  Extents,
  KeyKind,
  type Loci,
  Mappings,
  type PathStep,
  Paths,
  type SeqSet,
  Systems,
} from '@molstar/seq-core/model';
import { Coordinates, SeqSystemPrefix as PREFIX } from './coordinates.js';
import { type SystemDescriptor, Translator } from './translator.js';

/**
 * An own system's descriptor, valid while `owner` (the table or entity it reads) and `molecule` (a unit's molecule,
 * whose numberings it lists; replaced by an upsert that leaves the units' table shared) are the set's.
 */
interface Own {
  readonly set: SeqSet;
  readonly owner: unknown;
  readonly molecule: unknown;
  readonly system: CoordinateSystem;
  readonly descriptor: SystemDescriptor;
}
const refOf = Systems.refOf;
const domain = new Int32Array(2);
const addRun = (s: number, t: number, k: { out: Extents; to: CoordinateSystem }) => Extents.add(k.out, k.to, s, t);

export interface SystemRegistryOptions {
  /** Called after every registration, unregistration and added map (not on `update`). */
  readonly onChange?: () => void;
}

export class SystemRegistry extends Translator {
  /** Called after every registration, unregistration and added map: what holds `Coordinates` resolves them again. */
  onChange: (() => void) | null;
  /** Own descriptors built so far (a counter for tests). */
  built = 0;
  #set: SeqSet | null = null;
  readonly #own = new Map<string, Own>();
  readonly #key = { kind: 0, id: 0, row: -1 };
  readonly #one = Extents.create(1);

  constructor(options: SystemRegistryOptions = {}) {
    super();
    this.onChange = options.onChange ?? null;
  }

  /** The set its own keys name; null: none yet. */
  get set(): SeqSet | null {
    return this.#set;
  }

  /** The set (a viewer's shown one): own descriptors are kept within its lineage while what they read is kept. */
  update(set: SeqSet | null): void {
    if (!set || set === this.#set) return;
    if (this.#set?.lineage !== set.lineage) this.#own.clear();
    this.#set = set;
  }

  override register(
    d: SystemDescriptor,
    maps: readonly CoordinateMap[] | ((s: CoordinateSystem) => readonly CoordinateMap[]) = [],
    system?: CoordinateSystem,
  ): CoordinateSystem {
    const s = super.register(d, maps, system);
    this.onChange?.();
    return s;
  }
  override unregister(key: string): void {
    super.unregister(key);
    this.onChange?.();
  }
  override addMap(map: CoordinateMap): void {
    super.addMap(map);
    this.onChange?.();
  }

  /** A host's key, else the set's own (`seqstar:<ref>`). */
  override system(key: string): CoordinateSystem | null {
    return super.system(key) ?? this.ownOf(key)?.system ?? null;
  }
  override descriptor(key: string): SystemDescriptor | null {
    return super.descriptor(key) ?? this.ownOf(key)?.descriptor ?? null;
  }
  override descriptorOf(system: CoordinateSystem): SystemDescriptor | null {
    const d = super.descriptorOf(system);
    const ref = d || !this.#set ? null : refOf(this.#set, system.kind, system.id, system.row);
    return d ?? (ref === null ? null : (this.ownOf(PREFIX + ref)?.descriptor ?? null));
  }

  /** The key naming a system: a host's, else the set's own; null: none. */
  keyOf(system: CoordinateSystem): string | null {
    const d = super.descriptorOf(system);
    if (d) return d.key;
    const ref = this.#set ? refOf(this.#set, system.kind, system.id, system.row) : null;
    return ref === null ? null : PREFIX + ref;
  }

  /** Loci of the set as `Coordinates` (edge rate, allocates; `Coordinates.ofLoci`). */
  coordinatesOf(loci: Loci): Coordinates {
    return Coordinates.ofLoci(loci);
  }

  /** `Coordinates` appended to `out` as entries in the set's own systems (`Coordinates.toExtents`); the count added. */
  extentsOf(c: Coordinates, out: Extents): number {
    return Coordinates.toExtents(this, c, out, this.#one);
  }

  /** `Coordinates` as loci of the set, one per kind (edge rate, allocates); unknown keys give nothing. */
  lociOf(c: Coordinates): Loci[] {
    const set = this.#set;
    const ext = Extents.create(8);
    if (!set || this.extentsOf(c, ext) === 0) return [];
    return Coordinates.toLoci(set, ext);
  }

  /**
   * A target mapped into the system `to` names (seqspec:libraries#the-system-registry), with no viewer: a host's
   * system's runs first onto the set's systems through the maps registered with it, then through the set's own maps
   * (the model's shortest paths: mappings, row bindings, an alignment's columns, units onto their molecules), then, for
   * a host's `to`, through the maps registered with it (a composite `to`: each of its members). Edge rate, allocates.
   * Nothing: `to` unknown, or no route.
   */
  map(target: Coordinates | Loci, to: string): Coordinates {
    const set = this.#set;
    const dest = this.system(to);
    const none: Coordinates = { systems: [], runs: new Int32Array(0) };
    if (!set || !dest) return none;
    const from = Extents.create(8);
    if (Coordinates.is(target)) this.extentsOf(target, from);
    else if (target.set === set) Extents.addLoci(from, target);
    const out = Extents.create(8);
    if (dest.kind !== K.external) this.within(set, from, dest, out);
    else {
      const members = this.descriptorOf(dest)?.members;
      for (const key of members ?? [to]) {
        const m = this.system(key);
        if (m) this.onto(set, from, m, out);
      }
    }
    return Coordinates.ofExtents(set, out, (s) => this.keyOf(s));
  }

  /** Entries in the set's systems onto a host's system `to`, through each Seq* system a map registered with it reaches. */
  private onto(set: SeqSet, from: Extents, to: CoordinateSystem, out: Extents) {
    const mid = Extents.create(8);
    for (const n of this.neighbors(to)) {
      if (n.kind === K.external || n.kind === K.complex) continue;
      mid.count = 0;
      this.within(set, from, n, mid);
      if (mid.count) this.translate(mid, n, to, out);
    }
  }

  /** Entries onto one of the set's systems through the model's paths (the same system: as they are). */
  private within(set: SeqSet, from: Extents, to: CoordinateSystem, out: Extents) {
    const steps: PathStep[] = [];
    const sink = { out, to };
    for (let i = 0; i < from.count; i++) {
      const s = { kind: from.kind[i], id: from.id[i], row: from.row[i] };
      const [start, stop] = [from.start[i], from.stop[i]];
      if (Coords.sameSystem(s, to)) {
        Extents.add(out, to, start, stop);
        continue;
      }
      // -2: ambiguous, `steps` holds the first candidate (the model's rule: a status, never a drop).
      if (Mappings.path(set, s, to, steps) === -1) continue;
      const p = Mappings.compose(set, s, steps);
      if (Coords.sameSystem(p.from, s)) Paths.forEachTo(p, start, stop, addRun, sink);
      else Paths.forEachFrom(p, start, stop, addRun, sink);
    }
  }

  /** An own key's system and descriptor, cached; null: not `seqstar:`, or naming nothing in the set. */
  private ownOf(key: string): Own | null {
    const set = this.#set;
    const hit = this.#own.get(key);
    if (!set || (hit && hit.set === set)) return hit ?? null;
    const t = this.#key;
    if (!key.startsWith(PREFIX) || !set.byName(key.slice(PREFIX.length), t)) return null;
    // By `KeyKind` (molecule, unit, alignment, row, mapping, annotation, complex): the coordinate system it names.
    const kind = [K.molecule, K.unit, K.columns, K.row, -1, -1, K.complex][t.kind];
    const id = t.id;
    const owner =
      t.kind === KeyKind.unit
        ? set.units
        : t.kind === KeyKind.molecule
          ? set.molecules[id]
          : t.kind === KeyKind.complex
            ? set.complexes[id]
            : set.alignments[id];
    if (kind < 0 || !owner) return null;
    const molecule = t.kind === KeyKind.unit ? set.molecules[set.units.molecule[id]] : null;
    if (hit?.owner === owner && hit.molecule === molecule) {
      this.#own.set(key, { ...hit, set });
      return this.#own.get(key) ?? null;
    }
    const own = {
      set,
      owner,
      molecule,
      system: { kind, id, row: t.row },
      descriptor: describe(set, key, kind, id, t.row),
    };
    this.built++;
    this.#own.set(key, own);
    return own;
  }
}

function describe(set: SeqSet, key: string, kind: number, id: number, row: number): SystemDescriptor {
  const u = set.units;
  const str = (k: number) => set.strings.get(k);
  Systems.domainOf(set, kind, id, row, domain);
  const [start, stop] = [domain[0], Math.max(domain[0], domain[1])];
  const base = { key, start, stop, refs: {}, members: null, structure: null, numberings: [], owner: 'seqstar' };
  if (kind === K.unit) {
    const m = set.molecules[u.molecule[id]];
    const cx = set.complexes[u.complex[id]];
    const mvs = (u.custom[id]?.doc as { mvs?: SystemDescriptor['structure'] } | undefined)?.mvs ?? null;
    return {
      ...base,
      kind: 'unit',
      refs: { unit: str(u.key[id]), molecule: m ? str(m.key) : undefined, complex: cx ? str(cx.key) : undefined },
      structure: mvs,
      numberings: Systems.numberingNames(set, kind, id, row),
    };
  }
  if (kind === K.molecule) return { ...base, kind: 'molecule', refs: { molecule: key.slice(PREFIX.length) } };
  if (kind === K.complex) {
    const cx = set.complexes[id];
    const members = Array.from(cx?.units ?? [], (x) => PREFIX + str(u.key[x]));
    return { ...base, kind: 'complex', refs: { complex: key.slice(PREFIX.length) }, members };
  }
  const a = set.alignments[id];
  const alignment = a ? str(a.key) : undefined;
  if (kind === K.columns) return { ...base, kind: 'columns', refs: { alignment } };
  return { ...base, kind: 'row', refs: { alignment, row } };
}
