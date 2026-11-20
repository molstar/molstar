/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the structure translation reads of a model (_spec/README.md "Structures"): the static model and, when the
// complexes are a source's, its overlay, read as one: units, complexes, molecules and atoms systems by address, paths
// with `via`, a unit's numberings. Kept per model and overlay, so a hover a frame finds its indices made.

import type { QuerySystem } from '@molstar/seqx-core/query';
import {
  type AtomsEntry,
  type ComplexEntry,
  type MoleculeEntry,
  type Overlay,
  overlayQueryModel,
  type PathResult,
  type StaticModel,
  type SystemEntry,
  type UnitEntry,
} from '@molstar/seqx-core/resolve';

/** A model the translation reads: the static model, or it with the overlay a source's complexes are in. */
export type StructureModel = StaticModel | { readonly model: StaticModel; readonly overlay: Overlay };

/** One complex a translation is bound to, and the maps its paths take (`via`: their rivals left out). */
export interface StructureBinding {
  /** The complex's address (`/data/1fin`, `/data/1fin/(node :key "complex")`). */
  readonly complex: string;
  readonly via?: readonly string[];
}

/** @internal The reads, one object a model and overlay. */
export interface StructureReads {
  readonly model: StaticModel;
  unit(address: string): UnitEntry | undefined;
  complex(address: string): ComplexEntry | undefined;
  /** Every complex, the model's then the overlay's, each in its address order. */
  complexes(): readonly ComplexEntry[];
  molecule(address: string): MoleculeEntry | undefined;
  atoms(address: string): AtomsEntry | undefined;
  system(address: string): SystemEntry | undefined;
  /** The query forms' system (its numberings: a unit's own, else its molecule's). */
  querySystem(address: string): QuerySystem | undefined;
  path(from: string, to: string, via?: readonly string[]): PathResult;
  /** The units of a complex that instantiate a molecule (EC-48), in the complex's order. */
  unitsOf(complex: ComplexEntry, molecule: string): readonly UnitEntry[];
  /** The molecules a complex's units instantiate, each once. */
  moleculesOf(complex: ComplexEntry): Iterable<string>;
  /** The units of a complex with a map of their own (a step that is no unit step: the spec's or the overlay's). */
  ownMapped(complex: ComplexEntry): ReadonlySet<string>;
}

const Kept = new WeakMap<StaticModel, WeakMap<Overlay | StaticModel, StructureReads>>();

/** @internal The reads of a model (kept: one a model and overlay). */
export function structureReads(m: StructureModel): StructureReads {
  const model = 'overlay' in m ? m.model : m;
  const overlay = 'overlay' in m ? m.overlay : undefined;
  let byOverlay = Kept.get(model);
  if (!byOverlay) Kept.set(model, (byOverlay = new WeakMap()));
  const key = overlay ?? model;
  let r = byOverlay.get(key);
  if (!r) byOverlay.set(key, (r = makeReads(model, overlay)));
  return r;
}

function makeReads(model: StaticModel, overlay: Overlay | undefined): StructureReads {
  const own = model.molecular;
  const dyn = overlay?.molecular;
  const query = overlay ? overlayQueryModel(model, overlay) : model;
  let all: ComplexEntry[] | undefined;
  const byMolecule = new WeakMap<ComplexEntry, Map<string, UnitEntry[]>>();
  const mapped = new WeakMap<ComplexEntry, Set<string>>();
  const indexOf = (complex: ComplexEntry) => {
    let index = byMolecule.get(complex);
    if (!index) {
      byMolecule.set(complex, (index = new Map()));
      for (const u of complex.units) {
        let list = index.get(u.of);
        if (!list) index.set(u.of, (list = []));
        list.push(u);
      }
    }
    return index;
  };
  const edgesFrom = (s: string) => overlay?.edgesFrom?.(s) ?? model.edgesFrom(s);
  return {
    model,
    unit: (a) => own.units.get(a) ?? dyn?.units.get(a),
    complex: (a) => own.complexes.get(a) ?? dyn?.complexes.get(a),
    complexes: () => {
      if (all) return all;
      const byAddress = (m: ReadonlyMap<string, ComplexEntry>) =>
        [...m.values()].sort((a, b) => (a.address < b.address ? -1 : a.address > b.address ? 1 : 0));
      all = [...byAddress(own.complexes), ...(dyn ? byAddress(dyn.complexes) : [])];
      return all;
    },
    molecule: (a) => own.molecules.get(a) ?? dyn?.molecules.get(a),
    atoms: (a) => own.atoms.get(a) ?? dyn?.atoms.get(a),
    system: (a) => model.systems.get(a) ?? overlay?.system(a),
    querySystem: (a) => query.system(a),
    path: (from, to, via) => (overlay ? overlay.path(from, to, via) : model.path(from, to, via)),
    unitsOf: (complex, molecule) => indexOf(complex).get(molecule) ?? [],
    moleculesOf: (complex) => indexOf(complex).keys(),
    ownMapped(complex) {
      let set = mapped.get(complex);
      if (!set) {
        mapped.set(complex, (set = new Set()));
        for (const u of complex.units) if (edgesFrom(u.address).some((e) => !e.unit)) set.add(u.address);
      }
      return set;
    },
  };
}

/** @internal The bindings a call reads: those given, else every complex with no `via`. */
export function bindingsOf(
  r: StructureReads,
  complexes: readonly StructureBinding[] | undefined,
): readonly StructureBinding[] {
  return complexes ?? r.complexes().map((c): StructureBinding => ({ complex: c.address }));
}
