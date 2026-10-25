/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import {
  Structure,
  StructureElement,
  StructureProperties as P,
  StructureSelection,
} from '@molstar/model/model/structure';
import { Script } from '@molstar/model/script/script';
import { MolScriptBuilder as B } from '@molstar/query-language/language/builder';
import type { Expression } from '@molstar/query-language/language/expression';
import type { LabelPart, Sel } from '../language/ir.js';
import { isEmptyPart, selectionOn } from '../language/selection.js';

// Atoms for the host (_spec §5.2): a selection evaluated on one structure, an atom set written back as short MolQL,
// atom properties for labels and reads, and geometry. Selections are evaluated when a command runs, which is PyMOL's
// meaning (a selection is the atoms it selects then).

type Loci = StructureElement.Loci;

/** The atoms of a selection on a structure (its part for that object). */
export function lociOf(sel: Sel, object: string, structure: Structure): Loci {
  const molql = selectionOn(sel, object);
  if (isEmptyPart(molql)) return StructureElement.Loci.none(structure);
  const selection = Script.getStructureSelection(molql, structure);
  return StructureElement.Loci.remap(StructureSelection.toLociWithSourceUnits(selection), structure);
}

const PolymerQuery = B.struct.generator.atomGroups({
  'entity-test': B.core.rel.eq([B.ammp('entityType'), 'polymer']),
});

/** The atoms of `loci` in polymer entities: what a cartoon, a putty or a backbone draws (PyMOL draws nothing else). */
export function polymerAtoms(loci: Loci): Loci {
  if (StructureElement.Loci.isEmpty(loci)) return loci;
  const polymer = Script.getStructureSelection(PolymerQuery, loci.structure);
  return StructureElement.Loci.intersect(
    loci,
    StructureElement.Loci.remap(StructureSelection.toLociWithSourceUnits(polymer), loci.structure),
  );
}

/** The atoms of a component's structure on its parent structure. */
export function componentLoci(component: Structure | undefined, parent: Structure): Loci {
  if (!component) return StructureElement.Loci.none(parent);
  return StructureElement.Loci.remap(Structure.toStructureElementLoci(component), parent);
}

export const atomCount = (l: Loci) => StructureElement.Loci.size(l);

/**
 * An atom set as MolQL that selects the same atoms again: everything as `all`, else the atoms by their place in the
 * file (`Loci.toExpression`, ranges of source indices). Short whatever made the set (fifty show and hide lines).
 */
export function compactExpression(l: Loci): Expression {
  if (StructureElement.Loci.isEmpty(l)) return B.struct.generator.empty();
  if (StructureElement.Loci.isWholeStructure(l)) return B.struct.generator.all();
  return StructureElement.Loci.toExpression(l);
}

/** One atom's properties, for labels and `$sel get`. */
export function atomProperty(loc: StructureElement.Location, property: string): string {
  switch (property) {
    case 'resn':
    case 'resname':
      return P.atom.auth_comp_id(loc);
    case 'resi':
    case 'resid':
      return `${P.residue.auth_seq_id(loc)}${P.residue.pdbx_PDB_ins_code(loc) || ''}`;
    case 'chain':
    case 'segname':
      return P.chain.auth_asym_id(loc);
    case 'name':
      return P.atom.label_atom_id(loc);
    case 'elem':
    case 'element':
    case 'type':
      return P.atom.type_symbol(loc);
    case 'b':
    case 'beta':
      return String(P.atom.B_iso_or_equiv(loc));
    case 'occupancy':
      return String(P.atom.occupancy(loc));
    case 'index':
      return String(P.atom.sourceIndex(loc));
    case 'x':
      return P.atom.x(loc).toFixed(3);
    case 'y':
      return P.atom.y(loc).toFixed(3);
    case 'z':
      return P.atom.z(loc).toFixed(3);
    default:
      throw new Error(`The atom property '${property}' is not read here.`);
  }
}

export const ReadableProperties = [
  'resname',
  'resid',
  'chain',
  'name',
  'element',
  'type',
  'beta',
  'occupancy',
  'index',
  'x',
  'y',
  'z',
];

/** A label's text for one atom. */
export function labelText(loc: StructureElement.Location, parts: readonly LabelPart[]): string {
  return parts.map((p) => ('text' in p ? p.text : atomProperty(loc, p.property))).join('');
}

/** Each atom of a set, up to `max`: its location and its own loci. */
export function eachAtom(l: Loci, max: number): { location: StructureElement.Location; loci: Loci }[] {
  const out: { location: StructureElement.Location; loci: Loci }[] = [];
  for (const e of l.elements) {
    for (let i = 0, n = OrderedSet.size(e.indices); i < n && out.length < max; i++) {
      const index = OrderedSet.getAt(e.indices, i);
      const location = StructureElement.Location.create(l.structure, e.unit, e.unit.elements[index]);
      out.push({
        location,
        loci: StructureElement.Loci(l.structure, [{ unit: e.unit, indices: OrderedSet.ofSingleton(index) }]),
      });
    }
  }
  return out;
}

/** The position of the first atom of a set. */
export function positionOf(l: Loci): Vec3 | undefined {
  const loc = StructureElement.Loci.getFirstLocation(l);
  if (!loc) return undefined;
  return Vec3.create(P.atom.x(loc), P.atom.y(loc), P.atom.z(loc));
}

const Degrees = 180 / Math.PI;

/** A distance (Å), an angle or a dihedral (degrees) of points. */
export function geometryOf(kind: 'distance' | 'angle' | 'dihedral', p: readonly Vec3[]): number {
  if (kind === 'distance') return Vec3.distance(p[0], p[1]);
  if (kind === 'angle') {
    const a = Vec3.sub(Vec3(), p[0], p[1]);
    const b = Vec3.sub(Vec3(), p[2], p[1]);
    return Vec3.angle(a, b) * Degrees;
  }
  return Vec3.dihedralAngle(p[0], p[1], p[2], p[3]) * Degrees;
}

/** The centroid, the radius of gyration and the bounding box of a set. */
export function shapeOf(l: Loci): { center: Vec3; rgyr: number; min: Vec3; max: Vec3 } {
  const center = Vec3();
  const min = Vec3.create(Infinity, Infinity, Infinity);
  const max = Vec3.create(-Infinity, -Infinity, -Infinity);
  const ps: Vec3[] = [];
  StructureElement.Loci.forEachLocation(l, (loc) => {
    const v = Vec3.create(P.atom.x(loc), P.atom.y(loc), P.atom.z(loc));
    ps.push(v);
    Vec3.add(center, center, v);
    Vec3.min(min, min, v);
    Vec3.max(max, max, v);
  });
  if (ps.length) Vec3.scale(center, center, 1 / ps.length);
  let sq = 0;
  for (const v of ps) sq += Vec3.squaredDistance(v, center);
  return { center, rgyr: ps.length ? Math.sqrt(sq / ps.length) : 0, min, max };
}
