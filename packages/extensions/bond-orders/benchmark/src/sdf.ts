/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import type { Structure, StructureElement } from '@molstar/model/model/structure';
import { Unit } from '@molstar/model/model/structure';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { BondType } from '@molstar/model/model/structure/model/types';
import { calculateHydrogensCharge, explicitValence } from '@molstar/model/props/computed/chemistry/valence-model';
import { isHydrogen } from './input.ts';

/**
 * A V2000 mol block of every atom of `structure` with its covalent bonds (plan §4.1). Charges from Mol*'s valence model
 * with `assignCharge: 'always'` (`atomCharge`): the input's charges never reach the perceived side. `isCap`: tier C's
 * cap atoms, written neutral.
 */
export function toMolBlock(structure: Structure, name: string, isCap?: (element: ElementIndex) => boolean): string {
  const atomLines: string[] = [];
  const bondLines: string[] = [];
  const charges: [number, number][] = [];
  const lineOf = new Map<number, Int32Array>();

  let hasExplicitH = false;
  for (const unit of structure.units) {
    const { type_symbol } = unit.model.atomicHierarchy.atoms;
    for (let i = 0; i < unit.elements.length; i++) {
      if (isHydrogen(type_symbol.value(unit.elements[i]))) hasExplicitH = true;
    }
  }

  const props = { assignCharge: 'always', assignH: 'auto' } as const;
  for (const unit of structure.units) {
    if (!Unit.isAtomic(unit)) continue;
    const { type_symbol } = unit.model.atomicHierarchy.atoms;
    const c = unit.conformation;
    const lines = new Int32Array(unit.elements.length);
    for (let i = 0; i < unit.elements.length; i++) {
      const e = unit.elements[i];
      const index = i as StructureElement.UnitIndex;
      // a cap is neutral, its H filled by RDKit (plan §4.2)
      const charge = isCap?.(e) ? 0 : atomCharge(structure, unit, index, props, hasExplicitH);
      atomLines.push(
        `${f(c.x(e))}${f(c.y(e))}${f(c.z(e))} ${symbol(type_symbol.value(e)).padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`,
      );
      lines[i] = atomLines.length;
      if (charge !== 0) charges.push([atomLines.length, charge]);
    }
    lineOf.set(unit.id, lines);

    const { a, b, edgeProps } = unit.bonds;
    for (let t = 0; t < a.length; t++) {
      if (a[t] < b[t] && BondType.isCovalent(edgeProps.flags[t])) {
        bondLines.push(bond(lines[a[t]], lines[b[t]], edgeProps.order[t]));
      }
    }
  }
  for (const e of structure.interUnitBonds.edges) {
    if (e.unitA > e.unitB || !BondType.isCovalent(e.props.flag)) continue;
    const la = lineOf.get(e.unitA)?.[e.indexA];
    const lb = lineOf.get(e.unitB)?.[e.indexB];
    if (la && lb) bondLines.push(bond(la, lb, e.props.order));
  }

  const out = [name.slice(0, 80), `  ${'Molstar'.padEnd(8)}${''.padEnd(10)}3D`, ''];
  out.push(`${p3(atomLines.length)}${p3(bondLines.length)}  0  0  0  0  0  0  0  0999 V2000`);
  out.push(...atomLines, ...bondLines);
  for (let i = 0; i < charges.length; i += 8) {
    const chunk = charges.slice(i, i + 8);
    out.push(`M  CHG${p3(chunk.length)}${chunk.map(([l, c]) => `${p4(l)}${p4(c)}`).join('')}`);
  }
  out.push('M  END');
  return `${out.join('\n')}\n`;
}

/**
 * Bond-order sum past which a neutral atom is charged, and the charge that completes its octet (as in the extension's
 * valence table): a tetrahedral boronate B-, an ammonium or iminium N+, an oxonium O+ (pyrylium).
 */
const OctetCharge: Record<string, [number, number]> = { B: [4, -1], N: [4, 1], O: [3, 1] };

/**
 * Mol*'s valence model, except for the atoms it leaves neutral past their valence: boron, which it never charges, and
 * N and O without H (the model then reads the missing H as implicit).
 */
function atomCharge(
  structure: Structure,
  unit: Unit.Atomic,
  index: StructureElement.UnitIndex,
  props: Parameters<typeof calculateHydrogensCharge>[3],
  hasExplicitH: boolean,
) {
  const element = unit.model.atomicHierarchy.atoms.type_symbol.value(unit.elements[index]);
  const octet = OctetCharge[element];
  if (octet && explicitValence(structure, unit, index) === octet[0]) return octet[1];
  if (element === 'B') return 0;
  return calculateHydrogensCharge(structure, unit, index, props, hasExplicitH)[0];
}

const f = (v: number) => v.toFixed(4).padStart(10);
const p3 = (v: number) => String(v).padStart(3);
const p4 = (v: number) => String(v).padStart(4);
const bond = (a: number, b: number, order: number) => `${p3(a)}${p3(b)}${p3(order)}  0`;
/** `CL` → `Cl` */
const symbol = (s: string) => s.charAt(0).toUpperCase() + s.slice(1).toLowerCase();
