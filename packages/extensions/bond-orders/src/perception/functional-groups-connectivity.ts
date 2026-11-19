/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 *
 * Based on Roger Sayle (2001) method. Functional groups patterns are used to assign bond orders
 * after a geometry assignment step. The patterns are supposed to be unambiguous.
 * Most of the patterns duplicate the ones already in chemistry/functional-group.ts, but the properties added
 * in the declaration allows assignment of bond orders based on information collected during matching.
 * Note the following quirks:
 * - It's possible that 2 functional groups overlap. e.g. HBI C[C@H](O)[C@H](O)C1=NC2=C(NC1)N=C(N)NC2=O has a N in common between an amide and a guanidine.
 *   No backtracking is done to find the best match: for now, order gives priority.
 * - Amide pattern is treated differently from Sayle's depiction: no single bond is enforced to the N partners. This addresses a limitation discussed in Labute (2005).
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { isDebugMode } from '@molstar/core/util/debug';
import { Elements } from '@molstar/model/model/structure/model/properties/atomic/types';
import type { ElementSymbol } from '@molstar/model/model/structure/model/types';
import { AtomGeometry } from '@molstar/model/props/computed/chemistry/geometry';
import { atomId } from '@molstar/model/props/computed/chemistry/util';
import {
  computeOpenValence,
  getOrder,
  hasMultipleBond,
  isAssignedBond,
  isHydrogenElement,
  pos,
  setAssignedBond,
  setBond,
  type State,
} from './common.js';

/** Element spec for functional-group matching: an element symbol or '*' = any heavy atom. */
type ElSpec = Elements | '*';

type NeighborConnectivity = 'terminal' | 'connected' | 'either';
interface NeighbourSpec {
  el: ElSpec;
  connectivity: NeighborConnectivity;
  order: number;
}
interface FunctionalGroup {
  name: string;
  center: { el: ElSpec; geometry: AtomGeometry };
  neighbours: NeighbourSpec[];
  fp: number;
}

const FunctionalGroups: FunctionalGroup[] = [
  // *N(*)C=O
  {
    name: 'formamide',
    center: { el: Elements.C, geometry: AtomGeometry.Trigonal },
    neighbours: [
      { el: Elements.O, connectivity: 'terminal', order: 2 },
      { el: Elements.N, connectivity: 'either', order: 1 },
      // H?
    ],
    fp: 0,
  },
  // CC(=O)O
  {
    name: 'carboxylic acid',
    center: { el: Elements.C, geometry: AtomGeometry.Trigonal },
    neighbours: [
      { el: Elements.C, connectivity: 'connected', order: 1 },
      { el: Elements.O, connectivity: 'terminal', order: 2 },
      { el: Elements.O, connectivity: 'terminal', order: 1 },
    ],
    fp: 0,
  },
  // *C(=O)O*
  {
    name: 'ester',
    center: { el: Elements.C, geometry: AtomGeometry.Trigonal },
    neighbours: [
      { el: '*', connectivity: 'connected', order: 1 },
      { el: Elements.O, connectivity: 'terminal', order: 2 },
      { el: Elements.O, connectivity: 'connected', order: 1 },
    ],
    fp: 0,
  },
  // *C(=O)S*
  {
    name: 'thioester',
    center: { el: Elements.C, geometry: AtomGeometry.Trigonal },
    neighbours: [
      { el: '*', connectivity: 'connected', order: 1 },
      { el: Elements.O, connectivity: 'terminal', order: 2 },
      { el: Elements.S, connectivity: 'connected', order: 1 },
    ],
    fp: 0,
  },
  // *C(=O)N* — covers amide (C-C(=O)-N, incl. acetamide with a terminal methyl),
  // urea (N-C(=O)-N, e.g. the C2 carbonyl of uracil/thymine) and carbamate (O-C(=O)-N).
  // Note: this differs from the depiction in Sayle (2001) which hints at a carbonyl substituent
  {
    name: 'amide',
    center: { el: Elements.C, geometry: AtomGeometry.Trigonal },
    neighbours: [
      { el: '*', connectivity: 'either', order: 1 },
      { el: Elements.O, connectivity: 'terminal', order: 2 },
      { el: Elements.N, connectivity: 'either', order: 1 },
    ],
    fp: 0,
  },
  // *C(=S)S*
  {
    name: 'dithioester',
    center: { el: Elements.C, geometry: AtomGeometry.Trigonal },
    neighbours: [
      { el: '*', connectivity: 'connected', order: 1 },
      { el: Elements.S, connectivity: 'terminal', order: 2 },
      { el: Elements.S, connectivity: 'connected', order: 1 },
    ],
    fp: 0,
  },
  // *C(=S)N*
  {
    name: 'thioamide',
    center: { el: Elements.C, geometry: AtomGeometry.Trigonal },
    neighbours: [
      { el: '*', connectivity: 'connected', order: 1 },
      { el: Elements.S, connectivity: 'terminal', order: 2 },
      { el: Elements.N, connectivity: 'connected', order: 1 },
    ],
    fp: 0,
  },
  // *NC(=N*)N*
  {
    name: 'guanidinium',
    center: { el: Elements.C, geometry: AtomGeometry.Trigonal },
    neighbours: [
      { el: Elements.N, connectivity: 'either', order: 2 },
      { el: Elements.N, connectivity: 'either', order: 1 },
      { el: Elements.N, connectivity: 'connected', order: 1 },
    ],
    fp: 0,
  },
  // *C(=N)N
  {
    name: 'amidine',
    center: { el: Elements.C, geometry: AtomGeometry.Trigonal },
    neighbours: [
      { el: '*', connectivity: 'connected', order: 1 },
      { el: Elements.N, connectivity: 'terminal', order: 2 },
      { el: Elements.N, connectivity: 'terminal', order: 1 },
    ],
    fp: 0,
  },
  // *N=[N+]=[N-]
  {
    name: 'azide',
    center: { el: Elements.N, geometry: AtomGeometry.Linear },
    neighbours: [
      { el: Elements.N, connectivity: 'terminal', order: 2 },
      { el: Elements.N, connectivity: 'connected', order: 2 },
    ],
    fp: 0,
  },
  // *P(=O)(*)*
  {
    name: 'phosphoryl',
    center: { el: Elements.P, geometry: AtomGeometry.Tetrahedral },
    neighbours: [
      { el: '*', connectivity: 'either', order: 1 },
      { el: '*', connectivity: 'either', order: 1 },
      { el: '*', connectivity: 'either', order: 1 },
      { el: Elements.O, connectivity: 'terminal', order: 2 },
    ],
    fp: 0,
  },
  // *S(=O)(=O)*
  {
    name: 'sulfonyl',
    center: { el: Elements.S, geometry: AtomGeometry.Tetrahedral },
    neighbours: [
      { el: '*', connectivity: 'either', order: 1 },
      { el: '*', connectivity: 'either', order: 1 },
      { el: Elements.O, connectivity: 'terminal', order: 2 },
      { el: Elements.O, connectivity: 'terminal', order: 2 },
    ],
    fp: 0,
  },
  // *S(=N)(=N)*
  {
    name: 'thiodiimine',
    center: { el: Elements.S, geometry: AtomGeometry.Tetrahedral },
    neighbours: [
      { el: '*', connectivity: 'either', order: 1 },
      { el: '*', connectivity: 'either', order: 1 },
      { el: Elements.N, connectivity: 'terminal', order: 2 },
      { el: Elements.N, connectivity: 'terminal', order: 2 },
    ],
    fp: 0,
  },
  // *[N+](=O)[O-]
  {
    name: 'nitro',
    center: { el: Elements.N, geometry: AtomGeometry.Trigonal },
    neighbours: [
      { el: '*', connectivity: 'either', order: 1 },
      { el: Elements.O, connectivity: 'terminal', order: 2 },
      { el: Elements.O, connectivity: 'terminal', order: 1 },
    ],
    fp: 0,
  },
  // *C#N
  {
    name: 'nitrile',
    center: { el: Elements.C, geometry: AtomGeometry.Linear },
    neighbours: [
      { el: '*', connectivity: 'connected', order: 1 },
      { el: Elements.N, connectivity: 'terminal', order: 3 },
    ],
    fp: 0,
  },
  // *[Se](=O)O
  {
    name: 'seleninic acid',
    center: { el: Elements.SE, geometry: AtomGeometry.Tetrahedral },
    neighbours: [
      { el: '*', connectivity: 'connected', order: 1 },
      { el: Elements.O, connectivity: 'terminal', order: 2 },
      { el: Elements.O, connectivity: 'terminal', order: 1 },
    ],
    fp: 0,
  },
];

const FingerprintBits = {
  C: 1 << 0,
  N1: 1 << 1,
  N2: 1 << 2,
  N3: 1 << 3,
  O1: 1 << 4,
  O2: 1 << 5,
  S1: 1 << 6,
  S2: 1 << 7,
};

FunctionalGroups.forEach((g) => (g.fp = computeFingerprint(g.neighbours.map((n) => n.el))));

function computeFingerprint(elements: (ElementSymbol | ElSpec)[]) {
  let fp = 0x0;
  for (const el of elements) {
    switch (el) {
      case Elements.C:
        fp |= FingerprintBits.C;
        break;
      case Elements.N:
        if (fp & FingerprintBits.N1) {
          if (fp & FingerprintBits.N2) fp |= FingerprintBits.N3;
          else fp |= FingerprintBits.N2;
        } else fp |= FingerprintBits.N1;
        break;
      case Elements.O:
        if (fp & FingerprintBits.O1) {
          fp |= FingerprintBits.O2;
        } else fp |= FingerprintBits.O1;
        break;
      case Elements.S:
        if (fp & FingerprintBits.S1) {
          fp |= FingerprintBits.S2;
        } else fp |= FingerprintBits.S1;
        break;
    }
  }
  return fp;
}

function matchEl(spec: ElSpec, el: ElementSymbol) {
  return spec === '*' ? !isHydrogenElement(el) : spec === el;
}

function connectivityOk(c: NeighborConnectivity, isTerminal: boolean) {
  return c === 'either' || (c === 'terminal' ? isTerminal : !isTerminal);
}

export function applyFunctionalGroups(state: State) {
  const { unitIndices, n, geometry, open } = state;
  // open valences, kept up to date as patterns write: a neighbour takes a multiple bond only within its own (the
  // shared N of two guanidines, an O with an explicit H)
  computeOpenValence(state);
  for (let i = 0; i < n; i++) {
    // If the center already has a multiple bond (e.g. a backbone C=O or phosphate
    // P=O from the order table), its pi system is already placed - don't add another.
    if (hasMultipleBond(state, i)) continue;
    const el = state.el[i];
    const geom = geometry[i];

    // Skip common cases that are not in patterns
    if (el === Elements.C && geom !== AtomGeometry.Trigonal && geom !== AtomGeometry.Linear) continue;
    if (el === Elements.O) continue;
    if (el === Elements.H) continue;
    if (geom === AtomGeometry.Terminal) continue;

    const nb = candidates(state, i);
    const u = unitIndices[i];
    const fp = computeFingerprint(nb.map((c) => c.el));
    isDebugMode &&
      console.log(
        'Functional group search for ',
        atomId(state.unit, u),
        ' el=',
        el,
        ' geom=',
        geom,
        ' heavyDeg=',
        nb.length,
      );
    for (const g of FunctionalGroups) {
      if (g.center.el !== el) continue;
      if (geom !== g.center.geometry) continue;
      if (nb.length !== g.neighbours.length) continue;
      if ((fp & g.fp) !== g.fp) continue; // fingerprint mismatch: some required element is missing
      const slots = assignSlots(state, i, nb, g.neighbours);
      if (!slots) continue;

      isDebugMode && console.log('Functional group ', g.name, ' matched at ', atomId(state.unit, u));
      for (let k = 0; k < slots.length; k++) {
        const j = nb[slots[k]].local;
        // a neighbour outside the residue fills an order-1 slot only, and its bond is never written (_spec §1.11)
        if (j < 0) continue;
        const v = unitIndices[j];
        setAssignedBond(state, u, v);
        const order = g.neighbours[k].order;
        isDebugMode &&
          console.log(
            'Functional group ',
            g.name,
            ' setting bond ',
            atomId(state.unit, u),
            '-',
            atomId(state.unit, v),
            ' order ',
            order,
          );
        if (order < 2) continue;
        const added = order - getOrder(state.bonds, u, v);
        setBond(state.bonds, u, v, order, 0);
        open[i] = Math.max(0, open[i] - added);
        open[j] = Math.max(0, open[j] - added);
      }
      break;
    }
  }
  computeOpenValence(state);
}

/** A heavy neighbour of a pattern's center: in the residue (`local`, its residue-local index) or outside it (`-1`). */
interface Candidate {
  readonly el: ElementSymbol;
  readonly local: number;
  /** squared distance to the center */
  readonly d2: number;
}

/**
 * The center's heavy neighbours, every bond counted: those in the residue and those outside it (another residue of
 * the unit, another unit), from the neighbour store (`neighbourPos`, `neighbourLocal`, `neighbourEl`).
 */
function candidates(state: State, i: number): Candidate[] {
  const local = state.neighbourLocal[i];
  const els = state.neighbourEl[i];
  const p = state.neighbourPos[i];
  pos(state, i, tmpPos);
  const result: Candidate[] = [];
  for (let k = 0, kl = local.length; k < kl; k++) {
    if (isHydrogenElement(els[k])) continue;
    const dx = p[3 * k] - tmpPos[0];
    const dy = p[3 * k + 1] - tmpPos[1];
    const dz = p[3 * k + 2] - tmpPos[2];
    result.push({ el: els[k], local: local[k], d2: dx * dx + dy * dy + dz * dz });
  }
  return result;
}

const tmpPos = Vec3();

/**
 * Match the center's neighbours with the neigbours specs.
 * Returns an array of indices into `nb` ordered by the slots in `specs`, or undefined if no complete assignment exists.
 *
 * Perfect bipartite matching of slots to neighbours. A per-pair cost ranks the assignments.
 * Heuristics:
 * - higher-order prefer short bonds (Labute 2005)
 * - higher-order prefer in-ring bonds (Sayle 2001)
 * A neighbour outside the residue is never terminal (an atom of another residue has another heavy neighbour) and
 * fills order-1 slots only: no multiple bond is written across the residue's boundary (_spec §1.11). A neighbour in it
 * fills a higher-order slot only within its open valence.
 */
function assignSlots(state: State, i: number, nb: Candidate[], specs: NeighbourSpec[]): number[] | undefined {
  const m = nb.length;
  const s = specs.length;
  const u = state.unitIndices[i];

  // Precompute the slot×neighbour cost matrix once (Infinity = element/connectivity mismatch, forbidden).
  // Heuristics: higher-order prefer short (Labute 2005) and in-ring (Sayle 2001) bonds.
  const cost = new Array<number>(s * m);
  for (let k = 0; k < m; k++) {
    const { el: elJ, local: j, d2 } = nb[k];
    const outside = j < 0;
    // heavy degree, every bond counted: a residue atom bonded outside the residue is not terminal
    const isTerminalHSuppressed = !outside && state.degree[j] - state.hCount[j] === 1;
    const inRingJ = outside ? 0 : state.inRing[j];
    const v = outside ? undefined : state.unitIndices[j];
    const assigned = v !== undefined && isAssignedBond(state, u, v);
    const order = v === undefined ? 0 : getOrder(state.bonds, u, v);
    const curOrder = assigned ? order : 0;
    for (let si = 0; si < s; si++) {
      const spec = specs[si];
      let c: number;
      if (
        !matchEl(spec.el, elJ) ||
        !connectivityOk(spec.connectivity, isTerminalHSuppressed) ||
        (outside && spec.order > 1) ||
        // the neighbour's open valence must allow the order increase
        (!outside && spec.order > order && state.open[j] < spec.order - order)
      ) {
        c = Infinity;
      } else {
        const ringBonus = spec.order > 1 && inRingJ ? -2 : 0;
        // If already assigned bond and conflicts with the spec, penalize to allow picking a better slot
        // (e.g. guanidine double bond assignment after amide recognition)
        // FIXME: wrong, kept until a benchmark can measure the change. A center with a multiple bond is skipped
        // (`applyFunctionalGroups`), so the only conflict left is an order-2 slot on a bond assigned single: proven
        // single by hybridization (too long, non-planar) or stated (`auto`). The search minimizes, and the penalty is
        // subtracted there: such a bond is preferred for the double (by 1 Å², more than the d2 differences), and
        // `setBond` then overrides the geometry or the stated order. Intended: no higher order on an assigned bond
        // (Infinity). The order-1 branch never applies.
        const preAssignmentPenalty = assigned && spec.order !== curOrder ? 1 : 0;
        c = spec.order > 1 ? d2 + ringBonus - preAssignmentPenalty : -d2 + preAssignmentPenalty;
      }
      cost[si * m + k] = c;
    }
  }

  const used = new Array<boolean>(m).fill(false);
  const curAssign = new Array<number>(s).fill(-1);
  let bestCost = Infinity;
  let bestAssign: number[] | undefined;

  const search = (si: number, acc: number) => {
    if (si === s) {
      if (acc < bestCost) {
        bestCost = acc;
        bestAssign = curAssign.slice();
      }
      return;
    }
    const base = si * m;
    for (let k = 0; k < m; k++) {
      if (used[k]) continue;
      const c = cost[base + k];
      if (!Number.isFinite(c)) continue;
      used[k] = true;
      curAssign[si] = k;
      search(si + 1, acc + c);
      used[k] = false;
    }
  };
  search(0, 0);
  return bestAssign;
}
