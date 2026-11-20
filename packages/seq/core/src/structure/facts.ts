/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A complex's units as plain facts (_spec/README.md "Structures", Facts): what a host's adapter joins against its own
// residue list once, so per message it looks up and never walks the model. Its `auth` numbering is the unit's own,
// else its molecule's of that id shifted to the unit (EC-8), labels as written (EC-1 to EC-7: read, never re-derived).

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { findNumbering, type Numbering } from '@molstar/seqx-core/query';
import type { UnitEntry } from '@molstar/seqx-core/resolve';
import { type StructureModel, type StructureReads, structureReads } from './model.js';

/** A unit's molecule's kind: a sequence's unit `polymer`, a `branched` molecule's, else `other` (a ligand, an ion). */
export type UnitKind = 'polymer' | 'branched' | 'other';

/** A unit's author numbering, both ways. */
export interface AuthNumbering {
  /** The label at a unit position (`94`, `100A`), undefined where none (an unmodeled residue, past the unit). */
  label(position: number): string | undefined;
  /** The unit position of a label, or -1 (the first one with it, EC-4). */
  position(label: string): number;
}

/** A unit's residues as columns, a row a unit position (from 0): the table an adapter joins once. */
export interface ResidueRows {
  readonly length: number;
  /** `auth_seq_id`: the label's number; 0 where `labeled` is 0. */
  readonly authSeqId: Int32Array;
  /** `pdbx_PDB_ins_code`: what follows the number (`''`: none). */
  readonly insCode: readonly string[];
  /** 1 where the residue has an author label. */
  readonly labeled: Uint8Array;
  /** 1 where the residue was modeled (atoms to light). */
  readonly observed: Uint8Array;
}

export interface UnitFacts {
  /** The unit's system address: what a Location on it names. */
  readonly address: string;
  /** Its id in its complex (a chain and an operator: `A-1_555`). */
  readonly id: string;
  /** `label_asym_id`. */
  readonly chain?: string;
  readonly authChain?: string;
  /** The assembly operator (`instance_id`). */
  readonly operator?: string;
  readonly label?: string;
  /** The sequence or molecule it instantiates, and where it starts on it. */
  readonly molecule: string;
  readonly start: number;
  /** Undefined while its molecule's length is unknown. */
  readonly length?: number;
  /** `polymer` for a sequence's unit, else its molecule's class (`branched`, `ligand`, `ion`, any word). */
  readonly class: string;
  readonly kind: UnitKind;
  /** Its modeled runs on its own positions, `[start, end, …]` (every position when all were). */
  readonly observed: Int32Array;
  /** Its copy group (EC-65): units of one complex with equal molecule and range share it. */
  readonly copy: string;
  readonly auth?: AuthNumbering;
  /** Its residues, when its length is known. */
  readonly residues?: ResidueRows;
}

export interface StructureFacts {
  readonly complex: string;
  readonly name?: string;
  /** In the complex's order. */
  readonly units: readonly UnitFacts[];
}

const Facts = new WeakMap<UnitEntry, UnitFacts>();

/** The facts of a complex's units, or undefined when the address is no complex. */
export function structureFacts(model: StructureModel, complex: string): StructureFacts | undefined {
  const r = structureReads(model);
  const c = r.complex(complex);
  if (!c) return undefined;
  return {
    complex: c.address,
    ...(c.name !== undefined ? { name: c.name } : {}),
    units: c.units.map((u) => unitFacts(r, u)),
  };
}

/** @internal A unit's facts (kept a unit entry). */
export function unitFacts(r: StructureReads, u: UnitEntry): UnitFacts {
  let f = Facts.get(u);
  if (f) return f;
  const length = u.end === undefined ? undefined : u.end - u.start;
  const auth = authOf(r, u.address, length);
  const observed = observedOf(u);
  f = {
    address: u.address,
    id: u.id,
    ...(u.chain !== undefined ? { chain: u.chain } : {}),
    ...(u.authChain !== undefined ? { authChain: u.authChain } : {}),
    ...(u.operator !== undefined ? { operator: u.operator } : {}),
    ...(u.label !== undefined ? { label: u.label } : {}),
    molecule: u.of,
    start: u.start,
    ...(length !== undefined ? { length } : {}),
    class: u.class,
    kind: u.class === 'polymer' ? 'polymer' : u.class === 'branched' ? 'branched' : 'other',
    observed,
    copy: u.copy,
    ...(auth ? { auth } : {}),
    ...(length !== undefined ? { residues: residueRows(length, observed, auth) } : {}),
  };
  Facts.set(u, f);
  return f;
}

function authOf(r: StructureReads, unit: string, length: number | undefined): AuthNumbering | undefined {
  const system = r.querySystem(unit);
  const of = system && findNumbering(system, 'auth');
  if (!of) return undefined;
  const n: Numbering = of.numbering;
  const shift = of.shift;
  const max = length ?? system.length;
  return {
    label: (p) => (Number.isInteger(p) && p >= 0 && p < max ? n.labelAt(p - shift) : undefined),
    position: (label) => {
      const i = n.indexOf(label);
      if (i < 0) return -1;
      const p = i + shift;
      return p >= 0 && p < max ? p : -1;
    },
  };
}

const Observed = new WeakMap<UnitEntry, Int32Array>();

/**
 * @internal A unit's observed runs, sorted, merged and within it (every position when it names none): made once a
 * unit (a hover a frame reads it for every unit it reaches; the rest of the facts are made when asked).
 */
export function observedOf(u: UnitEntry): Int32Array {
  let o = Observed.get(u);
  if (o) return o;
  const length = u.end === undefined ? undefined : u.end - u.start;
  if (length === undefined) o = new Int32Array(0);
  else if (!u.observed) o = Int32Array.of(0, length);
  else o = HalfOpenRanges.intersect(HalfOpenRanges.normalize(u.observed), Int32Array.of(0, length));
  Observed.set(u, o);
  return o;
}

const LabelPattern = /^(-?\d+)(.*)$/;

function residueRows(length: number, observed: Int32Array, auth: AuthNumbering | undefined): ResidueRows {
  const authSeqId = new Int32Array(length);
  const insCode: string[] = new Array(length).fill('');
  const labeled = new Uint8Array(length);
  const obs = new Uint8Array(length);
  for (let i = 0; i + 1 < observed.length; i += 2)
    obs.fill(1, Math.max(0, observed[i]), Math.min(length, observed[i + 1]));
  if (auth)
    for (let p = 0; p < length; p++) {
      const m = LabelPattern.exec(auth.label(p) ?? '');
      if (!m) continue;
      authSeqId[p] = Number(m[1]);
      insCode[p] = m[2];
      labeled[p] = 1;
    }
  return { length, authSeqId, insCode, labeled, observed: obs };
}
