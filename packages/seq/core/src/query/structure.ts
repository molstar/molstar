/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the structure forms walk (_spec/README.md "Queries", Language and Properties): a prefix's frames, a frame one
// unit of a complex (or a lone molecule or sequence) with its residues system and, when it has them, its atoms; and how
// each property is read for an element of a frame. A prefix is nodes (a complex, a unit, a molecule, a source holding
// complexes), a Location (positions on residues or atoms systems) or Keys (read as their extents): every form is
// monadic over it, each frame a part, one with nothing giving nothing (locations-and-queries.md, The language).

import { type Atom, type Named, parseAddress, printAddress } from '@molstar/core/repr/tree';
import { type Ctx, fail, prefixNodes, stage, tick } from './form-kit.js';
import { keysExtent } from './key-extent.js';
import { parseLabel } from './labels.js';
import { isEmptyLocation, type Location, rangesLocation } from './location.js';
import type { PropertyRead, QueryAtoms, QueryComponents, QuerySystem, QueryUnit } from './model.js';
import { findNumbering, numberingPending } from './numbering.js';
import type { PropertyDef } from './properties.js';
import { fromNamed } from './result.js';

/** One unit of a complex, or a lone molecule or sequence: its residues, its atoms, what the prefix holds of them. */
export interface Frame {
  /** Its residues system: the unit's, else the molecule's or the sequence's own. */
  readonly system: QuerySystem;
  readonly unit?: QueryUnit;
  /** The complex the unit is of. */
  readonly complex?: string;
  /** The molecule or sequence it is (a unit's `of`). */
  readonly of: string;
  /** The residues the prefix holds (sorted half-open ranges); undefined: all. */
  readonly residues?: Int32Array;
  /** The atoms the prefix holds, on its atoms system; undefined: all its residues' atoms. */
  readonly atoms?: Int32Array;
}

/** A value a property gives an element: an atom, a set of words; undefined where the element lacks it. */
export type PropValue = Atom | readonly string[] | undefined;

/** A property read on a frame: a value per element, by its residue and (at the atom level) its atom. */
export type Reader = (residue: number, atom: number) => PropValue;

/** The parent address of a unit: its complex. */
function parentOf(address: string): string | undefined {
  const s = parseAddress(address).pattern?.segments;
  return s && s.length > 1 ? printAddress(s.slice(0, -1)) : undefined;
}

/** The unit at an address, with its complex, or undefined. */
function unitAt(ctx: Ctx, address: string): { unit: QueryUnit; complex: string } | undefined {
  const complex = parentOf(address);
  const unit = complex ? ctx.model.units?.(complex)?.find((u) => u.address === address) : undefined;
  return unit && complex ? { unit, complex } : undefined;
}

/** The frames of a complex's units, each on its system (a unit the model has no system for is left out). */
function complexFrames(ctx: Ctx, complex: string, units: readonly QueryUnit[], out: Frame[]) {
  for (const unit of units) {
    tick(ctx);
    const system = ctx.model.system(unit.address);
    if (system) out.push({ system, unit, complex, of: unit.of });
  }
}

/** The frame of a residues system's address: a unit, a molecule or a sequence; undefined when it is none. */
function systemFrame(ctx: Ctx, address: string): Frame | undefined {
  const system = ctx.model.system(address);
  if (!system) return undefined;
  const u = unitAt(ctx, address);
  if (u) return { system, unit: u.unit, complex: u.complex, of: u.unit.of };
  return ctx.model.components?.(address) || ctx.model.atoms?.(address) ? { system, of: address } : undefined;
}

/** The frames of a Location: each part on a residues system (its positions) or an atoms system (its atoms). */
function locationFrames(ctx: Ctx, l: Location, kind: string): Frame[] {
  if (l.rects.length) fail('seq.query.form.prefix', `A form of kind ${kind} applies to positions, not a rect`);
  const out: Frame[] = [];
  for (const p of l.parts) {
    const atoms = ctx.model.atoms?.(p.system);
    if (atoms && atoms.address === p.system) {
      const f = systemFrame(ctx, atoms.of);
      if (f) out.push({ ...f, residues: atomsResidues(atoms, p.ranges), atoms: p.ranges });
      continue;
    }
    const f = systemFrame(ctx, p.system);
    if (f) out.push({ ...f, residues: p.ranges });
  }
  return out;
}

/** The residues holding atoms of the ranges. */
function atomsResidues(atoms: QueryAtoms, ranges: Int32Array): Int32Array {
  const out: number[] = [];
  for (let k = 0; k < ranges.length; k += 2) {
    if (ranges[k] >= ranges[k + 1]) continue;
    const a = atoms.residue(ranges[k]);
    const b = atoms.residue(ranges[k + 1] - 1) + 1;
    if (out.length && out[out.length - 1] >= a) out[out.length - 1] = Math.max(out[out.length - 1], b);
    else out.push(a, b);
  }
  return Int32Array.from(out);
}

/**
 * The frames a prefix names: over nodes, a complex's units, a unit, a molecule or a sequence, else the complexes held
 * under the node (a source's); over a Location, its parts' systems with the positions it holds; over Keys, their
 * extents. A node with no frame gives nothing (monadic), never an error.
 */
export function framesOf(ctx: Ctx, prefix: Named, kind: string): Frame[] {
  const r = prefix.kind === 'nodes' ? undefined : fromNamed(prefix);
  if (r?.type === 'location') return locationFrames(ctx, r.location, kind);
  if (r?.type === 'keys') return locationFrames(ctx, keysExtent(ctx, r.keys), kind);
  const out: Frame[] = [];
  for (const m of prefixNodes(prefix, kind)) {
    const address = printAddress(m.address);
    const units = ctx.model.units?.(address);
    if (units) {
      complexFrames(ctx, address, units, out);
      continue;
    }
    const f = systemFrame(ctx, address);
    if (f) {
      out.push(f);
      continue;
    }
    for (const c of ctx.model.complexes?.(address) ?? []) complexFrames(ctx, c, ctx.model.units?.(c) ?? [], out);
  }
  return out;
}

// ---- reading properties ----

/** A frame's reads, made once a frame and a property; what lacks and what is pending counted for one issue each. */
export class FrameReads {
  private readonly readers = new Map<PropertyDef, Reader | 'pending' | undefined>();
  private components?: QueryComponents | null;
  private atomsSystem?: QueryAtoms | null;

  constructor(
    private readonly ctx: Ctx,
    readonly frame: Frame,
  ) {}

  /** The frame's atoms, or undefined (a polymer's units have none in this version). */
  atoms(): QueryAtoms | undefined {
    if (this.atomsSystem === undefined) this.atomsSystem = this.ctx.model.atoms?.(this.frame.system.address) ?? null;
    return this.atomsSystem ?? undefined;
  }

  private comps(): QueryComponents | undefined {
    if (this.components === undefined) this.components = this.ctx.model.components?.(this.frame.system.address) ?? null;
    return this.components ?? undefined;
  }

  /** A property's reader on the frame: undefined when the frame lacks it, `'pending'` while its values have not arrived. */
  reader(d: PropertyDef): Reader | 'pending' | undefined {
    if (this.readers.has(d)) return this.readers.get(d);
    const r = this.make(d);
    this.readers.set(d, r);
    return r;
  }

  private hooked(name: string, system: string | undefined, field?: string): Reader | 'pending' | undefined {
    if (!system) return undefined;
    return byColumn(this.ctx.model.property?.(system, name, field), (r, a) => (a >= 0 ? a : r));
  }

  private make(d: PropertyDef): Reader | 'pending' | undefined {
    const { ctx, frame } = this;
    const { model } = ctx;
    if (d.read) return this.extension(d);
    const scalar = (v: PropValue): Reader | undefined => (v === undefined ? undefined : () => v);
    const unit = frame.unit;
    switch (d.name) {
      case 'entity':
        return scalar(unit?.entity);
      case 'molecule-class': {
        const m = model.molecule?.(frame.of);
        return scalar(m ? m.class : model.components?.(frame.of)?.polymer ? 'polymer' : undefined);
      }
      case 'identifier': {
        const ids = model.molecule?.(frame.of)?.identifiers;
        return scalar(ids?.length ? ids : undefined);
      }
      case 'chain':
        return scalar(unit?.chain);
      case 'auth-chain':
        return scalar(unit?.authChain);
      case 'operator':
        return scalar(unit?.operator);
      case 'seq': {
        const start = frame.system.parent?.start ?? 0;
        return (r) => start + r + 1;
      }
      case 'index':
        return (r) => r;
      case 'auth-seq':
      case 'ins': {
        const n = findNumbering(frame.system, 'auth');
        if (!n) return numberingPending(frame.system, 'auth') ? 'pending' : undefined;
        const ins = d.name === 'ins';
        return (r) => {
          const l = n.numbering.labelAt(r - n.shift);
          if (l === undefined) return undefined;
          const p = parseLabel(l);
          if (p.int === undefined) return undefined;
          return ins ? (p.kind === 1 ? p.suffix : '') : p.int;
        };
      }
      case 'component':
      case 'component-class':
      case 'standard':
      case 'modified':
        return this.componentReader(d.name);
      case 'atom':
      case 'element':
      case 'atom-index': {
        const atoms = this.atoms();
        if (!atoms) return undefined;
        if (d.name === 'atom-index') return (_, a) => a;
        const col = d.name === 'atom' ? atoms.names : atoms.element;
        return (_, a) => col.get(a) ?? undefined;
      }
      default:
        // `b-factor`, `occupancy`: values channels on the atoms the hook answers, lacking otherwise
        return this.hooked(d.name, this.atoms()?.address);
    }
  }

  private componentReader(name: string): Reader | 'pending' | undefined {
    const c = this.comps();
    if (!c) return undefined;
    const polymer = c.polymer;
    if (c.pending && name !== 'standard' && name !== 'modified') return 'pending';
    if (name === 'component') return (r) => c.at(r);
    if (name === 'component-class') {
      const { model } = this.ctx;
      return (r) => {
        const words: string[] = [];
        for (const code of c.at(r))
          for (const w of (model.compoundClass?.(code) ?? '').toLowerCase().split(/[^a-z0-9]+/)) if (w) words.push(w);
        return words.length ? words : undefined;
      };
    }
    const non = new Set<number>();
    const m = c.modified;
    if (polymer && m)
      for (let i = 0; i < m.position.length; i++) for (let k = 0; k < m.length[i]; k++) non.add(m.position[i] + k);
    if (name === 'standard') return (r) => polymer && !non.has(r);
    const items = modifiedItems(this.ctx, this.frame.system.address);
    return (r) => (polymer && non.has(r)) || items.has(r);
  }

  private extension(d: PropertyDef): Reader | 'pending' | undefined {
    const { model } = this.ctx;
    const system =
      d.level === 'atom'
        ? this.atoms()?.address
        : d.level === 'residue'
          ? this.frame.system.address
          : d.level === 'unit'
            ? this.frame.unit?.address
            : this.frame.of;
    if (!system) return undefined;
    const read = d.read!;
    const col = typeof read === 'function' ? read(model, system) : model.property?.(system, read.values, read.field);
    return byColumn(col, d.level === 'atom' ? (_, a) => a : d.level === 'residue' ? (r) => r : () => 0);
  }
}

/** A hook's column as a reader, indexed by `at`. */
function byColumn(col: PropertyRead, at: (r: number, a: number) => number): Reader | 'pending' | undefined {
  if (col === undefined || col === 'pending') return col;
  return (r, a) => {
    const i = at(r, a);
    return i >= 0 && i < col.length ? col.get(i) : undefined;
  };
}

/** Annotation types that say a position is modified (a database's PTM known by name only). */
const ModifiedTypes = new Set(['modified', 'glycosylation', 'lipidation']);

/** The positions under the `modified`-type items of the annotations on a system (by `type`, else `category`). */
function modifiedItems(ctx: Ctx, system: string): Set<number> {
  const out = new Set<number>();
  for (const address of ctx.model.annotationsOn?.(system) ?? []) {
    const a = ctx.model.annotation(address);
    if (!a) continue;
    const read = a.items();
    stage(ctx, read.pending);
    for (const it of read.items) {
      tick(ctx);
      const t = it.fields?.type ?? it.fields?.category;
      if (typeof t === 'string' && ModifiedTypes.has(t.toLowerCase()) && (it.system ?? a.target) === system)
        for (let p = it.start; p < it.end; p++) out.add(p);
    }
  }
  return out;
}

/** Stages a frame whose property is pending: its whole residues system waits (`seq.query.property.pending`, once). */
export function stageFrame(ctx: Ctx, f: Frame) {
  const l = rangesLocation(f.system.address, [0, f.system.length]);
  if (!isEmptyLocation(l)) stage(ctx, l);
}
