/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { fail, type Span, warning } from '../diagnostics.js';
import type { RepType, Sel } from '../ir.js';
import { pymolSelection } from '../selection.js';
import type { Compiler, ParamSpec } from '../table.js';
import { PymolColorNames, pymolColor } from '../colors.js';
import { nearest } from '../diagnostics.js';
import { type PymolArg, type PymolCall, unquote } from './tokenizer.js';

// A PyMOL command's arguments bound to its signature (_spec §3.1): positionals in order, then `key=value` in any
// order by the signature's names. An argument PyMOL takes that Mol* ignores is a warning (`ignored`); one PyMOL does
// not take is an error.

export class PymolArgs {
  readonly #bound = new Map<string, PymolArg>();

  constructor(
    readonly call: PymolCall,
    readonly params: readonly ParamSpec[],
    readonly c: Compiler,
    /** Arguments PyMOL takes that change nothing here: their name and what to say. */
    ignored: Readonly<Record<string, string>> = {},
  ) {
    let i = 0;
    for (const a of call.args) {
      if (a.key) {
        const key = a.key.toLowerCase();
        const p = params.find((x) => x.name === key);
        if (p) {
          this.#bound.set(p.name, a);
        } else if (key in ignored) {
          c.report(warning('unsupported-argument', `${a.key} is ignored: ${ignored[key]}`, a.keyAt ?? a.at));
        } else {
          fail(
            'argument',
            `${call.name} takes no argument '${a.key}'.`,
            a.keyAt ?? a.at,
            `It takes: ${params.map((x) => x.name).join(', ') || 'none'}.`,
          );
        }
        continue;
      }
      while (i < params.length && this.#bound.has(params[i].name)) i++;
      if (i >= params.length) {
        // A positional past the signature: PyMOL's own extra arguments (state, quiet) are the ignored ones by order.
        const extra = Object.keys(ignored)[i - params.length];
        if (extra) {
          c.report(warning('unsupported-argument', `${extra} is ignored: ${ignored[extra]}`, a.at));
          i++;
          continue;
        }
        fail(
          'argument',
          `${call.name} takes at most ${params.length} argument(s).`,
          a.at,
          usageHint(call.name, params),
        );
      }
      this.#bound.set(params[i].name, a);
      i++;
    }
    for (const p of params) {
      if (!p.optional && !this.#bound.has(p.name)) {
        fail('argument', `${call.name} needs ${p.name}.`, call.at, usageHint(call.name, params));
      }
    }
  }

  has(name: string): boolean {
    const a = this.#bound.get(name);
    return !!a && a.text.length > 0;
  }

  at(name: string): Span {
    return this.#bound.get(name)?.at ?? this.call.at;
  }

  /** The argument's text, quotes removed; `def` when absent. */
  text(name: string, def = ''): string {
    const a = this.#bound.get(name);
    return a && a.text ? unquote(a.text) : def;
  }

  sel(name: string, def = 'all'): Sel {
    const a = this.#bound.get(name);
    const text = a && a.text ? unquote(a.text) : def;
    return pymolSelection(text, this.c.state, a?.at ?? this.call.at);
  }

  number(name: string, def?: number): number {
    if (!this.has(name)) {
      if (def === undefined) return fail('argument', `${this.call.name} needs ${name}.`, this.call.at);
      return def;
    }
    const t = this.text(name);
    const v = Number(t);
    if (!Number.isFinite(v)) fail('argument', `${name} is a number, not '${t}'.`, this.at(name));
    return v;
  }

  int(name: string, def?: number): number {
    const v = this.number(name, def);
    if (!Number.isInteger(v)) fail('argument', `${name} is a whole number, not '${this.text(name)}'.`, this.at(name));
    return v;
  }

  bool(name: string, def: boolean): boolean {
    if (!this.has(name)) return def;
    return readBool(this.text(name), this.at(name));
  }

  color(name: string): number {
    const t = this.text(name);
    const v = pymolColor(t);
    if (v === undefined) {
      const near = nearest(t, PymolColorNames);
      fail(
        'argument',
        `'${t}' is not a PyMOL color.`,
        this.at(name),
        near ? `Did you mean ${near}?` : 'Give a name (red, slate, grey70), 0xRRGGBB or [r, g, b].',
      );
    }
    return v;
  }
}

export function readBool(t: string, at: Span): boolean {
  const s = t.trim().toLowerCase();
  if (s === '1' || s === 'on' || s === 'true' || s === 'yes') return true;
  if (s === '0' || s === 'off' || s === 'false' || s === 'no') return false;
  return fail('argument', `'${t}' is not on or off.`, at, 'Give 1 or 0, on or off.');
}

export function usageHint(name: string, params: readonly ParamSpec[]): string {
  const ps = params.map((p) => (p.optional ? `[${p.name}]` : p.name));
  return `Usage: ${name}${ps.length ? ` ${ps.join(', ')}` : ''}.`;
}

/** PyMOL's representation names, as the IR's (`everything`: all). */
export const PymolReps: Readonly<Record<string, readonly RepType[] | 'all'>> = {
  cartoon: ['cartoon'],
  ribbon: ['backbone'],
  sticks: ['sticks'],
  lines: ['lines'],
  spheres: ['spheres'],
  surface: ['surface'],
  mesh: ['mesh'],
  dots: ['points'],
  ellipsoids: ['ellipsoids'],
  nb_spheres: ['spheres'],
  nonbonded: ['points'],
  putty: ['putty'],
  everything: 'all',
};

/** Representation names PyMOL has that Mol* does not show, with what to do instead. */
export const PymolRepsUnsupported: Readonly<Record<string, string>> = {
  labels: 'Labels come from `label selection, expression`.',
  slice: 'Map slices are not shown; maps are phase 2.',
  extent: 'Map extents are not shown.',
  cgo: 'CGO objects are not read.',
  cell: 'Show the unit cell from the Structure panel.',
  callback: 'Callbacks are Python.',
  dashes: 'Measurements are made with `distance`.',
  angles: 'Measurements are made with `angle`.',
  dihedrals: 'Measurements are made with `dihedral`.',
};

/** The representation of an argument (singular forms read too: `stick`, `line`, `sphere`). */
export function pymolRep(text: string, at: Span): readonly RepType[] | 'all' {
  const t = unquote(text).toLowerCase();
  const plural = t.endsWith('s') ? t : `${t}s`;
  const reps = PymolReps[t] ?? PymolReps[plural];
  if (reps) return reps;
  const unsupported = PymolRepsUnsupported[t] ?? PymolRepsUnsupported[plural];
  if (unsupported) fail('unsupported', `The representation '${t}' is not supported.`, at, unsupported);
  const near = nearest(t, Object.keys(PymolReps));
  return fail(
    'argument',
    `'${t}' is not a PyMOL representation.`,
    at,
    near ? `Did you mean ${near}?` : `Representations: ${Object.keys(PymolReps).join(', ')}.`,
  );
}
