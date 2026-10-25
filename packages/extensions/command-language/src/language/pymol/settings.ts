/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { pymolColor } from '../colors.js';
import { fail, nearest, type Span, warning } from '../diagnostics.js';
import type { RepType } from '../ir.js';
import type { Compiled, Compiler } from '../table.js';
import { type PymolArgs, readBool } from './args.js';

// PyMOL's settings that Mol* has (pymol.md §4): `set name, value [, selection]`. Any other fails, listing these.
// PyMOL accepts a unique prefix of a setting's name; so does this.

export interface PymolSetting {
  readonly summary: string;
  /** The value `unset` puts back, as text. */
  readonly default: string;
  set(value: string, at: Span, a: PymolArgs, c: Compiler): Compiled;
}

function number(value: string, at: Span, min = -Infinity, max = Infinity): number {
  const v = Number(value);
  if (!Number.isFinite(v) || v < min || v > max) {
    fail('argument', `'${value}' is not a number${Number.isFinite(min) ? ` from ${min} to ${max}` : ''}.`, at);
  }
  return v;
}

const transparency =
  (reps: readonly RepType[] | undefined, what: string): PymolSetting['set'] =>
  (value, at, a, c) => {
    const alpha = number(value, at, 0, 1);
    if (reps) {
      c.report(
        warning(
          'unsupported-argument',
          `Mol* makes the atoms transparent in every representation, not only ${what}.`,
          at,
        ),
      );
    }
    return [{ kind: 'transparency', sel: a.sel('selection'), alpha, ...(reps ? { reps } : {}) }];
  };

const size =
  (reps: readonly RepType[], pymolDefault: number): PymolSetting['set'] =>
  (value, at, a) => [{ kind: 'size', sel: a.sel('selection'), factor: number(value, at, 0) / pymolDefault, reps }];

const color: PymolSetting['set'] = (value, at) => {
  const rgb = pymolColor(value);
  if (rgb === undefined) fail('argument', `'${value}' is not a PyMOL color.`, at);
  return [{ kind: 'canvas', background: rgb }];
};

export const PymolSettings: Readonly<Record<string, PymolSetting>> = {
  bg_rgb: { summary: 'The background color.', default: 'black', set: color },
  bg_color: { summary: 'The background color.', default: 'black', set: color },
  orthoscopic: {
    summary: 'Orthographic projection (1) or perspective (0).',
    default: '0',
    set: (v, at) => [{ kind: 'canvas', projection: readBool(v, at) ? 'orthographic' : 'perspective' }],
  },
  field_of_view: {
    summary: 'The camera’s field of view, in degrees.',
    default: '20',
    set: (v, at) => [{ kind: 'canvas', fov: number(v, at, 1, 179) }],
  },
  depth_cue: {
    summary: 'Fog: far atoms fade (1) or not (0).',
    default: '1',
    set: (v, at) => [{ kind: 'canvas', fog: readBool(v, at) }],
  },
  fog: {
    summary: 'Fog: far atoms fade (1) or not (0).',
    default: '1',
    set: (v, at) => [{ kind: 'canvas', fog: number(v, at, 0) > 0 }],
  },
  transparency: {
    summary: 'Transparency of surfaces, 0 (opaque) to 1.',
    default: '0',
    set: transparency(['surface', 'mesh', 'gaussian-surface'], 'surfaces'),
  },
  cartoon_transparency: {
    summary: 'Transparency of cartoons, 0 to 1.',
    default: '0',
    set: transparency(['cartoon', 'putty', 'tube'], 'cartoons'),
  },
  stick_transparency: {
    summary: 'Transparency of sticks, 0 to 1.',
    default: '0',
    set: transparency(['sticks'], 'sticks'),
  },
  sphere_transparency: {
    summary: 'Transparency of spheres, 0 to 1.',
    default: '0',
    set: transparency(['spheres'], 'spheres'),
  },
  sphere_scale: { summary: 'The size of spheres (PyMOL’s 1 is Mol*’s).', default: '1', set: size(['spheres'], 1) },
  stick_radius: {
    summary: 'The radius of sticks (PyMOL’s 0.25 is Mol*’s).',
    default: '0.25',
    set: size(['sticks'], 0.25),
  },
  ray_opaque_background: {
    summary: 'Whether png keeps the background (1) or makes it transparent (0).',
    default: '1',
    set(v, at, _a, c) {
      c.state.session.vars.set('ray_opaque_background', readBool(v, at) ? '1' : '0');
      return [];
    },
  },
};

/** PyMOL settings Mol* has no counterpart for, but scripts set often: what to say. */
const PymolSettingsUnsupported: Readonly<Record<string, string>> = {
  cartoon_color: 'Color the atoms with `color`; a color per representation is not supported yet.',
  stick_color: 'Color the atoms with `color`; a color per representation is not supported yet.',
  surface_color: 'Color the atoms with `color`; a color per representation is not supported yet.',
  sphere_color: 'Color the atoms with `color`; a color per representation is not supported yet.',
  ray_trace_mode: 'Mol* renders live; its outlines are in the canvas settings.',
  ray_shadows: 'Mol*’s shadows are in the canvas settings.',
  antialias: 'Mol* antialiases by itself.',
  label_size: 'Label sizes are set from the Measurements panel.',
  label_color: 'Label colors are set from the Measurements panel.',
  dash_color: 'Measurement colors are set from the Measurements panel.',
  cartoon_fancy_helices: 'Mol*’s cartoon is its own.',
  cartoon_side_chain_helper: 'Show side chains with `show sticks, …`.',
};

/** A setting by name, or a unique prefix of one; else fails naming the ones Mol* has. */
export function settingNamed(name: string, at: Span): PymolSetting {
  const n = name.trim().toLowerCase();
  const exact = PymolSettings[n];
  if (exact) return exact;
  const prefixed = Object.keys(PymolSettings).filter((k) => k.startsWith(n));
  if (prefixed.length === 1) return PymolSettings[prefixed[0]];
  const unsupported = PymolSettingsUnsupported[n];
  if (unsupported) fail('unsupported', `The setting '${n}' is not supported.`, at, unsupported);
  const near = nearest(n, Object.keys(PymolSettings));
  return fail(
    'unsupported',
    `The setting '${n}' is not supported.`,
    at,
    `${near ? `Did you mean ${near}? ` : ''}Settings: ${Object.keys(PymolSettings).join(', ')}.`,
  );
}
