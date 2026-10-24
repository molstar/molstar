/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Unsupported } from '../table.js';

// VMD's commands (and Tcl's) that are not done here (vmd.md §4): each fails saying what it is and what to do instead.

const tcl = 'The console runs VMD commands one by one; Tcl’s control flow and procedures are not run.';
const coordinates = 'Mol* shows structures as they are loaded; it does not change coordinates.';

export const VmdUnsupported: readonly Unsupported[] = [
  ...['proc', 'if', 'for', 'foreach', 'while', 'expr', 'switch', 'eval', 'return', 'break', 'continue'].map((name) => ({
    name,
    why: 'Tcl’s control flow is not supported.',
    hint: tcl,
  })),
  ...['package', 'namespace', 'lindex', 'llength', 'lappend', 'list', 'incr', 'format', 'string', 'array', 'dict'].map(
    (name) => ({ name, why: 'Tcl’s commands are not supported (set, unset, puts and $variables are).', hint: tcl }),
  ),
  {
    name: 'measure fit',
    why: 'The matrix is for $sel move, which changes coordinates.',
    hint: 'Superpose with pymol align mobile, target.',
  },
  {
    name: 'measure rmsd',
    why: 'RMSD between selections is not supported yet.',
    hint: 'pymol align prints the RMSD of a superposition.',
  },
  {
    name: 'measure contacts',
    why: 'Contacts are not supported.',
    hint: 'Interactions are shown from the Structure panel.',
  },
  {
    name: 'measure hbonds',
    why: 'Hydrogen bonds are not supported here.',
    hint: 'Interactions are shown from the Structure panel.',
  },
  { name: 'measure sasa', why: 'Surface areas are not supported.' },
  {
    name: 'mol addfile',
    why: 'Adding coordinates to a molecule is not supported yet.',
    hint: 'Load a trajectory from the Load menu.',
  },
  { name: 'mol rename', why: 'Renaming molecules is not supported yet.' },
  { name: 'mol fix', why: 'Fixing molecules is not supported.' },
  { name: 'mol free', why: 'Fixing molecules is not supported.' },
  {
    name: 'mol active',
    why: 'Active and inactive molecules are not supported.',
    hint: 'mol on and mol off show and hide one.',
  },
  {
    name: 'mol inactive',
    why: 'Active and inactive molecules are not supported.',
    hint: 'mol on and mol off show and hide one.',
  },
  { name: 'color change', why: 'Changing VMD’s color definitions is not supported.' },
  { name: 'color scale', why: 'Color scales are not supported.' },
  { name: 'color Name', why: 'Changing the colors of names is not supported.' },
  { name: 'color Element', why: 'Changing the colors of elements is not supported.' },
  {
    name: 'pbc',
    why: 'Periodic boundaries are not supported here.',
    hint: 'Show a simulation box from the trajectory panel.',
  },
  { name: 'graphics', why: 'Drawing graphics is not supported.' },
  { name: 'draw', why: 'Drawing graphics is not supported.' },
  { name: 'material', why: 'Changing materials is not supported.', hint: 'mol material Opaque or Transparent.' },
  { name: 'menu', why: 'VMD’s menus are not here.' },
  { name: 'plugin', why: 'VMD’s plugins are not here.' },
  { name: 'light', why: 'Lights are set in the canvas settings.' },
  { name: 'logfile', why: 'Log files are not written.' },
  { name: 'topo', why: 'Changing the topology is not supported.', hint: coordinates },
  { name: 'quit', why: 'There is nothing to quit.' },
  { name: 'exit', why: 'There is nothing to quit.' },
];

export const VmdUnsupportedByName: ReadonlyMap<string, Unsupported> = new Map(
  VmdUnsupported.map((u) => [u.name.toLowerCase(), u]),
);
