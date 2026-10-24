/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Unsupported } from '../table.js';

// PyMOL's commands that Mol* does not do (pymol.md §5): each fails saying what it is and what to do instead, so a
// script never silently skips one. A name in neither this list nor the table is an unknown command.

const python = 'Scripts here are PyMOL commands; Python is not run.';
const objects = 'Mol* shows parts of a structure as components: `show sticks, …`, `select name, …`.';
const editing = 'Mol* shows structures as they are loaded; it does not edit them.';
const movies = 'Use the snapshots and animation panels.';

export const PymolUnsupported: readonly Unsupported[] = [
  { name: 'python', why: 'Python blocks are not run.', hint: python },
  { name: 'run', why: 'Python scripts are not run.', hint: python },
  { name: 'cmd', why: 'PyMOL’s Python API is not run.', hint: python },
  { name: 'alter', why: 'Changing atom properties with Python expressions is not supported.', hint: editing },
  { name: 'alter_state', why: 'Changing coordinates is not supported.', hint: editing },
  {
    name: 'iterate',
    why: 'Reading properties with Python is not supported yet.',
    hint: '`count_atoms` counts a selection.',
  },
  { name: 'iterate_state', why: 'Reading coordinates with Python is not supported.', hint: python },
  { name: 'create', why: 'New objects from selections are not supported.', hint: objects },
  { name: 'extract', why: 'New objects from selections are not supported.', hint: objects },
  { name: 'copy', why: 'Copying objects is not supported.', hint: 'Load it again with `fetch` or `load`.' },
  { name: 'split_states', why: 'Splitting states into objects is not supported.', hint: '`frame n` shows one state.' },
  { name: 'fab', why: 'Building peptides is not supported.', hint: editing },
  { name: 'fuse', why: 'Joining objects is not supported.', hint: editing },
  { name: 'h_add', why: 'Adding hydrogens is not supported.', hint: editing },
  { name: 'h_fill', why: 'Adding hydrogens is not supported.', hint: editing },
  { name: 'mutagenesis', why: 'The mutagenesis wizard is not supported.', hint: editing },
  { name: 'wizard', why: 'Wizards are not supported.' },
  { name: 'edit', why: 'Editing is not supported.', hint: editing },
  { name: 'bond', why: 'Adding bonds is not supported.', hint: editing },
  { name: 'unbond', why: 'Removing bonds is not supported.', hint: editing },
  { name: 'sort', why: 'Sorting atoms is not supported.', hint: editing },
  { name: 'sculpt_activate', why: 'Sculpting is not supported.', hint: editing },
  { name: 'mset', why: 'Movies are not supported.', hint: movies },
  { name: 'mview', why: 'Movies are not supported.', hint: movies },
  { name: 'mdo', why: 'Movies are not supported.', hint: movies },
  { name: 'save', why: 'Saving files from a script is not supported.', hint: 'Save the session from the menu (⌘S).' },
  {
    name: 'symexp',
    why: 'Symmetry mates as objects are not supported.',
    hint: 'Show symmetry from the Structure panel.',
  },
  { name: 'ramp_new', why: 'Color ramps are not supported.', hint: '`spectrum` colors along a property.' },
  { name: 'set_name', why: 'Renaming objects is not supported yet.' },
  { name: 'order', why: 'Ordering objects is not supported.' },
  { name: 'group', why: 'Grouping objects is not supported.' },
  { name: 'ungroup', why: 'Grouping objects is not supported.' },
  { name: 'dss', why: 'Mol* assigns secondary structure when it loads; there is nothing to recompute.' },
  {
    name: 'set_view',
    why: 'Setting PyMOL’s view matrix is not supported yet.',
    hint: '`zoom`, `orient` and `turn` move the camera.',
  },
  { name: 'scene', why: 'Scenes are not supported yet.', hint: 'Use the snapshots panel.' },
  { name: 'isomesh', why: 'Maps are not supported yet.', hint: 'Load density from the Volume panel.' },
  { name: 'isosurface', why: 'Maps are not supported yet.', hint: 'Load density from the Volume panel.' },
  { name: 'map_new', why: 'Maps are not supported.' },
  { name: 'quit', why: 'There is nothing to quit.' },
];

export const PymolUnsupportedByName: ReadonlyMap<string, Unsupported> = new Map(
  PymolUnsupported.map((u) => [u.name, u]),
);
