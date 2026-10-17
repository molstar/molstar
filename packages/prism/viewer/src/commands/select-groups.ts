/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandGroup } from '@molstar/prism-ui/core';

// The `select` group's tree (commands.md §5; the commands: `select.ts`): Select with its shared mode, By query (the
// registry's categories below it, made with the commands) and In the loaded structures. No plugin here: the stories'
// fixture shows the same group.

/** Select's modes (the `select.set` operation's), shared by the group: set once in the command menu's header. */
const Modes = [
  ['set', 'Replace'],
  ['add', 'Add'],
  ['remove', 'Remove'],
  ['intersect', 'Intersect'],
] as const;

export const SelectGroupId = 'select';

/** The groups that do not depend on what is loaded. */
export const SelectCommandGroups: readonly CommandGroup[] = [
  {
    id: SelectGroupId,
    title: 'Select',
    icon: 'mouse-pointer-click',
    order: 0,
    shared: { name: 'mode', title: 'Mode', options: Modes, default: 'set' },
  },
  { id: 'select.query', title: 'By query', icon: 'list-tree', parent: SelectGroupId, order: 0 },
  { id: 'select.structure', title: 'In the loaded structures', icon: 'structure', parent: SelectGroupId, order: 1 },
  { id: 'select.structure.entities', title: 'Entities', icon: 'pick-entity', parent: 'select.structure', order: 0 },
  { id: 'select.structure.chains', title: 'Chains', icon: 'pick-chain', parent: 'select.structure', order: 1 },
  {
    id: 'select.structure.residues',
    title: 'Ligands and other residues',
    icon: 'ligand',
    parent: 'select.structure',
    order: 2,
  },
  { id: 'select.structure.elements', title: 'Elements', icon: 'pick-atom', parent: 'select.structure', order: 3 },
];

const StructureQueryId = /^select\.(entity|chain|residue|element)\.[^\s]+$/;

/** The id of a command made from the loaded structures (an entity, chain, residue or element): it exists once one is loaded. */
export function isStructureQueryId(id: string): boolean {
  return StructureQueryId.test(id);
}
