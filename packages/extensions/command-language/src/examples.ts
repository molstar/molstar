/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Dialect } from './language/ir.js';

// Scripts as people write them (_spec §9, Import): each imported into the console with the structures it loads, then
// replayed. The stories show them (apps/prism-viewer, "PyMOL and VMD commands"); a headless test runs every one and
// fails on any error (src/prism-viewer/_test/examples.test.ts).

/** A structure a script loads, by the name it uses and the example file it is (`data/examples`). */
export interface ExampleStructure {
  readonly name: string;
  readonly example: '1cbs_updated.cif' | '1tup.bcif' | '1crn.cif';
}

export interface ExampleScript {
  readonly title: string;
  readonly dialect: Dialect;
  /** The script files, the first the one imported (the others it includes). */
  readonly scripts: readonly { readonly name: string; readonly text: string }[];
  readonly structures: readonly ExampleStructure[];
  /** The warnings it gives, known (where Mol* does less than PyMOL or VMD): a part of each. */
  readonly warns?: readonly string[];
}

/** A binding-site figure in PyMOL: a quiet cartoon, the ligand and the residues around it, a distance, labels. */
export const BindingSitePml: ExampleScript = {
  title: 'A binding-site figure (PyMOL)',
  dialect: 'pymol',
  scripts: [
    {
      name: 'binding-site.pml',
      text: `# Retinoic acid in CRABP-II (1CBS): the pocket around the ligand
load 1cbs.cif, cbs
hide everything
show cartoon
color grey80, cbs
set cartoon_transparency, 0.5

select lig, resn REA
select pocket, byres (polymer within 5 of lig)
show sticks, lig or pocket
hide sticks, hydro
util.cbag lig
util.cbac pocket

# The carboxylate's polar contacts: every pair within 4 Å
distance polar, lig and elem O, pocket and (elem N or elem O), 4.0
label pocket and name CA, resn

deselect
orient lig
bg_color white
count_atoms pocket
`,
    },
  ],
  structures: [{ name: '1cbs.cif', example: '1cbs_updated.cif' }],
  // The camera's: in the browser only (orient needs a canvas).
  warns: ['transparent in every representation', 'orients along the structure'],
};

/** The p53 core domain on DNA in PyMOL: chains, a B-factor spectrum, the zinc. */
export const ProteinDnaPml: ExampleScript = {
  title: 'Protein and DNA (PyMOL)',
  dialect: 'pymol',
  scripts: [
    {
      name: 'p53-dna.pml',
      text: `# p53 core domain bound to DNA (1TUP)
load 1tup.bcif, p53
as cartoon
util.cbc p53
show sticks, resn DA+DC+DG+DT
util.cnc resn DA+DC+DG+DT
spectrum b, blue_white_red, chain A
show spheres, resn ZN
set sphere_scale, 0.6
zoom chain A
get_names
`,
    },
  ],
  structures: [{ name: '1tup.bcif', example: '1tup.bcif' }],
  warns: ['colors whole components by a theme', 'sizes the atoms in every representation'],
};

/** A script that includes another (`@style.pml`): crambin's secondary structure and its disulfides. */
export const IncludePml: ExampleScript = {
  title: 'A script and the one it includes (PyMOL)',
  dialect: 'pymol',
  scripts: [
    {
      name: 'crambin.pml',
      text: `load 1crn.cif, crn
@style.pml
count_atoms resn CYS
`,
    },
    {
      name: 'style.pml',
      text: `# The house style: secondary structure, cysteines as sticks, sulfurs yellow
as cartoon
util.ss
show sticks, resn CYS and not name N+C+O
color yellow, resn CYS and elem S
`,
    },
  ],
  structures: [{ name: '1crn.cif', example: '1crn.cif' }],
};

/** VMD's way, from a .tcl: reps one by one, a transparent surface around the ligand, the display set. */
export const RepsTcl: ExampleScript = {
  title: 'Representations (VMD)',
  dialect: 'vmd',
  scripts: [
    {
      name: 'reps.tcl',
      text: `# 1CBS in VMD: cartoon, the ligand as licorice, a transparent surface around it
mol new 1cbs.cif
mol delrep 0 top

mol representation NewCartoon
mol color Structure
mol selection {protein}
mol addrep top

mol representation Licorice
mol color Name
mol selection {resname REA}
mol addrep top

mol representation QuickSurf
mol color ColorID 6
mol material Transparent
mol selection {protein and same residue as within 5 of resname REA}
mol addrep top

color Display Background white
axes location off
display projection Orthographic
puts "reps added"
`,
    },
  ],
  structures: [{ name: '1cbs.cif', example: '1cbs_updated.cif' }],
};

export const ExampleScripts: readonly ExampleScript[] = [BindingSitePml, ProteinDnaPml, IncludePml, RepsTcl];
