/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewerSkill } from '@molstar/prism-viewer/spec';

// The agent's way to the command language (_spec §9): a skill, not a core tool, for "run this .pml" and for people who
// paste PyMOL or VMD lines. The agent's own vocabulary stays the viewer's operations; this runs the person's commands
// as they wrote them, one undoable step, under the agent's policy (agent.md §4: `script.run` says what it reaches, and
// refuses what would need the person's yes inside it).

export const CommandLanguageSkill: ViewerSkill = {
  id: 'command-language',
  title: 'PyMOL and VMD commands',
  summary: 'run PyMOL or VMD commands or a script the person gives, as they wrote them',
  guidance: [
    'Skill command-language: run_commands runs PyMOL (dialect pymol) or VMD (dialect vmd) commands as one undoable',
    'step: what the person pasted or asked to run, as written. One command a line or separated by ;. Prefer the',
    "viewer's own tools for what you decide yourself; use this when the person speaks PyMOL or VMD. A line that does",
    'not read fails with the line and a hint: fix it, or tell the person. Each call starts without Tcl variables. Removing',
    'structures and saving images are refused here: use the viewer’s own tools, which ask the person.',
  ].join('\n'),
  tools: [
    {
      id: 'run_commands',
      operation: 'script.run',
      description: 'Runs PyMOL or VMD commands (the person’s, as written) as one undoable step.',
    },
  ],
};
