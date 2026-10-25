/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The agent's tool and skill declarations (agent.md §3): plain data, no imports, so a spec reads or extends them without
// the toolkit (`toolkit.ts`, on `@molstar/ai` and zod) in the page's startup bundle. The toolkit turns them into
// tools when the agent is made, in the assistant's chunk.

/** A tool over an operation (agent.md §2). */
export interface OperationToolSpec {
  /** The tool's id, in the model's words (`load_entry`, `color`). */
  readonly id: string;
  /** The operation it runs. */
  readonly operation: string;
  /** What the model reads: when to use it, what it does; the operation's title when not given. */
  readonly description?: string;
  /** The operation's params offered, by name (default: all that JSON can carry). */
  readonly fields?: readonly string[];
  /** Values set for the model (`mode: 'add'`). */
  readonly fixed?: Readonly<Record<string, unknown>>;
  /** Field formats said in their descriptions (`target`: a target string, targets.md). */
  readonly formats?: Readonly<Record<string, string>>;
}

/** A skill: its guidance and its tools over operations; `feature` names the extension that brings the operations. */
export interface SkillDeclaration {
  readonly id: string;
  readonly title: string;
  readonly summary: string;
  readonly guidance: string;
  readonly tools: readonly OperationToolSpec[];
  /** The extension whose operations they are (lazy ones load when the skill is disclosed). */
  readonly feature?: string;
}

/**
 * A change to a list of defaults (the spec's `agent.core`, `agent.skills`): `remove` takes entries out by id, `add`
 * appends, and an added entry whose id is a default's takes its place. The defaults are never restated.
 */
export interface AgentListChange<T extends { readonly id: string }> {
  readonly add?: readonly T[];
  readonly remove?: readonly string[];
}

/** The defaults with a change applied (`AgentListChange`); none: the defaults as they are. */
export function withChange<T extends { readonly id: string }>(
  defaults: readonly T[],
  change: AgentListChange<T> | undefined,
): readonly T[] {
  if (!change) return defaults;
  const removed = new Set(change.remove ?? []);
  const added = new Map((change.add ?? []).map((t) => [t.id, t] as const));
  const out: T[] = [];
  for (const t of defaults) {
    if (removed.has(t.id)) continue;
    const replaced = added.get(t.id);
    out.push(replaced ?? t);
    added.delete(t.id);
  }
  for (const t of added.values()) out.push(t);
  return out;
}

/** The core's operation tools (agent.md §3, the first scope: loading, selections, focus, styles, canvas, camera). */
export const CoreTools: readonly OperationToolSpec[] = [
  { id: 'get_selection', operation: 'select.get', description: 'What is selected now, per structure.' },
  {
    id: 'load_entry',
    operation: 'load.entry',
    description:
      'Loads an entry from a database (PDB, PDB-IHM, AlphaFold DB, EMDB, ModelArchive) by its id. An entry already loaded is not loaded again.',
  },
  { id: 'load_url', operation: 'load.url', description: 'Loads a structure from a URL (https).' },
  {
    id: 'select',
    operation: 'select.set',
    description:
      'Selects a target (mode set, add, remove, intersect). Targets: "ligand", "A:10-20", "resn HEM", MolQL.',
  },
  { id: 'clear_selection', operation: 'select.clear', description: 'Clears the selection.' },
  {
    id: 'focus',
    operation: 'focus.set',
    description: 'Focuses a target: the camera goes to it and its surroundings show.',
  },
  { id: 'apply_preset', operation: 'structure.preset', description: 'Applies a representation preset to structures.' },
  {
    id: 'apply_style',
    operation: 'styles.apply',
    description:
      'Applies a quick style: default, cartoon, spacefill, surface (representations), plain, illustrative (looks).',
  },
  {
    id: 'color',
    operation: 'representation.color',
    description: 'Colors a target: a CSS color (red, #ff8800) or a theme by name (chain-id, secondary-structure, …).',
  },
  { id: 'set_background', operation: 'canvas.background', description: 'Sets the background color (a CSS color).' },
  { id: 'reset_camera', operation: 'camera.reset', description: 'Resets the camera to show everything.' },
];

/** The skills (agent.md §3); each tool is left out where the viewer lacks its operation. */
export const Skills: readonly SkillDeclaration[] = [
  {
    id: 'selections',
    title: 'Selections',
    summary: 'select around, expand, invert; the target forms',
    guidance: [
      'Skill selections: targets are strings.',
      '- Keywords: "ligand", "polymer", "protein", "nucleic", "water", "ion", "all", "selection", "focus".',
      '- Addresses: "A" (chain), "A:10-20" (residues), "A:10@CA" (an atom), "resn HEM", "resi 10-20".',
      '- PyMOL-style text and MolQL also work; give `structure` when several are loaded.',
      'Use select_around for what is near a target, expand_selection to whole residues or chains.',
    ].join('\n'),
    tools: [
      { id: 'select_around', operation: 'select.around', description: 'Selects what is within a radius of a target.' },
      {
        id: 'expand_selection',
        operation: 'select.expand',
        description: 'Expands the selection to whole residues, chains or molecules.',
      },
      { id: 'invert_selection', operation: 'select.invert', description: 'Inverts the selection.' },
      { id: 'select_all', operation: 'select.all', description: 'Selects everything.' },
    ],
  },
  {
    id: 'styles',
    title: 'Styles',
    summary: 'representations, components, transparency, hiding parts, the unit cell, moving structures',
    guidance: [
      'Skill styles: a component is a part of a structure (ligand, polymer, a selection) with representations.',
      "add_component makes one from a target; set_representation sets a component's representation type;",
      'set_transparency and color act on targets; hide_contents hides entities or chains everywhere.',
      "show_unit_cell shows a crystal's unit cell; transform_structure moves or rotates a structure (about its",
      'center; compose: true adds to the current placement), reset_transform puts it back where it was loaded.',
      '- set_look: how representations are drawn (hydrogens, quality, ignoreLight, material, solidSurface, interior);',
      '  without node it sets the scene defaults, with node (a structure, component, representation, volume, shape: #ref)',
      "  those representations' own values, which a later default change keeps; reset_look takes them back to the",
      '  defaults, get_look lists the defaults and who has values of its own.',
    ].join('\n'),
    tools: [
      {
        id: 'add_component',
        operation: 'component.add',
        description: 'Makes a component of a target, with a representation.',
      },
      { id: 'remove_component', operation: 'component.remove', description: 'Removes a component.' },
      {
        id: 'set_look',
        operation: 'look.set',
        description: 'Sets the scene defaults of how representations are drawn, or the representations under a node.',
      },
      { id: 'reset_look', operation: 'look.reset', description: 'Takes representations back to the scene defaults.' },
      {
        id: 'get_look',
        operation: 'look.get',
        description: 'The scene defaults, and the representations with their own.',
      },
      {
        id: 'set_representation',
        operation: 'representation.set',
        description: "Sets a component's representation type.",
      },
      { id: 'set_transparency', operation: 'representation.transparency', description: 'Makes a target transparent.' },
      { id: 'hide_contents', operation: 'contents.visibility', description: 'Shows or hides an entity or a chain.' },
      {
        id: 'list_annotations',
        operation: 'annotations.list',
        description: 'What the structures carry and can show (quality scores, validation, B-factors, charges).',
      },
      {
        id: 'apply_annotation',
        operation: 'annotations.apply',
        description:
          'Shows what a structure carries as color: plddt, qmean, b-factor, a validation report, partial charges (list_annotations names them).',
      },
      {
        id: 'show_unit_cell',
        operation: 'structure.unit-cell',
        description: "Shows or hides a crystal's unit cell or a simulation's box, with its color and opacity.",
      },
      {
        id: 'transform_structure',
        operation: 'structure.transform',
        description:
          'Moves or rotates a structure: a translation in Å and a rotation (axis, angle in degrees), or a 4×4 matrix.',
      },
      {
        id: 'reset_transform',
        operation: 'structure.transform-reset',
        description: 'Puts a moved structure back where it was loaded.',
      },
      {
        id: 'set_palette',
        operation: 'representation.defaultPalette',
        description:
          'Colorblind-safe colors (Okabe–Ito) for chains, entities and other categories, or the default ones.',
      },
    ],
  },
  {
    id: 'canvas',
    title: 'Canvas',
    summary: 'lighting, illumination, camera mode, clipping',
    feature: 'canvas-settings',
    guidance:
      "Skill canvas: set_lighting turns shadows, ambient occlusion, outlines, fog and illumination on ({ name: 'on', params }) or off ({ name: 'off' }); what it leaves out stays.",
    tools: [
      { id: 'set_lighting', operation: 'canvas.lighting' },
      { id: 'set_camera_mode', operation: 'canvas.camera' },
      { id: 'set_clipping', operation: 'canvas.clipping' },
    ],
  },
  {
    id: 'camera',
    title: 'Camera',
    summary: 'where the camera is and what it looks at, in numbers; the recent foci',
    guidance: [
      'Skill camera: get_camera gives the position, target and up vectors (Å), the distance, field of view and mode.',
      'set_camera takes any of them ("x, y, z"; empty keeps); prefer focus or reset_camera to show a part.',
      'recent_foci lists what was focused, newest first, each with a target and structure that focus takes back.',
    ].join('\n'),
    tools: [
      { id: 'get_camera', operation: 'camera.get' },
      { id: 'set_camera', operation: 'camera.set' },
      { id: 'recent_foci', operation: 'focus.history' },
    ],
  },
  {
    id: 'measure',
    title: 'Measurements',
    summary: 'distances, angles, dihedrals, labels',
    feature: 'measurements',
    guidance: 'Skill measure: each point is a target that selects one atom or residue ("A:10@CA", "ligand").',
    tools: [
      { id: 'measure_distance', operation: 'measure.distance' },
      { id: 'measure_angle', operation: 'measure.angle' },
      { id: 'measure_dihedral', operation: 'measure.dihedral' },
      { id: 'add_label', operation: 'measure.label' },
      { id: 'list_measurements', operation: 'measure.list' },
      { id: 'remove_measurement', operation: 'measure.remove' },
    ],
  },
  {
    id: 'superpose',
    title: 'Superposition',
    summary: 'superpose structures by chains, TM-align, ligands',
    feature: 'superposition',
    guidance: 'Skill superpose: the first target (or structure) is the reference; results give the RMSD.',
    tools: [
      { id: 'superpose_chains', operation: 'superpose.chains' },
      { id: 'superpose_tm_align', operation: 'superpose.tm-align' },
      { id: 'superpose_ligands', operation: 'superpose.ligands' },
    ],
  },
  {
    id: 'snapshots',
    title: 'Snapshots',
    summary: 'save views of the scene and step through them',
    feature: 'snapshots',
    guidance: 'Skill snapshots: add_snapshot saves the scene as it is; apply_snapshot goes back to one.',
    tools: [
      { id: 'add_snapshot', operation: 'snapshot.add' },
      { id: 'apply_snapshot', operation: 'snapshot.apply' },
      { id: 'list_snapshots', operation: 'snapshot.list' },
      { id: 'play_snapshots', operation: 'snapshot.play' },
      { id: 'stop_snapshots', operation: 'snapshot.stop' },
    ],
  },
  {
    id: 'animation',
    title: 'Animation and trajectories',
    summary: 'spin, rock, wiggle; trajectory frames',
    guidance:
      'Skill animation: frames of a trajectory or an NMR ensemble are shown by index (from 0). A simulation in a periodic box shows its molecules whole by default; periodic_box moves them into the box (centered on the protein with center: polymers), keeps them from jumping between frames (nojump), or shows the file as written (off); periodic_images shows the neighbouring boxes (mode box tiles them). wiggle moves atoms by noise (dynamics) or by their uncertainty; type, speed, frequency shape it, tumble turns whole copies; what is not given stays.',
    tools: [
      { id: 'play_animation', operation: 'animation.play' },
      { id: 'stop_animation', operation: 'animation.stop' },
      { id: 'wiggle', operation: 'animation.wiggle' },
      { id: 'show_frame', operation: 'trajectory.frame' },
      { id: 'get_trajectory', operation: 'trajectory.get' },
      { id: 'periodic_box', operation: 'trajectory.periodic' },
      {
        id: 'periodic_images',
        operation: 'structure.symmetry',
        description:
          "A simulation box's periodic images (kind: cells, cells: 1 is the box and its 26 neighbours; off: the box alone), or a crystal's symmetry.",
      },
    ],
  },
  {
    id: 'export',
    title: 'Export',
    summary: 'images and session files',
    guidance: 'Skill export: export_image renders the view; save_session writes the scene with its data.',
    tools: [
      { id: 'export_image', operation: 'export.image' },
      { id: 'save_session', operation: 'session.save' },
    ],
  },
];
