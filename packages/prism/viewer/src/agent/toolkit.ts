/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  createToolkit,
  type ToolGroup,
  type Toolkit,
  type ToolkitPart,
  type ToolRegistration,
} from '@molstar/ai-tools';
import { type AgentPolicy, mayAsk } from './policies.js';
import { failed } from '@molstar/ai-tools';
import {
  type AgentViewer,
  type OperationToolSpec,
  operationTool,
  tool,
  toolResult,
  type ViewerAgentApi,
} from './tools.js';

// The core and the skills (agent.md §3): a small core from the first turn (Stories' is about 24 tools); the rest in
// skills, disclosed on request. A skill over a lazy feature is a lazy skill (plan track A6): in the catalog by its
// summary; disclosing it loads the feature, then its tools. A tool whose operation the viewer does not have (a feature
// left out of the spec) is left out of its skill, and a skill with none is left out.

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
    ].join('\n'),
    tools: [
      {
        id: 'add_component',
        operation: 'component.add',
        description: 'Makes a component of a target, with a representation.',
      },
      { id: 'remove_component', operation: 'component.remove', description: 'Removes a component.' },
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
        description: "Shows or hides a crystal structure's unit cell, with its color and opacity.",
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
        operation: 'canvas.palette',
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
    guidance: 'Skill canvas: set_lighting turns shadows, ambient occlusion, outlines and fog on or off.',
    tools: [
      { id: 'set_lighting', operation: 'canvas.lighting' },
      { id: 'set_illumination', operation: 'canvas.illumination' },
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
      { id: 'list_measurements', operation: 'measure.get' },
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
    guidance: 'Skill animation: frames of a trajectory or an NMR ensemble are shown by index (from 0).',
    tools: [
      { id: 'play_animation', operation: 'animation.play' },
      { id: 'stop_animation', operation: 'animation.stop' },
      { id: 'wiggle', operation: 'animation.wiggle' },
      { id: 'show_frame', operation: 'trajectory.frame' },
      { id: 'get_trajectory', operation: 'trajectory.get' },
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

/** The agent's own tools that are not operations: the scene's outline, asking the person, undoing its last turn. */
function ownTools(): ToolRegistration<ViewerAgentApi>[] {
  return [
    tool({
      id: 'undo_turn',
      description:
        'Undoes your previous turn at once (what it changed comes back as it was). Fails when the person changed something after it: then ask them.',
      input: { kind: 'object', fields: {} },
      modules: ['viewer'],
      run: async ({ api }) => {
        const turn = api.viewer.lastTurn;
        if (!turn) return failed('not-found', 'There is no earlier turn of yours to undo.');
        // Never the person's later steps (`discardLater`): that is theirs to choose.
        return toolResult(await api.viewer.api.run('history.undo-turn', { turn, discardLater: false }));
      },
    }),
    tool({
      id: 'get_scene',
      description: "The scene's outline: structures with their components and representations, volumes, measurements.",
      input: { kind: 'object', fields: {} },
      modules: ['account'],
      read: true,
      run: ({ api }) => ({ ok: true, did: 'The scene:', body: api.account.outline(400) }),
    }),
    tool({
      id: 'ask',
      description:
        'Asks the person a question when the request is unclear (which chain? load a second copy?); end the turn after it: their answer comes as the next request.',
      input: {
        kind: 'object',
        fields: {
          question: { kind: 'string' },
          options: { kind: 'array', item: { kind: 'string' }, optional: true, description: 'Choices to offer.' },
        },
      },
      modules: [],
      run: ({ args }) => ({
        ok: true,
        did: `Asked: ${args.question}${args.options?.length ? ` (${args.options.join(' / ')})` : ''}. End your turn now.`,
      }),
    }),
  ];
}

/** The tools of a skill the viewer has now (a tool whose operation is missing is left out). */
function toolsOf(
  viewer: AgentViewer,
  specs: readonly OperationToolSpec[],
  policy?: AgentPolicy,
): ToolRegistration<ViewerAgentApi>[] {
  return specs.filter((s) => viewer.describe(s.operation)).map((s) => operationTool(viewer, s, policy));
}

export interface ViewerToolkitOptions {
  /** Loads a lazy extension (a lazy skill's load); true when it is loaded. */
  readonly loadExtension: (id: string) => Promise<boolean>;
  /** Whether an extension is lazy and not loaded yet. */
  readonly isStub: (id: string) => boolean;
  /** What goes without asking (agent.md §4). */
  readonly policy?: AgentPolicy;
  /** The core's operation tools (default `CoreTools`; the spec's `agent.core`). */
  readonly core?: readonly OperationToolSpec[];
  /** The skills: the spec's (default `Skills`) and the extensions' (viewer-spec.md §3). */
  readonly skills?: readonly SkillDeclaration[];
}

/** The viewer's toolkit, and its refresh between turns. */
export interface ViewerToolkit extends Toolkit<ViewerAgentApi> {
  /**
   * The tools over operations described again as the viewer is now (a registry's choices, a feature added or gone:
   * agent.md §2); a part whose definitions changed is replaced, a tool gone is the eviction hint. Between turns only (a
   * turn's tool set is fixed at its start). True when anything changed.
   */
  refresh(): boolean;
}

interface RefreshablePart {
  readonly build: () => ToolkitPart<ViewerAgentApi>;
  signature: string;
  remove: () => void;
}

/** What a part tells the model: its tools' descriptions and schemas, its groups' tools. */
function signatureOf(part: ToolkitPart<ViewerAgentApi>): string {
  return JSON.stringify([
    (part.tools ?? []).map((t) => [t.definition.tool.id, t.definition.description, t.definition.inputSchema]),
    (part.groups ?? []).map((g) => [g.id, g.tools]),
  ]);
}

/** The viewer's toolkit: the core, and the skills (lazy ones over lazy features). */
export function createViewerToolkit(viewer: AgentViewer, options: ViewerToolkitOptions): ViewerToolkit {
  const skills = options.skills ?? Skills;
  const own = ownTools();
  const eager = skills.filter((s) => !s.feature || !options.isStub(s.feature));
  const lazy = skills.filter((s) => s.feature && options.isStub(s.feature));
  const toolkit = createToolkit<ViewerAgentApi>({ tools: own, core: own.map((t) => t.definition.tool.id), groups: [] });
  // Every part over operations is added on its own, so a refresh can replace it.
  const parts: RefreshablePart[] = [];
  const addPart = (build: () => ToolkitPart<ViewerAgentApi>) => {
    const part = build();
    parts.push({ build, signature: signatureOf(part), remove: toolkit.add(part) });
  };
  addPart(() => {
    const core = toolsOf(viewer, options.core ?? CoreTools, options.policy);
    return { tools: core, core: core.map((t) => t.definition.tool.id) };
  });
  const skillPart = (s: SkillDeclaration) => (): ToolkitPart<ViewerAgentApi> => {
    const tools = toolsOf(viewer, s.tools, options.policy);
    // A skill with no tool (its feature left out) is left out.
    return tools.length ? { tools, groups: [{ ...groupOf(s), tools: tools.map((t) => t.definition.tool.id) }] } : {};
  };
  for (const s of eager) addPart(skillPart(s));
  for (const s of lazy) {
    let removeStub = () => {};
    let loading: Promise<void> | undefined;
    removeStub = toolkit.add({
      groups: [
        {
          ...groupOf(s),
          tools: [],
          // The stubs declare their operations' flags (viewer-spec.md §3): a turn knows from its start which may ask,
          // by the rule the real tools use (`mayAsk`).
          pending: s.tools
            .filter((t) => viewer.describe(t.operation))
            .map((t) => ({ id: t.id, approval: mayAsk(viewer.describe(t.operation)!.definition) })),
          // Once (two disclosures share it); a failed load leaves the stub, and a later disclosure tries again.
          load: () =>
            (loading ??= (async () => {
              if (!(await options.loadExtension(s.feature!))) throw new Error(`${s.title} could not be loaded.`);
              removeStub();
              addPart(skillPart(s));
            })().catch((e) => {
              loading = undefined;
              throw e;
            })),
        },
      ],
    });
  }
  return Object.assign(toolkit, {
    refresh(): boolean {
      let changed = false;
      for (const part of parts) {
        const next = part.build();
        const signature = signatureOf(next);
        if (signature === part.signature) continue;
        // Out, then in: the same ids are taken by the new definitions.
        part.remove();
        part.remove = toolkit.add(next);
        part.signature = signature;
        changed = true;
      }
      return changed;
    },
  });
}

function groupOf(s: SkillDeclaration): Omit<ToolGroup<ViewerAgentApi>, 'tools'> {
  return { id: s.id, title: s.title, summary: s.summary, guidance: s.guidance };
}
