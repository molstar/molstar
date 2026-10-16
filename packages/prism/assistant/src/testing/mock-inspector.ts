/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { JsonObject } from '@molstar/ai-agent';
import type {
  AttemptView,
  LimitRowView,
  LimitsData,
  ModelStepView,
  StepView,
  TokensView,
  ToolDefinitionView,
  ToolStepView,
  ToolsData,
  ToolView,
  TraceData,
  TurnView,
  WindowEntryView,
} from '@molstar/ai-trace';
import { observableState } from '@molstar/prism-ui';
import { InspectorModel } from '../inspector/model.js';

// A three-turn session in the viewer mode, for the Inspector's stories and tests: turn 1 loads 1TQN and shows its heme;
// turn 2 has a rate-limited model call (429, a wait, the retry), loads the selections skill (the next request does not
// extend the last: a prefix break), a tool call that fails on a wrong keyword before the model fixes it, and a tool that
// runs a helper model (a nested run); turn 3 shows the active site, or is still running.

const tok = (text: string) => Math.max(1, Math.round(text.length / 4));

const Model = 'claude-sonnet-4-5';
const Helper = 'claude-haiku-4-5';

/** Steps numbered as they are made: the same ids on every build of the session. */
function steps() {
  let n = 0;
  const model = (
    start: number,
    end: number | undefined,
    out: {
      toolCalls?: string[];
      reply?: string;
      tokens?: [number, number];
      firstOutput?: number;
      under?: string;
      model?: string;
      cached?: boolean;
    },
    attempts?: AttemptView[],
  ): ModelStepView => {
    n += 1;
    const tokens: TokensView | undefined = out.tokens && {
      input: out.tokens[0],
      output: out.tokens[1],
      // After a turn's first call, the provider reads the earlier prompt from its cache.
      ...(out.cached !== false && out.tokens[0] > 3200 ? { cached: Math.round(out.tokens[0] * 0.82) } : {}),
    };
    return {
      kind: 'model',
      id: `s${n}`,
      callId: `call-${n}`,
      model: out.model ?? Model,
      attempts: attempts ?? [
        {
          ordinal: 1,
          start,
          ...(out.firstOutput ? { firstOutput: out.firstOutput } : {}),
          ...(end !== undefined ? { end } : {}),
          status: end === undefined ? 'running' : 'succeeded',
        },
      ],
      ...(tokens ? { tokens: { ...tokens, credits: creditsOf(tokens) } } : {}),
      ...(out.reply ? { reply: out.reply, finish: 'complete' as const } : {}),
      ...(out.toolCalls ? { toolCalls: out.toolCalls, finish: 'tool-calls' as const } : {}),
      ...(out.under ? { under: out.under } : {}),
    };
  };
  const tool = (
    start: number,
    end: number | undefined,
    name: string,
    args: JsonObject,
    result: string | undefined,
    opts: { title?: string; read?: boolean; changed?: [string, string][]; failed?: boolean } = {},
  ): ToolStepView => {
    n += 1;
    const status = end === undefined ? 'running' : opts.failed ? 'failed' : 'succeeded';
    return {
      kind: 'tool',
      id: `s${n}`,
      callId: `call-${n}`,
      toolId: name,
      name,
      ...(opts.title ? { title: opts.title } : {}),
      ...(opts.read ? { read: true } : {}),
      args,
      status,
      ...(result !== undefined ? { result } : {}),
      attempts: [
        {
          ordinal: 1,
          start,
          ...(end !== undefined ? { end } : {}),
          status,
          ...(opts.failed && result ? { failure: result.split('\n')[0] } : {}),
        },
      ],
      ...(opts.changed ? { changed: opts.changed.map(([id, label]) => ({ id, label })) } : {}),
    };
  };
  return { model, tool };
}

// Rates of the mock's price table: 1 credit per 1k uncached input, 0.1 per 1k cached, 5 per 1k output.
function creditsOf(t: TokensView): number {
  const cached = t.cached ?? 0;
  return ((t.input - cached) * 1 + cached * 0.1 + t.output * 5) / 1000;
}

const Guidance = `You work in Mol*, a molecular viewer, for the person at the screen. You act only through tools: the viewer's commands.

Each request starts with the state (what is loaded, selected and shown) and the focus (what the person points at).
Rules:
1. Refer to structures and components by the ids the state lists. Never invent ids.
2. Read before you change: scene.get when the state is not enough.
3. Make changes with tools; do not just describe them. Several small calls are fine.
4. "this", "here", "the selection" mean the focus.
5. Load a skill (disclose_skills) before you use its tools.
…`;

const SelectionsGuidance = `Skill selections: atoms, residues and chains by distance, neighborhood and logic.

1. selection.set takes an expression in PyMOL's syntax ("within 5 of resname HEM and polymer") or MolQL.
2. query_reference explains a keyword or symbol: its arguments, types and defaults.
3. bounds.get measures what a selection spans, before you frame it.`;

const VolumesGuidance = `Skill volumes: density maps drawn with the structure.

1. volume.load brings a map in: an EMDB entry, a PDB entry's X-ray density, or a file.
2. volume.set_look draws it: an isosurface at a level (relative or absolute), or a slice.`;

const target: JsonObject = {
  description: 'What it acts on: "selection", a structure id ("1tqn"), or a component id ("1tqn.ligand").',
  anyOf: [{ type: 'string' }, { type: 'object', properties: { ligand: { type: 'string' } } }],
};

function definition(
  name: string,
  description: string,
  input: JsonObject,
  opts: Partial<
    Pick<ToolView, 'read' | 'approval' | 'enabled' | 'state' | 'loadedAtTurn' | 'extra' | 'unavailable'>
  > = {},
): ToolView {
  return {
    id: name,
    name,
    description,
    input,
    read: opts.read ?? false,
    approval: opts.approval ?? false,
    enabled: opts.enabled ?? true,
    state: opts.state ?? 'offered',
    tokens: tok(JSON.stringify({ name, description, input })),
    ...(opts.loadedAtTurn !== undefined ? { loadedAtTurn: opts.loadedAtTurn } : {}),
    ...(opts.extra ? { extra: true } : {}),
    ...(opts.unavailable ? { unavailable: opts.unavailable } : {}),
  };
}

const Core: ToolView[] = [
  definition(
    'scene.get',
    'The scene in detail: every structure with its components, representations and colors, the selection, the camera. Read-only.',
    { type: 'object', properties: { structure: { type: 'string', description: 'Only this structure (its id).' } } },
    { read: true, state: 'loaded' },
  ),
  definition(
    'structure.load',
    'Load a structure: a PDB or AlphaFold DB entry by id, or a URL. Asks first: it fetches from a server.',
    {
      type: 'object',
      properties: {
        source: { type: 'string', enum: ['pdb', 'afdb', 'url'], description: 'Where it comes from.' },
        id: { type: 'string', description: 'The entry id ("1tqn", "P69905"), or the URL.' },
        assembly: { type: 'string', description: 'An assembly id; the first by default.' },
      },
      required: ['id'],
    },
    { approval: true, state: 'loaded' },
  ),
  definition(
    'selection.set',
    'Set the selection: a component, a ligand, residues, or an expression. Replaces what is selected unless mode says otherwise.',
    {
      type: 'object',
      properties: {
        target,
        expression: {
          type: 'string',
          description: 'An expression, in PyMOL syntax or MolQL (with the selections skill).',
        },
        mode: { type: 'string', enum: ['set', 'add', 'remove'], description: 'How it combines with the selection.' },
      },
    },
    { state: 'loaded' },
  ),
  definition(
    'representation.set',
    'Show a target as a representation, or restyle what is shown of it: its type, color and transparency.',
    {
      type: 'object',
      properties: {
        target,
        type: {
          type: 'string',
          enum: ['cartoon', 'ball-and-stick', 'spacefill', 'surface', 'putty'],
          description: 'The representation.',
        },
        color: { type: 'string', description: 'A color ("red", "#d2491f") or a theme in words ("by residue type").' },
        transparency: { type: 'number', description: '0 (opaque) to 1.' },
      },
      required: ['target'],
    },
    { state: 'loaded' },
  ),
  definition(
    'camera.focus',
    'Turn the camera to a target and frame it, optionally with its surroundings.',
    {
      type: 'object',
      properties: { target, radius: { type: 'number', description: 'Extra room around it, in Å.' } },
      required: ['target'],
    },
    { state: 'loaded' },
  ),
  definition(
    'disclose_skills',
    'Load skills: their guidance and their tools, all a request needs in one call. Skills: "selections" (by distance and logic), "volumes" (density maps).',
    {
      type: 'object',
      properties: {
        skills: {
          type: 'array',
          items: { type: 'string', enum: ['selections', 'volumes'] },
          description: 'Skill names.',
        },
      },
      required: ['skills'],
    },
    { state: 'loaded' },
  ),
];

const Selections: ToolView[] = [
  definition(
    'query_reference',
    'Selection keywords and MolQL symbols in detail: arguments, types, defaults, an example. Read-only.',
    {
      type: 'object',
      properties: { symbols: { type: 'array', items: { type: 'string' }, description: 'Keywords or symbols.' } },
      required: ['symbols'],
    },
    { read: true, state: 'loaded', loadedAtTurn: 2 },
  ),
  definition(
    'bounds.get',
    'What a selection spans: its center, radius and bounding box, in Å. Read-only.',
    { type: 'object', properties: { target } },
    { read: true, state: 'loaded', loadedAtTurn: 2 },
  ),
];

const Volumes: ToolView[] = [
  definition(
    'volume.load',
    'Load a density map: an EMDB entry, the X-ray density of a PDB entry, or a URL. Asks first: it fetches from a server.',
    {
      type: 'object',
      properties: {
        emdb: { type: 'string', description: 'An EMDB id ("EMD-1181").' },
        pdb_density: { type: 'string', description: 'A PDB id whose X-ray density to load.' },
        url: { type: 'string', description: 'A map file (CCP4, MRC, DSN6).' },
      },
    },
    { approval: true },
  ),
  definition('volume.set_look', 'Draw a map: an isosurface at a level, or a slice; its color and transparency.', {
    type: 'object',
    properties: {
      volume: { type: 'string', description: "The map's id." },
      kind: { type: 'string', enum: ['isosurface', 'slice'] },
      level: { type: 'number', description: 'The isovalue, relative (σ) unless absolute is set.' },
      absolute: { type: 'boolean' },
    },
    required: ['volume'],
  }),
];

const Code: ToolView[] = [
  definition(
    'code.run',
    'Run a short script against the viewer: for what no command does in one step.',
    { type: 'object', properties: { source: { type: 'string', description: 'TypeScript.' } }, required: ['source'] },
    { unavailable: 'Needs the code module, which this app does not include.' },
  ),
];

const Single: ToolView[] = [
  definition(
    'snapshot.add',
    'Save the view as a snapshot, with a title and a description.',
    {
      type: 'object',
      properties: {
        title: { type: 'string' },
        description: { type: 'string', description: 'Markdown.' },
      },
      required: ['title'],
    },
    { state: 'deferred' },
  ),
  definition(
    'structure.remove',
    'Remove a structure and everything shown of it. Asks first: it cannot be undone after the next action.',
    { type: 'object', properties: { structure: { type: 'string', description: 'Its id.' } }, required: ['structure'] },
    { approval: true, state: 'deferred' },
  ),
  definition(
    'export.image',
    'Save an image of the canvas: its size and transparency.',
    {
      type: 'object',
      properties: {
        width: { type: 'integer', description: 'Pixels.' },
        height: { type: 'integer', description: 'Pixels.' },
        transparent: { type: 'boolean', description: 'No background.' },
      },
    },
    { enabled: false, state: 'deferred' },
  ),
];

export const MockTools: ToolsData = {
  groups: [
    { id: 'core', kind: 'core', title: 'Core tools', tools: Core },
    {
      id: 'selections',
      kind: 'skill',
      title: 'Selections',
      summary: 'Select by distance, neighborhood and logic (PyMOL selection syntax, MolQL).',
      guidance: SelectionsGuidance,
      guidanceTokens: tok(SelectionsGuidance),
      tools: Selections,
    },
    {
      id: 'volumes',
      kind: 'skill',
      title: 'Volumes',
      summary: 'Density maps: EMDB entries, X-ray density, map files; drawn as isosurfaces or slices.',
      guidance: VolumesGuidance,
      guidanceTokens: tok(VolumesGuidance),
      tools: Volumes,
    },
    {
      id: 'code',
      kind: 'skill',
      title: 'Code',
      summary: 'Scripts against the viewer, for what no command does in one step.',
      guidance: 'Skill code: write short TypeScript against the viewer’s API and run it with code.run.',
      guidanceTokens: 24,
      unavailable: 'Needs the code module, which this app does not include.',
      tools: Code,
    },
    { id: 'single', kind: 'single', title: 'Tools loaded by name', tools: Single },
  ],
};

const definitions = (tools: readonly ToolView[]): ToolDefinitionView[] =>
  tools.map((t) => ({ id: t.id, name: t.name, description: t.description }));
const toolsText = (tools: readonly ToolView[]) => tools.map((t) => `${t.name}: ${t.description}`).join('\n');

const State1 = `Structures: none loaded.
Selection: empty.`;

const State2 = `Structures:
- 1tqn "1TQN" Cytochrome P450 3A4, chain A (486 residues)
  - 1tqn.polymer: cartoon, colored by chain
  - 1tqn.ligand HEM 508: ball-and-stick, colored by element
Selection: HEM 508 (A), 43 atoms
Camera: on HEM 508`;

const State3 = `Structures:
- 1tqn "1TQN" Cytochrome P450 3A4, chain A (486 residues)
  - 1tqn.polymer: cartoon, colored by chain
  - 1tqn.ligand HEM 508: ball-and-stick, colored by element
  - 1tqn.site: ball-and-stick, colored by residue name (24 residues)
Selection: 24 residues (A): PHE 57, ARG 105, ILE 118, ALA 297, THR 309, CYS 442 …
Camera: on HEM 508`;

const SceneResult = `1tqn "1TQN" Cytochrome P450 3A4 (PDB), 1 model, 3,898 atoms
  polymer: chain A, residues 27-499, cartoon, colored by chain
  ligand HEM 508: 43 atoms, ball-and-stick, colored by element
  site: 24 residues, ball-and-stick, colored by residue name
selection: 24 residues
camera: target HEM 508, radius 9.2 Å`;

function turnsOf(running: boolean): TurnView[] {
  const { model, tool } = steps();

  const t1: StepView[] = [
    model(0, 1650, { toolCalls: ['structure.load(source="pdb", id="1tqn")'], tokens: [3120, 24], firstOutput: 1300 }),
    tool(
      1660,
      2900,
      'structure.load',
      { source: 'pdb', id: '1tqn' },
      'Loaded 1TQN: Cytochrome P450 3A4 with heme and metyrapone. 1 model, chain A (residues 27-499); ligand HEM 508.\nAdded: 1tqn, 1tqn.polymer, 1tqn.ligand.',
      { title: 'Load 1TQN from PDBe', changed: [['1tqn', '1TQN']] },
    ),
    model(2910, 4700, {
      toolCalls: [
        'selection.set(target={"ligand":"HEM"})',
        'representation.set(target="selection", type="ball-and-stick")',
      ],
      tokens: [3390, 86],
      firstOutput: 4100,
    }),
    tool(4710, 4760, 'selection.set', { target: { ligand: 'HEM' } }, 'Selected HEM 508 (A): 43 atoms.', {
      title: 'Select HEM 508',
    }),
    tool(
      4761,
      5200,
      'representation.set',
      { target: 'selection', type: 'ball-and-stick' },
      'Ball-and-stick on HEM 508 (A).\nChanged: 1tqn.ligand.',
      { title: 'Show HEM 508 as ball-and-stick', changed: [['1tqn.ligand', 'HEM 508']] },
    ),
    model(5210, 6500, { toolCalls: ['camera.focus(target="selection")'], tokens: [3610, 31], firstOutput: 6000 }),
    tool(6510, 7300, 'camera.focus', { target: 'selection' }, 'The camera is on HEM 508 (radius 9.2 Å).', {
      title: 'Focus the camera on HEM 508',
    }),
    model(7310, 8400, {
      reply: 'Loaded 1TQN and showed its heme (HEM 508) as ball-and-stick, with the camera on it.',
      tokens: [3700, 28],
      firstOutput: 7900,
    }),
  ];

  const t2: StepView[] = [
    model(60000, 63100, { toolCalls: ['disclose_skills(skills=["selections"])'], tokens: [3980, 18] }, [
      {
        ordinal: 1,
        start: 60000,
        end: 60380,
        status: 'failed',
        failure: '429 Too many requests: the rate limit for claude-sonnet-4-5',
      },
      {
        ordinal: 2,
        start: 61880,
        firstOutput: 62600,
        end: 63100,
        status: 'succeeded',
        waitBefore: { reason: 'rate limited', ms: 1500 },
      },
    ]),
    tool(
      63110,
      63130,
      'disclose_skills',
      { skills: ['selections'] },
      'Loaded skill selections: its guidance and 2 tools (query_reference, bounds.get).',
      { title: 'Load the selections skill' },
    ),
    {
      ...model(63140, 65600, {
        toolCalls: ['query_reference(symbols=["within"])'],
        tokens: [5240, 22],
        firstOutput: 64900,
        cached: false,
      }),
      prefixBreak: { reason: 'the tool definitions changed (the selections skill loaded)', kept: 1, of: 9 },
    },
    tool(
      65610,
      65640,
      'query_reference',
      { symbols: ['within'] },
      'within D of S: atoms within D Å of selection S (D a number; S any selection). Example: within 5 of resname HEM.\nSee also: around, beyond, byres.',
      { title: 'Look up “within”', read: true },
    ),
    model(65650, 67700, {
      toolCalls: ['selection.set(expression="byres (within 5 of resn HEM) and polymer")'],
      tokens: [5520, 44],
      firstOutput: 67100,
    }),
    tool(
      67710,
      67760,
      'selection.set',
      { expression: 'byres (within 5 of resn HEM) and polymer' },
      'Could not read the expression: unknown keyword "resn" at 17.\nDid you mean resname? Keywords: resname, resi, chain, name, elem.',
      { title: 'Select the residues near HEM', failed: true },
    ),
    model(67770, 69500, {
      toolCalls: [
        'selection.set(expression="byres (within 5 of resname HEM) and polymer")',
        'representation.set(target="selection", type="ball-and-stick", color="by residue type")',
      ],
      tokens: [5700, 61],
      firstOutput: 68900,
    }),
    tool(
      69510,
      69600,
      'selection.set',
      { expression: 'byres (within 5 of resname HEM) and polymer' },
      'Selected 24 residues in chain A: PHE 57, ARG 105, ILE 118, ALA 297, THR 309, CYS 442 and 18 more.',
      { title: 'Select the residues near HEM' },
    ),
  ];
  const colors = tool(
    69601,
    71900,
    'representation.set',
    { target: 'selection', type: 'ball-and-stick', color: 'by residue type' },
    'Colored 24 residues by residue name (theme residue-name), as ball-and-stick.\nChanged: 1tqn.site.',
    { title: 'Color the selection by residue type', changed: [['1tqn.site', '24 residues']] },
  );
  t2.push(
    colors,
    // The tool asks a small model which color theme the words mean: a nested run under its call.
    model(69700, 71500, {
      reply: '{ "theme": "residue-name" }',
      tokens: [640, 12],
      firstOutput: 71300,
      under: colors.callId,
      model: Helper,
    }),
    model(71910, 73300, {
      reply: 'Selected the 24 residues within 5 Å of the heme and colored them by residue type.',
      tokens: [5980, 34],
      firstOutput: 72800,
    }),
  );

  const t3: StepView[] = [
    model(120000, 121900, { toolCalls: ['scene.get()'], tokens: [6200, 9], firstOutput: 121700 }),
    tool(121910, 121950, 'scene.get', {}, SceneResult, { title: 'Read the scene', read: true }),
    ...(running
      ? [model(121960, undefined, {})]
      : [
          model(121960, 123800, {
            toolCalls: [
              'representation.set(target="1tqn.site", type="ball-and-stick", transparency=0)',
              'camera.focus(target="1tqn.site", radius=4)',
            ],
            tokens: [6950, 70],
            firstOutput: 123100,
          }),
          tool(
            123810,
            124300,
            'representation.set',
            { target: '1tqn.site', type: 'ball-and-stick', transparency: 0 },
            'Ball-and-stick on 1tqn.site (24 residues).\nChanged: 1tqn.site.',
            { title: 'Show the site as ball-and-stick', changed: [['1tqn.site', '24 residues']] },
          ),
          tool(124310, 125100, 'camera.focus', { target: '1tqn.site', radius: 4 }, 'The camera is on 1tqn.site.', {
            title: 'Focus the camera on the site',
          }),
          model(125110, 126400, {
            reply: 'The active site: the heme with the 24 residues around it, the camera on them.',
            tokens: [7180, 26],
            firstOutput: 126000,
          }),
        ]),
  ];

  return [
    turn('turn-1', 'Load 1TQN and show the heme as ball-and-stick', undefined, t1),
    turn(
      'turn-2',
      'Select the residues within 5 Å of the heme and color them by residue type',
      'selection: HEM 508 (A), 43 atoms',
      t2,
    ),
    turn(
      'turn-3',
      'Show the active site around the heme and focus on it',
      'selection: 24 residues (A)',
      t3,
      running ? 'thinking' : 'done',
    ),
  ];
}

function turn(
  id: string,
  prompt: string,
  focus: string | undefined,
  steps: StepView[],
  status: TurnView['status'] = 'done',
): TurnView {
  const models = steps.filter((s): s is ModelStepView => s.kind === 'model');
  const first = steps[0]?.attempts[0]?.start ?? 0;
  const last = steps.at(-1)?.attempts.at(-1)?.end;
  const sum = (f: (t: TokensView) => number) => models.reduce((a, m) => a + (m.tokens ? f(m.tokens) : 0), 0);
  const ended = status === 'done' || status === 'stopped' || status === 'failed';
  return {
    id,
    prompt,
    ...(focus ? { focus } : {}),
    status,
    start: first,
    ...(ended && last !== undefined ? { end: last } : {}),
    tokens: {
      input: sum((t) => t.input),
      output: sum((t) => t.output),
      cached: sum((t) => t.cached ?? 0),
      credits: sum((t) => t.credits ?? 0),
    },
    steps,
    milestones: [
      { position: 0, at: first, kind: 'start', detail: 'request appended' },
      ...steps
        .filter((s) => s.under === undefined)
        .flatMap((s, i) => {
          const end = s.attempts.at(-1)?.end;
          if (end === undefined) return [];
          return [{ position: i + 1, at: end, kind: `${s.kind === 'model' ? 'model' : 'tools'} completed` }];
        }),
      ...(ended && last !== undefined ? [{ position: steps.length + 1, at: last, kind: 'end', detail: status }] : []),
    ],
  };
}

/** What each model call saw: the guidance, the tool definitions, then every turn's state, request, calls and results. */
function withWindows(turns: TurnView[]): TurnView[] {
  const history: WindowEntryView[] = [];
  let seen = 0;
  let skills = false;
  const states = [State1, State2, State3];
  const base = (): WindowEntryView[] => {
    const tools = skills ? [...Core, ...Selections] : Core;
    return [
      { role: 'system', part: 'guidance', title: 'Guidance', text: Guidance, tokens: tok(Guidance) * 4 },
      ...(skills
        ? [
            {
              role: 'system' as const,
              part: 'guidance' as const,
              title: 'Skill: selections',
              text: SelectionsGuidance,
              tokens: tok(SelectionsGuidance),
            },
          ]
        : []),
      {
        role: 'tools',
        part: 'tools',
        title: `${tools.length} tool definitions`,
        text: toolsText(tools),
        tokens: tools.reduce((a, t) => a + t.tokens, 0),
      },
    ];
  };
  return turns.map((t, ti) => {
    const state = states[ti] ?? '';
    history.push({ role: 'user', part: 'state', title: 'State', text: state, tokens: tok(state) * 3 });
    const request = `${t.focus ? `Focus: ${t.focus}\n\n` : ''}Request: ${t.prompt}`;
    history.push({ role: 'user', part: 'conversation', title: 'Request', text: request, tokens: tok(request) });
    const steps = t.steps.map((s): StepView => {
      if (s.kind === 'model') {
        if (s.under !== undefined) {
          // The helper sees only what its tool gives it.
          const ask = 'Which Mol* color theme do these words mean? "by residue type". Answer as JSON: { "theme": id }.';
          return {
            ...s,
            window: {
              entries: [
                { role: 'user', part: 'conversation', title: 'Request', text: ask, tokens: tok(ask), fresh: true },
              ],
              limit: 200_000,
            },
          };
        }
        const head = base();
        const entries = [...head, ...history.map((e, i) => (i >= seen ? { ...e, fresh: true } : e))];
        seen = history.length;
        const text = s.reply ?? (s.toolCalls ?? []).join('\n');
        if (s.attempts.at(-1)?.end !== undefined)
          history.push({
            role: 'assistant',
            part: 'conversation',
            title: s.reply ? 'Reply' : 'Tool calls',
            text,
            tokens: s.tokens?.output ?? tok(text),
          });
        return { ...s, window: { entries, limit: 200_000 } };
      }
      if (s.status !== 'running')
        history.push({
          role: 'tool',
          part: 'conversation',
          title: `${s.name} result`,
          text: s.result ?? '',
          tokens: tok(s.result ?? ''),
          ...(s.status === 'failed' ? { failed: true } : {}),
        });
      if (s.name === 'disclose_skills') skills = true;
      return s;
    });
    return { ...t, steps };
  });
}

/** The viewer-mode session: three turns, the last still running when `running`. */
export function mockTrace(options: { running?: boolean } = {}): TraceData {
  const running = !!options.running;
  const turns = withWindows(turnsOf(running));
  const latest = turns.at(-1)!;
  const uncached = latest.tokens.input - (latest.tokens.cached ?? 0);
  return {
    harnessId: 'h-viewer-7c1e2a',
    execution: {
      seq: running ? 118 : 131,
      headId: running ? '5dcdfe502d9e7e2f003fd34f12345509' : 'a41f0c9e8b27d6e35f1c0a9b8e7d6c5b',
      ...(running ? { pendingStep: 'model' as const } : {}),
    },
    budget: {
      consumed: uncached,
      ceiling: 400_000,
      kind: 'uncachedInput',
      dimension: 'uncached input',
      scope: 'latest turn',
      ...(latest.tokens.credits !== undefined ? { credits: latest.tokens.credits } : {}),
    },
    preamble: {
      guidance: Guidance,
      guidanceTokens: tok(Guidance) * 4,
      tools: definitions(Core),
      toolsTokens: Core.reduce((a, t) => a + t.tokens, 0),
    },
    turns,
    loaded: [
      ...Core.map((t) => ({ id: t.id, name: t.name })),
      ...Selections.map((t) => ({ id: t.id, name: t.name, turn: 2 })),
    ],
    declared: Single.map((t) => t.id),
  };
}

/** The turn budget, the model's limits, and credits as information (never a limit). */
export function mockLimits(trace: TraceData | undefined): LimitsData {
  const latest = trace?.turns.at(-1);
  const steps = latest?.steps ?? [];
  const tokens = latest?.tokens;
  const all = trace?.turns.reduce((a, t) => a + (t.tokens.credits ?? 0), 0);
  const credits = (n: number | undefined) => Math.round((n ?? 0) * 100) / 100;
  return {
    groups: [
      {
        id: 'budget',
        title: 'Turn budget',
        note: 'What one turn may use: it stops when one runs out. The bars are the latest turn; each limit goes up to its ceiling, from the next turn.',
        rows: [
          {
            key: 'modelCalls',
            label: 'Model calls',
            note: 'requests to the model in one turn',
            unit: 'calls',
            value: 24,
            default: 24,
            used: steps.filter((s) => s.kind === 'model').length,
            ceiling: 100,
            editable: true,
          },
          {
            key: 'toolCalls',
            label: 'Tool calls',
            note: 'commands the agent runs in one turn',
            unit: 'calls',
            value: 64,
            default: 64,
            used: steps.filter((s) => s.kind === 'tool').length,
            ceiling: 256,
            editable: true,
          },
          {
            key: 'uncachedInput',
            label: 'Uncached input',
            note: 'tokens the provider reads anew',
            unit: 'tokens',
            value: 400_000,
            default: 400_000,
            used: tokens ? tokens.input - (tokens.cached ?? 0) : 0,
            ceiling: 3_200_000,
            step: 1000,
            editable: true,
          },
          {
            key: 'cacheRead',
            label: 'Cache reads',
            note: 'input served from the provider’s cache',
            unit: 'tokens',
            value: 4_000_000,
            default: 4_000_000,
            used: tokens?.cached ?? 0,
            ceiling: 16_000_000,
            step: 1000,
            editable: true,
          },
          {
            key: 'output',
            label: 'Output',
            note: 'what the model writes, reasoning included',
            unit: 'tokens',
            value: 32_000,
            default: 32_000,
            used: tokens?.output ?? 0,
            ceiling: 128_000,
            step: 1000,
            editable: true,
          },
        ],
      },
      {
        id: 'model',
        title: 'Model',
        note: `For ${Model}; applies from the next request.`,
        rows: [
          {
            key: 'contextWindow',
            label: 'Context window',
            note: 'the most the model reads at once',
            unit: 'tokens',
            value: 200_000,
            default: 200_000,
            ceiling: 1_000_000,
            step: 1000,
            editable: true,
          },
          {
            key: 'compaction',
            label: 'Compaction trigger',
            note: 'older turns are summarized past this share of the window',
            unit: '%',
            value: 80,
            default: 80,
            ceiling: 95,
            editable: true,
          },
          {
            key: 'thinking',
            label: 'Thinking',
            note: 'how long it may reason before it answers',
            unit: '',
            value: 1,
            default: 1,
            choices: [
              { value: 0, label: 'Off' },
              { value: 1, label: 'Low' },
              { value: 2, label: 'Medium' },
              { value: 3, label: 'High' },
            ],
            editable: true,
          },
        ],
      },
      {
        id: 'credits',
        title: 'Credits',
        note: 'What the agent costs, in the app’s credits (1 credit is 1,000 uncached input tokens on the default model). Shown, never a limit.',
        rows: [
          {
            key: 'credits.turn',
            label: 'Latest turn',
            unit: 'credits',
            value: credits(tokens?.credits),
            default: credits(tokens?.credits),
            editable: false,
          },
          {
            key: 'credits.session',
            label: 'This session',
            unit: 'credits',
            value: credits(all),
            default: credits(all),
            editable: false,
          },
        ],
      },
    ],
  };
}

/**
 * The Inspector over the mock session (`empty`: before the agent has run): fixed data, and limits kept in memory (a
 * change, or a reset, rewrites them).
 */
export function mockInspector(options: { running?: boolean; empty?: boolean } = {}): InspectorModel {
  const trace = options.empty ? undefined : mockTrace(options);
  const limits = observableState<LimitsData>(mockLimits(trace));
  return new InspectorModel({
    trace: observableState<TraceData | undefined>(trace).asReadonly,
    tools: observableState<ToolsData>(MockTools).asReadonly,
    limits: limits.asReadonly,
    setLimit: (key, value) =>
      limits.set({
        groups: limits.value.groups.map((g) => ({
          ...g,
          rows: g.rows.map((r): LimitRowView => (r.key === key ? { ...r, value: value ?? r.default } : r)),
        })),
      }),
  });
}
