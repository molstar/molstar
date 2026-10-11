/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CommandManager } from '../../commands/manager.js';
import { rankBy } from '../../commands/match.js';
import type { CommandDefinition, CommandFormats, CommandResult } from '../../commands/types.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ShellModel } from '../../shell/shell.js';

// biome-ignore lint/complexity/noBannedTypes: the stories' commands see nothing of the plugin
type Ctx = {};

const Colors = ['red', 'orange', 'yellow', 'green', 'teal', 'blue', 'violet', 'magenta', 'white', 'grey', 'black'];
const Ligands = [
  { value: 'HEM 508', label: 'HEM 508', detail: 'Protoporphyrin IX containing Fe' },
  { value: 'SO4 1', label: 'SO4 1', detail: 'Sulfate ion' },
  { value: 'SO4 2', label: 'SO4 2', detail: 'Sulfate ion' },
];

/** A selection language: its examples, offered as completions of what is typed. */
interface Language {
  title: string;
  examples: readonly { value: string; label: string }[];
}

// What each language's expression looks like (Mol*'s transpilers read them for real; here, a light check).
const Languages: Record<string, Language> = {
  pymol: {
    title: 'PyMOL',
    examples: [
      { value: 'resn HEM', label: 'resn HEM' },
      { value: 'chain A and resi 100-120', label: 'chain A and resi 100-120' },
      { value: 'byres (resn HEM around 5)', label: 'byres (resn HEM around 5)' },
      { value: 'polymer and not hydro', label: 'polymer and not hydro' },
      { value: 'organic', label: 'organic' },
    ],
  },
  vmd: {
    title: 'VMD',
    examples: [
      { value: 'resname HEM', label: 'resname HEM' },
      { value: 'chain A and resid 100 to 120', label: 'chain A and resid 100 to 120' },
      { value: 'same residue as within 5 of resname HEM', label: 'same residue as within 5 of resname HEM' },
      { value: 'protein and backbone', label: 'protein and backbone' },
      { value: 'water', label: 'water' },
    ],
  },
  jmol: {
    title: 'Jmol',
    examples: [
      { value: 'HEM', label: 'HEM' },
      { value: ':A and 100-120', label: ':A and 100-120' },
      { value: 'within(5.0, HEM)', label: 'within(5.0, HEM)' },
      { value: 'protein and backbone', label: 'protein and backbone' },
      { value: 'water', label: 'water' },
    ],
  },
};

/** Balanced parentheses and something to read: the mocks' stand-in for a transpiler's parse. */
function checkExpression(language: Language, text: string) {
  const v = text.trim();
  if (!v) return { ok: false as const, message: `Type a ${language.title} expression.` };
  let depth = 0;
  for (const ch of v) {
    depth += ch === '(' ? 1 : ch === ')' ? -1 : 0;
    if (depth < 0) return { ok: false as const, message: 'A ")" without its "(".' };
  }
  if (depth > 0)
    return { ok: false as const, message: 'A "(" is not closed.', hint: `For example: ${language.examples[2].value}` };
  if (/\b(and|or|not)\s*$/.test(v))
    return { ok: false as const, message: `The expression ends in "${v.split(/\s+/).pop()}".` };
  return { ok: true as const, value: v };
}

/** Named argument types the stories' commands use: colors (names or hex), ligands of 1TQN, and expressions. */
const formats: CommandFormats<Ctx> = {
  get: (name) => {
    const language = name.startsWith('expression.') ? Languages[name.slice('expression.'.length)] : undefined;
    if (language)
      return {
        hint: () => `a ${language.title} expression`,
        check: (v) => checkExpression(language, v),
        // Examples to adapt: those that continue what is typed; taking one fills the expression for editing, and Enter
        // on the edited text runs it.
        complete: (text) => {
          const t = text.trim().toLowerCase();
          return language.examples
            .filter((e) => e.value.toLowerCase().startsWith(t) && e.value.toLowerCase() !== t)
            .map((e) => ({ ...e, detail: 'example', fill: true }));
        },
      };
    if (name === 'residues')
      return {
        hint: () => 'chains and residue ranges: A 10-20, B 5',
        check: (v) =>
          /^\s*[A-Za-z0-9]+(\s+\d+(-\d+)?)?(\s*,\s*[A-Za-z0-9]+(\s+\d+(-\d+)?)?)*\s*$/.test(v)
            ? { ok: true, value: v.trim() }
            : { ok: false, message: `"${v}" is not a residue list.`, hint: 'Chain, then residues: A 10-20, B 5' },
        complete: () => [
          { value: 'A 100-120', label: 'A 100-120', detail: 'a range in chain A' },
          { value: 'A 442, A 508', label: 'A 442, A 508', detail: 'Cys442 and the heme' },
        ],
      };
    if (name === 'color')
      return {
        hint: () => 'a color name or #hex',
        check: (v) =>
          Colors.includes(v.toLowerCase()) || /^#?[0-9a-f]{6}$/i.test(v)
            ? { ok: true, value: v.toLowerCase() }
            : { ok: false, message: `"${v}" is not a color.`, hint: 'Try a name (red, teal) or #hex.' },
        complete: (text) => rankBy(Colors, text, (c) => c).map((c) => ({ value: c, label: c })),
      };
    if (name === 'ligand')
      return {
        hint: () => 'a ligand of the loaded structures',
        check: (v) => {
          const hit = Ligands.find((l) => l.value.toLowerCase() === v.toLowerCase().trim());
          return hit ? { ok: true, value: hit.value } : { ok: false, message: `No ligand "${v}".` };
        },
        complete: (text) => rankBy(Ligands, text, (l) => `${l.label} ${l.detail}`),
      };
    return undefined;
  },
};

const done = (message: string): CommandResult => ({ message });

/** What a select command does to the current selection. */
const SelectMode = PD.Select(
  'add',
  [
    ['add', 'add to the selection'],
    ['remove', 'remove from the selection'],
    ['intersect', 'intersect with the selection'],
    ['set', 'replace the selection'],
  ],
  { label: 'mode' },
);
const modeLabel = (mode: string) =>
  ({ add: 'Added to the selection', remove: 'Removed from the selection', intersect: 'Intersected', set: 'Selected' })[
    mode
  ] ?? 'Selected';

/**
 * The stories' commands: what the plugin's would be (camera, view, panels, export, data, presets, selection,
 * focus), running against the mock shell. Each reports what it would do; panels and cards open for real.
 */
export function mockCommands(shell: ShellModel): CommandManager<Ctx> {
  const manager = new CommandManager<Ctx>({ context: observableState<Ctx>({}).asReadonly, formats });
  const open = (overlay: string, focus?: string) => () => {
    shell.openOverlay(overlay, undefined, focus);
    return done(`Opened ${overlay}`);
  };
  const show = (tab: string) => () => {
    shell.show(tab);
    return done(`Showing ${tab}`);
  };
  const commands: CommandDefinition<Ctx>[] = [
    { id: 'camera.reset', title: 'Reset camera', category: 'View', icon: 'focus', run: () => done('Camera reset') },
    { id: 'camera.reset-zoom', title: 'Reset zoom', category: 'View', icon: 'zoom-out', run: () => done('Zoom reset') },
    {
      id: 'camera.orient-axes',
      title: 'Orient axes',
      category: 'View',
      icon: 'axis-3d',
      run: () => done('Axes oriented'),
    },
    { id: 'view.expand', title: 'Expand viewport', category: 'View', icon: 'maximize-2', run: () => done('Expanded') },
    {
      id: 'view.illumination',
      title: 'Toggle illumination',
      category: 'View',
      icon: 'lightbulb',
      shortcut: { key: 'i', alt: true },
      run: () => done('Illumination toggled'),
    },
    {
      id: 'view.xr',
      title: 'Enter XR',
      category: 'View',
      icon: 'glasses',
      when: () => 'XR is not available in this browser',
      run: () => done(''),
    },
    {
      id: 'view.sequence',
      title: 'Show or hide sequence',
      category: 'View',
      icon: 'sequence',
      run: () => (shell.toggle('top'), done('Sequence toggled')),
    },
    { id: 'panel.structure', title: 'Show Structure', category: 'Panels', icon: 'layers', run: show('structure') },
    { id: 'panel.home', title: 'Show Home', category: 'Panels', icon: 'molstar', run: show('home') },
    { id: 'panel.log', title: 'Show Log', category: 'Panels', icon: 'scroll-text', run: show('log') },
    { id: 'panel.settings', title: 'Open settings', category: 'Panels', icon: 'settings', run: show('settings') },
    {
      id: 'canvas.settings',
      title: 'Canvas settings…',
      category: 'View',
      icon: 'aperture',
      run: open('canvas-settings'),
    },
    {
      id: 'canvas.background',
      title: 'Canvas background…',
      category: 'View',
      icon: 'palette',
      run: open('canvas-settings', 'background'),
    },
    { id: 'animation.open', title: 'Animation…', category: 'Animation', icon: 'clapperboard', run: open('animation') },
    {
      id: 'animation.wiggle',
      title: 'Wiggle by dynamics',
      category: 'Animation',
      icon: 'orbit',
      run: () => done('Wiggle: dynamics'),
    },
    { id: 'export.image', title: 'Export image…', category: 'Export', icon: 'image', run: open('export', 'image') },
    {
      id: 'export.animation',
      title: 'Export animation (MP4)…',
      category: 'Export',
      icon: 'film',
      run: open('export', 'animation'),
    },
    {
      id: 'export.geometry',
      title: 'Export geometry…',
      category: 'Export',
      icon: 'box',
      run: open('export', 'geometry'),
    },
    {
      id: 'state.save',
      title: 'Save state…',
      category: 'Snapshots',
      icon: 'save',
      shortcut: { key: 's', mod: true },
      run: open('save-state'),
    },
    {
      id: 'snapshot.next',
      title: 'Next snapshot',
      category: 'Snapshots',
      icon: 'chevron-right',
      run: () => done('Next snapshot'),
    },
    {
      id: 'snapshot.previous',
      title: 'Previous snapshot',
      category: 'Snapshots',
      icon: 'chevron-left',
      run: () => done('Previous snapshot'),
    },
    {
      id: 'data.download-structure',
      title: 'Download structure…',
      category: 'Data',
      icon: 'download',
      alias: 'download',
      params: () => ({
        id: PD.Text('', { label: 'ID' }),
        source: PD.Select(
          'pdb',
          [
            ['pdb', 'PDB'],
            ['alphafold', 'AlphaFold DB'],
            ['pdb-ihm', 'PDB-IHM'],
          ],
          { label: 'source' },
        ),
      }),
      fields: { source: { before: 'from' } },
      defaults: () => ({ source: 'pdb' }),
      action: 'Download',
      presets: () => [
        { id: 'alphafold', title: 'from AlphaFold DB', values: { source: 'alphafold' }, keywords: ['afdb', 'uniprot'] },
        { id: 'pdb-ihm', title: 'from PDB-IHM', values: { source: 'pdb-ihm' }, keywords: ['integrative'] },
      ],
      run: (a) => done(`Downloading ${a.id} from ${a.source}`),
    },
    {
      id: 'structure.preset',
      title: 'Apply preset…',
      category: 'Structure',
      icon: 'sparkles',
      alias: 'preset',
      params: () => ({
        preset: PD.Select(
          'default',
          [
            ['default', 'Default'],
            ['cartoon', 'Cartoon'],
            ['spacefill', 'Spacefill'],
            ['surface', 'Surface'],
            ['illustrative', 'Illustrative'],
          ],
          { label: 'preset' },
        ),
      }),
      action: 'Apply',
      run: (a) => done(`Preset: ${a.preset}`),
      defaults: () => ({ preset: 'default' }),
    },
    {
      id: 'selection.color',
      title: 'Color selection…',
      category: 'Selection',
      icon: 'palette',
      alias: ['color', 'colour'],
      params: () => ({
        color: PD.Text('', { label: 'color' }),
        on: PD.Select(
          'all',
          [
            ['all', 'all representations'],
            ['cartoon', 'cartoon'],
            ['ball-and-stick', 'ball & stick'],
          ],
          { label: 'representations' },
        ),
      }),
      fields: { color: { format: 'color' }, on: { before: 'on' } },
      defaults: () => ({ on: 'all' }),
      action: 'Color',
      run: (a) => done(`Colored the selection ${a.color} on ${a.on}`),
    },
    {
      id: 'selection.transparency',
      title: 'Set transparency…',
      category: 'Selection',
      icon: 'paintbrush',
      alias: 'transparency',
      params: () => ({ value: PD.Numeric(0.5, { min: 0, max: 1, step: 0.05 }, { label: 'value' }) }),
      fields: { value: { before: 'to' } },
      defaults: () => ({ value: 0.5 }),
      action: 'Apply',
      run: (a) => done(`Transparency ${a.value} on the selection`),
    },
    {
      id: 'focus.target',
      title: 'Focus on…',
      category: 'Focus',
      icon: 'locate',
      alias: 'focus',
      params: () => ({ target: PD.Text('', { label: 'ligand' }) }),
      fields: { target: { format: 'ligand' } },
      action: 'Focus',
      run: (a) => done(`Focused on ${a.target}`),
    },
    ...Object.entries(Languages).map(
      ([id, language]): CommandDefinition<Ctx> => ({
        id: `select.${id}`,
        title: `Select with ${language.title}…`,
        category: 'Select',
        icon: 'text-cursor-input',
        alias: id,
        params: () => ({
          expression: PD.Text('', { label: `${language.title} expression` }),
          mode: SelectMode,
        }),
        fields: { expression: { format: `expression.${id}` }, mode: { before: 'and' } },
        defaults: () => ({ mode: 'add' }),
        action: 'Select',
        run: (a) => done(`${modeLabel(a.mode as string)}: ${language.title} "${a.expression}"`),
      }),
    ),
    {
      id: 'select.residues',
      title: 'Select residues…',
      category: 'Select',
      icon: 'text-cursor-input',
      alias: 'residues',
      params: () => ({ residues: PD.Text('', { label: 'residues' }), mode: SelectMode }),
      fields: { residues: { format: 'residues' }, mode: { before: 'and' } },
      defaults: () => ({ mode: 'add' }),
      action: 'Select',
      run: (a) => done(`${modeLabel(a.mode as string)}: residues ${a.residues}`),
    },
    {
      id: 'select.ligands',
      title: 'Select ligands',
      category: 'Select',
      icon: 'mouse-pointer-click',
      run: () => done('Selected 3 ligands'),
    },
    {
      id: 'select.helices',
      title: 'Select helices',
      category: 'Select',
      icon: 'mouse-pointer-click',
      run: () => done('Selected 412 residues'),
    },
    {
      id: 'select.water',
      title: 'Select water',
      category: 'Select',
      icon: 'mouse-pointer-click',
      run: () => done('Selected 186 waters'),
    },
    {
      id: 'selection.clear',
      title: 'Clear selection',
      category: 'Selection',
      icon: 'x',
      run: () => done('Selection cleared'),
    },
    {
      id: 'measure.distance',
      title: 'Measure distance',
      category: 'Measure',
      icon: 'ruler',
      when: () => 'Select two atoms first',
      run: () => done(''),
    },
  ];
  // Categories in the order people reach for them (the command menu lists them so); the sort keeps each one's own order.
  const order = [
    'Data',
    'Structure',
    'Select',
    'Selection',
    'Focus',
    'Measure',
    'View',
    'Animation',
    'Snapshots',
    'Export',
    'Panels',
  ];
  commands.sort((a, b) => order.indexOf(a.category) - order.indexOf(b.category));
  manager.register(commands);
  return manager;
}
