/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CommandManager } from '@molstar/prism-ui/commands/manager';
import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import { rankBy } from '@molstar/prism-ui/commands/match';
import type { CommandDefinition, CommandFormats, CommandResult } from '@molstar/prism-ui/commands/types';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { SnapshotStepperModel } from '@molstar/prism-ui/model/playback';
import type { SceneModel, SceneNode } from '../model/scene.js';
import type { ShellModel } from '@molstar/prism-ui/shell/shell';
import type { CommandMenuPrefix } from '@molstar/prism-ui';

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
export function mockCommands(
  shell: ShellModel,
  /** How many structures are loaded (superposing needs two). */
  structureCount: () => number = () => 1,
): CommandManager<Ctx> {
  const manager = new CommandManager<Ctx>({ context: observableState<Ctx>({}).asReadonly, formats });
  const open = (overlay: string, focus?: string) => () => {
    shell.openOverlay(overlay, undefined, focus);
    return done(`Opened ${overlay}`);
  };
  const show = (tab: string, section?: string) => () => {
    shell.show(tab, section);
    return done(`Showing ${section ?? tab}`);
  };
  const commands: CommandDefinition<Ctx>[] = [
    {
      id: 'camera.reset',
      group: 'camera',
      title: 'Reset camera',
      alias: ['reset view', 'center', 'home view'],
      category: 'View',
      icon: 'focus',
      run: () => done('Camera reset'),
    },
    {
      id: 'camera.reset-zoom',
      group: 'camera',
      title: 'Reset zoom',
      alias: ['zoom to fit', 'fit'],
      category: 'View',
      icon: 'zoom-out',
      run: () => done('Zoom reset'),
    },
    {
      id: 'camera.orient-axes',
      group: 'camera',
      title: 'Orient axes',
      alias: ['orient', 'align axes', 'principal axes'],
      category: 'View',
      icon: 'axis-3d',
      run: () => done('Axes oriented'),
    },
    { id: 'view.expand', title: 'Expand viewport', category: 'View', icon: 'maximize-2', run: () => done('Expanded') },
    {
      id: 'view.illumination',
      title: 'Toggle illumination',
      alias: ['ray tracing', 'path tracing', 'global illumination', 'lighting'],
      category: 'View',
      icon: 'lightbulb',
      shortcut: { key: 'i', alt: true },
      run: () => done('Illumination toggled'),
    },
    {
      id: 'view.xr',
      title: 'Enter XR',
      alias: ['vr', 'ar', 'webxr', 'headset'],
      category: 'View',
      icon: 'glasses',
      when: () => 'XR is not available in this browser',
      run: () => done(''),
    },
    {
      id: 'view.sequence',
      title: 'Show or hide sequence',
      alias: ['sequence', 'sequence viewer'],
      category: 'View',
      icon: 'sequence',
      run: () => (shell.toggle('top'), done('Sequence toggled')),
    },
    { id: 'panel.structure', title: 'Show Structure', category: 'Panels', icon: 'structure', run: show('structure') },
    { id: 'panel.home', title: 'Show Home', category: 'Panels', icon: 'prism', run: show('home') },
    { id: 'panel.log', title: 'Show Log', category: 'Panels', icon: 'scroll-text', run: show('log') },
    { id: 'panel.settings', title: 'Open settings', category: 'Panels', icon: 'settings', run: show('settings') },
    {
      id: 'canvas.settings',
      title: 'Canvas settings…',
      alias: ['viewport settings', 'rendering', 'quality'],
      category: 'View',
      icon: 'monitor-cog',
      run: open('canvas-settings'),
    },
    {
      id: 'canvas.background',
      title: 'Canvas background…',
      alias: ['background', 'bg color'],
      category: 'View',
      icon: 'palette',
      run: open('canvas-settings', 'background'),
    },
    { id: 'animation.open', title: 'Animation…', category: 'Animation', icon: 'clapperboard', run: open('animation') },
    {
      id: 'animation.wiggle',
      title: 'Wiggle by dynamics',
      alias: ['wiggle', 'b-factor', 'dynamics'],
      category: 'Animation',
      icon: 'orbit',
      run: () => done('Wiggle: dynamics'),
    },
    {
      id: 'export.image',
      group: 'export',
      title: 'Export image…',
      alias: ['screenshot', 'save image', 'png', 'picture'],
      category: 'Export',
      icon: 'image',
      run: open('export', 'image'),
    },
    {
      id: 'export.animation',
      group: 'export',
      title: 'Export animation…',
      alias: ['video', 'movie', 'mp4', 'record'],
      category: 'Export',
      icon: 'film',
      run: open('export', 'animation'),
    },
    {
      id: 'export.geometry',
      group: 'export',
      title: 'Export geometry…',
      alias: ['glb', 'stl', 'obj', 'usdz', '3d print', 'mesh'],
      category: 'Export',
      icon: 'box',
      run: open('export', 'geometry'),
    },
    {
      id: 'snapshots.show',
      title: 'Snapshots',
      alias: ['states', 'scenes', 'save state', 'bookmarks'],
      category: 'Snapshots',
      icon: 'snapshots',
      shortcut: { key: 's', mod: true },
      run: show('snapshots'),
    },
    {
      id: 'snapshot.next',
      title: 'Next snapshot',
      alias: ['next scene', 'forward'],
      category: 'Snapshots',
      icon: 'chevron-right',
      run: () => done('Next snapshot'),
    },
    {
      id: 'snapshot.previous',
      title: 'Previous snapshot',
      alias: ['previous scene', 'prev', 'back'],
      category: 'Snapshots',
      icon: 'chevron-left',
      run: () => done('Previous snapshot'),
    },
    {
      id: 'data.download-structure',
      title: 'Download structure…',
      category: 'Data',
      icon: 'download',
      alias: ['download', 'fetch', 'load', 'open pdb'],
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
      id: 'selection.color',
      group: 'selection',
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
      group: 'selection',
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
        group: 'select.expression',
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
      group: 'select',
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
      group: 'select.query',
      title: 'Select ligands',
      alias: ['ligands', 'het', 'hetero'],
      category: 'Select',
      icon: 'mouse-pointer-click',
      run: () => done('Selected 3 ligands'),
    },
    {
      id: 'select.helices',
      group: 'select.query',
      title: 'Select helices',
      alias: ['helix', 'alpha helices'],
      category: 'Select',
      icon: 'mouse-pointer-click',
      run: () => done('Selected 412 residues'),
    },
    {
      id: 'select.water',
      group: 'select.query',
      title: 'Select water',
      alias: ['waters', 'solvent', 'hoh'],
      category: 'Select',
      icon: 'mouse-pointer-click',
      run: () => done('Selected 186 waters'),
    },
    {
      id: 'selection.clear',
      group: 'selection',
      title: 'Clear selection',
      alias: ['deselect', 'select none', 'unselect'],
      category: 'Selection',
      icon: 'x',
      run: () => done('Selection cleared'),
    },
    // Done with the selection: open the Selection tab at their section (only it open).
    {
      id: 'measure.open',
      title: 'Measure…',
      alias: ['measure', 'distance', 'angle', 'dihedral', 'ruler'],
      category: 'Measure',
      icon: 'ruler',
      run: show('selection', 'measure'),
    },
    {
      id: 'superpose.open',
      title: 'Superpose…',
      alias: ['align', 'superimpose', 'rmsd', 'tm-align', 'fit structures'],
      category: 'Structure',
      icon: 'superpose',
      when: () => (structureCount() > 1 ? true : 'Load two or more structures'),
      run: show('selection', 'superpose'),
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
  // Representation presets and quick styles: one command each, in their groups (the Structure tab's Preset and Styles
  // buttons open these groups in the command menu); a search anywhere finds them too.
  const preset = (id: string, title: string, category: string): CommandDefinition<Ctx> => ({
    id: `preset.${id}`,
    title,
    category,
    group: 'preset',
    icon: 'layers',
    run: () => done(`Preset: ${title}`),
  });
  commands.push(
    preset('empty', 'Empty', 'Presets'),
    preset('auto', 'Automatic', 'Presets'),
    preset('atomic-detail', 'Atomic detail', 'Basic'),
    preset('polymer-cartoon', 'Polymer cartoon', 'Basic'),
    preset('polymer-and-ligand', 'Polymer & ligand', 'Basic'),
    preset('protein-and-nucleic', 'Protein & nucleic', 'Basic'),
    preset('coarse-surface', 'Coarse surface', 'Basic'),
    preset('illustrative', 'Illustrative', 'Miscellaneous'),
    preset('molecular-surface', 'Molecular surface', 'Miscellaneous'),
    preset('auto-detail', 'Automatic detail', 'Miscellaneous'),
  );
  const style = (id: string, title: string, category: string, icon: IconName): CommandDefinition<Ctx> => ({
    id: `styles.${id}`,
    title,
    category,
    group: 'styles',
    icon,
    run: () => done(`${category}: ${title}`),
  });
  commands.push(
    style('repr-default', 'Default representation', 'Representation', 'layers'),
    style('cartoon', 'Cartoon', 'Representation', 'cartoon'),
    style('spacefill', 'Spacefill', 'Representation', 'spacefill'),
    style('surface', 'Surface', 'Representation', 'surface'),
    style('style-default', 'Default style', 'Style', 'paintbrush'),
    style('illustrative', 'Illustrative style', 'Style', 'paintbrush'),
  );

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
  // The command tree for browsing (search finds everything either way).
  manager.registerGroups([
    { id: 'preset', title: 'Preset', icon: 'shapes', order: -2 },
    { id: 'styles', title: 'Styles', icon: 'paintbrush', order: -1 },
    { id: 'select', title: 'Select', icon: 'mouse-pointer-click', order: 0 },
    { id: 'select.query', title: 'By query', icon: 'list-tree', parent: 'select', order: 0 },
    { id: 'select.expression', title: 'Expression', icon: 'text-cursor-input', parent: 'select', order: 1 },
    { id: 'selection', title: 'Selection', icon: 'palette', order: 1 },
    { id: 'camera', title: 'Camera', icon: 'focus', order: 2 },
    { id: 'export', title: 'Export', icon: 'download', order: 3 },
  ]);
  return manager;
}

/**
 * The command menu's prefixes in the stories: `@` goes to a scene node (picks it, shows the Structure tab), `#` to a
 * snapshot, `:` selects by residues or an expression (`:A 10-20`, `:pymol resn HEM`), opening that command's form.
 */
export function mockPrefixes(options: {
  manager: CommandManager<Ctx>;
  shell: ShellModel;
  scene: SceneModel;
  snapshots: SnapshotStepperModel;
}): CommandMenuPrefix[] {
  const { manager, shell, scene, snapshots } = options;
  const nodes = (list: readonly SceneNode[], path: string[] = []): { node: SceneNode; path: string[] }[] =>
    list.flatMap((n) => [{ node: n, path }, ...nodes(n.children.value, [...path, n.label.value])]);
  return [
    {
      char: '@',
      title: 'Go to',
      placeholder: 'A structure, component, volume…',
      items: (text) =>
        rankBy(nodes(scene.roots.value), text, (x) => `${x.node.label.value} ${x.path.join(' ')}`).map((x) => ({
          id: x.node.id,
          label: x.node.label.value,
          detail: x.path.join(' › ') || x.node.kind,
          icon: x.node.icon,
          run: () => {
            shell.show('structure');
            if (scene.picked.value !== x.node) scene.pick(x.node);
          },
        })),
    },
    {
      char: '#',
      title: 'Snapshots',
      placeholder: 'A snapshot by title',
      items: (text) =>
        rankBy(
          snapshots.titles.value.map((title, index) => ({ title, index })),
          text,
          (x) => x.title,
        ).map((x) => ({
          id: `snapshot:${x.index}`,
          label: x.title,
          detail: `${x.index + 1} / ${snapshots.count.value}`,
          icon: 'image',
          run: () => snapshots.goTo(x.index),
        })),
    },
    {
      char: ':',
      title: 'Select',
      placeholder: 'Residues (A 10-20) or an expression (pymol resn HEM)',
      items: (text) => {
        const language = /^(pymol|vmd|jmol)\b\s*(.*)$/i.exec(text);
        const open = (id: string, tail: string) => () => {
          manager.open(id, tail || undefined);
          return true;
        };
        if (language) {
          const id = language[1].toLowerCase();
          return [
            {
              id: `select.${id}`,
              label: `Select with ${Languages[id].title}${language[2] ? `: ${language[2]}` : '…'}`,
              icon: 'text-cursor-input',
              run: open(`select.${id}`, language[2]),
            },
          ];
        }
        return [
          {
            id: 'select.residues',
            label: text ? `Select residues ${text}` : 'Select residues…',
            detail: 'A 10-20, B 5',
            icon: 'text-cursor-input',
            run: open('select.residues', text),
          },
          ...Object.entries(Languages).map(([id, l]) => ({
            id: `select.${id}`,
            label: `Select with ${l.title}…`,
            detail: `:${id} …`,
            icon: 'text-cursor-input' as const,
            run: open(`select.${id}`, ''),
          })),
        ];
      },
    },
  ];
}
