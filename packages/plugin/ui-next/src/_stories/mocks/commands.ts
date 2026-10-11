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

/** Named argument types the stories' commands use: colors (names or hex) and ligands of 1TQN. */
const formats: CommandFormats<Ctx> = {
  get: (name) => {
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
    { id: 'view.expand', title: 'Expand viewport', category: 'View', icon: 'maximize', run: () => done('Expanded') },
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
  manager.register(commands);
  return manager;
}
