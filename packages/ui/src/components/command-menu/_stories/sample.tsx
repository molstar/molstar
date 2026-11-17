/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A sample command set over a fake context (no viewer), for the command menu's stories: a structure of chains and
// residues, a few looks, a Select group with a shared mode, formats with completions, an own form and a field editor.

import { State } from '@molstar/core/reactive';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CommandManager } from '../../../commands/manager.js';
import type { CommandDefinition, CommandFormat, CommandMenuPrefix } from '../../../commands/types.js';
import { setAttr, watch } from '../../../dom/index.js';
import type { CommandFieldView, CommandFormView, CommandMenuAssistant } from '../command-menu.js';

export interface SampleContext {
  loaded: boolean;
  chains: readonly string[];
  selection: string;
}

const Colors: readonly (readonly [name: string, hex: string])[] = [
  ['red', '#d1495b'],
  ['orange', '#ed8a3f'],
  ['yellow', '#edae49'],
  ['green', '#3a9e6b'],
  ['teal', '#2a9d8f'],
  ['blue', '#3f7cc0'],
  ['purple', '#7b5ea7'],
  ['grey', '#8a8f98'],
];

/** "A 10-20", "B", "ligands": what a selection reads. */
const selection: CommandFormat<SampleContext> = {
  hint: (ctx) => `A chain (${ctx.chains.join(', ')}), a range (A 10-20), or ligands, water`,
  check: (text, ctx) => {
    const t = text.trim();
    if (/^(ligands|water|all)$/i.test(t)) return { ok: true, value: t.toLowerCase() };
    const m = /^([A-Za-z])(?:\s+(\d+)(?:-(\d+))?)?$/.exec(t);
    if (m && ctx.chains.includes(m[1].toUpperCase())) return { ok: true, value: t.toUpperCase() };
    return { ok: false, message: `“${t}” is not a selection.`, hint: 'Try “A 10-20”.', unreadable: true };
  },
  complete: (text, ctx) => {
    const all = [
      ...ctx.chains.map((c) => ({ value: c, label: `Chain ${c}`, detail: 'chain' })),
      { value: 'ligands', label: 'Ligands', detail: '3 residues' },
      { value: 'water', label: 'Water', detail: '186 residues' },
      { value: 'A 10-20', label: 'A 10-20', detail: 'an example to edit', fill: true },
    ];
    const q = text.trim().toLowerCase();
    return all.filter((x) => !q || x.label.toLowerCase().includes(q) || x.value.toLowerCase().startsWith(q));
  },
};

const color: CommandFormat<SampleContext> = {
  hint: () => 'A color name or #rrggbb',
  check: (text) => {
    const t = text.trim().toLowerCase();
    if (Colors.some(([n]) => n === t) || /^#[0-9a-f]{6}$/.test(t)) return { ok: true, value: t };
    return { ok: false, message: `“${text}” is not a color.` };
  },
  complete: (text) =>
    Colors.filter(([n]) => n.startsWith(text.trim().toLowerCase())).map(([n, hex]) => ({
      value: n,
      label: n,
      detail: hex,
    })),
};

export const formats = {
  get: (name: string) => ({ selection, color })[name] as CommandFormat<SampleContext> | undefined,
};

function commands(log: (s: string) => void): CommandDefinition<SampleContext>[] {
  const done = (message: string) => {
    log(message);
    return { message };
  };
  const loaded = (ctx: SampleContext) => (ctx.loaded ? true : 'Load a structure first.');
  return [
    {
      id: 'file.open',
      title: 'Open a file…',
      category: 'File',
      icon: 'folder-open',
      shortcut: { key: 'o', mod: true },
      run: () => done('Opened.'),
    },
    {
      id: 'file.download',
      title: 'Download structure…',
      category: 'File',
      icon: 'download',
      alias: 'fetch',
      params: () => ({
        id: PD.Text('', { label: 'PDB id' }),
        source: PD.Select('pdbe', [
          ['pdbe', 'PDBe'],
          ['rcsb', 'RCSB PDB'],
          ['afdb', 'AlphaFold DB'],
        ]),
        assembly: PD.Optional(PD.Text('1')) as PD.Any,
      }),
      fields: { source: { before: 'from' }, assembly: { before: 'assembly' } },
      presets: () => [{ id: 'afdb', title: 'from AlphaFold DB', values: { source: 'afdb' }, keywords: ['alphafold'] }],
      action: 'Download',
      run: (args) => done(`Downloaded ${args.id} from ${args.source}.`),
    },
    {
      id: 'view.reset',
      title: 'Reset the camera',
      category: 'View',
      icon: 'rotate-ccw',
      shortcut: { key: 'r' },
      run: () => done('Camera reset.'),
    },
    {
      id: 'view.focus',
      title: 'Focus on…',
      category: 'View',
      icon: 'focus',
      shortcut: { key: 'f' },
      when: loaded,
      params: () => ({ target: PD.Text('') }),
      fields: { target: { format: 'selection' } },
      action: 'Focus',
      run: (args) => done(`Focused on ${args.target}.`),
    },
    {
      id: 'look.show',
      title: 'Show as…',
      category: 'Look',
      icon: 'layers',
      alias: 'show',
      when: loaded,
      params: () => ({
        what: PD.Text(''),
        type: PD.Select('cartoon', PD.arrayToOptions(['cartoon', 'surface', 'ball-and-stick', 'spacefill'])),
        size: PD.Optional(PD.Numeric(1, { min: 0.1, max: 4, step: 0.1 })) as PD.Any,
        quality: PD.Select('auto', PD.arrayToOptions(['auto', 'high', 'medium', 'low'])),
        opacity: PD.Numeric(1, { min: 0, max: 1, step: 0.05 }, { isHidden: false }),
        ignoreHydrogens: PD.Boolean(false, { isHidden: true }),
        material: PD.Group(
          {
            metalness: PD.Numeric(0, { min: 0, max: 1, step: 0.05 }),
            roughness: PD.Numeric(1, { min: 0, max: 1, step: 0.05 }),
          },
          { isFlat: true },
        ),
        tint: PD.Color(0x3f7cc0 as any),
      }),
      fields: {
        what: { format: 'selection' },
        type: { before: 'as' },
        size: { before: 'size', presets: [0.5, 1, [2, '2 (large)']] },
        quality: { before: 'quality' },
        opacity: {
          before: 'opacity',
          presets: [
            [0.4, '40%'],
            [0.6, '60%'],
            [1, 'opaque'],
          ],
        },
      },
      action: 'Show',
      run: (args) => done(`Shown ${args.what} as ${args.type}.`),
    },
    {
      id: 'look.color',
      title: 'Color…',
      category: 'Look',
      icon: 'palette',
      alias: 'color',
      when: loaded,
      params: () => ({ what: PD.Text(''), color: PD.Text('') }),
      fields: { what: { format: 'selection' }, color: { before: 'in', format: 'color' } },
      assist: { instruction: 'What to color and with which color.' },
      action: 'Color',
      run: (args) => done(`Colored ${args.what} ${args.color}.`),
    },
    {
      id: 'look.failing',
      title: 'Export image…',
      category: 'File',
      icon: 'camera',
      params: () => ({ width: PD.Numeric(1920, { min: 64, max: 8192, step: 1 }) }),
      fields: { width: { presets: [1280, 1920, 3840] } },
      view: { title: 'Export › Image', open: () => log('Opened Export › Image.') },
      action: 'Export',
      run: () => {
        throw new Error('The canvas is too small to render at that width.');
      },
    },
    {
      id: 'edit.remove',
      title: 'Remove everything',
      category: 'Edit',
      icon: 'trash-2',
      when: (ctx) => (ctx.loaded ? true : 'Nothing is loaded.'),
      run: () => done('Removed.'),
    },
    ...['Ligands', 'Water', 'Polymer', 'Chain A', 'Chain B'].map(
      (title, i): CommandDefinition<SampleContext> => ({
        id: `select.${title.toLowerCase().replace(' ', '-')}`,
        title,
        category: i < 3 ? 'Components' : 'Chains',
        icon: 'crosshair',
        group: 'select',
        preview: () => title,
        count: async () => `${[3, 186, 1204, 640, 564][i]} residues`,
        run: (args) => done(`Selected ${title} (${args.mode}).`),
      }),
    ),
  ];
}

/** A manager over the sample context, the log of what ran. */
export function sampleManager(options: { loaded?: boolean } = {}) {
  const log = new State<readonly string[]>([]);
  const context = new State<SampleContext>({ loaded: options.loaded ?? true, chains: ['A', 'B'], selection: '' });
  const manager = new CommandManager<SampleContext>({
    context: context.asReadonly,
    formats,
    preview: { set: () => {}, clear: () => {} },
  });
  manager.registerGroups([
    {
      id: 'select',
      title: 'Select',
      icon: 'crosshair',
      shared: {
        name: 'mode',
        title: 'Mode',
        options: [
          ['set', 'Replace'],
          ['add', 'Add'],
          ['remove', 'Subtract'],
        ],
        default: 'set',
      },
    },
  ]);
  manager.register(commands((s) => log.set([s, ...log.value].slice(0, 4))));
  return { manager, context, log };
}

export const prefixes: readonly CommandMenuPrefix[] = [
  {
    char: '@',
    title: 'Scene',
    placeholder: 'A node of the scene',
    items: (text) =>
      ['1cbs › Polymer', '1cbs › Ligands', '1cbs › Water', 'Measurements']
        .filter((x) => x.toLowerCase().includes(text.toLowerCase()))
        .map((label) => ({ id: label, label, icon: 'box' as const, detail: 'scene', run: () => {} })),
  },
  {
    char: '#',
    title: 'Snapshots',
    placeholder: 'A saved snapshot',
    items: () => [{ id: 's1', label: 'Overview', detail: 'snapshot 1', icon: 'camera', run: () => {} }],
  },
];

/** A fake assistant: words naming a color propose Color…; others are "asked". */
export function fakeAssistant(ready = true): CommandMenuAssistant {
  return {
    ready: () => ready,
    propose: (text) => {
      const c = Colors.find(([n]) => text.toLowerCase().includes(n));
      return c ? { command: 'look.color', tail: `ligands in ${c[0]}` } : undefined;
    },
    ask: () => {},
    setup: () => {},
  };
}

/** A swatch grid for every field of the `color` format: the arrows move, a click takes. */
export const colorSwatches: CommandFieldView = {
  render: (p) => {
    const buttons: HTMLElement[] = [];
    const current = () => String(p.invocation.state.value.values[p.field] ?? '');
    const el = (
      <div
        role="group"
        aria-label="Colors"
        style="display: grid; grid-template-columns: repeat(8, 24px); gap: var(--msui-space-1); padding-inline-start: calc(var(--msui-icon-slot) + var(--msui-icon-gap));"
      >
        {Colors.map(([name, hex]) => (
          <button
            type="button"
            tabIndex={-1}
            data-tip={name}
            aria-label={name}
            ref={(b: HTMLElement) => buttons.push(b)}
            style={`width: 24px; height: 24px; border-radius: var(--msui-radius-md); background: ${hex};`}
            onClick={() => p.set(name)}
          />
        ))}
      </div>
    );
    watch([p.invocation.state], () => {
      const v = current();
      for (let i = 0; i < buttons.length; i++)
        setAttr(
          buttons[i],
          'style',
          `width: 24px; height: 24px; border-radius: var(--msui-radius-md); background: ${Colors[i][1]};${
            Colors[i][0] === v ? ' box-shadow: 0 0 0 2px var(--msui-surface-raised), 0 0 0 3px var(--msui-text);' : ''
          }`,
        );
    });
    p.keys((e) => {
      if (e.key !== 'ArrowRight' && e.key !== 'ArrowLeft') return false;
      if (p.invocation.state.value.text) return false;
      const i = Colors.findIndex(([n]) => n === current());
      const next = Colors[(i + (e.key === 'ArrowRight' ? 1 : Colors.length - 1) + (i < 0 ? 1 : 0)) % Colors.length];
      p.invocation.set(p.field, next[0], { stay: true });
      return true;
    });
    return el;
  },
};

/** Show as…'s own form: the representations as tiles (its `type` field), the size under them. */
export const showAsForm: CommandFormView = {
  fields: ['type'],
  render: (p) => {
    const types = ['cartoon', 'surface', 'ball-and-stick', 'spacefill'];
    const tiles: HTMLElement[] = [];
    const el = (
      <div style="display: grid; grid-template-columns: repeat(4, 1fr); gap: var(--msui-space-2); padding: var(--msui-space-2) var(--msui-inset);">
        {types.map((t) => (
          <button
            type="button"
            tabIndex={-1}
            ref={(b: HTMLElement) => tiles.push(b)}
            onClick={() => p.invocation.set('type', t)}
            style="display: grid; place-items: center; height: calc(14 * var(--msui-unit)); border-radius: var(--msui-radius-md); background: var(--msui-control); color: var(--msui-text-muted);"
          >
            {t}
          </button>
        ))}
      </div>
    );
    watch([p.invocation.state], () => {
      const s = p.invocation.state.value;
      for (let i = 0; i < tiles.length; i++) {
        const on = s.values.type === types[i];
        tiles[i].style.boxShadow = on ? 'inset 0 0 0 1px var(--msui-accent)' : '';
        tiles[i].style.color = on ? 'var(--msui-accent)' : 'var(--msui-text-muted)';
      }
    });
    return el;
  },
};
