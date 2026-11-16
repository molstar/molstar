/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Material } from '@molstar/core/util/material';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { derive } from '@molstar/core/reactive';
import { observeResize } from '../../../dom/index.js';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../_stories/support.js';
import { ParamsFormModel } from '../../../forms/params-form.js';
import { state } from '../../../state.js';
import { TextInput } from '../../text-input/text-input.js';
import { installTooltips } from '../../tooltip/tooltip.js';
import { ParamsForm } from '../params-form.js';

/**
 * Every PD control, alone, in each state: on a surface and in a well, at a panel's width, with the values it produces
 * below. For checking spacing and alignment in isolation (styles §5): labels at the inset, fields from the label column,
 * their text at the control inset, every row 24px.
 */
const meta: Meta = { title: 'Components/Params form' };
export default meta;

const surfaces =
  'display: grid; grid-template-columns: repeat(auto-fill, 320px); gap: var(--msui-space-6); padding: var(--msui-space-5);';
const caption = 'margin: 0 0 var(--msui-space-2); color: var(--msui-text-subtle); font-size: var(--msui-text-xs);';
const values =
  'margin: var(--msui-space-3) 0 0; padding: var(--msui-space-2); border-radius: var(--msui-radius-md); background: var(--msui-surface); color: var(--msui-text-muted); font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); white-space: pre-wrap; word-break: break-all;';

/** The values as JSON, files by name. */
function show(v: unknown) {
  return JSON.stringify(
    v,
    (_, x) => (x && typeof x === 'object' && 'kind' in x && x.kind === 'file' ? `<file ${x.name}>` : x),
    1,
  );
}

/**
 * A form on the two grounds it lives on: a surface (a panel, an inspector) and a well (a section sunk into a panel);
 * `width` narrower to check long labels (they wrap; values keep their column).
 */
function Surfaces(props: { params: PD.Params; width?: number; filter?: boolean }) {
  const width = props.width ?? 320;
  const block = (title: string, ground: string) => {
    const form = new ParamsFormModel(props.params);
    const query = state('');
    return (
      <div style={`width: ${width}px`}>
        <p style={caption}>{title}</p>
        <div class="msui-surface" style={`background: ${ground}`} data-story-form="">
          {props.filter ? (
            // A full-width row (24px, the md field): its box from the surface's padding to the other side's, so its
            // text (at the control inset) lands on the label column, the inset; the rows follow on the 4px grid.
            <div style="display: grid; grid-template-columns: minmax(0, 1fr);">
              <TextInput search placeholder="Search" value={query} onInput={(v) => query.set(v)} />
            </div>
          ) : null}
          <ParamsForm form={form} filter={props.filter ? query : undefined} />
        </div>
        <pre style={values}>{derive(form.values, show)}</pre>
      </div>
    );
  };
  return (
    // The tooltips (a value's full text when it is cut off, an icon button's name), as an application installs them.
    <div style={surfaces} ref={(el: HTMLElement) => installTooltips(el)}>
      {block('On a surface', 'var(--msui-surface)')}
      {block('In a well', 'var(--msui-surface-sunken)')}
    </div>
  );
}

const render = (params: PD.Params, options: { width?: number; filter?: boolean } = {}): StoryObj => ({
  render: () => story(() => <Surfaces params={params} width={options.width} filter={options.filter} />),
});

/** A curve's points (a line graph's value). */
const curve = (points: readonly (readonly [number, number])[]) => points.map(([x, y]) => Vec2.create(x, y));

const onOff = (on: boolean, params: PD.Params = {}, info: PD.Info = {}) =>
  PD.MappedStatic(on ? 'on' : 'off', { on: PD.Group(params), off: PD.Group({}) }, info);

/** One of each type the renderer draws, as a panel would show them. */
export const AllTypes = render({
  boolean: PD.Boolean(true),
  number: PD.Numeric(5, { min: 0, max: 20, step: 0.5 }),
  integer: PD.Numeric(3),
  text: PD.Text('1tqn', { placeholder: 'PDB ID' }),
  select: PD.Select('cartoon', [
    ['cartoon', 'Cartoon'],
    ['ball-and-stick', 'Ball & Stick'],
    ['spacefill', 'Spacefill'],
  ]),
  color: PD.Color(0x3362b2 as any),
  file: PD.File({ accept: '.cif,.pdb' }),
  files: PD.FileList({ multiple: true }),
  outline: onOff(true, { scale: PD.Numeric(1, { min: 1, max: 5, step: 1 }) }, { label: 'Outline' }),
  mode: PD.MappedStatic(
    'auto',
    {
      auto: PD.Group({}, { isFlat: true }),
      custom: PD.Group({ width: PD.Numeric(1920), height: PD.Numeric(1080) }, { isFlat: true }),
    },
    {
      options: [
        ['auto', 'Auto'],
        ['custom', 'Custom'],
      ],
    },
  ),
  group: PD.Group(
    { samples: PD.Numeric(32, { min: 1, max: 256, step: 1 }), bias: PD.Numeric(0.8) },
    { label: 'Group' },
  ),
});

/** Switches: plain, described (hover the row), and on/off mapped params with and without their own settings. */
export const Booleans = render({
  visuals: PD.Boolean(true),
  transparent: PD.Boolean(false, { description: 'Render the background transparent.' }),
  occlusion: onOff(
    true,
    {
      samples: PD.Numeric(32, { min: 1, max: 256, step: 1 }),
      radius: PD.Numeric(5, { min: 0, max: 20, step: 0.1 }),
    },
    { label: 'Ambient occlusion' },
  ),
  shadows: onOff(false, { steps: PD.Numeric(1, { min: 1, max: 64, step: 1 }) }, { label: 'Shadows' }),
  handle: onOff(false, {}, { label: 'Handle (no settings)' }),
});

/** Numbers: a slider with a range (integer, fractional, negative), and a plain field without one. */
export const Numbers = render({
  samples: PD.Numeric(32, { min: 1, max: 256, step: 1 }),
  opacity: PD.Numeric(0.35, { min: 0, max: 1, step: 0.01 }),
  isoValue: PD.Numeric(-3, { min: -5, max: 5, step: 0.1 }, { label: 'Iso value (σ)' }),
  durationMs: PD.Numeric(1500, { min: 100, max: 15000, step: 100 }, { label: 'Duration (ms)' }),
  modelIndex: PD.Numeric(0, undefined, { label: 'Model index (no range)' }),
  width: PD.Numeric(1920, undefined, { label: 'Width (px)' }),
});

/** Text: a value, empty with a placeholder, a long value, multiline. */
export const Text = render({
  id: PD.Text('1tqn', { label: 'PDB ID' }),
  url: PD.Text('', { label: 'URL', placeholder: 'https://…' }),
  label: PD.Text('Cytochrome P450 3A4 with the heme cofactor and a long name', { label: 'Label' }),
  description: PD.Text('A note\non two lines', { multiline: true, label: 'Description (multiline)' }),
});

/** Selects: a few options, many (the list scrolls), long labels; opened, the options over the field. */
export const Selects = render({
  representation: PD.Select('cartoon', [
    ['cartoon', 'Cartoon'],
    ['ball-and-stick', 'Ball & Stick'],
    ['spacefill', 'Spacefill'],
    ['molecular-surface', 'Molecular Surface'],
  ]),
  colorTheme: PD.Select(
    'chain-id',
    [
      'atom-id',
      'carbohydrate-symbol',
      'chain-id',
      'element-index',
      'element-symbol',
      'entity-id',
      'entity-source',
      'hydrophobicity',
      'illumination',
      'model-index',
      'molecule-type',
      'occupancy',
      'operator-name',
      'partial-charge',
      'polymer-id',
      'residue-name',
      'secondary-structure',
      'sequence-id',
      'shape-group',
      'uncertainty',
      'uniform',
    ].map((v) => [v, v.replace(/-/g, ' ').replace(/^./, (c) => c.toUpperCase())] as [string, string]),
    { label: 'Color theme (many)' },
  ),
  server: PD.Select(
    'pdbe',
    [
      ['pdbe', 'PDBe · Protein Data Bank in Europe (EMBL-EBI)'],
      ['rcsb', 'RCSB PDB · Research Collaboratory for Structural Bioinformatics'],
    ],
    { label: 'Server (long labels)' },
  ),
});

/** Colors: the field; opened, the picker in a panel under it (the form does not move). */
export const Colors = render({
  color: PD.Color(0x3362b2 as any),
  background: PD.Color(0xffffff as any, { label: 'Background (white)' }),
  edge: PD.Color(0x000000 as any, { label: 'Edge (black)' }),
});

/** Files: one, many (drop files on the field, or pick); cleared with ✕. */
export const Files = render({
  file: PD.File({ label: 'Map file', accept: '.ccp4,.mrc,.map' }),
  files: PD.FileList({ label: 'Files', multiple: true }),
});

/** Groups: flat (more rows at the same level), folded, open by default, and a group in a group. */
export const Groups = render({
  flat: PD.Group({ id: PD.Text('1tqn', { label: 'PDB ID' }), assembly: PD.Text('1') }, { isFlat: true }),
  folded: PD.Group(
    { samples: PD.Numeric(32, { min: 1, max: 256, step: 1 }), bias: PD.Numeric(0.8) },
    { label: 'Folded' },
  ),
  expanded: PD.Group({ exposure: PD.Numeric(1, { min: 0, max: 3, step: 0.01 }) }, { label: 'Open', isExpanded: true }),
  nested: PD.Group(
    {
      inner: PD.Group({ depth: PD.Numeric(2, { min: 1, max: 5, step: 1 }) }, { label: 'Inner', isExpanded: true }),
      after: PD.Boolean(true),
    },
    { label: 'Group in a group', isExpanded: true },
  ),
});

/**
 * Presets (PD's): named values to start from. A group's in its header (Mol*'s Material Style: Matte, Plastic, Glossy,
 * Metallic), the matching one named on the value column and marked in the list (none: the values are custom); a flat
 * group's as a row before its fields; an object list's beside Add; a curve's, asked for when opened (from the data).
 */
export const Presets = render({
  material: Material.getParam(),
  materialOpen: { ...Material.getParam({ isExpanded: true }), label: 'Material (open)' },
  size: PD.Group(
    {
      width: PD.Numeric(1920, { min: 1, max: 8192, step: 1 }),
      height: PD.Numeric(1080, { min: 1, max: 8192, step: 1 }),
    },
    {
      isFlat: true,
      presets: [
        [{ width: 1920, height: 1080 }, 'Full HD'],
        [{ width: 3840, height: 2160 }, '4K'],
        [{ width: 1080, height: 1080 }, 'Square'],
      ],
    },
  ),
  channels: PD.ObjectList(
    { name: PD.Text('channel'), level: PD.Numeric(1.5, { min: 0, max: 5, step: 0.1 }) },
    (e) => `${e.name} at ${e.level}σ`,
    {
      defaultValue: [{ name: '2fo-fc', level: 1.5 }],
      presets: [
        [[{ name: '2fo-fc', level: 1.5 }], 'Electron density'],
        [
          [
            { name: 'fo-fc(+)', level: 3 },
            { name: 'fo-fc(-)', level: 3 },
          ],
          'Difference map',
        ],
      ],
    },
  ),
  transfer: PD.LineGraph(
    curve([
      [0, 0],
      [0.5, 0.2],
      [1, 1],
    ]),
    {
      label: 'Transfer curve',
      getPresets: () => [
        [
          curve([
            [0, 0],
            [1, 1],
          ]),
          'Linear',
        ],
        [
          curve([
            [0, 0],
            [0.5, 0.2],
            [1, 1],
          ]),
          'Soft',
        ],
        [
          curve([
            [0, 0],
            [0.4, 0],
            [0.6, 1],
            [1, 1],
          ]),
          'Step',
        ],
      ],
    },
  ),
});

/** Mapped: a choice with its own params under it, flat or nested, or a single value. */
export const Mapped = render({
  source: PD.MappedStatic(
    'pdb',
    {
      pdb: PD.Group({ id: PD.Text('1tqn', { label: 'PDB ID' }) }, { isFlat: true }),
      url: PD.Group({ url: PD.Text('', { label: 'URL' }), binary: PD.Boolean(false) }, { isFlat: true }),
    },
    {
      label: 'Source (flat)',
      options: [
        ['pdb', 'PDB'],
        ['url', 'URL'],
      ],
    },
  ),
  isoValue: PD.MappedStatic(
    'relative',
    {
      relative: PD.Group({ value: PD.Numeric(1.5, { min: -5, max: 5, step: 0.1 }, { label: 'σ' }) }),
      absolute: PD.Group({ value: PD.Numeric(0.035, { min: 0, max: 1, step: 0.001 }) }),
    },
    {
      label: 'Iso value (nested)',
      options: [
        ['relative', 'Relative'],
        ['absolute', 'Absolute'],
      ],
    },
  ),
  quality: PD.MappedStatic(
    'auto',
    {
      auto: PD.Group({}),
      custom: PD.Numeric(3, { min: 0, max: 6, step: 1 }, { label: 'Detail' }),
    },
    {
      label: 'Quality (single value)',
      options: [
        ['auto', 'Auto'],
        ['custom', 'Custom'],
      ],
    },
  ),
});

/** Visibility: rows hidden by another's value (hideIf), hidden params, and essentials first with the rest under Advanced. */
export const Visibility = render({
  mode: PD.Select(
    'perspective',
    [
      ['perspective', 'Perspective'],
      ['orthographic', 'Orthographic'],
    ],
    { isEssential: true },
  ),
  fov: PD.Numeric(
    45,
    { min: 10, max: 130, step: 1 },
    {
      label: 'Field of view (perspective only)',
      isEssential: true,
      hideIf: (p: any) => p.mode !== 'perspective',
    },
  ),
  internal: PD.Numeric(1, undefined, { isHidden: true }),
  resetDuration: PD.Numeric(250, { min: 0, max: 1000, step: 1 }, { label: 'Reset duration (ms)' }),
  manualReset: PD.Boolean(false, { label: 'Manual reset' }),
});

/** Search (Settings): only matching rows; groups open; a group whose name matches shows all its rows. Try "zoom", "trackball". */
export const Search = render(
  {
    pickPadding: PD.Numeric(3, { min: 0, max: 10, step: 1 }),
    trackball: PD.Group(
      {
        rotateSpeed: PD.Numeric(5, { min: 1, max: 10, step: 1 }),
        zoomSpeed: PD.Numeric(7, { min: 1, max: 15, step: 1 }),
      },
      { label: 'Trackball' },
    ),
    tooltips: PD.Boolean(true),
  },
  { filter: true },
);

/**
 * Long labels at a narrow width: the label column fits the longest label up to half the form; past that a label stays on
 * one line and fades (its full text in its tooltip); rows stay 24px and values keep one column.
 */
export const LongLabels = render(
  {
    checkeredTransparentBackground: PD.Boolean(true),
    userInteractionReleaseMs: PD.Numeric(250, { min: 0, max: 1000, step: 1 }),
    sceneRadiusFactor: PD.Numeric(1, { min: 1, max: 10, step: 0.1 }),
    highlightEdgeColor: PD.Color(0x999999 as any),
    transparencyMode: PD.Select('blended', [
      ['blended', 'Weighted, blended'],
      ['wboit', 'Weighted, order-independent'],
    ]),
    // With help: the ? after the last word.
    backgroundOcclusionRadiusScale: PD.Numeric(
      0.5,
      { min: 0, max: 2, step: 0.05 },
      {
        description: 'How far around each point occlusion is sampled, relative to the scene radius.',
      },
    ),
  },
  { width: 240 },
);

/** A density level, σ or absolute; a map with a σ of 0.1 around 0 converts between them. */
type Level = { kind: 'relative' | 'absolute'; value: number };
const toAbsolute = (sigmas: number) => Math.round(sigmas * 0.1 * 1000) / 1000;
const toRelative = (value: number) => Math.round((value / 0.1) * 10) / 10;

/**
 * The less common types: multi-select, vec3, mat4 (opens a grid), interval, a line graph (opens an editor), a color
 * list (opens its presets), an object list (items with their own forms), url, converted, conditioned, script; a plain
 * value has no row (as in the classic UI).
 */
export const MoreTypes = render({
  chains: PD.MultiSelect(
    ['A', 'C'],
    [
      ['A', 'Chain A'],
      ['B', 'Chain B'],
      ['C', 'Chain C'],
      ['D', 'Chain D'],
    ],
    { label: 'Chains (multi-select)' },
  ),
  axis: PD.Vec3(Vec3.create(0, 0, 1), undefined, { label: 'Axis (vec3)' }),
  transform: PD.Mat4(Mat4.identity(), { label: 'Transform (mat4)' }),
  range: PD.Interval([2, 7.5], { min: 0, max: 10, step: 0.5 }, { label: 'Range (interval)' }),
  transfer: PD.LineGraph(
    [
      [0, 0],
      [0.3, 0.05],
      [0.55, 0.7],
      [0.8, 0.9],
      [1, 1],
    ] as any,
    { label: 'Transfer (line graph)' },
  ),
  colors: PD.ColorList('viridis' as any, { label: 'Colors (color list)' }),
  channels: PD.ObjectList(
    {
      name: PD.Text('Channel'),
      isoValue: PD.Numeric(1.5, { min: -5, max: 5, step: 0.1 }),
      color: PD.Color(0x3362b2 as any),
    },
    (e: any) => `${e.name} · ${e.isoValue}σ`,
    {
      label: 'Channels (object list)',
      defaultValue: [
        { name: '2Fo-Fc', isoValue: 1.5, color: 0x3362b2 },
        { name: 'Fo-Fc', isoValue: 3, color: 0x33bb33 },
      ] as any,
      ctor: () => ({ name: 'Channel', isoValue: 1, color: 0x808080 }),
    },
  ),
  source: PD.Url('https://www.ebi.ac.uk/pdbe/entry-files/download/1tqn.bcif', { label: 'Source (url)' }),
  opacity: PD.Converted(
    (v: number) => v * 100,
    (v: number) => v / 100,
    PD.Numeric(35, { min: 0, max: 100, step: 1 }, { label: 'Opacity % (converted)' }),
  ),
  // The value carries its kind, as Mol*'s volume iso values do: the branch never follows the number's size.
  level: PD.Conditioned(
    { kind: 'relative', value: 1.5 } as Level,
    {
      relative: PD.Converted(
        (v: Level) => v.value,
        (value: number): Level => ({ kind: 'relative', value }),
        PD.Numeric(1.5, { min: -5, max: 5, step: 0.1 }, { label: 'σ' }),
      ),
      absolute: PD.Converted(
        (v: Level) => v.value,
        (value: number): Level => ({ kind: 'absolute', value }),
        PD.Numeric(0.15, { min: -0.5, max: 0.5, step: 0.001 }, { label: 'Value' }),
      ),
    },
    (v: Level) => v.kind,
    (v: Level, c: string): Level =>
      c === v.kind
        ? v
        : c === 'absolute'
          ? { kind: 'absolute', value: toAbsolute(v.value) }
          : { kind: 'relative', value: toRelative(v.value) },
    { label: 'Level (conditioned)' },
  ) as any,
  query: PD.Script({ language: 'pymol', expression: 'resn HEM' } as any, { label: 'Query (script)' }),
  data: PD.Value({ a: 1 }),
});

const ChainColors: [string, number][] = [
  ['A', 0x1b9e77],
  ['B', 0xd95f02],
  ['C', 0x7570b3],
  ['D', 0xe7298a],
];

/**
 * Help: a ? by a label (on hover, or while open) opens the param's description and legend under its row: a table of
 * colors, a scale; a theme choice shows the chosen theme's (`help(value)`), as representation themes do.
 */
export const HelpAndLegends = render({
  quality: PD.Select(
    'auto',
    [
      ['auto', 'Auto'],
      ['high', 'High'],
      ['low', 'Low'],
    ],
    { description: 'How finely surfaces and cartoons are tessellated; `auto` follows the structure size.' },
  ),
  chainColors: PD.Boolean(true, {
    label: 'Chain colors',
    description: 'Gives every chain a color based on its `asym_id` value.',
    legend: { kind: 'table-legend', table: ChainColors as any },
  }),
  bFactor: PD.Boolean(false, {
    label: 'B-factor',
    description: 'Colors by the atoms’ B-factor, low to high.',
    legend: { kind: 'scale-legend', minLabel: '12.4', maxLabel: '86.1', colors: [0x0000ff, 0xffffff, 0xff0000] as any },
  }),
  colorTheme: PD.MappedStatic(
    'chain-id',
    {
      'chain-id': PD.Group({}),
      'element-symbol': PD.Group({ carbonColor: PD.Color(0x808080 as any, { label: 'Carbon color' }) }),
      uncertainty: PD.Group({}),
    },
    {
      label: 'Color theme',
      options: [
        ['chain-id', 'Chain ID'],
        ['element-symbol', 'Element'],
        ['uncertainty', 'Uncertainty'],
      ],
      help: (v: { name: string }) =>
        v.name === 'chain-id'
          ? {
              description: 'Gives every chain a color based on its `asym_id` value.',
              legend: { kind: 'table-legend', table: ChainColors as any },
            }
          : v.name === 'element-symbol'
            ? {
                description: 'Colors atoms by their element; carbons in the chosen carbon color.',
                legend: {
                  kind: 'table-legend',
                  table: [
                    ['C', 0x808080],
                    ['N', 0x3050f8],
                    ['O', 0xff0d0d],
                    ['S', 0xffff30],
                    ['Fe', 0xe06633],
                  ] as any,
                },
              }
            : {
                description: 'Colors by the atoms’ uncertainty (B-factor or pLDDT).',
                legend: {
                  kind: 'scale-legend',
                  minLabel: 'certain',
                  maxLabel: 'uncertain',
                  colors: [0x0053d6, 0x65cbf3, 0xffdb13, 0xff7d45] as any,
                },
              },
    },
  ),
});

/** The form's columns (styles §5), drawn over a form with nested levels: every row's text and field on them. */
const line = (left: string, color = 'var(--msui-accent)') =>
  `position: absolute; top: 0; bottom: 0; left: ${left}; width: 1px; background: ${color}; opacity: 0.5; pointer-events: none;`;

export const Columns: StoryObj = {
  render: () =>
    story(() => {
      const form = new ParamsFormModel({
        name: PD.Text('1tqn', { label: 'Name' }),
        opacity: PD.Numeric(0.6, { min: 0, max: 1, step: 0.01 }),
        visible: PD.Boolean(true),
        type: PD.MappedStatic(
          'cartoon',
          {
            cartoon: PD.Group({
              quality: PD.Select('auto', PD.arrayToOptions(['auto', 'high', 'low'])),
              detail: PD.Group({ level: PD.Numeric(2, { min: 0, max: 4, step: 1 }) }, { isExpanded: true }),
            }),
            spacefill: PD.Group({}),
          },
          { label: 'Representation' },
        ),
        color: PD.Color(0x3362b2 as any),
        group: PD.Group(
          { samples: PD.Numeric(32), inner: PD.Group({ depth: PD.Boolean(false) }, { isExpanded: true }) },
          { isExpanded: true, label: 'A group' },
        ),
      });
      form.fold('/type#params').set(true);
      const inset = 'var(--msui-inset)';
      // The field column: where the form's label column (fitted to its labels) ends, plus the gap; read once shown.
      const field = 'var(--story-field, calc(var(--msui-inset) + var(--msui-label) + var(--msui-space-3)))';
      return (
        <div style="padding: var(--msui-space-5);" ref={(el: HTMLElement) => installTooltips(el)}>
          <div
            class="msui-surface"
            style="position: relative; width: 320px;"
            ref={(el: HTMLElement) =>
              observeResize(el, () => {
                const control = el.querySelector('[data-field] > :last-child');
                if (!control) return;
                const x = control.getBoundingClientRect().left - el.getBoundingClientRect().left;
                el.style.setProperty('--story-field', `${x}px`);
              })
            }
          >
            <ParamsForm form={form} />
            <span style={line(inset)} />
            <span style={line(`calc(${inset} + var(--msui-form-indent))`, 'var(--msui-text-subtle)')} />
            <span style={line(`calc(${inset} + 2 * var(--msui-form-indent))`, 'var(--msui-text-subtle)')} />
            <span style={line(field)} />
            <span style={line(`calc(${field} + var(--msui-control-inset))`, 'var(--msui-success)')} />
            <span style={line(`calc(100% - ${inset})`)} />
          </div>
        </div>
      );
    }),
};
