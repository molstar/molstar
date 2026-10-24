/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type {
  PresetsModel,
  QuickStylesModel,
  SceneInspection,
  SceneMenuItem,
  SceneModel,
  SceneNode,
  SceneNodeKind,
} from '../model/scene.js';
import { type ObservableState, observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { DefaultRows } from './structure-rows.js';

/** Representation params shaped like Mol*'s (type, color theme as a mapped param, size, opacity, quality, material). */
const RepresentationTypes = [
  ['cartoon', 'Cartoon'],
  ['ball-and-stick', 'Ball & Stick'],
  ['spacefill', 'Spacefill'],
  ['molecular-surface', 'Molecular Surface'],
  ['gaussian-surface', 'Gaussian Surface'],
] as const;

const RepresentationParams = {
  type: PD.Select<string>('cartoon', RepresentationTypes, PD.Essential),
  colorTheme: PD.MappedStatic(
    'chain-id',
    {
      'chain-id': PD.Group({
        palette: PD.Select('many-distinct', [
          ['many-distinct', 'Many distinct'],
          ['set-1', 'Set 1'],
          ['pastel', 'Pastel'],
        ]),
      }),
      'secondary-structure': PD.EmptyGroup(),
      'element-symbol': PD.Group({
        carbonColor: PD.Select('chain-id', [
          ['chain-id', 'Chain ID'],
          ['uniform', 'Uniform'],
        ]),
      }),
      uniform: PD.Group({ value: PD.Color(0x2e86c1 as any) }),
    },
    { label: 'Color', isEssential: true },
  ),
  sizeFactor: PD.Numeric(1, { min: 0.1, max: 3, step: 0.05 }, { label: 'Size', isEssential: true }),
  alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Opacity', isEssential: true }),
  quality: PD.Select('auto', [
    ['auto', 'Auto'],
    ['medium', 'Medium'],
    ['high', 'High'],
    ['highest', 'Highest'],
  ]),
  flatShaded: PD.Boolean(false),
  ignoreLight: PD.Boolean(false),
  material: PD.Group(
    {
      metalness: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
      roughness: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
      bumpiness: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
    },
    { label: 'Material' },
  ),
};

/** Volume representation params shaped like Mol*'s: the iso value relative to σ or absolute, first. */
const VolumeRepresentationParams = {
  isoValue: PD.MappedStatic(
    'relative',
    {
      relative: PD.Group({ value: PD.Numeric(1.5, { min: -5, max: 10, step: 0.05 }, { label: 'σ' }) }),
      absolute: PD.Group({ value: PD.Numeric(0.035, { min: -0.2, max: 0.5, step: 0.001 }, { label: 'Value' }) }),
    },
    {
      label: 'Iso value',
      isEssential: true,
      options: [
        ['relative', 'Relative (σ)'],
        ['absolute', 'Absolute'],
      ],
    },
  ),
  type: PD.Select(
    'isosurface',
    [
      ['isosurface', 'Isosurface'],
      ['slice', 'Slice'],
      ['direct-volume', 'Direct volume'],
    ],
    PD.Essential,
  ),
  color: PD.Color(0x3362b2 as any, { isEssential: true }),
  alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Opacity', isEssential: true }),
  wireframe: PD.Boolean(false),
  flatShaded: PD.Boolean(false),
};

// Mol*'s structure hierarchy presets (how the structure is built from its model), then the structure transform's own
// params under Advanced.
const HierarchyPresets: [string, string][] = [
  ['default', 'Default (Assembly)'],
  ['unit-cell', 'Unit cell'],
  ['supercell', 'Super cell'],
  ['crystal-contacts', 'Crystal contacts'],
];

const StructureParams = {
  hierarchy: PD.Select('default', HierarchyPresets, {
    label: 'Hierarchy',
    isEssential: true,
    description: 'How the structure is built from its model (Mol*’s hierarchy presets).',
  }),
  type: PD.Select('assembly', [
    ['model', 'Model'],
    ['assembly', 'Assembly'],
    ['symmetry', 'Symmetry mates'],
  ]),
  assembly: PD.Select('1', [
    ['1', 'Assembly 1'],
    ['2', 'Assembly 2'],
  ]),
  model: PD.Numeric(1, { min: 1, max: 20, step: 1 }, { label: 'Model index' }),
};

/** The crystallographic unit cell (Mol*'s unit cell representation): drawn as a box, hidden by default. */
const UnitCellParams = {
  color: PD.Color(0x808080 as any, { isEssential: true }),
  size: PD.Numeric(0.33, { min: 0.1, max: 2, step: 0.01 }, { label: 'Line size', isEssential: true }),
  alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Opacity' }),
};

/** A measurement's look (the measurement manager's line and label params), adjustable after it is added. */
const MeasurementStyleParams = {
  lineColor: PD.Color(0xffee88 as any, { label: 'Line color', isEssential: true }),
  lineSize: PD.Numeric(0.075, { min: 0.01, max: 0.5, step: 0.005 }, { label: 'Line size', isEssential: true }),
  dashed: PD.Boolean(true, { isEssential: true }),
  textColor: PD.Color(0x000000 as any, { label: 'Text color', isEssential: true }),
  textSize: PD.Numeric(0.5, { min: 0.1, max: 5, step: 0.1 }, { label: 'Text size', isEssential: true }),
};

const ComponentParams = {
  label: PD.Text(''),
};

const ColorLabels: Record<string, string> = {
  'chain-id': 'Chain ID',
  'secondary-structure': 'Secondary structure',
  'element-symbol': 'Element',
  uniform: 'Uniform',
};

const VolumeTypeLabels: Record<string, string> = {
  isosurface: 'Isosurface',
  slice: 'Slice',
  'direct-volume': 'Direct volume',
};

const ReprIcons: Record<string, IconName> = {
  cartoon: 'cartoon',
  'ball-and-stick': 'ball-and-stick',
  spacefill: 'spacefill',
  'molecular-surface': 'surface',
  'gaussian-surface': 'surface',
};

class MockNode implements SceneNode {
  readonly #label: ObservableState<string>;
  readonly #description: ObservableState<string | undefined>;
  readonly #hidden = observableState(false);
  readonly #expanded: ObservableState<boolean>;
  readonly #children = observableState<readonly SceneNode[]>([]);
  parent: MockNode | undefined;
  menu?: readonly SceneMenuItem[];
  /** False: no eye (Entities, Layers). */
  hasEye?: boolean;
  removable?: boolean;

  constructor(
    readonly id: string,
    readonly kind: SceneNodeKind,
    public icon: IconName,
    label: string,
    description: string | undefined,
    expanded: boolean,
    readonly onHighlight: ((node: MockNode, on: boolean) => void) | undefined,
    readonly onRemove: (node: MockNode) => void,
  ) {
    this.#label = observableState(label);
    this.#description = observableState(description);
    this.#expanded = observableState(expanded);
  }

  get label() {
    return this.#label.asReadonly;
  }
  get description() {
    return this.#description.asReadonly;
  }
  get isHidden() {
    return this.#hidden.asReadonly;
  }
  get isExpanded() {
    return this.#expanded.asReadonly;
  }
  get children() {
    return this.#children.asReadonly;
  }

  /** Builds the children on first expand (a structure's contents); undefined once loaded. */
  loader?: () => MockNode[];
  get lazy() {
    return !!this.loader;
  }
  load() {
    const loader = this.loader;
    if (!loader) return;
    this.loader = undefined;
    // As the plugin would compute entities and chains: a moment later, not on the way in.
    setTimeout(() => this.add(...loader()), 250);
  }

  /** Puts children first (a structure's contents, above its components). */
  addFirst(...children: MockNode[]): this {
    for (const c of children) c.parent = this;
    this.#children.set([...children, ...this.#children.value]);
    return this;
  }

  add(...children: MockNode[]): this {
    for (const c of children) c.parent = this;
    this.#children.set([...this.#children.value, ...children]);
    return this;
  }

  drop(child: SceneNode) {
    this.#children.set(this.#children.value.filter((c) => c !== child));
  }

  setDescription(d: string | undefined) {
    this.#description.set(d);
  }
  setLabel(l: string) {
    this.#label.set(l);
  }

  /** What the eye does: hide or show (a streaming row turns streaming on and off instead). */
  onToggle?: () => void;
  toggleVisibility() {
    if (this.onToggle) this.onToggle();
    else this.#hidden.set(!this.#hidden.value);
  }
  setHidden(hidden: boolean) {
    this.#hidden.set(hidden);
  }
  toggleExpanded() {
    this.#expanded.set(!this.#expanded.value);
  }
  highlight(on: boolean) {
    this.onHighlight?.(this, on);
  }
  focus() {}
  remove() {
    if (this.parent) this.parent.drop(this);
    else this.onRemove(this);
  }
}

/** What the stories know about a mock structure (what the plugin reads from the model). */
export interface MockStructureInfo {
  method: 'x-ray' | 'em';
  /** The assembly's global symmetry, when it has one (C1 has none). */
  symmetry?: string;
}

export interface MockScene extends SceneModel {
  /** The structure behind a node, for sections' `when`. */
  structureInfo(node: SceneNode): MockStructureInfo | undefined;
  /** An extension row's options for a structure (its `enabled` switch turns the row on). */
  rowForm(structure: SceneNode, rowId: string): ParamsFormModel<any> | undefined;
  /** The volume streaming row's options for a structure. */
  streaming(structure: SceneNode): ParamsFormModel<any>;
  /** The measurements made: the rows under Measurements (each picked shows its style). */
  measurements: ReadonlyState<readonly SceneNode[]>;
  /** A measurement's style (shared by its row's inspector and the Selection tab's list). */
  measurementStyle(node: SceneNode): ParamsFormModel<any> | undefined;
  /** Adds a measurement row under Measurements (what the measurement manager would create). */
  addMeasurement(label: string, description: string, icon: IconName): void;
  /** Removes everything (the plugin's clear). */
  clear(): void;
}

/** What an extension's structure row can put under itself once on: volume channels, a shape (symmetry axes). */
export interface MockRowChildren {
  volume(id: string, label: string, iso: { relative?: number; absolute?: number }, color: number): MockNode;
  shape(id: string, label: string, description: string, icon: IconName): MockNode;
}

/**
 * A row an extension adds under each structure it applies to (volume streaming, assembly symmetry), as its UI
 * contribution would register it: off (hidden) until its eye turns it on; on, what it shows appears under it; picked,
 * its options. The tree knows no extension by name.
 */
export interface MockRowSpec {
  id: string;
  label: string;
  icon: IconName;
  when(info: MockStructureInfo): boolean;
  /** Its options; `enabled` (a boolean) is the row's on/off. */
  params(info: MockStructureInfo): PD.Params;
  describe(values: any, info: MockStructureInfo): string;
  children(values: any, info: MockStructureInfo, make: MockRowChildren, id: string): MockNode[];
}

export type MockScenePreset = 'single' | 'many' | 'em' | 'particles';

/**
 * Scenes as the plugin would show them: `single` is 1TQN with the default preset (components, representations, a
 * measurement); `many` is four CYP3A4 structures superposed; `em` is an EMDB map with its fitted model; `particles` is
 * cryo-ET ribosome picks (spacefill, orientation, the unit cell) beside a Simularium trajectory. Forms change
 * the nodes as the plugin would.
 */
export function mockScene(
  options: {
    preset?: MockScenePreset;
    onHighlight?: (label: string | undefined) => void;
    /** Reports a menu action (what the plugin would do). */
    onAction?: (text: string) => void;
    /** Extensions' structure rows (volume streaming, assembly symmetry by default). */
    rows?: readonly MockRowSpec[];
    /** Adding a measurement: the Selection tab's Measure and Label section. */
    onMeasure?: () => void;
  } = {},
): MockScene {
  const highlight = (node: MockNode, on: boolean) => options.onHighlight?.(on ? node.label.value : undefined);
  const forms = new Map<string, ParamsFormModel<any>>();
  const facts = new Map<string, readonly (readonly [string, string])[]>();
  const infos = new Map<string, MockStructureInfo>();
  const rowForms = new Map<string, ParamsFormModel<any>>();
  const rows = options.rows ?? DefaultRows;
  const roots = observableState<readonly SceneNode[]>([]);
  const removeRoot = (node: MockNode) => roots.set(roots.value.filter((n) => n !== node));
  const node = (
    id: string,
    kind: SceneNodeKind,
    icon: IconName,
    label: string,
    description?: string,
    expanded = true,
  ) => new MockNode(id, kind, icon, label, description, expanded, highlight, removeRoot);

  const repr = (id: string, type: string, color: string, alpha = 1) => {
    const n = node(id, 'representation', ReprIcons[type] ?? 'cartoon', '', undefined, false);
    const values = {
      ...PD.getDefaultValues(RepresentationParams),
      type,
      alpha,
      colorTheme: {
        name: color,
        params:
          color === 'uniform'
            ? { value: 0xe74c3c }
            : PD.getDefaultValues((RepresentationParams.colorTheme.map(color) as PD.Group<any>).params),
      },
    } as PD.Values<typeof RepresentationParams>;
    const sync = (v: PD.Values<typeof RepresentationParams>) => {
      n.setLabel(RepresentationTypes.find(([t]) => t === v.type)?.[1] ?? v.type);
      n.icon = ReprIcons[v.type] ?? 'cartoon';
      n.setDescription(
        `${ColorLabels[v.colorTheme.name] ?? v.colorTheme.name}${v.alpha < 1 ? ` · ${Math.round(v.alpha * 100)}%` : ''}`,
      );
    };
    sync(values);
    forms.set(id, new ParamsFormModel(RepresentationParams, values, { onChange: sync }));
    return n;
  };

  const volumeRepr = (
    id: string,
    label: string | undefined,
    iso: { relative?: number; absolute?: number },
    color: number,
  ) => {
    const n = node(id, 'volume-representation', 'volume', '', undefined, false);
    const values = {
      ...PD.getDefaultValues(VolumeRepresentationParams),
      isoValue:
        iso.absolute !== undefined
          ? { name: 'absolute', params: { value: iso.absolute } }
          : { name: 'relative', params: { value: iso.relative ?? 1.5 } },
      color,
    } as PD.Values<typeof VolumeRepresentationParams>;
    const sync = (v: PD.Values<typeof VolumeRepresentationParams>) => {
      const type = VolumeTypeLabels[v.type] ?? v.type;
      n.setLabel(label ?? type);
      const level = v.isoValue.name === 'relative' ? `${v.isoValue.params.value}σ` : `${v.isoValue.params.value}`;
      n.setDescription(label ? `${type} · ${level}` : level);
    };
    sync(values);
    forms.set(id, new ParamsFormModel(VolumeRepresentationParams, values, { onChange: sync }));
    return n;
  };

  const act = (text: string) => () => options.onAction?.(text);
  // Mol*'s component menu.
  const componentMenu = (label: string): SceneMenuItem[] => [
    {
      id: 'add-repr',
      label: 'Add representation',
      icon: 'plus',
      items: RepresentationTypes.map(([t, l]) => ({
        id: t,
        label: l,
        icon: ReprIcons[t],
        run: act(`${label}: add ${l}`),
      })),
    },
    {
      id: 'coloring',
      label: 'Set coloring',
      icon: 'palette',
      items: Object.entries(ColorLabels).map(([c, l]) => ({ id: c, label: l, run: act(`${label}: color by ${l}`) })),
    },
    {
      id: 'modify',
      label: 'Modify by selection',
      icon: 'mouse-pointer-click',
      items: [
        { id: 'add', label: 'Include', icon: 'squares-unite', run: act(`${label}: include the selection`) },
        { id: 'remove', label: 'Subtract', icon: 'squares-subtract', run: act(`${label}: subtract the selection`) },
        {
          id: 'intersect',
          label: 'Intersect',
          icon: 'squares-intersect',
          run: act(`${label}: intersect with the selection`),
        },
      ],
    },
    { id: 'select', label: 'Select this', icon: 'mouse-pointer-2', run: act(`Selected ${label}`) },
  ];

  // A structure's contents: entities, then their chains, built when first expanded (mesoscale structures have
  // thousands of chains); each row can be selected, focused, colored, or made a component.
  const contentsMenu = (label: string): SceneMenuItem[] => [
    { id: 'select', label: 'Select', icon: 'mouse-pointer-2', run: act(`Selected ${label}`) },
    { id: 'focus', label: 'Focus', icon: 'focus', run: act(`Focused on ${label}`) },
    { id: 'color', label: 'Color…', icon: 'palette', run: act(`Coloring ${label}…`) },
    { id: 'component', label: 'Create component', icon: 'layers', run: act(`Component from ${label}`) },
  ];
  interface EntitySpec {
    id: string;
    label: string;
    type: 'polymer' | 'non-polymer' | 'water';
    chains: readonly { id: string; description: string }[];
  }
  const EntityIcons: Record<EntitySpec['type'], IconName> = {
    polymer: 'pick-chain',
    'non-polymer': 'ligand',
    water: 'pick-atom',
  };
  const contents = (sid: string, entities: readonly EntitySpec[]) => {
    const chains = new Set(entities.flatMap((e) => e.chains.map((c) => c.id))).size;
    const count = (k: number, one: string, many: string) => `${k} ${k === 1 ? one : many}`;
    // What the structure is made of (entities, then chains): first under the structure, named for what it shows.
    const n = node(
      `${sid}-contents`,
      'contents',
      'list-tree',
      'Entities',
      `${count(entities.length, 'entity', 'entities')} · ${count(chains, 'chain', 'chains')}`,
      false,
    );
    n.removable = false;
    n.hasEye = false;
    n.loader = () =>
      entities.map((e) => {
        const en = node(
          `${sid}-entity-${e.id}`,
          'entity',
          EntityIcons[e.type],
          e.label,
          `Entity ${e.id} · ${e.chains.length} chain${e.chains.length > 1 ? 's' : ''}`,
          false,
        );
        en.menu = contentsMenu(e.label);
        en.removable = false;
        facts.set(en.id, [
          ['Type', e.type],
          ['Chains', e.chains.map((c) => c.id).join(', ')],
        ]);
        en.loader = () =>
          e.chains.map((c) => {
            const cn = node(
              `${sid}-chain-${e.id}-${c.id}`,
              'chain',
              'pick-chain',
              `Chain ${c.id}`,
              c.description,
              false,
            );
            cn.menu = contentsMenu(`chain ${c.id}`);
            cn.removable = false;
            facts.set(cn.id, [
              ['Chain', c.id],
              ['Entity', e.label],
              ['Contains', c.description],
            ]);
            return cn;
          });
        return en;
      });
    return n;
  };

  const component = (id: string, label: string, description: string, ...reprs: MockNode[]) => {
    const n = node(id, 'component', 'layers', label, description).add(...reprs);
    n.menu = componentMenu(label);
    forms.set(id, new ParamsFormModel(ComponentParams, { label }, { onChange: (v) => n.setLabel(v.label || label) }));
    return n;
  };

  const structure = (
    id: string,
    label: string,
    description: string,
    info: MockStructureInfo & { facts: readonly (readonly [string, string])[]; expanded?: boolean },
    ...components: MockNode[]
  ) => {
    // Its rows grouped (viewer-spec.md §3): its components, then the others (the unit cell, the extensions' rows).
    const comps = node(
      `${id}#components`,
      'components',
      'shapes',
      'Components',
      `${components.length} component${components.length === 1 ? '' : 's'}`,
    ).add(...components);
    comps.removable = false;
    const layers = node(`${id}#other`, 'other', 'layers', 'Other');
    layers.removable = false;
    layers.hasEye = false;
    const n = node(id, 'structure', 'structure', label, description, info.expanded ?? true).add(comps, layers);
    forms.set(id, new ParamsFormModel(StructureParams));
    facts.set(id, info.facts);
    infos.set(id, info);
    n.menu = [
      {
        id: 'hierarchy',
        label: 'Hierarchy preset',
        icon: 'layers',
        items: HierarchyPresets.map(([p, l]) => ({
          id: p,
          label: l,
          run: () => {
            forms.get(id)?.set('hierarchy', p);
            options.onAction?.(`${label}: ${l}`);
          },
        })),
      },
    ];
    for (const spec of rows) if (spec.when(info)) layers.add(extensionRow(n, info, spec));
    return n;
  };

  // An extension's row under a structure: off (hidden) until its eye turns it on; on, what it shows appears under it
  // (expanded); picked, its options. Its description follows its options.
  const make: MockRowChildren = {
    volume: (id, label, iso, color) => volumeRepr(id, label, iso, color),
    shape: (id, label, description, icon) => node(id, 'representation', icon, label, description, false),
  };
  const extensionRow = (owner: MockNode, info: MockStructureInfo, spec: MockRowSpec) => {
    const id = `${owner.id}-${spec.id}`;
    const n = node(id, 'extension', spec.icon, spec.label, undefined, false);
    const form = new ParamsFormModel(spec.params(info), undefined, {
      onChange: (v) => {
        n.setHidden(!v.enabled);
        for (const c of n.children.value) n.drop(c);
        if (v.enabled) {
          n.add(...spec.children(v, info, make, id));
          if (!n.isExpanded.value && n.children.value.length) n.toggleExpanded();
        }
        n.setDescription(spec.describe(v, info));
      },
    });
    n.setHidden(true);
    n.setDescription(spec.describe(form.values.value, info));
    n.onToggle = () => form.set('enabled', !form.values.value.enabled);
    n.removable = false;
    forms.set(id, form);
    rowForms.set(id, form);
    return n;
  };

  const volume = (id: string, label: string, description: string, f: readonly (readonly [string, string])[]) => {
    const n = node(id, 'volume', 'volume', label, description);
    facts.set(id, f);
    return n;
  };

  let measurementsNode: MockNode | undefined;
  // The measurements made (the group's rows), followed wherever the group is (it appears with the first).
  const measurementList = observableState<readonly SceneNode[]>([]);
  const followMeasurements = (n: MockNode) => n.children.subscribe((c) => measurementList.set(c));
  const measurementNode = (id: string, icon: IconName, label: string, description: string) => {
    const m = node(id, 'measurement', icon, label, description, false);
    forms.set(id, new ParamsFormModel(MeasurementStyleParams));
    return m;
  };
  const measurementsDescription = (n: MockNode) => {
    const count = n.children.value.length;
    n.setDescription(count ? `${count} measurement${count === 1 ? '' : 's'}` : 'None yet');
  };
  const measurements = () => {
    const n = node('measurements', 'measurements', 'ruler', 'Measurements', '1 distance').add(
      measurementNode('d1', 'measure-distance', 'CYS 442 SG – HEM FE', '2.31 Å'),
    );
    // The group stays; its menu adds (in the Selection tab) or removes all the measurements, not the group.
    n.removable = false;
    n.menu = [
      { id: 'add', label: 'Add measurement…', icon: 'plus', run: () => options.onMeasure?.() },
      {
        id: 'remove-all',
        label: 'Remove all',
        icon: 'trash-2',
        danger: true,
        run: () => {
          for (const c of n.children.value) n.drop(c);
          measurementsDescription(n);
        },
      },
    ];
    measurementsNode = n;
    followMeasurements(n);
    return n;
  };
  let measurementCount = 1;

  const cyp = (id: string, label: string, ligand: string, description: string, expanded: boolean) =>
    structure(
      id,
      label,
      description,
      {
        method: 'x-ray',
        expanded,
        facts: [
          ['Method', 'X-ray'],
          ['Chains', '1'],
          ['Residues', '470'],
          ['Atoms', '3,624'],
        ],
      },
      component(`${id}-polymer`, 'Polymer', 'Chain A', repr(`${id}-polymer-cartoon`, 'cartoon', 'chain-id')),
      component(`${id}-ligand`, 'Ligand', ligand, repr(`${id}-ligand-bs`, 'ball-and-stick', 'element-symbol')),
    );

  const preset = options.preset ?? 'single';
  if (preset === 'single') {
    const s = structure(
      '1tqn',
      '1TQN',
      'Assembly 1 · 3,624 atoms',
      {
        method: 'x-ray',
        facts: [
          ['Method', 'X-ray · 2.05 Å'],
          ['Chains', '1'],
          ['Residues', '470'],
          ['Atoms', '3,624'],
        ],
      },
      component('polymer', 'Polymer', 'Chain A · 470 residues', repr('polymer-cartoon', 'cartoon', 'chain-id')),
      component(
        'ligand',
        'Ligand',
        'HEM 508',
        repr('ligand-bs', 'ball-and-stick', 'element-symbol'),
        repr('ligand-surface', 'molecular-surface', 'uniform', 0.35),
      ),
      component('water', 'Water', '186 molecules', repr('water-bs', 'ball-and-stick', 'element-symbol', 0.6)),
      component('ions', 'Ions', 'SO4 ×2', repr('ions-bs', 'ball-and-stick', 'element-symbol')),
    );
    const groups = s.children.value as readonly MockNode[];
    const water = groups[0].children.value[2] as MockNode;
    water.toggleVisibility();
    water.toggleExpanded();
    // The unit cell: part of the structure's record, hidden until asked for.
    const cell = node('1tqn-cell', 'unit-cell', 'box', 'Unit cell', 'I 2 2 2', false);
    cell.toggleVisibility();
    cell.menu = [{ id: 'select', label: 'Select this', icon: 'mouse-pointer-2', run: act('Selected the unit cell') }];
    forms.set('1tqn-cell', new ParamsFormModel(UnitCellParams));
    facts.set('1tqn-cell', [['Space group', 'I 2 2 2']]);
    groups[1].addFirst(cell);
    s.addFirst(
      contents('1tqn', [
        { id: '1', label: 'Cytochrome P450 3A4', type: 'polymer', chains: [{ id: 'A', description: '470 residues' }] },
        {
          id: '2',
          label: 'Protoporphyrin IX containing Fe',
          type: 'non-polymer',
          chains: [{ id: 'A', description: 'HEM 508' }],
        },
        { id: '3', label: 'Sulfate ion', type: 'non-polymer', chains: [{ id: 'A', description: 'SO4 ×2' }] },
        { id: '4', label: 'Water', type: 'water', chains: [{ id: 'A', description: '186 molecules' }] },
      ]),
    );
    roots.set([s, measurements()]);
  } else if (preset === 'many') {
    // Four CYP3A4 structures superposed on the first: collapsed but the first, so the list stays short.
    roots.set([
      cyp('1tqn', '1TQN', 'HEM 508', 'Reference', true),
      cyp('1w0e', '1W0E', 'HEM 1499', 'Superposed · RMSD 0.38 Å', false),
      cyp('1w0f', '1W0F', 'HEM · PROGESTERONE', 'Superposed · RMSD 0.41 Å', false),
      cyp('1w0g', '1W0G', 'HEM · METYRAPONE', 'Superposed · RMSD 0.35 Å', false),
      measurements(),
    ]);
  } else if (preset === 'particles') {
    // Particle lists (the lazy particles feature's rows): picks with their representations and unit cell, and a
    // trajectory's list.
    const picks = node('ribosomes', 'particles', 'shapes', 'cryoet-ribosomes.ndjson', '200 particles').add(
      node('ribosomes-spacefill', 'particles-representation', 'shapes', 'Spacefill', 'Particle Index', false),
      node('ribosomes-orientation', 'particles-representation', 'shapes', 'Orientation', 'Particle Index', false),
      node('ribosomes-cell', 'unit-cell', 'box', 'Unit cell', undefined, false),
    );
    facts.set('ribosomes', [
      ['Particles', '200'],
      ['Source', 'cryoet-data-portal-ndjson'],
      ['Orientations', 'yes'],
      ['Targets', '1: target 0: 1CRN'],
    ]);
    const trajectory = node('simularium', 'particles', 'shapes', 'small.simularium', '5 particles · 2 targets').add(
      node('simularium-spacefill', 'particles-representation', 'shapes', 'Spacefill', 'Particle Entity', false),
    );
    roots.set([picks, trajectory, measurements()]);
  } else {
    // An EMDB entry: the map (top level, like a structure) and the model fitted into it.
    const map = volume('emd-2984', 'EMD-2984', 'β-galactosidase · 2.2 Å', [
      ['Source', 'EMDB · CCP4'],
      ['Grid', '320 × 320 × 320'],
      ['Voxel', '0.64 Å'],
      ['Mean · σ', '0.000 · 0.012'],
      ['Range', '−0.081 … 0.214'],
    ]).add(volumeRepr('emd-2984-iso', undefined, { absolute: 0.035 }, 0x3362b2));
    const model = structure(
      '5a1a',
      '5A1A',
      'Assembly 1 · 4 chains',
      {
        method: 'em',
        symmetry: 'D2',
        facts: [
          ['Method', 'EM · 2.2 Å'],
          ['Chains', '4'],
          ['Residues', '4,084'],
        ],
      },
      contents('5a1a', [
        {
          id: '1',
          label: 'Beta-galactosidase',
          type: 'polymer',
          chains: ['A', 'B', 'C', 'D'].map((id) => ({ id, description: '1,021 residues' })),
        },
      ]),
      component('5a1a-polymer', 'Polymer', 'Chains A–D', repr('5a1a-polymer-cartoon', 'cartoon', 'chain-id')),
    );
    roots.set([model, map]);
  }

  const picked = observableState<SceneNode | undefined>(undefined);
  const inspections = new Map<SceneNode, SceneInspection>();
  return {
    roots: roots.asReadonly,
    picked: picked.asReadonly,
    pick: (n) => picked.set(n === picked.value ? undefined : n),
    inspect: (n) => {
      let inspection = inspections.get(n);
      if (!inspection) {
        inspection = {
          node: n,
          form: forms.get(n.id),
          slot: `inspector:${n.kind}`,
          facts: facts.get(n.id),
          // Measurements are added from the selection: a way there.
          actions:
            n.kind === 'measurements'
              ? [{ id: 'measure', label: 'Measure the selection…', icon: 'ruler', run: () => options.onMeasure?.() }]
              : undefined,
        };
        inspections.set(n, inspection);
      }
      return inspection;
    },
    structureInfo: (n) => infos.get(n.id),
    rowForm: (n, rowId) => rowForms.get(`${n.id}-${rowId}`),
    measurements: measurementList.asReadonly,
    measurementStyle: (n) => forms.get(n.id),
    addMeasurement: (label, description, icon) => {
      let owner = measurementsNode;
      if (!owner) {
        // No measurements yet: the row appears with the first.
        owner = measurements();
        owner.drop(owner.children.value[0]);
        measurementCount = 0;
        roots.set([...roots.value, owner]);
      }
      measurementCount++;
      owner.add(measurementNode(`m${measurementCount}-${Date.now()}`, icon, label, description));
      measurementsDescription(owner);
      if (!owner.isExpanded.value) owner.toggleExpanded();
    },
    streaming: (n) => rowForms.get(`${n.id}-volume-streaming`) as ParamsFormModel<any>,
    clear: () => {
      picked.set(undefined);
      roots.set([]);
      measurementsNode = undefined;
      measurementList.set([]);
    },
  };
}

/** The scene with a node picked (the inspector open). */
export function mockSceneWithPicked(id: string, options?: Parameters<typeof mockScene>[0]): MockScene {
  const scene = mockScene(options);
  const find = (nodes: readonly SceneNode[]): SceneNode | undefined => {
    for (const n of nodes) {
      if (n.id === id) return n;
      const c = find(n.children.value);
      if (c) return c;
    }
    return undefined;
  };
  const n = find(scene.roots.value);
  if (n) scene.pick(n);
  return scene;
}

/** Mol*'s quick styles: representation presets and visual styles; applying reports what would happen. */
export function mockQuickStyles(post: (text: string) => void = () => {}): QuickStylesModel {
  return {
    representations: [
      { id: 'default', label: 'Default' },
      { id: 'cartoon', label: 'Cartoon' },
      { id: 'spacefill', label: 'Spacefill' },
      { id: 'surface', label: 'Surface' },
    ],
    styles: [
      { id: 'default', label: 'Default' },
      { id: 'illustrative', label: 'Illustrative' },
    ],
    applyRepresentation: (id) => post(`Representation preset: ${id}`),
    applyStyle: (id) => post(`Style: ${id}`),
  };
}

/** Mol*'s representation presets, grouped as its registry groups them, with their common options. */
export function mockPresets(post: (text: string) => void = () => {}): PresetsModel {
  const current = observableState<string | undefined>('auto');
  return {
    groups: [
      {
        presets: [
          { id: 'empty', label: 'Empty', description: 'No representations' },
          { id: 'auto', label: 'Automatic', description: 'Picks a preset by the size of the structure' },
        ],
      },
      {
        label: 'Basic',
        presets: [
          { id: 'atomic-detail', label: 'Atomic detail', description: 'Ball & stick for everything' },
          { id: 'polymer-cartoon', label: 'Polymer cartoon' },
          { id: 'polymer-and-ligand', label: 'Polymer & ligand' },
          { id: 'protein-and-nucleic', label: 'Protein & nucleic' },
          { id: 'coarse-surface', label: 'Coarse surface' },
        ],
      },
      {
        label: 'Miscellaneous',
        presets: [
          { id: 'illustrative', label: 'Illustrative' },
          { id: 'molecular-surface', label: 'Molecular surface' },
          { id: 'auto-detail', label: 'Automatic detail' },
        ],
      },
    ],
    current: current.asReadonly,
    apply: (id) => {
      current.set(id);
      post(`Preset: ${id}`);
    },
  };
}

/** Mol*'s component options (StructureComponentManager.OptionsParams), shaped for the form. */
export function mockComponentOptions(): ParamsFormModel<any> {
  return new ParamsFormModel({
    hydrogens: PD.Select(
      'all',
      [
        ['all', 'Show all'],
        ['hide-all', 'Hide all'],
        ['polar', 'Only polar'],
      ],
      { label: 'Hydrogens', isEssential: true },
    ),
    visualQuality: PD.Select(
      'auto',
      [
        ['auto', 'Auto'],
        ['medium', 'Medium'],
        ['high', 'High'],
        ['highest', 'Highest'],
      ],
      { label: 'Visual quality', isEssential: true },
    ),
    ignoreLight: PD.Boolean(false, { label: 'Ignore light', isEssential: true }),
    materialStyle: PD.Group(
      {
        metalness: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
        roughness: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
        bumpiness: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
      },
      { label: 'Material style', isEssential: true },
    ),
    clipObjects: PD.Group(
      {
        variant: PD.Select('instance', [
          ['instance', 'Instance'],
          ['pixel', 'Pixel'],
        ]),
        objects: PD.Select('none', [
          ['none', 'None'],
          ['plane', 'Plane'],
          ['sphere', 'Sphere'],
          ['cube', 'Cube'],
        ]),
      },
      { label: 'Clip objects', isEssential: true },
    ),
    interactions: PD.Group(
      {
        hydrogenBonds: PD.Boolean(true, { label: 'Hydrogen bonds' }),
        ionic: PD.Boolean(true, { label: 'Ionic' }),
        piStacking: PD.Boolean(true, { label: 'π-stacking' }),
        halogenBonds: PD.Boolean(true, { label: 'Halogen bonds' }),
      },
      { label: 'Non-covalent interactions', isEssential: true },
    ),
    interior: PD.Group(
      {
        color: PD.Color(0x8c8c8c as any),
        colorStrength: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Color strength' }),
      },
      { label: 'Interior', isEssential: true },
    ),
    animation: PD.Group(
      {
        wiggle: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
        tumble: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
      },
      { label: 'Animation', isEssential: true },
    ),
  });
}
