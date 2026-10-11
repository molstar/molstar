/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { IconName } from '../../components/icon/icon.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import type { SceneInspection, SceneModel, SceneNode, SceneNodeKind } from '../../model/scene.js';
import { type ObservableState, observableState } from '../../reactive/observable-state.js';

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

const StructureParams = {
  type: PD.Select('assembly', [
    ['model', 'Model'],
    ['assembly', 'Assembly'],
    ['symmetry', 'Symmetry mates'],
  ]),
  assembly: PD.Select('1', [
    ['1', 'Assembly 1 (monomer)'],
    ['2', 'Assembly 2'],
  ]),
  model: PD.Numeric(1, { min: 1, max: 20, step: 1 }, { label: 'Model index' }),
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

  constructor(
    readonly id: string,
    readonly kind: SceneNodeKind,
    public icon: IconName,
    label: string,
    description?: string,
    expanded = true,
    readonly onHighlight?: (node: MockNode, on: boolean) => void,
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

  add(...children: MockNode[]): this {
    for (const c of children) c.parent = this;
    this.#children.set([...this.#children.value, ...children]);
    return this;
  }

  setDescription(d: string | undefined) {
    this.#description.set(d);
  }
  setLabel(l: string) {
    this.#label.set(l);
  }

  toggleVisibility() {
    this.#hidden.set(!this.#hidden.value);
  }
  toggleExpanded() {
    this.#expanded.set(!this.#expanded.value);
  }
  highlight(on: boolean) {
    this.onHighlight?.(this, on);
  }
  focus() {}
  remove() {
    const parent = this.parent;
    if (parent) parent.#children.set(parent.#children.value.filter((c) => c !== this));
    else removeRoot?.(this);
  }
}

let removeRoot: ((node: MockNode) => void) | undefined;

/**
 * A scene like 1TQN loaded with the default preset: the structure, its components (protein, heme, waters, ions) with
 * their representations, and a distance measurement. Forms change the nodes as the plugin would.
 */
export function mockScene(options: { onHighlight?: (label: string | undefined) => void } = {}): SceneModel {
  const highlight = (node: MockNode, on: boolean) => options.onHighlight?.(on ? node.label.value : undefined);
  const forms = new Map<string, ParamsFormModel<any>>();

  const repr = (id: string, type: string, color: string, alpha = 1) => {
    const node = new MockNode(id, 'representation', ReprIcons[type] ?? 'cartoon', '', undefined, false, highlight);
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
      node.setLabel(RepresentationTypes.find(([t]) => t === v.type)?.[1] ?? v.type);
      node.icon = ReprIcons[v.type] ?? 'cartoon';
      node.setDescription(
        `${ColorLabels[v.colorTheme.name] ?? v.colorTheme.name}${v.alpha < 1 ? ` · ${Math.round(v.alpha * 100)}%` : ''}`,
      );
    };
    sync(values);
    forms.set(id, new ParamsFormModel(RepresentationParams, values, { onChange: sync }));
    return node;
  };
  const component = (id: string, label: string, description: string, ...reprs: MockNode[]) => {
    const node = new MockNode(id, 'component', 'layers', label, description, true, highlight).add(...reprs);
    forms.set(
      id,
      new ParamsFormModel(ComponentParams, { label }, { onChange: (v) => node.setLabel(v.label || label) }),
    );
    return node;
  };

  const structure = new MockNode('1tqn', 'structure', 'atom', '1TQN', 'Assembly 1 · 3,624 atoms', true, highlight).add(
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
  forms.set('1tqn', new ParamsFormModel(StructureParams));
  const water = structure.children.value[2] as MockNode;
  water.toggleVisibility();
  water.toggleExpanded();

  const measurements = new MockNode(
    'measurements',
    'measurements',
    'ruler',
    'Measurements',
    '1 distance',
    true,
    highlight,
  ).add(new MockNode('d1', 'measurement', 'measure-distance', 'CYS 442 SG – HEM FE', '2.31 Å', false, highlight));

  const roots = observableState<readonly SceneNode[]>([structure, measurements]);
  removeRoot = (node) => roots.set(roots.value.filter((n) => n !== node));
  const picked = observableState<SceneNode | undefined>(undefined);
  const inspections = new Map<SceneNode, SceneInspection>();
  return {
    roots: roots.asReadonly,
    picked: picked.asReadonly,
    pick: (node) => picked.set(node === picked.value ? undefined : node),
    inspect: (node) => {
      let inspection = inspections.get(node);
      if (!inspection) {
        inspection = { node, form: forms.get(node.id), slot: `inspector:${node.kind}` };
        inspections.set(node, inspection);
      }
      return inspection;
    },
  };
}

/** The scene with a node picked (the inspector open). */
export function mockSceneWithPicked(id: string, options?: Parameters<typeof mockScene>[0]): SceneModel {
  const scene = mockScene(options);
  const find = (nodes: readonly SceneNode[]): SceneNode | undefined => {
    for (const n of nodes) {
      if (n.id === id) return n;
      const c = find(n.children.value);
      if (c) return c;
    }
    return undefined;
  };
  const node = find(scene.roots.value);
  if (node) scene.pick(node);
  return scene;
}
