/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ParamsFormModel } from '../../forms/params-form.js';
import type { HomeAction } from '../../model/home.js';
import type { StateTreeInspection, StateTreeModel, StateTreeNode } from '../../model/state-tree.js';
import { observableState } from '../../reactive/observable-state.js';

type Spec = Omit<StateTreeNode, 'children'> & { children?: Spec[] };

// 1TQN as Mol*'s default preset leaves it.
const Tree: Spec[] = [
  {
    ref: 'download',
    label: '1tqn.bcif',
    description: 'PDBe · BinaryCIF',
    icon: 'download',
    transform: 'Download',
    children: [
      {
        ref: 'cif',
        label: 'CIF File',
        description: '1 block',
        icon: 'file-code',
        transform: 'Parse CIF',
        children: [
          {
            ref: 'trajectory',
            label: 'Trajectory',
            description: '1 model',
            icon: 'trajectory',
            transform: 'Trajectory from mmCIF',
            children: [
              {
                ref: 'model',
                label: 'Model 1',
                description: '3,999 atoms',
                icon: 'box',
                transform: 'Model from Trajectory',
                decorators: ['Custom Model Properties'],
                children: [
                  {
                    ref: 'structure',
                    label: 'Assembly 1',
                    description: '3,999 elements',
                    icon: 'structure',
                    transform: 'Structure from Model',
                    decorators: ['Custom Structure Properties'],
                    visible: true,
                    children: [
                      {
                        ref: 'polymer',
                        label: 'Polymer',
                        description: '3,766 elements',
                        icon: 'layers',
                        transform: 'Structure Component',
                        visible: true,
                        children: [
                          {
                            ref: 'polymer-cartoon',
                            label: 'Cartoon',
                            description: 'Chain ID',
                            icon: 'cartoon',
                            transform: '3D Representation',
                            visible: true,
                          },
                        ],
                      },
                      {
                        ref: 'ligand',
                        label: 'Ligand',
                        description: '43 elements',
                        icon: 'layers',
                        transform: 'Structure Component',
                        visible: true,
                        children: [
                          {
                            ref: 'ligand-bs',
                            label: 'Ball & Stick',
                            description: 'Element',
                            icon: 'ball-and-stick',
                            transform: '3D Representation',
                            visible: true,
                          },
                        ],
                      },
                      {
                        ref: 'water',
                        label: 'Water',
                        description: '186 elements',
                        icon: 'layers',
                        transform: 'Structure Component',
                        visible: false,
                        children: [
                          {
                            ref: 'water-bs',
                            label: 'Ball & Stick',
                            description: 'Element',
                            icon: 'ball-and-stick',
                            transform: '3D Representation',
                            visible: false,
                          },
                        ],
                      },
                    ],
                  },
                ],
              },
            ],
          },
        ],
      },
    ],
  },
  {
    ref: 'density',
    label: 'emd_2984.map',
    icon: 'download',
    transform: 'Download',
    status: 'error',
    error: 'Not found (404)',
  },
];

// Update forms by transform, shaped like Mol*'s transformer params (a few each).
const UpdateParams: Record<string, () => PD.Params> = {
  Download: () => ({
    url: PD.Text('https://www.ebi.ac.uk/pdbe/entry-files/download/1tqn.bcif', { label: 'URL' }),
    isBinary: PD.Boolean(true, { label: 'Binary' }),
    label: PD.Text('1tqn.bcif'),
  }),
  'Model from Trajectory': () => ({ modelIndex: PD.Numeric(0, { min: 0, max: 0, step: 1 }, { label: 'Model index' }) }),
  'Structure from Model': () => ({
    type: PD.MappedStatic(
      'assembly',
      {
        auto: PD.Group({}, { isFlat: true }),
        model: PD.Group({}, { isFlat: true }),
        assembly: PD.Group({ id: PD.Select('1', [['1', 'Assembly 1']]) }, { isFlat: true }),
        symmetry: PD.Group({ ijkMin: PD.Numeric(-1), ijkMax: PD.Numeric(1) }, { isFlat: true }),
      },
      {
        label: 'Type',
        options: [
          ['auto', 'Auto'],
          ['model', 'Model'],
          ['assembly', 'Assembly'],
          ['symmetry', 'Symmetry'],
        ],
      },
    ),
  }),
  'Structure Component': () => ({
    label: PD.Text('Polymer'),
    type: PD.Select('static', [
      ['static', 'Static'],
      ['expression', 'Expression'],
      ['bundle', 'Bundle'],
    ]),
    nullIfEmpty: PD.Boolean(true, { label: 'Empty is none' }),
  }),
  '3D Representation': () => ({
    type: PD.Select('cartoon', [
      ['cartoon', 'Cartoon'],
      ['ball-and-stick', 'Ball & Stick'],
      ['spacefill', 'Spacefill'],
      ['molecular-surface', 'Molecular Surface'],
    ]),
    colorTheme: PD.Select('chain-id', [
      ['chain-id', 'Chain ID'],
      ['element-symbol', 'Element'],
      ['uniform', 'Uniform'],
      ['secondary-structure', 'Secondary structure'],
    ]),
    sizeTheme: PD.Select('uniform', [
      ['uniform', 'Uniform'],
      ['physical', 'Physical'],
    ]),
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Opacity' }),
  }),
};

// Decorators' params (A → A transforms folded into their parent cell).
const DecoratorParams: Record<string, () => PD.Params> = {
  'Custom Model Properties': () => ({
    autoAttach: PD.Boolean(true, { label: 'Auto attach' }),
    validation: PD.Boolean(false, { label: 'wwPDB validation' }),
    sifts: PD.Boolean(true, { label: 'SIFTS mapping' }),
  }),
  'Custom Structure Properties': () => ({
    accessibleSurfaceArea: PD.Boolean(false, { label: 'Accessible surface area' }),
    interactions: PD.Boolean(false, { label: 'Interactions' }),
    secondaryStructure: PD.Boolean(true, { label: 'Secondary structure' }),
  }),
};

// Actions applicable by transform (Mol*'s "Apply Action").
const Actions: Record<string, (post: (t: string) => void) => HomeAction[]> = {
  'Structure from Model': (post) => [
    action(
      'add-component',
      'Add Component',
      'layers',
      { label: PD.Text('Component'), selection: PD.Text('ligand') },
      'Add',
      post,
    ),
    action(
      'download-density',
      'Volume Streaming',
      'volume',
      {
        method: PD.Select('em', [
          ['em', 'EM'],
          ['x-ray', 'X-ray'],
        ]),
      },
      'Enable',
      post,
    ),
  ],
  'Structure Component': (post) => [
    action(
      'add-representation',
      'Add Representation',
      'plus',
      {
        type: PD.Select('ball-and-stick', [
          ['cartoon', 'Cartoon'],
          ['ball-and-stick', 'Ball & Stick'],
          ['spacefill', 'Spacefill'],
        ]),
      },
      'Add',
      post,
    ),
  ],
  'Trajectory from mmCIF': (post) => [
    action('load-coords', 'Load Coordinates', 'trajectory', { url: PD.Text('', { label: 'URL' }) }, 'Load', post),
  ],
};

function update(id: string, transform: string, params: PD.Params, what: string, post: (t: string) => void): HomeAction {
  return {
    id,
    title: `Update ${transform}`,
    icon: 'pencil',
    form: new ParamsFormModel(params),
    applyLabel: 'Update',
    apply: () => post(`Updated ${what}`),
  };
}

function action(
  id: string,
  title: string,
  icon: HomeAction['icon'],
  params: PD.Params,
  applyLabel: string,
  post: (t: string) => void,
): HomeAction {
  return { id, title, icon, form: new ParamsFormModel(params), applyLabel, apply: () => post(`${title}: applied`) };
}

/** The state tree for the stories: 1TQN as loaded, plus a download that failed. */
export function mockStateTree(post: (text: string) => void): StateTreeModel {
  const hidden = observableState<ReadonlySet<string>>(new Set(['water', 'water-bs']));
  const removed = observableState<ReadonlySet<string>>(new Set());
  const build = (specs: Spec[]): StateTreeNode[] =>
    specs
      .filter((s) => !removed.value.has(s.ref))
      .map((s) => ({
        ...s,
        visible: s.visible === undefined ? undefined : !hidden.value.has(s.ref),
        children: s.children ? build(s.children) : undefined,
      }));
  const roots = observableState<readonly StateTreeNode[]>(build(Tree));
  const rebuild = () => roots.set(build(Tree));
  const selected = observableState<string | undefined>('structure');
  const find = (ref: string, list: readonly StateTreeNode[] = roots.value): StateTreeNode | undefined => {
    for (const n of list) {
      if (n.ref === ref) return n;
      const c = n.children && find(ref, n.children);
      if (c) return c;
    }
  };
  const cache = new Map<string, StateTreeInspection>();
  return {
    roots: roots.asReadonly,
    selected: selected.asReadonly,
    select: (ref) => selected.set(ref),
    inspect: (ref) => {
      const node = find(ref);
      if (!node) return undefined;
      let inspection = cache.get(ref);
      if (!inspection) {
        const params = UpdateParams[node.transform]?.();
        inspection = {
          node,
          updates: [
            ...(params ? [update('update', node.transform, params, node.label, post)] : []),
            ...(node.decorators ?? []).flatMap((d) => {
              const p = DecoratorParams[d]?.();
              return p ? [update(`update:${d}`, d, p, `${node.label} · ${d}`, post)] : [];
            }),
          ],
          actions: Actions[node.transform]?.(post) ?? [],
        };
        cache.set(ref, inspection);
      }
      return { ...inspection, node };
    },
    toggleVisible: (ref) => {
      const next = new Set(hidden.value);
      if (next.has(ref)) next.delete(ref);
      else next.add(ref);
      hidden.set(next);
      rebuild();
    },
    remove: (ref) => {
      post(`Removed ${find(ref)?.label}`);
      removed.set(new Set([...removed.value, ref]));
      if (selected.value === ref) selected.set(undefined);
      rebuild();
    },
  };
}
