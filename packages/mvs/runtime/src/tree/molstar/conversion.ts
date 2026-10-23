/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { omitObjectKeys, pickObjectKeys } from '@molstar/core/util/object';
import { canonicalJsonString } from '@molstar/mvs-builder/json';
import { addParamDefaults } from '@molstar/mvs-builder/tree/generic/params-schema';
import {
  type ConversionRules,
  addDefaults,
  condenseTree,
  convertTree,
  dfs,
  resolveUri,
  resolveUris,
  windowUrl,
} from '@molstar/mvs-builder/tree/generic/tree-utils';
import { type FullMVSTree, type MVSTree, MVSTreeSchema } from '@molstar/mvs-builder/tree/mvs/mvs-tree';
import { MolstarParseFormatT, ParseFormatT } from '@molstar/mvs-builder/tree/mvs/param-types';
import type { MolstarKind, MolstarNode, MolstarSubtree, MolstarTree } from './molstar-tree.js';

/** The structure's keys its shown model reads (`molstar_periodic`, load.md §6); none: no custom. */
function modelCustom(custom: { [key: string]: unknown } | undefined): { [key: string]: unknown } | undefined {
  return custom && 'molstar_periodic' in custom ? { molstar_periodic: custom.molstar_periodic } : undefined;
}

/** Convert `format` parameter of `parse` node in `MolstarTree`
 * into `format` and `is_binary` parameters in `MolstarTree` */
export const ParseFormatMvsToMolstar = {
  // trajectory
  mmcif: { format: 'cif', is_binary: false },
  bcif: { format: 'cif', is_binary: true },
  pdb: { format: 'pdb', is_binary: false },
  pdbqt: { format: 'pdbqt', is_binary: false },
  gro: { format: 'gro', is_binary: false },
  xyz: { format: 'xyz', is_binary: false },
  mol: { format: 'mol', is_binary: false },
  sdf: { format: 'sdf', is_binary: false },
  mol2: { format: 'mol2', is_binary: false },
  lammpstrj: { format: 'lammpstrj', is_binary: false },
  // coordinates
  xtc: { format: 'xtc', is_binary: true },
  nctraj: { format: 'nctraj', is_binary: true },
  dcd: { format: 'dcd', is_binary: true },
  trr: { format: 'trr', is_binary: true },
  // topology
  psf: { format: 'psf', is_binary: false },
  prmtop: { format: 'prmtop', is_binary: false },
  top: { format: 'top', is_binary: false },
  // maps
  map: { format: 'map', is_binary: true },
  dx: { format: 'dx', is_binary: false },
  dxbin: { format: 'dxbin', is_binary: true },
  // shape
  vtp: { format: 'vtp', is_binary: true }, // XML, but ParseVtp takes Uint8Array and decodes the header itself
  ply: { format: 'ply', is_binary: true }, // parsePly auto-detects ascii vs binary from the header
  obj: { format: 'obj', is_binary: false },
} satisfies { [p in ParseFormatT]: { format: MolstarParseFormatT; is_binary: boolean } };

const TopologyFormats = new Set<ParseFormatT>(['psf', 'prmtop', 'top']);

/** Conversion rules for conversion from `MVSTree` (with all parameter values) to `MolstarTree` */
const mvsToMolstarConversionRules: ConversionRules<FullMVSTree, MolstarTree> = {
  download: (node) => ({ subtree: [] }),
  parse: (node, parent) => {
    const { format, is_binary } = ParseFormatMvsToMolstar[node.params.format];
    if (parent?.kind === 'download') {
      return {
        subtree: [
          { kind: 'download', params: { ...parent.params, is_binary }, custom: parent.custom, ref: parent.ref },
          { kind: 'parse', params: { ...node.params, format }, custom: node.custom, ref: node.ref },
        ] satisfies MolstarNode[],
      };
    } else {
      console.warn('"parse" node is not being converted, this is suspicious');
      return {
        subtree: [
          { kind: 'parse', params: { ...node.params, format }, custom: node.custom, ref: node.ref },
        ] satisfies MolstarNode[],
      };
    }
  },
  coordinates: (node, parent) => {
    if (parent?.kind !== 'parse') throw new Error(`Parent of "coordinates" must be "parse", not "${parent?.kind}".`);
    const { format } = ParseFormatMvsToMolstar[parent.params.format];
    return {
      subtree: [
        { kind: 'coordinates', params: { format }, custom: node.custom, ref: node.ref },
      ] satisfies MolstarNode[],
    };
  },
  structure: (node, parent) => {
    if (parent?.kind !== 'parse') throw new Error(`Parent of "structure" must be "parse", not "${parent?.kind}".`);
    const { format } = ParseFormatMvsToMolstar[parent.params.format];

    if (TopologyFormats.has(parent.params.format)) {
      if (!node.params.coordinates_ref) {
        throw new Error(
          `"structure" node with topology format "${parent.params.format}" must have "coordinates_ref" parameter.`,
        );
      }
      return {
        subtree: [
          { kind: 'topology_with_coordinates', params: { format, coordinates_ref: node.params.coordinates_ref } },
          { kind: 'model', params: pickObjectKeys(node.params, ['model_index']), custom: modelCustom(node.custom) },
          {
            kind: 'structure',
            params: omitObjectKeys(node.params, ['block_header', 'block_index', 'model_index', 'coordinates_ref']),
            custom: node.custom,
            ref: node.ref,
          },
        ] satisfies MolstarNode[],
      };
    } else if (node.params.coordinates_ref) {
      return {
        subtree: [
          { kind: 'trajectory', params: { format, ...pickObjectKeys(node.params, ['block_header', 'block_index']) } },
          { kind: 'model', params: { model_index: 0 } },
          { kind: 'trajectory_with_coordinates', params: { coordinates_ref: node.params.coordinates_ref } },
          { kind: 'model', params: pickObjectKeys(node.params, ['model_index']), custom: modelCustom(node.custom) },
          {
            kind: 'structure',
            params: omitObjectKeys(node.params, ['block_header', 'block_index', 'model_index', 'coordinates_ref']),
            custom: node.custom,
            ref: node.ref,
          },
        ] satisfies MolstarNode[],
      };
    } else {
      return {
        subtree: [
          { kind: 'trajectory', params: { format, ...pickObjectKeys(node.params, ['block_header', 'block_index']) } },
          { kind: 'model', params: pickObjectKeys(node.params, ['model_index']), custom: modelCustom(node.custom) },
          {
            kind: 'structure',
            params: omitObjectKeys(node.params, ['block_header', 'block_index', 'model_index', 'coordinates_ref']),
            custom: node.custom,
            ref: node.ref,
          },
        ] satisfies MolstarNode[],
      };
    }
  },
};

/** Node kinds in `MolstarTree` that it makes sense to condense */
const molstarNodesToCondense = new Set<MolstarKind>([
  'download',
  'parse',
  'trajectory',
  'model',
] satisfies MolstarKind[]);

/** Convert MolViewSpec tree into MolStar tree */
export function convertMvsToMolstar(mvsTree: MVSTree, sourceUrl: string | undefined): MolstarTree {
  const full = addDefaults<typeof MVSTreeSchema>(mvsTree, MVSTreeSchema) as FullMVSTree;
  if (sourceUrl) resolveUris(full, sourceUrl, ['uri', 'url']);
  if (sourceUrl) resolveCustomUrls(full, sourceUrl);
  const converted = convertTree<FullMVSTree, MolstarTree>(full, mvsToMolstarConversionRules);
  if (converted.kind !== 'root')
    throw new Error("Root's type is not 'root' after conversion from MVS tree to Molstar tree.");
  const condensed = condenseTree<MolstarTree>(converted, molstarNodesToCondense);
  return condensed;
}

/**
 * Mol*-specific keys with addresses resolved as `url` params are (load.md §6): `custom.molstar_boxes.url` on a
 * `coordinates` node. A new `custom` object: the caller's tree is not changed.
 */
function resolveCustomUrls(tree: MVSTree, sourceUrl: string) {
  dfs(tree, (node) => {
    const boxes = (node.custom as { molstar_boxes?: { url?: unknown } } | undefined)?.molstar_boxes;
    if (!boxes || typeof boxes.url !== 'string') return;
    node.custom = { ...node.custom, molstar_boxes: { ...boxes, url: resolveUri(boxes.url, sourceUrl, windowUrl()) } };
  });
}

/** Kinds whose children are converted by the rules (split, condensed); below them a node converts on its own, and so
 * do the root's children other than `download`. */
const TopParentKinds = new Set<string>(['root', 'download', 'parse']);

/**
 * Converted subtrees by MVS node, across the frames of an animation (`load.md §4`). A node whose parent is not `root`,
 * `download` or `parse`, or a child of `root` other than `download`, converts to a copy of itself with defaults and
 * resolved URIs, whatever its parent and siblings (no rule and no condensing applies there), so an MVS subtree that is
 * the same object converts to the same object.
 */
export interface MolstarConversionMemo {
  readonly sourceUrl: string | undefined;
  readonly subtrees: WeakMap<object, MolstarSubtree>;
  /** Converted params and custom by the MVS object: a node converted again because a descendant changed keeps them */
  readonly values: WeakMap<object, object>;
  /** Converted params of nodes without params, by kind */
  readonly defaultParams: Map<string, object>;
  /** Converted `download` subtrees under the root, by MVS node, with the units inside them */
  readonly downloads: WeakMap<object, { node: MolstarSubtree; units: MolstarSubtree[] }>;
}

export function createMolstarConversionMemo(sourceUrl: string | undefined): MolstarConversionMemo {
  return {
    sourceUrl,
    subtrees: new WeakMap(),
    values: new WeakMap(),
    defaultParams: new Map(),
    downloads: new WeakMap(),
  };
}

const UnitToken = '\u0000mvs-unit:';

/**
 * `convertMvsToMolstar`, with the subtrees below the top nodes (`root`, `download`, `parse`, the children of `download`
 * and `parse`) taken from `memo` when the MVS subtree is the same object. Returns the tree and its "units": the
 * converted subtrees right below the top nodes, whose identity says that their MVS source did not change.
 */
export function convertMvsToMolstarMemo(
  mvsTree: MVSTree,
  memo: MolstarConversionMemo,
): { tree: MolstarTree; units: Set<MolstarSubtree> } {
  const children = (mvsTree.children ?? []) as any[];
  // each `download` under the root converts on its own unless it has a twin, which the rules condense with it
  const downloadKeys = new Set<string>();
  let twins = false;
  for (const child of children) {
    if (child.kind !== 'download') continue;
    const key = canonicalJsonString(addParamDefaults(MVSTreeSchema.nodes.download.params, (child as any).params));
    if (downloadKeys.has(key)) twins = true;
    downloadKeys.add(key);
  }
  if (twins) return convertTopPart(mvsTree, memo);

  const units = new Set<MolstarSubtree>();
  const converted: MolstarSubtree[] = [];
  for (const child of children) {
    if (child.kind !== 'download') {
      const unit = convertUnit(child, memo);
      units.add(unit);
      converted.push(unit);
      continue;
    }
    let download = memo.downloads.get(child);
    if (!download) {
      const own = convertTopPart({ kind: 'root', children: [child] } as MVSTree, memo);
      // (a download that is not one Molstar download, e.g. none without a `parse`, two for two formats: the whole
      // top part, as `convertMvsToMolstar` makes it)
      if (own.tree.children?.length !== 1) return convertTopPart(mvsTree, memo);
      download = { node: own.tree.children[0] as MolstarSubtree, units: Array.from(own.units) };
      memo.downloads.set(child, download);
    }
    units.add(download.node);
    for (const u of download.units) units.add(u);
    converted.push(download.node);
  }
  const tree = convertMvsToMolstar({ ...mvsTree, children: [] } as MVSTree, memo.sourceUrl);
  tree.children = converted as MolstarTree['children'];
  return { tree, units };
}

/** The top part through the rules (`convertMvsToMolstar`), with each unit as a leaf placeholder (its ref is a token). */
function convertTopPart(
  mvsTree: MVSTree,
  memo: MolstarConversionMemo,
): { tree: MolstarTree; units: Set<MolstarSubtree> } {
  const units: MVSTree[] = [];
  const top = (node: MVSTree): MVSTree => {
    if (!node.children) return node;
    const isTopParent = TopParentKinds.has(node.kind);
    const children = node.children.map((child) => {
      // (a child of the root other than `download` has no rule and no twin to condense with: a unit too)
      if (isTopParent && (node.kind !== 'root' || child.kind === 'download')) return top(child as MVSTree);
      units.push(child as MVSTree);
      // (a kind with no schema, rule or condensing: the conversion only copies it)
      return { kind: UnitToken, ref: UnitToken + (units.length - 1) } as any;
    });
    return { ...node, children } as MVSTree;
  };
  const tree = convertMvsToMolstar(top(mvsTree), memo.sourceUrl);

  const converted = new Set<MolstarSubtree>();
  const replace = (node: MolstarSubtree) => {
    const children = node.children as MolstarSubtree[] | undefined;
    if (!children) return;
    for (let i = 0; i < children.length; i++) {
      const ref = children[i].ref;
      if (typeof ref === 'string' && ref.startsWith(UnitToken)) {
        const unit = convertUnit(units[+ref.substring(UnitToken.length)], memo);
        children[i] = unit;
        converted.add(unit);
      } else {
        replace(children[i]);
      }
    }
  };
  replace(tree);
  return { tree, units: converted };
}

function memoized<T extends object>(memo: MolstarConversionMemo, key: object, compute: () => T): T {
  let value = memo.values.get(key) as T | undefined;
  if (!value) memo.values.set(key, (value = compute()));
  return value;
}

/** The params of a unit node with defaults and resolved URIs; the same object for the same MVS params object. */
function convertUnitParams(kind: string, mvsParams: any, memo: MolstarConversionMemo): any {
  const schema = (MVSTreeSchema.nodes as any)[kind] as { params: any } | undefined;
  if (!schema && !mvsParams) return undefined;
  const compute = () => {
    const params: any = schema ? (addParamDefaults as any)(schema.params, mvsParams) : { ...mvsParams };
    if (memo.sourceUrl) {
      for (const name in params) {
        if (name !== 'uri' && name !== 'url') continue;
        const uri = params[name];
        if (typeof uri === 'string') params[name] = resolveUri(uri, memo.sourceUrl, windowUrl());
      }
    }
    return params;
  };
  if (mvsParams) return memoized(memo, mvsParams, compute);
  let params = memo.defaultParams.get(kind);
  if (!params) memo.defaultParams.set(kind, (params = compute()));
  return params;
}

/** What `convertMvsToMolstar` makes of a node below the top nodes: defaults, URIs, copies, children always an array. */
function convertUnit(node: MVSTree, memo: MolstarConversionMemo): MolstarSubtree {
  const cached = memo.subtrees.get(node);
  if (cached) return cached;
  const n = node as { kind: string; params?: any; custom?: any; ref?: string; children?: MVSTree[] };
  const result = {
    kind: n.kind,
    params: convertUnitParams(n.kind, n.params, memo),
    custom: n.custom ? memoized(memo, n.custom, () => ({ ...n.custom })) : undefined,
    ref: n.ref,
    children: n.children ? n.children.map((c) => convertUnit(c, memo)) : [],
  } as unknown as MolstarSubtree;
  memo.subtrees.set(node, result);
  return result;
}

type FileExtension = `.${Lowercase<string>}`;
function fileExtensionMatches(filename: string, extensions: (FileExtension | '*')[]): boolean {
  filename = filename.toLowerCase();
  return extensions.some((ext) => ext === '*' || filename.endsWith(ext));
}

const StructureFormatExtensions: Record<ParseFormatT, (FileExtension | '*')[]> = {
  // trajectory
  mmcif: ['.cif', '.mmif'],
  bcif: ['.bcif'],
  pdb: ['.pdb', '.ent'],
  pdbqt: ['.pdbqt'],
  gro: ['.gro'],
  xyz: ['.xyz'],
  mol: ['.mol'],
  sdf: ['.sdf'],
  mol2: ['.mol2'],
  lammpstrj: ['.lammpstrj'],
  // coordinates
  xtc: ['.xtc'],
  nctraj: ['.nc', '.nctraj'],
  dcd: ['.dcd'],
  trr: ['.trr'],
  // topology
  psf: ['.psf'],
  prmtop: ['.prmtop', '.parm7'],
  top: ['.top'],
  // volumes
  map: ['.map', '.ccp4', '.mrc', '.mrcs'],
  dx: ['.dx'],
  dxbin: ['.dxbin'],
  // shape
  vtp: ['.vtp'],
  ply: ['.ply'],
  obj: ['.obj'],
};

/** Run some sanity check on a MVSTree. Return a list of potential problems (`undefined` if there are none) */
export function mvsSanityCheckIssues(tree: MVSTree): string[] | undefined {
  const result: string[] = [];
  dfs(tree, (node, parent) => {
    if (node.kind === 'parse' && parent?.kind === 'download') {
      const source = parent.params.url;
      const extensions = StructureFormatExtensions[node.params.format];
      if (!fileExtensionMatches(source, extensions)) {
        result.push(
          `Parsing data from ${source} as ${node.params.format} format might be a mistake. The file extension doesn't match recommended file extensions (${extensions.join(', ')})`,
        );
      }
    }
  });
  return result.length > 0 ? result : undefined;
}

/** Run some sanity check on a MVSTree and print potential issues to the console. */
export function mvsSanityCheck(tree: MVSTree): void {
  const issues = mvsSanityCheckIssues(tree);
  if (issues) {
    console.warn('There are potential issues in the MVS tree:');
    for (const issue of issues) {
      console.warn(' ', issue);
    }
  }
}
