/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The overlay's molecules (_spec/README.md "Overlay", Molecules and complexes): a source's sequences and molecules
// as systems, its complexes' units as systems on them with their steps, links held on per-row units, compounds by
// code. A unit and its links pass the language's checks (`checkUnit`, `checkLinks`) against what they name, as a spec's
// do; one that fails is left out with its issues. A complex past `maxUnits` units is refused whole.

import {
  CategoryColumn,
  type Issue,
  issueOf,
  type Node,
  parseAddress,
  printAddress,
  type Value,
} from '@molstar/core/repr/tree';
import { checkLinks, checkUnit, LangIssues } from '@molstar/seqx-core/lang';
import { labelsNumbering, type Numbering } from '@molstar/seqx-core/query';
import { type AnnotationEntry, type ResolveConfig, resolveSequenceNode, type SystemEntry } from './model.js';
import { compoundEntry, linksEntry, unitsIssue } from './model-molecules.js';
import { ResolveIssues } from './issues.js';
import type { CompoundEntry, Molecular } from './molecule-query.js';
import {
  type AtomsEntry,
  type ComplexEntry,
  complexEntry,
  type ModifiedIndex,
  type MoleculeEntry,
  modifiedIndex,
  moleculeEntry,
  type UnitEntry,
  unitEdges,
  unitEntry,
} from './molecules.js';
import type { Edge } from './paths.js';

/** An address as its canonical spelling (forms too), or undefined. */
function canonical(v: Value | undefined): string | undefined {
  if (typeof v !== 'string') return undefined;
  const p = parseAddress(v);
  if (p.issues.length || !p.pattern || p.pattern.subtree || !p.pattern.segments.length) return undefined;
  return printAddress(p.pattern.segments);
}

/** A child's address under a node's (a unit under its complex). */
function childAddress(parent: string, id: string): string {
  return printAddress([...(parseAddress(parent).pattern?.segments ?? []), id]);
}

/** The molecule entries of an overlay's nodes, built as the nodes are checked (sequences and molecules first). */
export class DynamicMolecular implements Molecular {
  readonly molecules = new Map<string, MoleculeEntry>();
  readonly complexes = new Map<string, ComplexEntry>();
  readonly units = new Map<string, UnitEntry>();
  readonly atoms = new Map<string, AtomsEntry>();
  readonly modified = new Map<string, ModifiedIndex>();
  readonly nodes = new Map<string, Node>();
  readonly numberings = new Map<string, Map<string, Numbering>>();
  readonly compounds = new Map<string, CompoundEntry>();
  /** The systems they make: sequences, molecules' residues and atoms, units. */
  readonly systems = new Map<string, SystemEntry>();
  /** Complexes' and molecules' links, held whole. */
  readonly annotations = new Map<string, AnnotationEntry>();
  /** The units' steps up and down. */
  readonly edges: Edge[] = [];

  constructor(
    private readonly config: ResolveConfig,
    /** A system of the static model or the overlay, by address. */
    private readonly systemOf: (address: string) => SystemEntry | undefined,
    private readonly push: (issue: Issue) => void,
  ) {}

  /** A dynamic node of a kind that makes systems or is reference data; false when it is left out. */
  add(address: string, node: Node): boolean {
    switch (node.kind) {
      case 'sequence': {
        const source = canonical(node.params?.source);
        this.systems.set(address, resolveSequenceNode(address, node, source));
        const m = modifiedIndex(node);
        if (m) this.modified.set(address, m);
        return true;
      }
      case 'molecule': {
        const e = moleculeEntry(address, node, canonical(node.params?.source));
        this.molecules.set(address, e.molecule);
        this.nodes.set(address, node);
        for (const s of e.systems) this.systems.set(s.address, s);
        if (e.molecule.atoms) this.atoms.set(e.molecule.atoms.address, e.molecule.atoms);
        this.links(address, node, undefined);
        return true;
      }
      case 'compounds':
        for (const c of node.children ?? [])
          if (c.kind === 'compound' && c.id !== undefined) this.compounds.set(c.id, compoundEntry(c));
        return true;
      case 'complex':
        return this.complex(address, node);
      default:
        return true;
    }
  }

  /** A language issue at a node, under its params. */
  private say(address: string, ok: { value: boolean }) {
    return (code: Parameters<Parameters<typeof checkUnit>[2]>[0], message: string, ...tail: (string | number)[]) => {
      const issue = issueOf(LangIssues, code, message, { address, path: ['params', ...tail] });
      if (issue.severity === 'error') ok.value = false;
      this.push(issue);
    };
  }

  private complex(address: string, node: Node): boolean {
    const children = (node.children ?? []).filter((c) => c.kind === 'unit' && c.id !== undefined);
    if (children.length > this.config.maxUnits) {
      this.push(unitsIssue(address, children.length, this.config.maxUnits));
      return false;
    }
    const list: UnitEntry[] = [];
    for (const u of children) {
      const at = childAddress(address, u.id!);
      const of = canonical(u.params?.of);
      const system = of === undefined ? undefined : this.systemOf(of);
      if (!system || (system.kind !== 'sequence' && system.kind !== 'molecule')) {
        this.push(
          issueOf(ResolveIssues, 'seq.resolve.node.ref', 'The unit’s `of` names no sequence or molecule held', {
            address: at,
          }),
        );
        continue;
      }
      const ok = { value: true };
      checkUnit(u, system.length, this.say(at, ok));
      if (!ok.value) continue;
      const e = unitEntry(at, address, u, { system, molecule: this.molecules.get(system.address) });
      list.push(e.unit);
      this.units.set(at, e.unit);
      this.nodes.set(at, u);
      for (const s of e.systems) this.systems.set(s.address, s);
      if (e.atoms) this.atoms.set(e.atoms.address, e.atoms);
      this.edges.push(...unitEdges(e.unit));
      for (const c of u.children ?? []) {
        const labels = c.params?.labels;
        if (c.kind !== 'numbering' || c.id === undefined || !(labels instanceof CategoryColumn)) continue;
        let own = this.numberings.get(at);
        if (!own) this.numberings.set(at, (own = new Map()));
        own.set(c.id, labelsNumbering(c.id, labels));
      }
    }
    const links = this.links(address, node, list);
    this.complexes.set(address, complexEntry(address, node, list, links));
    return true;
  }

  /** A complex's or a molecule's links, checked against its units (or residues) and held; their addresses. */
  private links(address: string, node: Node, units: readonly UnitEntry[] | undefined): string[] {
    const out: string[] = [];
    for (const c of node.children ?? []) {
      if (c.kind !== 'pairs' || c.id === undefined) continue;
      const at = childAddress(address, c.id);
      const ok = { value: true };
      const byId = new Map(units?.map((u) => [u.id, u]));
      checkLinks(
        c,
        units
          ? {
              kind: 'complex',
              unit: (id) => {
                const u = byId.get(id);
                return u && { length: u.end === undefined ? undefined : u.end - u.start };
              },
            }
          : { kind: 'molecule', length: this.systems.get(address)?.length },
        this.say(at, ok),
      );
      if (!ok.value || (units && !units.length)) continue;
      this.annotations.set(at, linksEntry(at, c, address, units));
      out.push(at);
    }
    return out;
  }
}
