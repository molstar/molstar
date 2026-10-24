/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { DocNode, SeqViewSpecDoc } from '@molstar/seq-core/format';
import type { StructureFormat } from './molstar/source.js';
import { LettersTrack, SecondaryData, SecondaryTrack } from './tracks.js';

// The document the sequence model builds (sequence.md §3): one view with one section over a `molstar:structure`
// complex whose file the host already parsed; with no blocks the section takes every polymer unit in order, then the
// ligands, else the chain filter's blocks. Its tracks are the track providers' (`tracks.ts`), put in by the model. The
// theme (§10) is written into the letters' layer: Seq*'s themes draw the colors.

/** The complex's ref in the document, and the name the link maps the 3D structure onto. */
export const ComplexRef = 'structure';
/** The section's ref. */
export const SectionRef = 'main';
/** The data node's ref: behaviors add their data under it (sequence.md §3). */
export const DataRef = 'data';
/** The letters track's ref. */
export const LettersRef = LettersTrack;
/** The tracks' refs by what they show (the track menus are by these: `trackActions`). */
export const SecondaryRef = SecondaryData;
export { SecondaryTrack };

/** The host file the document names its structure by: the parsed one, which the loaders' `parsed` hook returns. */
export const hostFileOf = (key: string) => `structure:${key}`;

/** The color scheme the letters use, by its registered name (`KernelColorSchemes`); null: not colored. */
export interface SequenceTheme {
  readonly letters: string | null;
}

/**
 * The letters' color: a scheme, or for no colors one quiet color (`plain`, the surface a shade toward the text): a
 * neutral strip, the letters in the text's contrast. SVS has no "no color" for a layer (a null color is the type's default, which
 * colors by residue); text mode's plain letters are the render's (`wrappedColors` off).
 */
export function lettersColor(theme: SequenceTheme, plain: string): Record<string, unknown> {
  return theme.letters === null
    ? { theme: 'uniform', params: { color: plain } }
    : { theme: 'scheme', params: { name: theme.letters } };
}
export const DefaultSequenceTheme: SequenceTheme = { letters: 'prism:residues' };

/** A track and its data as a provider builds them, for `sequenceDocument` and `insertTracks`. */
export interface BuiltTrack {
  readonly data?: readonly DocNode[];
  readonly track: DocNode;
}

/** A block's ref for a unit: a ref's characters only (`seqspec:format#nodes-and-refs`). */
const blockRef = (unit: string) => `block.${unit.replace(/[^A-Za-z0-9_.-]/g, '_')}`;

/**
 * A section's blocks: one per unit, by its name in the complex (the chain filter, sequence.md §3), each with a ref of
 * its own so a patch can take it out.
 */
export const blocksOf = (units: readonly string[]): DocNode[] =>
  units.length
    ? units.map((u) => ({ kind: 'block', ref: blockRef(u), params: { unit: `${ComplexRef}/${u}` } }))
    : // No chain: one block of a unit there is not, left out quietly (no blocks at all would be every unit).
      [{ kind: 'block', ref: 'block.none', params: { unit: `${ComplexRef}/-` } }];

export function sequenceDocument(o: {
  title: string;
  key: string;
  format: StructureFormat;
  /** Only the ligands of these chemical components become units (the ones the 3D view draws); absent: all. */
  ligands?: readonly string[];
  /** The units the section shows (the chain filter); absent: Seq*'s default, every polymer unit then the ligands. */
  units?: readonly string[];
  /** Every copy a block of its own (an assembly's); absent: copies merged (Seq*'s `one`). */
  each?: boolean;
  /** The tracks in, in order (the providers', `tracks.ts`). One document for both styles: Text draws the ruler nowhere,
   * secondary structure as underlines and values (a pLDDT track) as strips (Seq*'s lanes); Tracks draws each as a track. */
  tracks: readonly BuiltTrack[];
}): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: o.title });
  b.data.complex({
    ref: ComplexRef,
    provider: 'molstar:structure',
    source: { host_file: hostFileOf(o.key), format: o.format },
    ...(o.ligands ? { options: { ligand_comps: [...o.ligands] } } : {}),
  });
  b.view({ ref: 'sequence', label: o.title }).section({
    ref: SectionRef,
    style: 'compact',
    ...(o.each ? { copies: 'each' } : {}),
  });
  // Behaviors add under the data node by its ref (the builder gives it none).
  const doc = b.getDocument();
  const data = doc.root.children?.find((n) => n.kind === 'data');
  if (data) data.ref = DataRef;
  const section = findRef(doc.root, SectionRef)!;
  section.children = [...(o.units ? blocksOf(o.units) : []), ...(section.children ?? [])];
  if (data) insertTracks(data, section, o.tracks);
  return doc;
}

/** Puts tracks into a document's data node and section (at its end, or after `after`). */
export function insertTracks(data: DocNode, section: DocNode, tracks: readonly BuiltTrack[], after?: string): void {
  for (const t of tracks) data.children = [...(data.children ?? []), ...(t.data ?? [])];
  const list = [...(section.children ?? [])];
  const at = after ? list.findIndex((n) => n.ref === after) : -1;
  list.splice(at < 0 ? list.length : at + 1, 0, ...tracks.map((t) => t.track));
  section.children = list;
}

/** The node with this ref under `node`. */
export function findRef(node: DocNode, ref: string): DocNode | undefined {
  if (node.ref === ref) return node;
  for (const c of node.children ?? []) {
    const n = findRef(c, ref);
    if (n) return n;
  }
  return undefined;
}

/** A document's first `molstar:structure` complex: its ref and its source's params (by `source_ref`, or a source child). */
export interface SvsStructureSource {
  readonly complex: string;
  readonly source: Readonly<Record<string, unknown>>;
}

function walk(node: DocNode, visit: (n: DocNode, parent: DocNode | undefined) => void, parent?: DocNode) {
  visit(node, parent);
  for (const c of node.children ?? []) walk(c, visit, node);
}

function sourceNodeOf(doc: SeqViewSpecDoc, complex: DocNode): DocNode | undefined {
  const ref = complex.params?.source_ref;
  let found: DocNode | undefined;
  if (typeof ref === 'string') walk(doc.root, (n) => (found ??= n.kind === 'source' && n.ref === ref ? n : undefined));
  return found ?? complex.children?.find((c) => c.kind === 'source');
}

function firstStructureComplex(doc: SeqViewSpecDoc): { complex: DocNode; parent: DocNode | undefined } | undefined {
  let found: { complex: DocNode; parent: DocNode | undefined } | undefined;
  walk(doc.root, (n, parent) => {
    if (!found && n.kind === 'complex' && n.params?.provider === 'molstar:structure') found = { complex: n, parent };
  });
  return found;
}

/** What a loaded SeqViewSpec maps onto a structure (sequence.md §6): its first `molstar:structure` complex; none: undefined. */
export function svsStructureSource(doc: SeqViewSpecDoc): SvsStructureSource | undefined {
  const found = firstStructureComplex(doc);
  const source = found && sourceNodeOf(doc, found.complex);
  if (!found?.complex.ref || !source?.params) return undefined;
  return { complex: found.complex.ref, source: source.params };
}

/** A loaded SeqViewSpec as the sequence shows it over the followed structure: the document and its refs. */
export interface SvsShown {
  readonly doc: SeqViewSpecDoc;
  readonly complex: string;
  /** The complex's data node (behaviors add theirs under it). */
  readonly data: string;
  /** Every section (the person's style applies to each); the first with tracks is where behaviors add theirs. */
  readonly sections: readonly string[];
  readonly section: string;
  readonly after: string;
}

/**
 * A loaded SeqViewSpec over the structure the host parsed (sequence.md §6, David 2026-10-15): a copy whose first
 * `molstar:structure` complex reads the host's parse (its source becomes the host file the loaders' `parsed` answers),
 * with the person's letter colors over its own; nodes the viewer needs by ref (the data node, sections, the last track)
 * get one. Undefined: the document has no such complex.
 */
export function svsDocument(
  svs: SeqViewSpecDoc,
  o: { key: string; format: StructureFormat; theme: SequenceTheme; plain: string },
): SvsShown | undefined {
  const doc = structuredClone(svs);
  const found = firstStructureComplex(doc);
  const source = found && sourceNodeOf(doc, found.complex);
  if (!found?.complex.ref || !source) return undefined;
  source.params = { host_file: hostFileOf(o.key), format: o.format };
  const taken = new Set<string>();
  walk(doc.root, (n) => n.ref && taken.add(n.ref));
  const refFor = (n: DocNode, base: string) => {
    if (n.ref) return n.ref;
    let ref = base;
    for (let i = 2; taken.has(ref); i++) ref = `${base}-${i}`;
    taken.add(ref);
    n.ref = ref;
    return ref;
  };
  const data = found.parent && found.parent.kind === 'data' ? refFor(found.parent, DataRef) : undefined;
  const sections: string[] = [];
  let section: string | undefined;
  let after: string | undefined;
  walk(doc.root, (n) => {
    if (n.kind === 'section') {
      sections.push(refFor(n, 'section'));
      const tracks = n.children?.filter((c) => c.kind === 'track') ?? [];
      if (!section && tracks.length) {
        section = n.ref;
        after = refFor(tracks[tracks.length - 1], 'track');
      }
    }
    // The person's colors over the document's (sequence.md §6).
    if (n.kind === 'layer' && n.params?.type === 'letters') {
      const { color: _, ...rest } = n.params;
      n.params = { ...rest, color: lettersColor(o.theme, o.plain) };
    }
  });
  if (!data || !section || !after) return undefined;
  return { doc, complex: found.complex.ref, data, sections, section, after };
}
