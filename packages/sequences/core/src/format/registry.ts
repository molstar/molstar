/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Registries (seqspec:format#extensions): entries are pure data, providers grouped by kind, nothing run at
// registration.

import type {
  AlignmentInput,
  Alphabet,
  Annotation,
  AnnotationInput,
  MappingInput,
  MoleculeInput,
  Severity,
  StringTable,
  TreeInput,
  UnitInput,
} from '@molstar/seq-core/model';
import type { AlignerInput, FormatRequest, FormatResult, StageContext } from './contract.js';
import type { Geometry, ThemeRef } from './doc.js';
import { KernelRegistry } from './kernel-registry.js';
import type { SchemeDef } from './palettes/schemes.js';
import type { CheckContext } from './schema/field-schema.js';
import { isIssueCode, isNamespace, isNamespaced, ReservedNamespaces } from './schema/field-schema.js';
import type { ParamsSchema } from './schema/params-schema.js';
import { MvsRegistry } from './tables/kernel.js';

export interface FormatEntry {
  readonly name: string;
  readonly extensions?: readonly string[];
  readonly gives: readonly FormatResult['kind'][];
  /** The request fields its parse depends on: one parse per source and distinct values of these. */
  readonly reads: readonly (keyof FormatRequest)[];
  read(bytes: Uint8Array, request: FormatRequest, ctx: StageContext): FormatResult | Promise<FormatResult>;
}
export interface ServiceEntry {
  readonly name: string;
  readonly resources: readonly string[];
  resolve(id: string, resource: string): { uri: string; format: string } | null;
  /** The read result made the service's; `ctx` takes its notes (a host may call it without). */
  after?(result: FormatResult, id: string, resource: string, ctx?: StageContext): FormatResult;
}
export interface ComplexProviderEntry {
  readonly name: string;
  readonly options: ParamsSchema;
  readonly formats: readonly string[];
  build(
    ctx: StageContext,
    source: FormatResult,
    options: unknown,
  ): Promise<{ molecules: MoleculeInput[]; units: UnitInput[] }>;
}
/** `align` gives a CIGAR from a to b; null: these cannot be aligned. */
export interface AlignerEntry {
  readonly name: string;
  readonly options: ParamsSchema;
  align(ctx: StageContext, a: AlignerInput, b: AlignerInput, options: unknown): Promise<string | null>;
  /** An MSA aligner (`realign`) has this too: the gapped rows of `inputs`, in their order; null: not aligned, the reason reported. */
  alignMany?(
    ctx: StageContext,
    inputs: readonly AlignerInput[],
    options: unknown,
    signal: AbortSignal,
  ): Promise<string[] | null>;
}
export interface ProviderEntry {
  readonly name: string;
  readonly options: ParamsSchema;
  readonly geometry: Geometry;
  /** A layer's color when all its data is this provider's and its `color` is null. */
  readonly defaultColor?: ThemeRef;
  compute(ctx: StageContext, options: unknown): Promise<AnnotationInput | null>;
}
/**
 * A node provider (seqspec:model#alignments-and-trees): computes a `tree`, `alignment` or `mapping` node from the nodes
 * it reads. `inputs` are those nodes' inputs: a tree's is its parent alignment's `AlignmentInput`; an alignment's is
 * its inputs' sequences (`RealignInput`); a mapping's is its two endpoint alignments (`AlignmentSummary`). `signal` is
 * the stage's: a long run honors it, and its throw is not an unavailable node. null: unavailable (the provider reports
 * why). A `mapping` provider may give several mappings (`codon_rows`: one per row pair).
 */
export interface NodeProviderEntry {
  readonly name: string;
  readonly kind: 'tree' | 'alignment' | 'mapping';
  readonly options: ParamsSchema;
  compute(
    ctx: StageContext,
    inputs: readonly unknown[],
    options: unknown,
    signal: AbortSignal,
  ): Promise<TreeInput | AlignmentInput | MappingInput | readonly MappingInput[] | null>;
}
/**
 * A field read from an annotation by a registered provider: its name (`<ns>:<name>`) is usable wherever a field is
 * (`field_palette`, `rules` measures). Read once per annotation, in the format half, so every viewer has it.
 */
export interface FieldProviderEntry {
  readonly name: string;
  readonly type: 'number' | 'string';
  applies(a: Annotation): boolean;
  /** One value per row into `out` (NaN, or -1 for a string: none); strings as codes interned in `strings`. */
  read(a: Annotation, out: Float64Array | Int32Array, strings: StringTable): void;
}
export interface LayerTypeEntry {
  readonly name: string;
  /** What it shows (the reference). */
  readonly description?: string;
  /** Its params: top-level layer params for a kernel type, the layer's `options` for an extension's. */
  readonly options: ParamsSchema;
  readonly accepts: readonly Geometry[] | 'none';
  readonly layouts: readonly ('tracks' | 'text')[];
  readonly axes: readonly ('blocks' | 'columns')[];
  readonly sizing: 'row' | 'free';
  /** A kernel type or `"omit"`; null: none (`letters`, `ranges`: every viewer draws them). */
  readonly fallback: string | null;
  /** What its layers take as channel params, `color` first (seqspec:format#layer-types). */
  readonly channels: readonly ChannelEntry[];
}
export type ChannelKind = 'color' | 'size' | 'glyph';
export interface ChannelEntry {
  readonly name: string;
  readonly kind: ChannelKind;
  /** A size channel's unit; null for the others. */
  readonly unit: 'px' | 'band' | null;
  /** null: the viewer's (the ink, a size). */
  readonly default: ThemeRef | null;
}
export interface ThemeEntry {
  readonly name: string;
  readonly description?: string;
  /** The channel kinds it serves, each its params table. */
  readonly params: { readonly [K in ChannelKind]?: ParamsSchema };
  readonly fallback: ThemeRef | null;
  /** A named table (seqspec:format#themes): the theme it stands for, so resolution hands viewers `rules` alone. */
  expand?(params: Readonly<Record<string, unknown>>, kind: ChannelKind): ThemeRef;
  /** Checks across its params, after each param's own (`scheme`: `name` or `custom`, `overrides` against it). */
  check?(params: Readonly<Record<string, unknown>>, path: string, ctx: CheckContext): void;
}
/** A named letter map: single letters to the letter shown, unlisted letters kept. */
export interface LetterMapEntry {
  readonly name: string;
  readonly map: Readonly<Record<string, string>>;
  /** The alphabets it is for: a layer showing letters of another keeps them as stored. Omitted: any. */
  readonly alphabets?: readonly Alphabet[];
}
/**
 * A named color scheme (seqspec:format#color-schemes): data a `scheme` theme names; what it fits is the scheme's own.
 */
export interface ColorSchemeEntry {
  readonly name: string;
  /** What a picker shows. */
  readonly label: string;
  readonly description: string;
  /** Where its colors come from. */
  readonly provenance: string;
  /** Designed or rated safe for the common color vision deficiencies. */
  readonly colorblind: boolean;
  readonly scheme: SchemeDef;
}
export interface CustomEntry {
  readonly namespace: string;
  readonly kind: string;
  readonly schema: ParamsSchema;
}
export interface IssueCodeEntry {
  readonly code: string;
  readonly severity: Severity;
}

export interface SvsRegistryEntry {
  readonly formats?: readonly FormatEntry[];
  readonly services?: readonly ServiceEntry[];
  readonly complexProviders?: readonly ComplexProviderEntry[];
  readonly providers?: readonly ProviderEntry[];
  readonly nodeProviders?: readonly NodeProviderEntry[];
  readonly fields?: readonly FieldProviderEntry[];
  readonly aligners?: readonly AlignerEntry[];
  readonly layerTypes?: readonly LayerTypeEntry[];
  readonly themes?: readonly ThemeEntry[];
  readonly letterMaps?: readonly LetterMapEntry[];
  readonly colorSchemes?: readonly ColorSchemeEntry[];
  readonly custom?: readonly CustomEntry[];
  readonly issueCodes?: readonly IssueCodeEntry[];
}

export interface RegistryEntries {
  formats: FormatEntry;
  services: ServiceEntry;
  complexProviders: ComplexProviderEntry;
  providers: ProviderEntry;
  nodeProviders: NodeProviderEntry;
  fields: FieldProviderEntry;
  aligners: AlignerEntry;
  layerTypes: LayerTypeEntry;
  themes: ThemeEntry;
  letterMaps: LetterMapEntry;
  colorSchemes: ColorSchemeEntry;
  custom: CustomEntry;
  issueCodes: IssueCodeEntry;
}
export type RegistryKind = keyof RegistryEntries;
export type EntryOf<K extends RegistryKind> = RegistryEntries[K];

export interface SvsRegistries {
  /** Bumped by every registration and undo. */
  readonly version: number;
  /** By name; a custom table by `<namespace>/<kind>`, an issue code by code. */
  get<K extends RegistryKind>(kind: K, name: string): EntryOf<K> | undefined;
  /** Registration order; cached until a change. */
  list<K extends RegistryKind>(kind: K): readonly EntryOf<K>[];
}

/** The interface's order: entries register kind by kind in it. */
const KINDS: readonly RegistryKind[] = [
  'formats',
  'services',
  'complexProviders',
  'providers',
  'nodeProviders',
  'fields',
  'aligners',
  'layerTypes',
  'themes',
  'letterMaps',
  'colorSchemes',
  'custom',
  'issueCodes',
];

function itemsOf(e: SvsRegistryEntry, kind: RegistryKind): readonly unknown[] {
  return e[kind] ?? [];
}

function nameOf(kind: RegistryKind, item: unknown): string {
  const e = item as Record<string, string>;
  if (kind === 'custom') return `${e.namespace}/${e.kind}`;
  return kind === 'issueCodes' ? e.code : e.name;
}

class Registries implements SvsRegistries {
  version = 0;
  readonly maps = new Map(KINDS.map((k) => [k, new Map<string, unknown>()]));
  private readonly cache = new Map<RegistryKind, readonly unknown[]>();
  get<K extends RegistryKind>(kind: K, name: string) {
    return this.maps.get(kind)?.get(name) as EntryOf<K> | undefined;
  }
  list<K extends RegistryKind>(kind: K) {
    let l = this.cache.get(kind);
    if (!l) {
      l = [...(this.maps.get(kind)?.values() ?? [])];
      this.cache.set(kind, l);
    }
    return l as readonly EntryOf<K>[];
  }
  changed() {
    this.version++;
    this.cache.clear();
  }
}

export function createRegistries(): SvsRegistries {
  return new Registries();
}

const FLAT: ReadonlySet<RegistryKind> = new Set(['formats', 'services']);
let kernelNames: Set<string> | undefined;

/** Why a name cannot be registered (namespaced, except the kernel's and formats' and services'); `owned`:
 * the reserved namespace this entry may use (the viewer's `seqv`). */
function invalidName(kind: RegistryKind, name: string, owned: string | null): string | null {
  kernelNames ??= new Set(
    [KernelRegistry, MvsRegistry].flatMap((e) =>
      KINDS.flatMap((k) => itemsOf(e, k).map((i) => `${k} ${nameOf(k, i)}`)),
    ),
  );
  if (FLAT.has(kind) || kernelNames.has(`${kind} ${name}`)) return null;
  const head = (kind === 'custom' ? name.split('/') : kind === 'issueCodes' ? name.split('.') : name.split(':'))[0];
  if (ReservedNamespaces.has(head) && head !== owned)
    return `${kind}: "${kind === 'issueCodes' ? name : head}" is reserved`;
  if (kind === 'custom') return isNamespace(head) ? null : `custom: "${head}" is not a namespace`;
  const ok = kind === 'issueCodes' ? isIssueCode(name) : isNamespaced(name);
  return ok ? null : `${kind}: "${name}" is not namespaced (${kind === 'issueCodes' ? '<ns>.<code>' : '<ns>:<name>'})`;
}

/** `register`'s error: every conflict, nothing registered; the message lists them one per line. */
export class SvsRegistrationError extends Error {
  readonly conflicts: readonly string[];
  constructor(conflicts: readonly string[]) {
    super(`register: ${conflicts.length} conflict(s), nothing registered:\n- ${conflicts.join('\n- ')}`);
    this.name = 'SvsRegistrationError';
    this.conflicts = conflicts;
  }
}

export interface RegisterOptions {
  /**
   * The viewer's own entries (the viewer owns `seqv`): names under `seqv` are allowed in these entry objects, and only
   * in them. Only the viewer passes it; an extension's entry naming `seqv` is still a conflict.
   */
  readonly seqv?: ReadonlySet<SvsRegistryEntry>;
}

/**
 * Registers entries in the order given, kinds in the interface's order. Every conflict (a name another registration
 * holds, a name listed twice with different entries, a name not namespaced) throws one error listing them all, and
 * nothing is registered; the same entry object twice is not a conflict. Returns an idempotent undo.
 */
export function register(
  reg: SvsRegistries,
  entry: SvsRegistryEntry | readonly SvsRegistryEntry[],
  opts: RegisterOptions = {},
): () => void {
  if (!(reg instanceof Registries)) throw new Error('register: registries not made by createRegistries()');
  const entries: readonly SvsRegistryEntry[] = Array.isArray(entry) ? entry : [entry as SvsRegistryEntry];
  const added: [Map<string, unknown>, string, unknown][] = [];
  const pending = new Map<string, unknown>();
  const conflicts: string[] = [];
  for (const e of entries) {
    const owned = opts.seqv?.has(e) ? 'seqv' : null;
    for (const kind of KINDS) {
      for (const item of itemsOf(e, kind)) {
        const name = nameOf(kind, item);
        const map = reg.maps.get(kind) as Map<string, unknown>;
        const held = map.get(name) ?? pending.get(`${kind} ${name}`);
        if (held === item) continue;
        if (held !== undefined) {
          conflicts.push(`${kind}: "${name}" is ${map.has(name) ? 'already registered' : 'listed twice'}`);
          continue;
        }
        const invalid = invalidName(kind, name, owned);
        if (invalid) conflicts.push(invalid);
        pending.set(`${kind} ${name}`, item);
        added.push([map, name, item]);
      }
    }
  }
  if (conflicts.length > 0) {
    throw new SvsRegistrationError(conflicts);
  }
  for (const [map, name, item] of added) map.set(name, item);
  reg.changed();
  let undone = false;
  return () => {
    if (undone) return;
    undone = true;
    for (const [map, name, item] of added.reverse()) if (map.get(name) === item) map.delete(name);
    reg.changed();
  };
}

/** The longest of `known` (and the reserved) namespaces an issue code starts with: `org.example.bad` is `org.example`'s. */
export function namespaceOfCode(known: ReadonlySet<string>, code: string): string | null {
  for (let i = code.lastIndexOf('.'); i > 0; i = code.lastIndexOf('.', i - 1)) {
    const ns = code.slice(0, i);
    if (known.has(ns) || ReservedNamespaces.has(ns)) return ns;
  }
  return null;
}

/** Every namespace some registration owns (`requires`); a code owns none by itself. */
export function registeredNamespaces(reg: SvsRegistries): Set<string> {
  const out = new Set<string>();
  for (const kind of KINDS) {
    if (kind === 'issueCodes') continue;
    for (const e of reg.list(kind)) {
      const name = nameOf(kind, e);
      if (kind === 'custom') out.add(name.split('/')[0]);
      else if (name.includes(':')) out.add(name.split(':')[0]);
    }
  }
  const named = new Set(out);
  for (const c of reg.list('issueCodes')) {
    const ns = namespaceOfCode(named, c.code);
    if (ns) out.add(ns);
  }
  return out;
}
