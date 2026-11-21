/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The genome forms (_spec/README.md "Queries", Forms; "Addressing"): an assembly's derived sequence by name
// (`(sequence :name chr17)`, any alias read, the assembly's name printed), a locus as readers type it
// (`(locus :text "chr17:7,668,402-7,687,550")`: 1-based, inclusive, commas, the name longest first) and a derived
// block of a section of assemblies by its sequence's name (`(block :name chrM)`). Names go through the assembly's
// names rule; one it lacks is `seq.query.locus.name`.

import { type Node, printAddress, type Selected, type Value } from '@molstar/core/repr/tree';
import { checkParams, type Ctx, fail, prefixNodes, type Resolve, required, tick } from './form-kit.js';
import { type Location, location } from './location.js';

/**
 * A form's name, a word or a string. A bare number is refused, never printed back: the address reader has made it a
 * number and its text is lost (`:name 01` reads as 1), so a name taken from it could be another sequence's.
 */
function nameOf(form: Node, key: string): string {
  const v = required(form, key);
  if (typeof v === 'number')
    fail('seq.query.form.param', `A name for :${key} reads as the number ${v}`, { hint: `Quote it: :${key} "…"` });
  if (typeof v !== 'string') fail('seq.query.form.param', `A form of kind ${form.kind} takes a name for :${key}`);
  return v as string;
}

/** The assemblies a prefix node stands for: itself (an assembly), or a section's `of` (one or a list). */
function assembliesOf(m: Selected): string[] {
  if (m.node.kind !== 'section') return [printAddress(m.address)];
  const of = m.node.params?.of as Value | undefined;
  return (Array.isArray(of) ? of : [of]).filter((x): x is string => typeof x === 'string');
}

/**
 * The derived sequence a name is under one of `assemblies`; `pending` while an assembly's source has not opened (the
 * form waits: incomplete, never unknown); else undefined.
 */
function find(ctx: Ctx, assemblies: readonly string[], name: string) {
  let pending = false;
  for (const a of assemblies) {
    const s = ctx.model.assemblySequence?.(a, name);
    if (s === 'pending') pending = true;
    else if (s) return s;
  }
  if (pending) ctx.complete = false;
  return pending ? 'pending' : undefined;
}

type Found = Exclude<ReturnType<typeof find>, 'pending' | undefined>;

/** A name found, or the form's issue; undefined while its assembly is pending. */
function named(ctx: Ctx, assemblies: readonly string[], name: string, at: string): Found | undefined {
  const s = find(ctx, assemblies, name);
  return s === 'pending' ? undefined : (s ?? unknown(name, at));
}

const unknown = (name: string, address: string): never =>
  fail('seq.query.locus.name', `No sequence of the assembly is named ${name}`, { address });

/** `(sequence :name n)` under an assembly: its derived sequence, addressed by the form with its printed name. */
const sequence: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['name']);
  const name = nameOf(form, 'name');
  const out: Selected[] = [];
  for (const m of prefixNodes(prefix, 'sequence')) {
    tick(ctx);
    const s = named(ctx, [printAddress(m.address)], name, printAddress(m.address));
    if (!s) continue;
    const node: Node = { kind: 'sequence', params: { name: s.name } };
    out.push({ node, address: [...m.address, node] });
  }
  return { type: 'nodes', nodes: out };
};

/**
 * A locus's text read (_spec/README.md "Addressing"): the whole text a name first, else a known name before the last
 * `:` and its positions (`7668402` or `7,668,402-7,687,550`, 1-based and inclusive) as `[start, end)`; undefined when the
 * positions do not parse.
 */
export function parseLocus(
  text: string,
  known: (name: string) => boolean,
): { name: string; range?: readonly [number, number] } | undefined {
  const t = text.trim();
  if (known(t)) return { name: t };
  const colon = t.lastIndexOf(':');
  // the name before the colon unknown too: the whole text is the name said unknown, never a shorter one at a position
  if (colon <= 0 || !known(t.slice(0, colon))) return { name: t };
  const m = /^(\d[\d,]*)(?:\s*[-–]\s*(\d[\d,]*))?$/.exec(t.slice(colon + 1).trim());
  if (!m) return undefined;
  const a = Number(m[1].replaceAll(',', ''));
  const b = m[2] === undefined ? a : Number(m[2].replaceAll(',', ''));
  if (!(a >= 1) || b < a) return undefined;
  return { name: t.slice(0, colon), range: [a - 1, b] };
}

/** `(locus :text t)` under an assembly or a section of assemblies: positions on the sequence it names. */
const locus: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['text']);
  const text = nameOf(form, 'text');
  const parts: { system: string; ranges: number[] }[] = [];
  for (const m of prefixNodes(prefix, 'locus')) {
    tick(ctx);
    const at = printAddress(m.address);
    const assemblies = assembliesOf(m);
    const read = parseLocus(text, (n) => typeof find(ctx, assemblies, n) === 'object');
    if (!read) return fail('seq.query.position.invalid', `Not a locus: ${text}`, { address: at });
    const s = named(ctx, assemblies, read.name, at);
    if (!s) continue;
    parts.push({ system: s.address, ranges: [...rangeOn(ctx, s.address, read.range)] });
  }
  return { type: 'location', location: location({ parts }) as Location };
};

/** `(block :name n)` under a section of assemblies: the derived block's positions (its sequence, whole). */
const block: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['name']);
  const name = nameOf(form, 'name');
  const parts: { system: string; ranges: number[] }[] = [];
  for (const m of prefixNodes(prefix, 'block')) {
    tick(ctx);
    const at = printAddress(m.address);
    const assemblies = m.node.kind === 'section' ? assembliesOf(m) : [];
    if (!assemblies.length || !assemblies.every((a) => ctx.model.isAssembly?.(a)))
      return fail('seq.query.form.prefix', 'A block form is under a section of assemblies', { address: at });
    const s = named(ctx, assemblies, name, at);
    if (!s) continue;
    parts.push({ system: s.address, ranges: [...rangeOn(ctx, s.address)] });
  }
  return { type: 'location', location: location({ parts }) as Location };
};

/** A range on a system, clamped to its length when known; the whole system without one. */
function rangeOn(ctx: Ctx, system: string, range?: readonly [number, number]): readonly [number, number] {
  const n = ctx.model.system(system)?.length;
  if (!range) return [0, n ?? 0];
  return n === undefined ? range : [Math.min(range[0], n), Math.min(range[1], n)];
}

export const GenomeForms: Readonly<Record<'sequence' | 'locus' | 'block', Resolve>> = { sequence, locus, block };
