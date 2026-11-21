/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { afterEach, describe, expect, it, vi } from 'vitest';

// Genome resolve loaded when used (_spec/README.md "Resolve", Assemblies, Loaded when used): a spec holding an
// assembly is pending until the genome chunk lands, then resolves; a spec without one never imports it; and no module
// the eager entries import reaches the chunk but by a dynamic import. Each case on a fresh module registry
// (`vi.resetModules`: this file runs isolated), so the chunk's state is the case's own.

const here = dirname(fileURLToPath(import.meta.url));
const src = resolve(here, '../..');

/** The chunk's modules: imported by `import()` only. */
const Chunk = [
  'resolve/assemblies.ts',
  'resolve/assembly-names.ts',
  'resolve/genome-model.ts',
  'resolve/genome-chunk.ts',
  'query/genome-forms.ts',
];

async function fresh() {
  vi.resetModules();
  const lang = await import('@molstar/seqx-core/lang');
  const resolveEntry = await import('@molstar/seqx-core/resolve');
  return { lang, resolve: resolveEntry };
}

afterEach(() => {
  vi.useRealTimers();
});

describe('the genome chunk', () => {
  it('a spec without an assembly never imports it', async () => {
    const { lang, resolve: r } = await fresh();
    const b = lang.createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.views.view('v', (v) => v.section('main', { of: s }));
    const m = r.resolveSpec(b.build());
    expect(m.loading).toBeUndefined();
    expect(r.genomeLoaded()).toBeUndefined();
  });

  it('a spec holding an assembly is pending, then resolved once the chunk lands, on a fake clock (EC-112, EC-169)', async () => {
    vi.useFakeTimers();
    const { lang, resolve: r } = await fresh();
    const b = lang.createSeqSpec();
    b.data.assembly(
      'g',
      lang.assembly('toy', [
        { name: 'chr1', length: 10 },
        { name: 'chr2', length: 4 },
      ]),
    );
    const root = b.build();
    const resolver = r.createResolver();
    const before = resolver.resolve(root);
    // pending, never wrong: no sequence derived, no issue, the model says what it waits for
    expect(before.assemblies.get('/data/g')?.pending).toBe(true);
    expect(before.assemblies.get('/data/g')?.systems).toEqual([]);
    expect(before.issues).toEqual([]);
    expect(before.loading).toBeDefined();
    // a resolve before it lands is pending again (never the cached pending model taken for the settled one)
    expect(resolver.resolve(root).loading).toBeDefined();
    await before.loading;
    const after = resolver.resolve(root);
    expect(after.loading).toBeUndefined();
    expect(after.assemblies.get('/data/g')?.systems).toEqual([
      '/data/g/(sequence :name chr1)',
      '/data/g/(sequence :name chr2)',
    ]);
    expect(after.hash).not.toBe(before.hash);
  });

  it('a query naming a locus before the chunk lands waits, then resolves', async () => {
    const { lang, resolve: r } = await fresh();
    const q = await import('@molstar/seqx-core/query');
    const b = lang.createSeqSpec();
    b.data.assembly('g', lang.assembly('toy', [{ name: 'chr1', length: 10 }]));
    const root = b.build();
    const at = (m: ReturnType<typeof r.resolveSpec>) =>
      q.resolveAddressIn(root, r.overlayQueryModel(m, r.EmptyOverlay), '/data/g/(locus :text "chr1:2-3")');
    const before = r.resolveSpec(root);
    expect(at(before).status).toBe('waiting');
    await before.loading;
    const after = at(r.resolveSpec(root));
    expect(after.status).toBe('resolved');
  });

  it('is reached by no static import of the eager entries', () => {
    const seen = new Set<string>();
    const visit = (file: string) => {
      if (seen.has(file)) return;
      seen.add(file);
      const text = readFileSync(file, 'utf8');
      // static imports and re-exports (a type-only import pulls no code)
      for (const m of text.matchAll(/^(?:import|export)\s+(?!type\b)[^;]*?from\s+'([^']+)'/gms)) {
        const spec = m[1];
        const target = spec.startsWith('.')
          ? join(dirname(file), spec.replace(/\.js$/, '.ts'))
          : spec.startsWith('@molstar/seqx-core/')
            ? entryOf(spec.slice('@molstar/seqx-core/'.length))
            : undefined;
        if (target) visit(target);
      }
    };
    for (const e of ['index.ts', 'lang/index.ts', 'resolve/index.ts', 'query/index.ts', 'view/index.ts'])
      visit(join(src, e));
    const reached = Chunk.filter((c) => seen.has(join(src, c)));
    expect(reached).toEqual([]);
    // and the walk did read the doors
    expect(seen.has(join(src, 'resolve/genome.ts'))).toBe(true);
    expect(seen.has(join(src, 'query/genome-door.ts'))).toBe(true);
  });
});

/** A package entry's source file (`resolve` → `resolve/index.ts`, `data/channels` → `data/channels.ts`). */
function entryOf(entry: string): string {
  try {
    readFileSync(join(src, entry, 'index.ts'));
    return join(src, entry, 'index.ts');
  } catch {
    return join(src, `${entry}.ts`);
  }
}
