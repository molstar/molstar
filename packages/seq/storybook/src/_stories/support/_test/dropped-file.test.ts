/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer } from '@molstar/seqx-viewer/viewer';
import { describe, expect, it } from 'vitest';
import { droppedUri, DropRuntimes, formatOfName, specOfDropped } from '../dropped-file.js';

// A file dropped on a view, for each format built: named, opened by its runtime, its spec shown by a headless viewer whose
// assets port answers the dropped bytes only (nothing is read from anywhere else), laid out with no issue and rendered.

const bytes = (file: string) => new Uint8Array(readFileSync(`${PackageFolder}${file}`));

async function show(name: string, data: Uint8Array) {
  const dropped = await specOfDropped(name, data);
  if ('refused' in dropped)
    return { dropped, reads: [] as string[], issues: [] as string[], sections: [] as string[], svg: '' };
  const base = DefaultSeqViewerSpec();
  const reads: string[] = [];
  const viewer = createSeqViewer({
    ...base,
    ui: undefined,
    registry: [
      ...(base.registry ?? []),
      ...DropRuntimes.map((runtime) => ({ kind: 'format-runtime' as const, runtime })),
    ],
    ports: {
      assets: () => ({
        reader: {
          readText: async (uri) => {
            reads.push(uri);
            return new TextDecoder().decode(data);
          },
          readBytes: async (uri) => {
            reads.push(uri);
            return data;
          },
        },
      }),
    },
  });
  viewer.resize(1000);
  await viewer.load(dropped.spec);
  await viewer.settle();
  const issues = [...viewer.modelIssues.value, ...viewer.layoutIssues.value].map((i) => `${i.code} ${i.address ?? ''}`);
  const sections = viewer.view.resolved.value?.sections.map((s) => s.address) ?? [];
  const svg = await viewer.renderSvg({});
  viewer.retire();
  return { dropped, reads, issues, sections, svg: svg.svg };
}

describe('a dropped file', () => {
  it.each([
    ['2OCJ.seq', 'fixtures/p53/2OCJ.seq', 'seq'],
    ['fab.fasta', 'fixtures/antibody/1N8Z.fasta', 'fasta'],
    ['PF00042.sto', 'fixtures/globins/PF00042-seed.sto', 'stockholm'],
    ['domains.tsv', 'fixtures/files/domains.tsv', 'tsv'],
    ['domains.csv', 'fixtures/files/domains.csv', 'csv'],
    ['domains.json', 'fixtures/files/domains-aos.json', 'json'],
    ['2ocj.cif.gz', 'fixtures/structures/2OCJ.cif.gz', 'mmcif'],
  ])('%s is read as %s: shown with no issue, read from the dropped bytes only', async (name, file, format) => {
    const r = await show(name, bytes(file));
    expect('refused' in r.dropped, JSON.stringify(r.dropped)).toBe(false);
    if ('refused' in r.dropped) return;
    expect(r.dropped.format).toBe(format);
    expect(r.issues).toEqual([]);
    expect(r.sections.length).toBeGreaterThan(0);
    expect(r.svg.startsWith('<svg')).toBe(true);
    // every location the viewer read is the dropped file's
    expect([...new Set(r.reads)]).toEqual([droppedUri(name)]);
  });

  it('GenBank is refused, saying so; an unknown extension names what is read; a file that does not parse is refused with the reason', async () => {
    expect(formatOfName('pUC19.gb')).toEqual({ refused: 'GenBank files are not read yet' });
    expect(
      'refused' in formatOfName('notes.txt') && (formatOfName('notes.txt') as { refused: string }).refused,
    ).toContain('.seq');
    const bad = await specOfDropped('broken.seq', new TextEncoder().encode('this is not a seq file\n\u0000'));
    expect('refused' in bad && bad.refused).toContain('broken.seq');
    const empty = await specOfDropped('empty.tsv', new TextEncoder().encode('a\tb\n'));
    expect('refused' in empty).toBe(true);
  });

  it('a name with characters of a URL is one location', () => {
    expect(droppedUri('my file #1.seq')).toBe('dropped/my%20file%20%231.seq');
  });
});
