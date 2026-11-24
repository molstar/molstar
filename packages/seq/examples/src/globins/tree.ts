/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Globins seed's tree: InterPro serves no tree of the Pfam seed, so it is computed here and committed as Newick
// text (`SeedTree`), compared with this generator over the seed's fixture by a test. The views write it inline: a tree
// node takes its Newick as text.

/** The distance of two aligned rows: 1 less their identity over the columns both have a residue in (1 when none). */
export function rowDistance(a: string, b: string): number {
  let both = 0;
  let same = 0;
  for (let j = 0; j < a.length; j++) {
    const x = a.charCodeAt(j) | 32;
    const y = b.charCodeAt(j) | 32;
    // `-` and `.` are gaps (45, 46; unchanged by the case bit)
    if (x === 45 || x === 46 || y === 45 || y === 46) continue;
    both++;
    if (x === y) same++;
  }
  return both ? 1 - same / both : 1;
}

const length = (x: number) => Math.max(0, x).toFixed(4);

/**
 * Saitou and Nei's neighbor joining over the rows' distances: each step joins the pair of least Q (the first such pair
 * in row order), the last three nodes the root's children. Labels as given (none may hold `(),:;` or a space).
 */
export function neighborJoining(labels: readonly string[], letters: readonly string[]): string {
  const n = labels.length;
  if (n < 3) return `(${labels.join(',')});`;
  const d: number[][] = labels.map((_, i) => labels.map((_, j) => (i === j ? 0 : rowDistance(letters[i], letters[j]))));
  let nodes = labels.map((l) => l);
  let dist = d;
  while (nodes.length > 3) {
    const m = nodes.length;
    const r = dist.map((row) => row.reduce((s, x) => s + x, 0));
    let bi = 0;
    let bj = 1;
    let best = Infinity;
    for (let i = 0; i < m; i++)
      for (let j = i + 1; j < m; j++) {
        const q = (m - 2) * dist[i][j] - r[i] - r[j];
        if (q < best - 1e-12) {
          best = q;
          bi = i;
          bj = j;
        }
      }
    const li = dist[bi][bj] / 2 + (r[bi] - r[bj]) / (2 * (m - 2));
    const lj = dist[bi][bj] - li;
    const joined = `(${nodes[bi]}:${length(li)},${nodes[bj]}:${length(lj)})`;
    const keep = nodes.map((_, k) => k).filter((k) => k !== bi && k !== bj);
    const du = keep.map((k) => (dist[bi][k] + dist[bj][k] - dist[bi][bj]) / 2);
    dist = [...keep.map((k, a) => [...keep.map((q) => dist[k][q]), du[a]]), [...du, 0]];
    nodes = [...keep.map((k) => nodes[k]), joined];
  }
  const [a, b, c] = [0, 1, 2];
  const la = (dist[a][b] + dist[a][c] - dist[b][c]) / 2;
  const lb = dist[a][b] - la;
  const lc = dist[a][c] - la;
  return `(${nodes[a]}:${length(la)},${nodes[b]}:${length(lb)},${nodes[c]}:${length(lc)});`;
}
