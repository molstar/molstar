/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Small meshes made in code (a few KB, no fixture files): a torus as OBJ, an octahedron with colored corners as PLY. For
// the Shapes group's stories and tests: what a mesh file loads as (`shape.*`, `look.set` by node).

/** A torus as OBJ: `rings` × `sides` quads (as triangles) around `center`, in Å. */
export function torusObj(
  options: {
    center?: readonly [number, number, number];
    radius?: number;
    tube?: number;
    rings?: number;
    sides?: number;
  } = {},
): string {
  const [cx, cy, cz] = options.center ?? [0, 0, 0];
  const R = options.radius ?? 20;
  const r = options.tube ?? 3;
  const rings = options.rings ?? 32;
  const sides = options.sides ?? 12;
  const lines = ['# a torus made in code (testing/meshes.ts)'];
  for (let i = 0; i < rings; i++) {
    const u = (2 * Math.PI * i) / rings;
    for (let j = 0; j < sides; j++) {
      const v = (2 * Math.PI * j) / sides;
      const x = (R + r * Math.cos(v)) * Math.cos(u);
      const y = (R + r * Math.cos(v)) * Math.sin(u);
      const z = r * Math.sin(v);
      lines.push(`v ${(cx + x).toFixed(3)} ${(cy + y).toFixed(3)} ${(cz + z).toFixed(3)}`);
    }
  }
  // OBJ indices start at 1
  const at = (i: number, j: number) => (i % rings) * sides + (j % sides) + 1;
  for (let i = 0; i < rings; i++) {
    for (let j = 0; j < sides; j++) {
      const a = at(i, j);
      const b = at(i + 1, j);
      const c = at(i + 1, j + 1);
      const d = at(i, j + 1);
      lines.push(`f ${a} ${b} ${c}`, `f ${a} ${c} ${d}`);
    }
  }
  return `${lines.join('\n')}\n`;
}

/** An octahedron as ASCII PLY, each corner its own color, `size` Å from `center` to a corner. */
export function octahedronPly(options: { center?: readonly [number, number, number]; size?: number } = {}): string {
  const [cx, cy, cz] = options.center ?? [0, 0, 0];
  const s = options.size ?? 8;
  const corners: [number, number, number, number, number, number][] = [
    [s, 0, 0, 230, 80, 60],
    [-s, 0, 0, 60, 120, 230],
    [0, s, 0, 240, 200, 60],
    [0, -s, 0, 70, 190, 120],
    [0, 0, s, 200, 90, 210],
    [0, 0, -s, 240, 240, 240],
  ];
  const faces = [
    [0, 2, 4],
    [2, 1, 4],
    [1, 3, 4],
    [3, 0, 4],
    [2, 0, 5],
    [1, 2, 5],
    [3, 1, 5],
    [0, 3, 5],
  ];
  return [
    'ply',
    'format ascii 1.0',
    'comment an octahedron made in code (testing/meshes.ts)',
    `element vertex ${corners.length}`,
    'property float x',
    'property float y',
    'property float z',
    'property uchar red',
    'property uchar green',
    'property uchar blue',
    `element face ${faces.length}`,
    'property list uchar int vertex_indices',
    'end_header',
    ...corners.map(([x, y, z, r, g, b]) => `${cx + x} ${cy + y} ${cz + z} ${r} ${g} ${b}`),
    ...faces.map((f) => `3 ${f.join(' ')}`),
    '',
  ].join('\n');
}

/** The meshes as files, as a person drops them. */
export const MeshFiles = {
  torus: (options?: Parameters<typeof torusObj>[0]) => new File([torusObj(options)], 'torus.obj'),
  octahedron: (options?: Parameters<typeof octahedronPly>[0]) => new File([octahedronPly(options)], 'octahedron.ply'),
};
