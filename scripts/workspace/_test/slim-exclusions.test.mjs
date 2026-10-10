import assert from 'node:assert/strict';
import { test } from 'node:test';
import { exclusionErrors, findExcluded, globToRegExp, loadExclusions } from '../slim-exclusions.mjs';

const matches = (glob, file) => globToRegExp(glob).test(file);

test('globs: *, ?, ** and {a,b} match whole posix paths', () => {
  assert.ok(matches('a/b/*.ts', 'a/b/c.ts'));
  assert.ok(!matches('a/b/*.ts', 'a/b/c/d.ts'));
  assert.ok(!matches('a/b/*.ts', 'a/b/c.tsx'));
  assert.ok(matches('a/mmcif*', 'a/mmcif.ts'));
  assert.ok(matches('a/mmcif*', 'a/mmcif-format.ts'));
  assert.ok(matches('a/?.ts', 'a/x.ts'));
  assert.ok(!matches('a/?.ts', 'a/xy.ts'));
  assert.ok(matches('a/**', 'a/b/c/d.ts'));
  assert.ok(!matches('a/**', 'b/a/c.ts'));
  assert.ok(matches('a/**/d.ts', 'a/d.ts'));
  assert.ok(matches('a/**/d.ts', 'a/b/c/d.ts'));
  assert.ok(matches('**/catalog.ts', 'catalog.ts'));
  assert.ok(matches('**/catalog.ts', 'packages/x/src/catalog.ts'));
  assert.ok(!matches('**/catalog.ts', 'packages/x/src/not-catalog.ts'));
  assert.ok(matches('a/{x,y}.ts', 'a/y.ts'));
  assert.ok(!matches('a/{x,y}.ts', 'a/z.ts'));
  assert.ok(matches('a/{x,y}/{p,q}.ts', 'a/y/q.ts'));
  assert.ok(matches('a/b.c', 'a/b.c'));
  assert.ok(!matches('a/b.c', 'a/bxc'));
  assert.throws(() => globToRegExp('a/{x,y'), /Unbalanced/);
});

test('findExcluded splits violations from known leaks and reports unused leaks', () => {
  const exclusions = {
    excluded: [
      { name: 'One', paths: ['lib/one/**'] },
      { name: 'Two', paths: ['lib/two.ts', 'lib/three.ts'] },
    ],
    knownLeaks: [
      { path: 'lib/one/leak.ts', reason: 'r', decision: 'd' },
      { path: 'lib/unused.ts', reason: 'r', decision: 'd' },
    ],
  };
  const result = findExcluded(['lib/one/leak.ts', 'lib/one/other.ts', 'lib/two.ts', 'lib/fine.ts'], exclusions);
  assert.deepEqual(
    [...result.violations],
    [
      ['lib/one/other.ts', 'One'],
      ['lib/two.ts', 'Two'],
    ],
  );
  assert.deepEqual([...result.known.keys()], ['lib/one/leak.ts']);
  assert.deepEqual(
    result.unusedLeaks.map((leak) => leak.path),
    ['lib/unused.ts'],
  );
});

test('exclusionErrors validates the shape', () => {
  assert.deepEqual(exclusionErrors({ entry: 'e', excluded: [{ name: 'n', paths: ['p'] }], knownLeaks: [] }), []);
  assert.equal(exclusionErrors({ excluded: [] }).length, 2);
  assert.equal(exclusionErrors({ entry: 'e', excluded: [{ name: 'n' }] }).length, 1);
  assert.equal(
    exclusionErrors({ entry: 'e', excluded: [{ name: 'n', paths: ['p'] }], knownLeaks: [{ path: 'p', reason: 'r' }] })
      .length,
    1,
  );
});

test('the checked-in exclusions pin the spec §12 table', () => {
  const exclusions = loadExclusions();
  assert.deepEqual(exclusionErrors(exclusions), []);
  assert.equal(exclusions.entry, 'examples/slim-plugin/src/index.ts');
  const { violations } = findExcluded(
    [
      'packages/model/src/formats/structure/mmcif.ts',
      'packages/plugin/core/src/state/formats/trajectory/mmcif.ts',
      'packages/io/src/reader/ccp4/parser.ts',
      'packages/plugin/core/src/state/formats/volume/ccp4.ts',
      'packages/graphics/src/repr/structure/representation/cartoon.ts',
      'packages/graphics/src/repr/volume/direct-volume.ts',
      'packages/graphics/src/repr/volume/isosurface.ts',
      'packages/graphics/src/repr/volume/slice.ts',
      'packages/graphics/src/repr/volume/dot.ts',
      'packages/graphics/src/repr/volume/segment.ts',
      'packages/plugin/core/src/state/builder/structure/representation-presets/auto.ts',
      'packages/plugin/core/src/state/builder/structure/representation-presets/catalog.ts',
      'packages/plugin/core/src/state/builder/structure/hierarchy-presets/supercell.ts',
      'packages/plugin/core/src/state/builder/structure/hierarchy-presets/catalog.ts',
      'packages/plugin/core/src/state/queries/structure/catalog.ts',
      'packages/model/src/script/transpilers/all.ts',
      'packages/model/src/script/transpilers/pymol/parser.ts',
      'extensions/mp4-export/src/index.ts',
      'packages/graphics/src/theme/color/catalog.ts',
      'packages/plugin/core/src/default-spec.ts',
      'packages/plugin/core/src/default-registry.ts',
      'packages/plugin/ui/src/default-spec.ts',
    ],
    exclusions,
  );
  assert.equal(violations.size, 22);

  // What the slim example builds and runs is not excluded.
  const kept = [
    'packages/plugin/core/src/state/builder/structure/representation-presets/ball-and-stick.ts',
    'packages/plugin/core/src/state/builder/structure/representation-presets/types.ts',
    'packages/plugin/core/src/state/builder/structure/hierarchy-presets/default.ts',
    'packages/plugin/core/src/state/builder/structure/hierarchy-presets/types.ts',
    'packages/plugin/core/src/state/formats/trajectory/sdf.ts',
    'packages/plugin/core/src/registry/structure/ball-and-stick.ts',
    'packages/model/src/formats/structure/mmcif-format.ts',
    'packages/plugin/core/src/state/queries/structure/type.ts',
    'packages/model/src/script/script.ts',
  ];
  assert.equal(findExcluded(kept, exclusions).violations.size, 0);
});
