# Zero-allocation hot loops: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## Allocate in setup, reuse in the loop

```ts
// BAD: allocates per frame (arrays, closures, objects); points of a fading pointer trail
function tick(points: TrailPoint[]) {
  const alive = points.filter((p) => p.life > 0);             // new array + closure
  return alive.map((p) => ({ x: p.x + p.vx, y: p.y + p.vy })); // N new objects
}

// GOOD: struct-of-arrays, in-place update, no per-frame allocation
export class TrailPoints {
  readonly x: Float32Array;
  readonly y: Float32Array;
  readonly vx: Float32Array;
  readonly vy: Float32Array;
  readonly life: Float32Array;
  count = 0;

  constructor(capacity: number) {
    this.x = new Float32Array(capacity);
    this.y = new Float32Array(capacity);
    this.vx = new Float32Array(capacity);
    this.vy = new Float32Array(capacity);
    this.life = new Float32Array(capacity);
  }

  step(dt: number) {
    let w = 0;
    for (let r = 0; r < this.count; r++) {                     // integrate and compact in one pass
      const life = this.life[r] - dt;
      if (life <= 0) continue;
      this.x[w] = this.x[r] + this.vx[r] * dt; this.y[w] = this.y[r] + this.vy[r] * dt;
      this.vx[w] = this.vx[r]; this.vy[w] = this.vy[r]; this.life[w] = life;
      w++;
    }
    this.count = w;
  }
}
```

## Out-parameters and pools

```ts
export function normalizeInto(out: Float32Array, x: number, y: number, z: number) {
  const inv = 1 / (Math.hypot(x, y, z) || 1);   // zero vector stays zero
  out[0] = x * inv; out[1] = y * inv; out[2] = z * inv;
  return out;                                    // the caller supplies `out`
}

export class Pool<T> {
  private free: T[] = [];
  constructor(private make: () => T, private reset: (t: T) => void) {}
  acquire(): T { return this.free.pop() ?? this.make(); }
  release(t: T) { this.reset(t); this.free.push(t); }
}
```

## Scratch arenas

What a scratch arena is and when to use it: see [README.md](README.md#shape).

In JS the arena hands out **offsets into one typed array**, not new views: `new Float32Array(buffer, off, n)` itself allocates an object per call.

```ts
export class F32Arena {
  readonly data: Float32Array;
  private top = 0;
  private max = 0;
  constructor(capacity: number, private readonly debug = false) { this.data = new Float32Array(capacity); }

  /** Reserve n floats; returns the offset into `data`. O(1), no allocation. */
  alloc(n: number): number {
    const off = this.top;
    if (off + n > this.data.length) throw new Error(`arena exhausted (${off + n} > ${this.data.length}); raise capacity`);
    this.top = off + n;
    if (this.top > this.max) this.max = this.top;
    return off;
  }
  mark(): number { return this.top; }
  release(mark: number) {
    if (this.debug) this.data.fill(NaN, mark, this.top);   // poison: use-after-release shows up as NaN in tests
    this.top = mark;
  }
  reset() { this.release(0); }
  /** Highest use since creation: the number a guard test budgets. */
  get peak() { return this.max; }
}
```

```ts
// The arena is created once by the runtime/composition and passed in: no module-level singleton (§dependencies-passed-in).
function smoothNormals(arena: F32Arena, positions: Float32Array, out: Float32Array) {
  const m = arena.mark();
  try {
    const tmp = arena.alloc(positions.length);
    for (let i = 0; i < positions.length; i++) arena.data[tmp + i] = positions[i] * 0.5;
    /* ... combine arena.data[tmp + i] into `out` ... */
  } finally {
    arena.release(m);                                       // temporaries freed, even on error
  }
}

// Per-frame arena: reset once per frame; passes allocate from it and free nothing individually.
function renderTick(ctx: { frameArena: F32Arena }) {
  ctx.frameArena.reset();
  /* ... */
}
```

Invariants (scope, nesting, one arena per lane, budgeted capacity): see [README.md](README.md#invariants).

- Other element types (`U32Arena`, ...) and typed views over WASM memory work the same way. For long-lived or irregular lifetimes use pools or normal allocation.

```ts
// frame-arena.guard.test.ts
it("frame arena stays within budget", () => {
  const ctx = { frameArena: new F32Arena(1 << 20, true) };
  for (let i = 0; i < 100; i++) renderTick(ctx);
  expect(ctx.frameArena.peak).toBeLessThan(ctx.frameArena.data.length * 0.5);
});
```

## Avoid accidentally quadratic

```ts
// BAD: O(n*m), includes inside a loop
const slowHits = a.filter((x) => b.includes(x));

// GOOD: O(n+m)
const bSet = new Set(b);
const hits = a.filter((x) => bSet.has(x));
```

If N is bounded and small, say so where it matters: `// N <= 8 (max open panes); linear scan is intentional.`

## Choose structures for the access pattern

```ts
// Membership / lookup by key, many queries      -> Set / Map           O(1)
// Ordered iteration, few lookups, small N       -> array               (N <= ~16: scan beats hashing)
// Ordered + range queries, frequent inserts     -> balanced tree / sorted array + binary search
// Numeric data iterated in bulk (hot loops)     -> typed arrays (SoA)   cache-friendly, no per-item objects
// FIFO with many pops                           -> ring buffer          (Array.shift() is O(n))
// Repeated "min/max next"                       -> binary heap          not re-sorting each time
// Incremental recompute                         -> dirty flags / keyed diff, not full rebuild

// N = open documents, expected < 50, cold path: array scan is fine.
const open = docs.find((d) => d.id === id);
```

Allocation-budget and scaling tests live in the guard tier: see [§test-tiers](../../../testing/bootstrap-patterns/test-tiers/typescript.md#guard-tier-allocation-and-complexity).
