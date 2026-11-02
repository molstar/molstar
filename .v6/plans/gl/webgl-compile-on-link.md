# Plan: WebGL shaders compiled on link

## Status

- **Done**: a WebGL program assembles its sources, compiles its shaders and creates its GL program when it links, not
  when its render item creates it, so that the variants no frame or interaction needs are never compiled. From an
  issue report (WebGL2, GTX 980) on merged render objects ([merged-render-objects.md](merged-render-objects.md)). Open:
  the first draw on the reporting machine (below).

## Question

A merged object's first frame waited for six programs where the unmerged one compiled two (cartoon on 4V6X, empty
shader cache: 741–773 against 342–362 ms first draw; ball-and-stick 18 programs and 1.63 s against 12 and 1.37 s).
Linking already waited for the variants a frame needs (`resources.linkPrograms`, `ShaderManager.getRequired`), and the
warm-up links the interaction's (pick, marking) after the frame's. But `createRenderItem` creates every variant's
program (`createProgramVariant` per `GraphicsRenderVariants`), and `createProgram` assembled the sources and called
`getShader`, which calls `compileShader`, at creation. With ANGLE's parallel compile the compile jobs started there and
shared the compiler's threads with the frame's.

## Done

- **Compile in `link()`** (`webgl/program.ts`): `createProgram` keeps what the cache and the render item need at
  creation (the id, the variant, the attribute layout and its limit check); `link()` assembles the sources
  (`addShaderDefines`, `getShaderDeclarations`), gets the shaders (`getShader`: compiled once per source in the
  resources' cache) and creates the GL program; `destroy` releases what `link` created. Compute programs link at
  creation as before. Not taken: creating the variants on first use (the resources' pending set holds created
  programs, `linkPrograms(variants)` reaches them through it).
- **Test** (`renderer.test.ts`, the GL shim): after a commit no shader is created; linking the color variant creates its
  two, linking all the rest.
- **Measured** (new cases `merging/4v6x-{cartoon,ball-and-stick}-{merged,unmerged}`, each alone, `--browser-arg
  --disable-gpu-program-cache --browser-arg --disable-gpu-shader-disk-cache`, WebGL2, RTX 4080 class; the setup note's
  first frame now the first that draws everything, its compile included): add and commit, cartoon merged 61 → 52 ms,
  unmerged 66-70 → 64-68 ms; ball-and-stick merged 60-61 → 39-41 ms, unmerged 91-94 → 66-72 ms; the first frame 3-6 ms
  either way. Visual suites (WebGL2, WebGL1) identical.

- **Without parallel compile** (`--no-parallel-compile` hides `KHR_parallel_shader_compile` from the page: compiles
  block where they are called, possibly the reporting setup): add and commit about doubles, and compiling on link saves
  more, merged most: ball-and-stick merged 142-147 → 114-117 ms, unmerged 163-167 → 146-152 ms; cartoon merged 119-120
  → 110-114 ms, unmerged 120-125 → 119-122 ms; the first frame 2-4 ms either way.

## Open questions

- The reported first draw (hundreds of ms of compile) does not reproduce here: compiles are fast on this machine even
  with the browser's shader caches and parallel compile off. Confirmed only on the reporting setup (GTX 980), or
  another machine with slow compiles.
- Whether ANGLE starts the native compile at `compileShader` on every backend (D3D11, Metal, Vulkan) or defers it to
  link on some: with compiles on link it no longer matters for the variants never linked.
