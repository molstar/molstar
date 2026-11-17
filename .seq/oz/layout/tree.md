# Codebase tree {#codebase-tree}

How a codebase built on these [rules](../README.md) is laid out at the top: the one layout rule and a minimal tree. A
project may keep its full tree, with the rule ids behind each folder, in its instance at
`.<instance>/principles/layout/tree.md` (the adopting procedure, step 1,
[../_core/patterns/adopting.md](../_core/patterns/adopting.md)).

## The top-level layout

**Top-level layout separates compositions, building blocks and tooling.** {#top-level-layout} `lib/` holds common modules (§common-vs-product); `packages/<product>/<module>` holds a product's building blocks (§family-folders); `apps/` (user-facing) and `services/` (deployable backends) hold only top-level compositions with their composition roots (§one-composition-root); `tooling/` holds everything that supports development but never ships: Storybook (§design-system), the eval runner (§evals), scripts (test runner, gate checks), code generators and local agent tool hosts. Tooling may depend on everything else; nothing depends on tooling. Process and knowledge folders (the package, the instance with its designs, sources, records and work folders, `_spec/`, local scratch) and `prototypes/` (§prototype-outside) sit at the root.

```
repo/
├── AGENTS.md                # root agent map: modules, tasks, where the rules and skills live     [§agents-md]
├── <entry point>            # ONE entry point: pinned toolchains + tasks (setup, test, gate)       [§one-entry-point]
├── .principles/             # the package; its checks in _checks/                                   [§rules-are-source]
├── .<instance>/             # this project's instance (tracked, except its work folder)
│   ├── principles/          # the instance's manifest (README.md) and its own modules
│   ├── design/              # agreed designs, a folder per area, one file per topic                [§design-first]
│   ├── sources/             # what the principles are drawn from and tested against
│   ├── records/             # one record per execution session                                     [§permanent-record]
│   ├── execution/           # GITIGNORED: work in flight (plan, loop, state)                       [§execution-unreferenced]
│   └── _skills/             # the instance skill
├── _spec/                   # global semantic mirror                                              [§spec-folders]
├── .tmp/                    # GITIGNORED: local scratch                                           [§machine-invariant]
├── lib/                     # COMMON: abstract, configurable, no product knowledge                [§common-vs-product]
│   └── core/                # the base layer every lib package may use; depends on nothing else
├── packages/<product>/      # PRODUCT BUILDING BLOCKS, one family folder per product              [§family-folders]
├── apps/<product>/<app>/    # USER-FACING COMPOSITIONS: a composition root, an entry, smoke tests [§one-composition-root]
├── services/<product>/<svc>/   # BACKEND COMPOSITIONS: same shape as apps
├── prototypes/<name>/       # PROTOTYPES: outside the build, laid out like this tree               [§prototype-outside]
└── tooling/                 # supports development, never ships
```

## Dependency direction

```
tooling/  ──▶  apps/, services/  ──▶  packages/<product>/*  ──▶  lib/ (domain modules)  ──▶  lib/core/*
```

- Arrows only point right; an import-boundary lint enforces it (§common-vs-product). Apps and services never import
  each other; what two of them share is a package.
- `apps/` and `services/` stay thin: a composition root, an entry point, smoke tests and deploy config. Anything worth
  unit-testing belongs in a package.
- Product packages of different products never import each other; shared code is promoted to `lib/`
  (§promote-on-reuse).
- Provider SDKs appear only in `*-<provider>-runtime` packages and services (§provider-runtime-packages); environment
  variables are read only in composition roots.
