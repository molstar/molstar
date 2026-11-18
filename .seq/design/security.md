# Security: loading specs and data you did not write

## Intent

The second review found no security policy: the assets port was said to have a policy `AssetManager` does not have;
specs come from URLs, links and agents; a headless renderer or server fetches whatever `uri` a spec names; regular
expressions, SVG text and `custom` payloads are unbounded. Specs and data are untrusted input, and work they size is
bounded before it starts (§bounded-work).

## Decision

- **A policy wraps the asset manager**: `policyAssets(manager, policy)` is the assets port's default adapter. The policy
  says which locations may be read (schemes, hosts; relative locations under the spec's own base), the largest size of
  one asset and of all in flight, whether credentials are sent (never by default), and how long a read may take. A read
  the policy refuses is an issue on its source, never an exception. Built: locations, an asset's size and a read's
  time; the bound on all in flight is not (TODO).
- **What the policy enforces**: locations (schemes, hosts, base, a headless root), and the size of an asset while its
  body streams and once it is read. Before the second slice only the first and the last; the asset manager's interface
  now carries the rest (the author, 2026-11-18, Gate 0 of the second slice; built by lane M):
  - **The interface carries what any app needs**: `AssetManager.ResolveOptions` gains `credentials`
    (`'omit' | 'same-origin' | 'include'`), `redirect` (`'follow' | 'error'`), `signal` and `maxBytes` (enforced while
    the body streams). A manager honors them or refuses the read. Why: Seq\* is composed into apps that need
    credentials (a lab's server behind a login), so the core must not fix `'omit'`.
  - **Each app composes its own manager**: a derived implementation, or `DefaultAssetManager`, whose defaults stay what
    they are today for every Mol\* app. The policy passes its choice on every read: in the browser viewer
    `credentials: 'omit'` unless the composition says otherwise, `redirect: 'error'` unless the policy says `'follow'`
    (an option of the policy for a server the composition trusts), the read's signal, the size bound. A size refusal
    is a typed error, `SizeLimitExceeded` (core's `util/data-source`), so the policy tells it from other failures.
  - **A conformance suite** every implementation satisfies (`checkAssetManager`, data and a function as the codec
    suite is, on a fake fetch): the four options honored or the read refused, references and release, the shared
    download, decompression, errors not cached. `DefaultAssetManager` runs it in core, and an app runs it on its own.
  - Landed in the second slice: the policy passes the four on every read, and also checks the size once the body is
    read and abandons a read the manager did not end in time, for a manager that falls short. Not covered by the
    suite: a manager that downloads by other means than `fetch` (TODO); `maxBytes` on `file://` reads in Node.
- **Defaults by host**:
  - in a browser viewer: `https:` and relative locations; no credentials; 256 MB an asset;
  - headless (the render CLI, a server): only files under a root the caller names and hosts it lists; no network by
    default (the previous CLI's rule, EC-126); this closes server-side request forgery through a spec's `uri`;
  - tests: no network at all (the test adapter).
- **Bounds on everything a spec sizes** (refused with a coded issue before allocating), each a configuration item with
  this default:

  | What | Default bound |
  |---|---|
  | spec size, decoded | 64 MB |
  | nesting depth (tree, queries) | 256, queries 64 (EC-76) |
  | issues kept | 100 a node, the first 10 places per code (EC-86) |
  | nested sources | depth 4 (EC-179) |
  | a column's length | 2^31 − 1, and the spec size above |
  | a grid or matrix block | 2^24 values (EC-58) |
  | regular expressions in search | linear-time matching only (no backreferences), a step budget |
  | `custom` payloads | counted in the spec size; never interpreted by a reader that does not own them |
  | an archive | entries, total size, no path above its root (EC-91) |

  Built in the second slice (lane Q): the spec's size, nesting, query depth, a column's length and the issues kept
  (`SpecBounds`, held by the check, the viewer's patches, previews and story steps), a block's values and an answer's
  dynamic nodes (at most 10,000 and 16 MB together, `AnswerBounds`), nested sources. The regular expressions' and the
  archive's bounds land with search and archives.

- **Output is escaped**: SVG text and attributes are escaped by the SVG backend; labels from data never become markup;
  links act only through commands on addresses ([addressing.md](addressing.md)).
- **Agents and links act as a user would**: through commands, within the same bounds; a patch is validated before it is
  applied.

## Consequences

- [ports.md](ports.md)'s assets port is the policy wrapper over `AssetManager`; a host passing its own manager passes
  its policy too.
- Tests: each bound refused before allocating (a 20,000-deep spec, a huge declared column, a backtracking regex), the
  headless default refusing the network and paths above its root, SVG escaping.
