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
  the policy refuses is an issue on its source, never an exception.
- **What the policy enforces today**: locations (schemes, hosts, base), and the size of an asset once its body is read.
  What waits for a core asset-manager change (proposed core change; the author's call): the manager's fetch sends
  cookies on relative reads, follows redirects to disallowed hosts, ignores the signal and checks size only after the
  whole body is read. The change is `credentials: 'omit'`, `redirect: 'error'`, the signal passed through and a streamed
  byte limit. Until it lands the policy does not claim those four, and says so in its documentation.
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

- **Output is escaped**: SVG text and attributes are escaped by the SVG backend; labels from data never become markup;
  links act only through commands on addresses ([addressing.md](addressing.md)).
- **Agents and links act as a user would**: through commands, within the same bounds; a patch is validated before it is
  applied.

## Consequences

- [ports.md](ports.md)'s assets port is the policy wrapper over `AssetManager`; a host passing its own manager passes
  its policy too.
- Tests: each bound refused before allocating (a 20,000-deep spec, a huge declared column, a backtracking regex), the
  headless default refusing the network and paths above its root, SVG escaping.
