# Representations

- **Kind:** rule set
- **Status:** stable
- **Scope:** declarative state as recursive, streamable trees; committed vs transient state; stable addresses, hashing, collaboration; transactions and app-wide optimistic state. A project builds its own representation rules on these (encodings, storage, references, history) in its instance.
- **Sources:** the principles as first written, hardened by a reference prototype in the package's home
- **Depends on:** [_core](../_core/README.md), [architecture](../architecture/README.md), [testing](../testing/README.md), [ui](../ui/README.md)

The rules are in [rules.md](rules.md). Code patterns: [bootstrap-patterns/](bootstrap-patterns/README.md) (§committed-transient-identity, §representations, §transactions).
