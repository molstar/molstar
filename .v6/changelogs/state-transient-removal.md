## Fixes from the perf audit

- Fix: a state builder that removes a node no longer writes into the state's tree before it commits. A node added
  under, updated or depended on by a node whose ref hashes like the removed one's (common for generated refs) appeared
  in `state.tree`, and the next update could fail with "Cannot read properties of undefined (reading 'ref')". Nodes
  added to the root after the root was cleared in the same builder are kept
