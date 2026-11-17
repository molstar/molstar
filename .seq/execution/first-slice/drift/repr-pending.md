# Drift: staged pending across a chain of resolved segments

Where: the tree core's staged selection keeps only the last resolver's `pending` (`unknown`), so earlier segments'
pending regions are lost (lane C). Options: repr keeps a list of pendings (recommended), or a merge hook. A tree-core
change: a follow-up for repr before stage 6.
