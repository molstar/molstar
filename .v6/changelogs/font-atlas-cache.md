## Graphics: labels

- Fix: label, measurement and MVS label texts share one font atlas per font (the 4 most recent are kept). Each label
  size step or MVS label set kept a new 43 MB atlas for the session: 20 steps of a label slider held 0.86 GB
  (perf-audit-2 08/B6, #2028)
