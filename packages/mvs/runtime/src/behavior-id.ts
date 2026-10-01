/** Stable transformer identity used by Mol* state snapshots. */
export const MVS_BEHAVIOR_NAME = 'molviewspec';

/**
 * Keep the behavior id in a leaf module so the loader can check registration
 * without importing the behavior module and creating a load/behavior cycle.
 */
export const MVS_BEHAVIOR_ID = `ms-plugin.${MVS_BEHAVIOR_NAME}`;
