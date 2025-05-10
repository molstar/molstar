import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    globals: true,
    environment: 'node',
    include: [
      'src/mol-util/__tests__/vitest/**/*.test.ts',
      'src/apps/mesoscale-explorer/data/__tests__/vitest/**/*.test.ts'
    ],
    exclude: [
      'node_modules',
      'dist',
      'build',
      'src/**/__tests__/jest/**',
      'src/**/_spec/**',
    ],
    coverage: {
      reporter: ['text', 'text-summary', 'html'],
    },
  },
});
