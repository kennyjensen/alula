import { defineConfig } from '@playwright/test';
export default defineConfig({
  testDir: './tests/browser', timeout: 45000,
  // Solver workers may share CPU with another numerical/browser test.
  expect: { timeout: 15000 },
  use: { baseURL: 'http://localhost:5173', headless: true, viewport: { width: 1440, height: 1100 } },
  webServer: { command: 'node scripts/serve.js', url: 'http://localhost:5173', reuseExistingServer: !process.env.CI },
});
