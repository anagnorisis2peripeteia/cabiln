const fs = require('node:fs');
const path = require('node:path');
const os = require('node:os');
const net = require('node:net');
const { spawn, execFileSync } = require('node:child_process');
const { once } = require('node:events');
const { test: base, expect } = require('@playwright/test');

const repo = path.resolve(__dirname, '../..');
const appRoot = path.resolve(process.env.CABILN_APP_ROOT || repo);
const localPython = path.join(repo, '.venv', process.platform === 'win32' ? 'Scripts/python.exe' : 'bin/python');
const python = process.env.CABILN_PYTHON || (fs.existsSync(localPython) ? localPython : 'python3');
const installed = process.env.CABILN_BROWSER_INSTALLED === '1';
const packageData = installed ? execFileSync(python, ['-I', '-c',
  'from importlib.resources import files; print(files("pyPept.data"))'],
{ encoding: 'utf8' }).trim() : path.join(appRoot, 'src/pyPept/data');

async function freePort() {
  const server = net.createServer();
  server.listen(0, '127.0.0.1');
  await once(server, 'listening');
  const port = server.address().port;
  await new Promise(resolve => server.close(resolve));
  return port;
}

async function startApp(workerInfo, registration) {
  const temporary = fs.mkdtempSync(path.join(os.tmpdir(), 'cabiln-browser-library-'));
  for (const name of ['monomers.sdf', 'monomers.csv']) {
    fs.copyFileSync(path.join(packageData, name), path.join(temporary, name));
  }
  const port = await freePort();
  const url = `http://127.0.0.1:${port}`;
  const logPath = path.join(workerInfo.project.outputDir, `server-${registration ? 'writable' : 'readonly'}.log`);
  fs.mkdirSync(path.dirname(logPath), { recursive: true });
  const log = fs.createWriteStream(logPath);
  fs.writeFileSync(path.join(workerInfo.project.outputDir, 'application.json'), JSON.stringify({
    appRoot, python, installed, browserChannel: process.env.CABILN_BROWSER_CHANNEL || 'chromium',
    viewport: workerInfo.project.use.viewport,
  }, null, 2));
  const server = spawn(python, [
    ...(installed ? ['-I'] : []),
    '-m', 'uvicorn', 'pyPept.web.app:app', '--host', '127.0.0.1', '--port', String(port),
    '--no-access-log',
  ], {
    cwd: appRoot,
    env: {
      ...process.env,
      // Explicitly select the app source even when python belongs to another editable checkout.
      ...(installed ? {} : { PYTHONPATH: path.join(appRoot, 'src') }),
      CABILN_MONOMER_LIBRARY: path.join(temporary, 'monomers.sdf'),
      CABILN_ENABLE_REGISTRATION: registration ? '1' : '0',
    },
    stdio: ['ignore', 'pipe', 'pipe'],
  });
  server.stdout.pipe(log);
  server.stderr.pipe(log);
  let spawnError;
  server.on('error', error => { spawnError = error; });
  async function stop() {
    if (server.exitCode === null && !spawnError) {
      const exited = once(server, 'exit');
      server.kill('SIGTERM');
      const force = setTimeout(() => server.kill('SIGKILL'), 5000);
      await exited;
      clearTimeout(force);
    }
    await new Promise(resolve => log.end(resolve));
    fs.rmSync(temporary, { recursive: true, force: true });
  }
  try {
    const deadline = Date.now() + 30_000;
    while (true) {
      if (spawnError) throw spawnError;
      if (server.exitCode !== null) throw new Error(`Application exited (${server.exitCode}); see ${logPath}`);
      try {
        if ((await fetch(`${url}/ready`, { signal: AbortSignal.timeout(1000) })).ok) break;
      } catch { /* Wait for this process to bind its local socket. */ }
      if (Date.now() > deadline) throw new Error(`Application did not become healthy; see ${logPath}`);
      await new Promise(resolve => setTimeout(resolve, 100));
    }
    return { url, appRoot, temporary, logPath, stop };
  } catch (error) {
    await stop();
    throw error;
  }
}

const test = base.extend({
  app: [async ({}, use, workerInfo) => {
    const app = await startApp(workerInfo, true);
    try { await use(app); } finally { await app.stop(); }
  }, { scope: 'worker' }],
  readonlyApp: [async ({}, use, workerInfo) => {
    const app = await startApp(workerInfo, false);
    try { await use(app); } finally { await app.stop(); }
  }, { scope: 'worker' }],
  page: async ({ page, context, app }, use) => {
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    const watchPage = child => child.on('pageerror', error => errors.push(error.message));
    context.on('page', watchPage);
    await page.goto(app.url);
    await use(page);
    context.off('page', watchPage);
    expect(errors, 'No unhandled browser exceptions').toEqual([]);
  },
});

function isCompletedResponse(response) {
  // Match the builder's bounded admission retry policy. Other errors must reach
  // the existing assertions; exhausting 503 retries still fails the wait.
  const retryAfter = Number(response.headers()['retry-after']);
  return response.status() !== 503 || !Number.isFinite(retryAfter) ||
    retryAfter <= 0 || retryAfter > 2;
}

async function render(page, source) {
  const received = page.waitForResponse(response =>
    isCompletedResponse(response) && new URL(response.url()).pathname === '/render' &&
    response.request().postDataJSON()?.cabiln === source);
  await page.locator('#cabiln-input').fill(source);
  const response = await received;
  expect(response.status(), await response.text()).toBe(200);
  const data = await response.json();
  expect(data.error).toBeUndefined();
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(page.locator('#render-inner svg')).toBeVisible();
  await expect(page.locator('#residue-chips [data-residue]')).toHaveCount(data.residues.length);
  return data;
}

async function tile(page, abbr, button = 'left') {
  await page.locator('#lib-search').fill(abbr);
  await page.locator('.lib-row').filter({ has: page.locator('.lib-abbr', { hasText: new RegExp(`^${abbr}$`) }) }).click({ button });
}

async function site(page, side, slot) {
  const button = page.locator(`#build-${side}-rgroups button`).filter({ hasText: new RegExp(`^R${slot} `) });
  await expect(button).not.toHaveClass(/used/);
  await button.click();
  await expect(button).toHaveClass(/selected/);
}

async function selectChip(page, idx, side, abbr) {
  await page.locator(`#residue-chips [data-residue="${idx}"]`).click();
  await expect(page.locator(`#build-${side}-abbr`)).toHaveText(abbr);
  await expect(page.locator(`#build-${side}-rgroups button`).first()).toBeVisible();
}

async function connect(page) {
  await expect(page.locator('#build-connect')).toBeEnabled();
  const result = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/insert_bond');
  await page.locator('#build-connect').click();
  const response = await result;
  expect(response.status(), await response.text()).toBe(200);
  const data = await response.json();
  expect(data.error).toBeUndefined();
  await expect(page.locator('#cabiln-input')).toHaveValue(data.result);
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  return { data, submitted: response.request().postDataJSON() };
}

async function capture(page, testInfo, name) {
  // Compare settled visual states; live interaction tests retain the shipped
  // transitions. Fast-forwarding only for capture avoids timing noise in PNGs.
  await page.screenshot({ path: testInfo.outputPath(`${name}.png`), animations: 'disabled' });
  const state = await page.evaluate(() => {
    const ids = ['lib-panel', 'examples-panel', 'build-panel', 'verify-pane', 'render-canvas'];
    return {
      source: document.querySelector('#cabiln-input').value,
      notation: document.querySelector('#notation-select').value,
      chips: [...document.querySelectorAll('#residue-chips .res-chip')].map(el => ({
        text: el.textContent, occurrence: el.dataset.residue, members: el.dataset.members,
        outline: el.style.outline, hover: el.classList.contains('hover'),
      })),
      panels: Object.fromEntries(ids.map(id => {
        const el = document.getElementById(id);
        const rect = el.getBoundingClientRect();
        return [id, { visible: !!(rect.width && rect.height), open: el.classList.contains('open'),
          x: rect.x, y: rect.y, width: rect.width, height: rect.height }];
      })),
      controls: [...document.querySelectorAll('button')].filter(el => el.id).map(el => ({
        id: el.id, text: el.textContent, disabled: el.disabled, active: el.classList.contains('active'),
      })),
      highlighted: [...document.querySelectorAll('#render-inner .res-hl')].map(el => el.getAttribute('class')),
      previewVisible: document.getElementById('lib-preview').getBoundingClientRect().width > 0,
    };
  });
  const statePath = testInfo.outputPath(`${name}.json`);
  fs.writeFileSync(statePath, JSON.stringify(state, null, 2));
  await testInfo.attach(`${name}.json`, { path: statePath, contentType: 'application/json' });
}

async function residuePoint(page, rendered, occurrence, label = false) {
  const point = await page.evaluate(({ atoms, label }) => {
    for (const path of document.querySelectorAll('#render-inner svg path[class]')) {
      const classes = path.getAttribute('class');
      const owners = [...classes.matchAll(/atom-(\d+)/g)].map(match => Number(match[1]));
      if (!owners.length || !owners.every(atom => atoms.includes(atom))) continue;
      if (classes.includes('bond-') === label) continue;
      const box = path.getBoundingClientRect();
      let x = box.x + box.width / 2, y = box.y + box.height / 2;
      if (!label) {
        const length = path.getTotalLength(), matrix = path.getScreenCTM();
        const start = path.getPointAtLength(0).matrixTransform(matrix);
        const end = path.getPointAtLength(length).matrixTransform(matrix);
        const mid = path.getPointAtLength(length / 2).matrixTransform(matrix);
        const span = Math.hypot(end.x - start.x, end.y - start.y);
        if (span < 20) continue;
        x = mid.x - (end.y - start.y) / span * 5;
        y = mid.y + (end.x - start.x) / span * 5;
      }
      const hit = document.elementFromPoint(x, y);
      if (hit?.closest('#render-inner') && !/atom-\d+/.test(hit.getAttribute('class') || '')) return { x, y };
    }
    return null;
  }, { atoms: rendered.residue_map[String(occurrence)], label });
  expect(point, label ? 'An unpainted atom-label interior' : 'A point 5px beside an owned bond').not.toBeNull();
  return point;
}

module.exports = { test, expect, render, tile, site, selectChip, connect, capture, isCompletedResponse, residuePoint };
