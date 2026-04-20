import { execFileSync, spawn } from 'node:child_process';
import { existsSync, rmSync, writeFileSync } from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const projectRoot = path.resolve(__dirname, '..');
const manifestPath = path.join(projectRoot, 'rust', 'ft-promoter-builder', 'Cargo.toml');

function run(command, args, options = {}) {
  return new Promise((resolve, reject) => {
    const child = spawn(command, args, {
      cwd: projectRoot,
      stdio: 'inherit',
      ...options,
    });

    child.on('error', reject);
    child.on('exit', (code) => {
      if (code === 0) {
        resolve();
        return;
      }

      reject(new Error(`${command} exited with code ${code ?? 'unknown'}`));
    });
  });
}

function locateVcvars() {
  const programFilesX86 = process.env['ProgramFiles(x86)'] ?? '';
  const vswhere = path.join(programFilesX86, 'Microsoft Visual Studio', 'Installer', 'vswhere.exe');

  if (existsSync(vswhere)) {
    try {
      const installationPath = execFileSync(
        vswhere,
        ['-latest', '-products', '*', '-requires', 'Microsoft.VisualStudio.Component.VC.Tools.x86.x64', '-property', 'installationPath'],
        { cwd: projectRoot, encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] },
      ).trim();

      if (installationPath) {
        const vcvars = path.join(installationPath, 'VC', 'Auxiliary', 'Build', 'vcvars64.bat');
        if (existsSync(vcvars)) {
          return vcvars;
        }
      }
    } catch {
      // Fall through to the conventional install paths below.
    }
  }

  const candidates = [
    path.join('C:', 'Program Files (x86)', 'Microsoft Visual Studio', '2022', 'BuildTools', 'VC', 'Auxiliary', 'Build', 'vcvars64.bat'),
    path.join('C:', 'Program Files', 'Microsoft Visual Studio', '2022', 'BuildTools', 'VC', 'Auxiliary', 'Build', 'vcvars64.bat'),
    path.join('C:', 'Program Files (x86)', 'Microsoft Visual Studio', '2022', 'Community', 'VC', 'Auxiliary', 'Build', 'vcvars64.bat'),
    path.join('C:', 'Program Files', 'Microsoft Visual Studio', '2022', 'Community', 'VC', 'Auxiliary', 'Build', 'vcvars64.bat'),
  ];

  return candidates.find((candidate) => existsSync(candidate)) ?? null;
}

function locateCargo() {
  try {
    const matches = execFileSync('where.exe', ['cargo'], {
      cwd: projectRoot,
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', 'ignore'],
    })
      .split(/\r?\n/)
      .map((line) => line.trim())
      .filter(Boolean);

    if (matches.length > 0) {
      return matches[0];
    }
  } catch {
    // Fall through to the default cargo home path.
  }

  const cargoHome = process.env.CARGO_HOME ?? path.join(process.env.USERPROFILE ?? '', '.cargo');
  const cargoExe = path.join(cargoHome, 'bin', 'cargo.exe');
  if (existsSync(cargoExe)) {
    return cargoExe;
  }

  return 'cargo';
}

function loadVcvarsEnvironment(vcvars) {
  const scriptPath = path.join(os.tmpdir(), `copilot-vcvars-${process.pid}.cmd`);
  writeFileSync(scriptPath, `@echo off\r\ncall "${vcvars}" >nul\r\nset\r\n`);

  let commandOutput;
  try {
    commandOutput = execFileSync(
      process.env.ComSpec ?? 'cmd.exe',
      ['/d', '/c', scriptPath],
      { cwd: projectRoot, encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] },
    );
  } finally {
    rmSync(scriptPath, { force: true });
  }

  const environment = { ...process.env };
  for (const line of commandOutput.split(/\r?\n/)) {
    const separatorIndex = line.indexOf('=');
    if (separatorIndex <= 0) {
      continue;
    }

    const key = line.slice(0, separatorIndex);
    const value = line.slice(separatorIndex + 1);
    environment[key] = value;
  }

  return environment;
}

async function main() {
  const cargoExecutable = locateCargo();
  const cargoArgs = ['run', '--release', '--manifest-path', manifestPath];

  if (process.platform !== 'win32') {
    await run(cargoExecutable, cargoArgs);
    return;
  }

  const vcvars = locateVcvars();
  if (!vcvars) {
    throw new Error('Could not find Visual Studio C++ build tools. Run npm run build:data:js or install the Visual Studio Build Tools VC workload.');
  }

  const env = loadVcvarsEnvironment(vcvars);
  await run(cargoExecutable, cargoArgs, { env });
}

main().catch((error) => {
  console.error(error.message || error);
  process.exitCode = 1;
});