<#
.SYNOPSIS
  Installs the RDKit + OPSIN third-party dependencies chemcad needs on Windows.

.DESCRIPTION
  Windows equivalent of scripts/setup_deps.sh. RDKit has no Windows package and
  building it from source needs a multi-hour C++ build, so this pulls a
  prebuilt conda-forge RDKit via micromamba instead -- no compiler required
  for this step.

  Idempotent: re-running with everything already present is a no-op and says
  so for each step. Does not require administrator rights -- everything lands
  under -Prefix (default $env:LOCALAPPDATA\chemcad-deps).

.PARAMETER Prefix
  Install root. Defaults to $env:CHEMCAD_DEPS_PREFIX, else
  "$env:LOCALAPPDATA\chemcad-deps".
#>
[CmdletBinding()]
param(
  [string]$Prefix = $(if ($env:CHEMCAD_DEPS_PREFIX) { $env:CHEMCAD_DEPS_PREFIX }
                      else { Join-Path $env:LOCALAPPDATA 'chemcad-deps' })
)

$ErrorActionPreference = 'Stop'
# The default Invoke-WebRequest progress bar renders one frame per byte over
# a remote session and makes multi-hundred-MB downloads look hung.
$ProgressPreference = 'SilentlyContinue'

function Write-Step {
  param([string]$Message)
  Write-Host "[deps] $Message"
}

$opsinVersion = if ($env:OPSIN_VERSION) { $env:OPSIN_VERSION } else { '2.9.0' }

$toolsDir = Join-Path $Prefix 'tools'
$rdkitPrefix = Join-Path $Prefix 'rdkit'
$opsinJar = Join-Path $Prefix 'share\opsin\opsin.jar'

New-Item -ItemType Directory -Force -Path $Prefix | Out-Null
New-Item -ItemType Directory -Force -Path $toolsDir | Out-Null

# ---------------------------------------------------------------- micromamba
# micromamba is only a means to an end here: it is the smallest tool that can
# pull a prebuilt conda-forge RDKit without needing a compiler on this box.
$micromambaExe = Join-Path $toolsDir 'micromamba.exe'
if (Test-Path $micromambaExe) {
  Write-Step "micromamba already installed at $micromambaExe"
} else {
  Write-Step 'downloading micromamba (win-64)'
  $archive = Join-Path $toolsDir 'micromamba.tar.bz2'
  Invoke-WebRequest -Uri 'https://micro.mamba.pm/api/micromamba/win-64/latest' -OutFile $archive
  # The archive is a .tar.bz2; the bsdtar shipped with Windows 10 1803+ and
  # Windows 11 handles bzip2 transparently, so plain "tar -xf" is enough.
  tar -xf $archive -C $toolsDir
  if ($LASTEXITCODE -ne 0) {
    throw "tar failed to extract $archive (exit code $LASTEXITCODE)"
  }
  Move-Item -Force (Join-Path $toolsDir 'Library\bin\micromamba.exe') $micromambaExe
  Remove-Item -Force $archive
  Remove-Item -Recurse -Force (Join-Path $toolsDir 'Library') -ErrorAction SilentlyContinue
  Remove-Item -Recurse -Force (Join-Path $toolsDir 'info') -ErrorAction SilentlyContinue
  Write-Step "micromamba installed -> $micromambaExe"
}

# micromamba needs a root prefix for its package cache even though the target
# environment lives at -p; keep it inside our own prefix so this script never
# touches a real conda/mamba install the user may already have.
$env:MAMBA_ROOT_PREFIX = Join-Path $toolsDir 'mamba_root'
New-Item -ItemType Directory -Force -Path $env:MAMBA_ROOT_PREFIX | Out-Null

# ---------------------------------------------------------------- RDKit
$rdkitConfig = Join-Path $rdkitPrefix 'Library\lib\cmake\rdkit\rdkit-config.cmake'
if ($env:CHEMCAD_SKIP_RDKIT -eq '1') {
  Write-Step 'CHEMCAD_SKIP_RDKIT set -- skipping the RDKit environment'
} elseif (Test-Path $rdkitConfig) {
  Write-Step "RDKit already installed at $rdkitPrefix"
} else {
  Write-Step 'creating conda-forge RDKit environment (a few minutes, mostly download)'
  & $micromambaExe create -y -p $rdkitPrefix -c conda-forge `
    librdkit-dev libboost-devel libcurl cairo freetype eigen zlib
  if ($LASTEXITCODE -ne 0) {
    throw "micromamba create failed (exit code $LASTEXITCODE)"
  }
  Write-Step "RDKit installed -> $rdkitPrefix"
}

# ---------------------------------------------------------------- OPSIN
if (Test-Path $opsinJar) {
  Write-Step "OPSIN already present at $opsinJar"
} else {
  Write-Step "downloading OPSIN $opsinVersion"
  New-Item -ItemType Directory -Force -Path (Split-Path $opsinJar) | Out-Null
  $opsinTmp = "$opsinJar.tmp"
  $opsinUrl = "https://github.com/dan2097/opsin/releases/download/$opsinVersion/" +
              "opsin-cli-$opsinVersion-jar-with-dependencies.jar"
  Invoke-WebRequest -Uri $opsinUrl -OutFile $opsinTmp
  Move-Item -Force $opsinTmp $opsinJar
  Write-Step "OPSIN installed -> $opsinJar"
}
Write-Step 'OPSIN also needs a JRE on PATH; without one, name -> structure falls back to PubChem.'

Write-Step 'done. Configure chemcad with:'
Write-Host "  cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=`"$rdkitPrefix\Library`""
