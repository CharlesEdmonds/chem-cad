#!/usr/bin/env bash
# Builds/install third-party deps that are not available as Fedora packages.
#   1. RDKit (C++ toolkit)  -> $PREFIX
#   2. OPSIN cli jar        -> $PREFIX/share/opsin/opsin.jar
#
# Distributions that package RDKit (Debian/Ubuntu: librdkit-dev) only need the
# OPSIN jar -- set CHEMCAD_SKIP_RDKIT=1 to skip the 20-minute source build.
#
# Idempotent: re-running with everything present is a no-op.
set -euo pipefail

PREFIX="${CHEMCAD_DEPS_PREFIX:-$HOME/.local/share/chemcad-deps}"
SRC="${CHEMCAD_DEPS_SRC:-$HOME/.cache/chemcad-deps}"
RDKIT_TAG="${RDKIT_TAG:-Release_2026_03_5}"
OPSIN_VERSION="${OPSIN_VERSION:-2.9.0}"
JOBS="${JOBS:-$(nproc)}"

mkdir -p "$PREFIX" "$SRC"

# ---------------------------------------------------------------- RDKit
if [ -n "${CHEMCAD_SKIP_RDKIT:-}" ]; then
  echo "[deps] CHEMCAD_SKIP_RDKIT set -- using the system RDKit"
elif [ -f "$PREFIX/rdkit/lib64/cmake/rdkit/rdkit-config.cmake" ] || \
     [ -f "$PREFIX/rdkit/lib/cmake/rdkit/rdkit-config.cmake" ]; then
  echo "[deps] RDKit already installed at $PREFIX/rdkit"
else
  echo "[deps] building RDKit $RDKIT_TAG (this takes 15-30 min)"
  if [ ! -d "$SRC/rdkit/.git" ]; then
    git clone --depth 1 --branch "$RDKIT_TAG" https://github.com/rdkit/rdkit.git "$SRC/rdkit"
  fi
  cmake -S "$SRC/rdkit" -B "$SRC/rdkit-build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PREFIX/rdkit" \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    -DCMAKE_INSTALL_RPATH="$PREFIX/rdkit/lib" \
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
    -DRDK_INSTALL_INTREE=OFF \
    -DRDK_BUILD_CPP_TESTS=OFF \
    -DRDK_BUILD_SLN_SUPPORT=OFF \
    -DRDK_TEST_MULTITHREADED=OFF \
    -DRDK_BUILD_DESCRIPTORS3D=OFF \
    -DRDK_BUILD_MAEPARSER_SUPPORT=OFF \
    -DRDK_BUILD_COORDGEN_SUPPORT=ON \
    -DRDK_BUILD_FREETYPE_SUPPORT=ON \
    -DRDK_BUILD_CAIRO_SUPPORT=OFF \
    -DRDK_BUILD_INCHI_SUPPORT=OFF \
    -DRDK_BUILD_AVALON_SUPPORT=OFF \
    -DRDK_BUILD_YAEHMOP_SUPPORT=OFF \
    -DRDK_BUILD_XYZ2MOL_SUPPORT=OFF \
    -DRDK_INSTALL_COMIC_FONTS=OFF \
    -DRDK_BUILD_TEST_GZIP=OFF
  cmake --build "$SRC/rdkit-build" -j "$JOBS"
  cmake --install "$SRC/rdkit-build"
  echo "[deps] RDKit installed -> $PREFIX/rdkit"
fi

# ---------------------------------------------------------------- OPSIN
OPSIN_JAR="$PREFIX/share/opsin/opsin.jar"
if [ -s "$OPSIN_JAR" ]; then
  echo "[deps] OPSIN already present at $OPSIN_JAR"
else
  echo "[deps] downloading OPSIN $OPSIN_VERSION"
  mkdir -p "$(dirname "$OPSIN_JAR")"
  curl -fL --retry 3 -o "$OPSIN_JAR.tmp" \
    "https://github.com/dan2097/opsin/releases/download/${OPSIN_VERSION}/opsin-cli-${OPSIN_VERSION}-jar-with-dependencies.jar"
  mv "$OPSIN_JAR.tmp" "$OPSIN_JAR"
  echo "[deps] OPSIN installed -> $OPSIN_JAR"
fi

echo "[deps] done. Configure chemcad with:"
echo "  cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release"
