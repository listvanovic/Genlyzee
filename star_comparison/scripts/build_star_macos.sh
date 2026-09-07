#!/usr/bin/env bash
# Gradi STAR 2.7.11b iz izvornog koda uz popravak za macOS/libc++.
#
# Homebrew "rna-star" (kao i obični `make STAR` na macOS-u) linka se na Appleov
# libc++, gdje je std::stringbuf::pubsetbuf() no-op (standard ga ostavlja
# implementation-defined; GNU libstdc++ ga implementira kao alias na zadani
# buffer). STAR-ovi bufferi za čitanje readova i pisanje SAM zapisa oslanjaju se
# upravo na to GNU ponašanje, pa na libc++-u svaki alignReads run tiho obradi
# nula reada i zapiše prazne SAM zapise, uz izlazni status 0.
# Popravak je u ../patches/star_2.7.11b_macos_libcxx_fix.patch (zamjenjuje
# pubsetbuf eksplicitnim stringbuf::str() kopiranjem, što je prenosivo i na
# Linuxu/libstdc++-u se ponaša jednako).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$STAR_DIR/build/STAR-src"
BIN_OUT="$STAR_DIR/bin/STAR"
PATCH_FILE="$STAR_DIR/patches/star_2.7.11b_macos_libcxx_fix.patch"
STAR_VERSION="2.7.11b"

LIBOMP_PREFIX="$(brew --prefix libomp)"

if [ ! -d "$BUILD_DIR" ]; then
  git clone --depth 1 --branch "$STAR_VERSION" https://github.com/alexdobin/STAR.git "$BUILD_DIR"
fi

cd "$BUILD_DIR"
git checkout -- . 2>/dev/null || true
git apply "$PATCH_FILE"

cd "$BUILD_DIR/source"
make clean
make STAR \
  CXXFLAGSextra="-D'COMPILE_FOR_MAC'" \
  CXXFLAGS_SIMD= \
  LDFLAGSextra="-L$LIBOMP_PREFIX/lib -lomp"

mkdir -p "$STAR_DIR/bin"
cp STAR "$BIN_OUT"
echo "== STAR binary izgrađen: $BIN_OUT =="
"$BIN_OUT" --version
