#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
INSTALL_ROOT=${STERICX_TOOL_ROOT:-"$PROJECT_ROOT/.stericx/tools"}
DOWNLOAD_ROOT="$PROJECT_ROOT/.stericx/downloads"
XTB_ARCHIVE="$DOWNLOAD_ROOT/xtb-210201.tar.xz"
CREST_ARCHIVE="$DOWNLOAD_ROOT/crest-2.12.zip"
XTB_URL="https://github.com/grimme-lab/xtb/releases/download/v6.4.0/xtb-210201.tar.xz"
CREST_URL="https://github.com/crest-lab/crest/releases/download/v2.12/crest.zip"
XTB_SHA256="c31f1c446a5a78a1e5e558b6e688904ae9b0398272b07f260f6e68a18fa27412"
CREST_SHA256="c55e0f075a6223317b33a5f0fae593ce0ad55c1229c382937b0a0c2dcaf72ef6"

require_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "error: required command not found: $1" >&2
    exit 2
  fi
}

verified_download() {
  local url=$1
  local destination=$2
  local expected=$3
  if [[ ! -f "$destination" ]]; then
    echo "downloading $(basename "$destination")"
    curl --fail --location --retry 3 --output "$destination.part" "$url"
    mv "$destination.part" "$destination"
  fi
  local actual
  actual=$(sha256sum "$destination" | awk '{print $1}')
  if [[ "$actual" != "$expected" ]]; then
    echo "error: checksum mismatch for $destination" >&2
    echo "expected=$expected" >&2
    echo "actual=$actual" >&2
    exit 2
  fi
}

require_command curl
require_command sha256sum
require_command tar
require_command unzip

mkdir -p "$DOWNLOAD_ROOT" "$INSTALL_ROOT" "$INSTALL_ROOT/bin"
verified_download "$XTB_URL" "$XTB_ARCHIVE" "$XTB_SHA256"
verified_download "$CREST_URL" "$CREST_ARCHIVE" "$CREST_SHA256"

if [[ ! -x "$INSTALL_ROOT/xtb-6.4.0/bin/xtb" ]]; then
  mkdir -p "$INSTALL_ROOT/xtb-6.4.0"
  tar -xJf "$XTB_ARCHIVE" --strip-components=1 -C "$INSTALL_ROOT/xtb-6.4.0"
fi

if [[ ! -x "$INSTALL_ROOT/crest-2.12/crest" ]]; then
  mkdir -p "$INSTALL_ROOT/crest-2.12"
  unzip -q -o "$CREST_ARCHIVE" -d "$INSTALL_ROOT/crest-2.12"
  crest_binary=$(find "$INSTALL_ROOT/crest-2.12" -type f -name crest -print -quit)
  if [[ -z "$crest_binary" ]]; then
    echo "error: CREST archive did not contain a crest executable" >&2
    exit 2
  fi
  chmod +x "$crest_binary"
fi

XTB_REAL="$INSTALL_ROOT/xtb-6.4.0/bin/xtb"
CREST_REAL=$(find "$INSTALL_ROOT/crest-2.12" -type f -name crest -perm -u+x -print -quit)
if [[ ! -x "$XTB_REAL" || ! -x "$CREST_REAL" ]]; then
  echo "error: extracted toolchain is incomplete" >&2
  exit 2
fi

# Relative symlinks so the tool tree survives being moved or copied between
# machines; an absolute link would dangle at a foreign checkout path.
ln -sfnr "$XTB_REAL" "$INSTALL_ROOT/bin/xtb"
ln -sfnr "$CREST_REAL" "$INSTALL_ROOT/bin/crest"

echo "quantum_toolchain_installed=true"
echo "tool_root=$INSTALL_ROOT"
echo "xtb=$INSTALL_ROOT/bin/xtb"
echo "crest=$INSTALL_ROOT/bin/crest"
"$INSTALL_ROOT/bin/xtb" --version | head -n 3
"$INSTALL_ROOT/bin/crest" --version | head -n 5
