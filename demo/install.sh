#!/usr/bin/env bash
set -euo pipefail

echo "===================================="
echo "IsoSeq + TAMA + SQANTI3 INSTALLER"
echo "HPC PRODUCTION MODE"
echo "===================================="


# ---------------------------------------------------------
# 0. Non-interactive / HPC-safe settings
# ---------------------------------------------------------
export GIT_TERMINAL_PROMPT=0
export GIT_ASKPASS=""
export PYTHONUNBUFFERED=1


# ---------------------------------------------------------
# 1. Detect platform
# ---------------------------------------------------------
OS=$(uname -s)

if [[ "$OS" != "Linux" ]]; then
    echo "ERROR: This pipeline is Linux-only (HPC required)."
    exit 1
fi


# ---------------------------------------------------------
# 2. Set project directories (IMPORTANT FOR QUOTA SAFETY)
# ---------------------------------------------------------
PROJECT_DIR=$(pwd)

export CONDA_PKGS_DIRS="${PROJECT_DIR}/.conda/pkgs"
export CONDA_ENVS_PATH="${PROJECT_DIR}/.conda/envs"

mkdir -p "$CONDA_PKGS_DIRS"
mkdir -p "$CONDA_ENVS_PATH"

# Redirect pip's cache too -- some conda envs (e.g. SQANTI3's)
# pull in pip packages that build wheels, and pip defaults to
# caching those under $HOME, which can hit small home-dir quotas
# even when the conda cache itself is safely on project storage.
export PIP_CACHE_DIR="${PROJECT_DIR}/.pip/cache"
export TMPDIR="${PROJECT_DIR}/.tmp"

mkdir -p "$PIP_CACHE_DIR"
mkdir -p "$TMPDIR"

echo "Using project conda cache:"
echo "  PKGS: $CONDA_PKGS_DIRS"
echo "  ENVS: $CONDA_ENVS_PATH"
echo "Using project pip cache:"
echo "  PIP_CACHE_DIR: $PIP_CACHE_DIR"
echo "Using project tmpdir:"
echo "  TMPDIR: $TMPDIR"


# ---------------------------------------------------------
# 3. Ensure conda exists
# ---------------------------------------------------------
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found. Load Miniconda first."
    exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"


# ---------------------------------------------------------
# 4. Retry helper (VERY IMPORTANT for HPC networks)
# ---------------------------------------------------------
retry () {
    local n=0
    local max=3
    local delay=5

    until [[ $n -ge $max ]]; do
        "$@" && break
        n=$((n+1))
        echo "Command failed. Retry $n/$max in $delay sec..."
        sleep $delay
    done

    if [[ $n -ge $max ]]; then
        echo "ERROR: Command failed after $max attempts: $*"
        exit 1
    fi
}


# ---------------------------------------------------------
# 5. Function defs: conda env creation + git cloning
#    (defined once, called below -- nothing calls these
#    until AFTER both functions exist)
# ---------------------------------------------------------
create_env () {
    ENV_NAME=$1
    YAML_FILE=$2

    if conda env list | awk '{print $1}' | grep -q "^${ENV_NAME}$"; then
        echo "[SKIP] Conda env exists: ${ENV_NAME}"
    else
        echo "[CREATE] ${ENV_NAME}"
        retry conda env create -f "$YAML_FILE" -n "$ENV_NAME"
    fi
}

clone_repo () {
    local url=$1
    local dir=$2

    if [[ -d "$dir" ]]; then
        echo "[SKIP] Exists: $dir"
    else
        echo "[CLONE] $url"
        retry git clone --depth 1 "$url" "$dir"
    fi
}


# ---------------------------------------------------------
# 6. Clone tools FIRST -- SQANTI3's own conda_env.yml lives
#    inside the repo we're about to clone, so cloning has
#    to happen before we try to create that env from it.
# ---------------------------------------------------------
mkdir -p tools

clone_repo "https://github.com/GenomeRIK/tama.git" tools/tama
clone_repo "https://github.com/ConesaLab/SQANTI3.git" tools/SQANTI3
clone_repo "https://github.com/Magdoll/cDNA_Cupcake.git" tools/cDNA_Cupcake


# ---------------------------------------------------------
# 7. Create conda environments (idempotent)
# ---------------------------------------------------------
echo "Creating conda environments..."

create_env "IsoSeq3" "envs/isoseq3.yml"
create_env "TAMA" "envs/tama.yml"
create_env "SQANTI3" "tools/SQANTI3/SQANTI3.conda_env.yml"


# ---------------------------------------------------------
# 8. cDNA_Cupcake note
#    SQANTI3 only imports from cupcake/sequence/ (pure Python,
#    no compilation needed). Do NOT run `python setup.py build`
#    here -- that also tries to Cythonize cupcake's legacy
#    ICE/TOFU clustering modules (cupcake/tofu/branch/c_branch.pyx,
#    cupcake/ice/find_ECE.pyx), which are unrelated to SQANTI3
#    and still contain Python 2 `print` statements that fail to
#    compile under Python 3 (known upstream issue, unresolved as
#    of this writing: https://github.com/ConesaLab/SQANTI3/issues/208).
#
#    Instead, cDNA_Cupcake just needs to be on PYTHONPATH at
#    runtime -- add this to the SQANTI3 Snakemake rule's shell
#    block later:
#
#      export PYTHONPATH=$PYTHONPATH:{PROJECT_DIR}/tools/cDNA_Cupcake/sequence:{PROJECT_DIR}/tools/cDNA_Cupcake
#
echo "[SKIP] cDNA_Cupcake build (not required -- see comment above; PYTHONPATH is set at rule runtime instead)"

# ---------------------------------------------------------
# 10. Install summary
# ---------------------------------------------------------
echo ""
echo "===================================="
echo "INSTALL COMPLETE"
echo "===================================="
echo ""
echo "Project-based conda cache enabled - OK"
echo "Non-interactive git enabled - OK"
echo "Retry logic enabled - OK"
echo "Idempotent installs enabled - OK"
echo ""
echo "Next:"
echo "  snakemake --cores 8"
echo ""

