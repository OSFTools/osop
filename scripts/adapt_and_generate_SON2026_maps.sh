#!/bin/bash
# =============================================================================
# ADAPTATION SCRIPTS CUSTOM JJA 2026 -> SON 2026
# =============================================================================
# Auteur     : Wafae Badi, DGM Climate Service / RCC-NorAf
# Date       : 5 septembre 2026
# Objectif   : Adapter tous les scripts Python custom de JJA a SON,
#              puis les lancer pour generer les cartes SON 2026 avec
#              la meme qualite visuelle que le rapport JJA 2026 draft2.
# Structure  : RCC_NA_maps/SON2026/{Precipitation, TMP}/
# =============================================================================

set -e
set -u

# =============================================================================
# 0. CONFIG
# =============================================================================
export SCRATCH="/c/Users/wafae/projets/osop-data"

OSOP_SCRIPTS="$HOME/projets/osop/scripts"
MAPS_BASE="$SCRATCH/seafoam/data/master/RCC_NA_maps"
SON_DIR="$MAPS_BASE/SON2026"

mkdir -p "$SON_DIR/Precipitation"
mkdir -p "$SON_DIR/TMP"

cd "$OSOP_SCRIPTS"

# Activer conda osop (avec set +u pour eviter bug _CONDA_SET_CPT_BIN_DIR)
set +u
source $(conda info --base)/etc/profile.d/conda.sh
conda activate osop
set -u

TS=$(date +%Y%m%d_%H%M%S)
LOG_DIR="$OSOP_SCRIPTS/logs/SON2026_maps_${TS}"
mkdir -p "$LOG_DIR"

echo "======================================================================"
echo "  ADAPTATION CARTES JJA 2026 -> SON 2026"
echo "======================================================================"
echo "  Date         : $(date)"
echo "  Scripts      : $OSOP_SCRIPTS"
echo "  Sortie       : $SON_DIR"
echo "  Logs         : $LOG_DIR"
echo "======================================================================"

# =============================================================================
# 1. FONCTION D'ADAPTATION JJA -> SON
# =============================================================================
adapt_script() {
    local src="$1"
    local dst="$2"

    if [ ! -f "$src" ]; then
        echo "  [SKIP] $src introuvable"
        return 1
    fi

    cp "$src" "$dst"

    # Remplacements dans l'ordre
    sed -i 's/JJA2026/SON2026/g' "$dst"
    sed -i "s/'JJA'/'SON'/g" "$dst"
    sed -i 's/"JJA"/"SON"/g' "$dst"
    sed -i 's/_JJA\./_SON\./g' "$dst"
    sed -i 's/_JJA_/_SON_/g' "$dst"
    sed -i 's/JJA Hindcast/SON Hindcast/g' "$dst"
    sed -i 's/JJA hindcast/SON hindcast/g' "$dst"
    sed -i 's/June-July-August/September-October-November/g' "$dst"
    sed -i 's/June - July - August/September - October - November/g' "$dst"
    sed -i 's/Initialization: May 2026/Initialization: August 2026/g' "$dst"
    sed -i 's/init.*May 2026/init August 2026/g' "$dst"
    # init=5 -> init=8 dans les patterns de fichiers
    sed -i 's|_5_234_|_8_234_|g' "$dst"

    echo "  [OK]   $src -> $dst"
    return 0
}

# =============================================================================
# 2. ADAPTATION DES SCRIPTS ESSENTIELS
# =============================================================================
echo ""
echo "======================================================================"
echo "  ETAPE 1 : ADAPTATION DES SCRIPTS JJA -> SON"
echo "======================================================================"

adapt_script "generate_RCC_NA_maps_v5.py"                     "generate_RCC_NA_maps_SON2026.py"
adapt_script "generate_OSOP_multisystem_TP_v3.py"             "generate_OSOP_multisystem_TP_SON2026.py"
adapt_script "generate_OSOP_multisystem_T2m_v1.py"            "generate_OSOP_multisystem_T2m_SON2026.py"
adapt_script "generate_OSOP_multisystem_skill_complete_v2.py" "generate_OSOP_multisystem_skill_SON2026.py"
adapt_script "generate_consistency_map_NA_v3.py"              "generate_consistency_map_SON2026.py"
adapt_script "apply_drymask_TP_JJA2026.py"                    "apply_drymask_TP_SON2026.py"

# =============================================================================
# 3. LANCEMENT DES SCRIPTS ADAPTES
# =============================================================================
echo ""
echo "======================================================================"
echo "  ETAPE 2 : GENERATION DES CARTES SON 2026"
echo "======================================================================"

run_script() {
    local script="$1"
    local logname="$2"
    if [ -f "$script" ]; then
        echo ""
        echo "  >>> Lancement : $script"
        python "$script" 2>&1 | tee "$LOG_DIR/${logname}.log"
        echo "  <<< Fin : $script"
    else
        echo "  [SKIP] $script introuvable"
    fi
}

run_script "generate_RCC_NA_maps_SON2026.py"             "01_rcc_maps"
run_script "generate_OSOP_multisystem_TP_SON2026.py"     "02_multi_TP"
run_script "generate_OSOP_multisystem_T2m_SON2026.py"    "03_multi_T2m"
run_script "generate_OSOP_multisystem_skill_SON2026.py"  "04_multi_skill"
run_script "generate_consistency_map_SON2026.py"         "05_consistency"
run_script "apply_drymask_TP_SON2026.py"                 "06_drymask"

# =============================================================================
# 4. INVENTAIRE FINAL
# =============================================================================
echo ""
echo "======================================================================"
echo "  INVENTAIRE FINAL - CARTES SON 2026"
echo "======================================================================"

echo ""
echo ">>> Racine SON2026/ :"
ls -la "$SON_DIR/"*.png 2>/dev/null || echo "  (aucun PNG racine)"

echo ""
echo ">>> Precipitation/ :"
ls -la "$SON_DIR/Precipitation/"*.png 2>/dev/null || echo "  (vide)"

echo ""
echo ">>> TMP/ :"
ls -la "$SON_DIR/TMP/"*.png 2>/dev/null || echo "  (vide)"

echo ""
echo "======================================================================"
echo "  TERMINE - $(date)"
echo "======================================================================"
echo ""
echo ">>> Ouvrir dans l'explorateur Windows :"
echo "    explorer \"\$(cygpath -w $SON_DIR)\""
echo ""

exit 0
