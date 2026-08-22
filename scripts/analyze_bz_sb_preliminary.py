from pathlib import Path

from molecule2fbx.comparison import run_preliminary_comparison


WORKSPACE = Path(__file__).resolve().parents[1]


if __name__ == "__main__":
    run_preliminary_comparison(
        WORKSPACE / "outputs" / "1Bz-LSD_RR",
        WORKSPACE / "outputs" / "1SB-LSD_RR",
        WORKSPACE / "outputs" / "Bz_vs_SB_preliminary_comparison",
        update_ensemble_json=True,
    )
