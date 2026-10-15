"""Export the main and supplementary figures as 600-dpi TIFF (LZW) for journal upload (figures/tiff/)."""
import subprocess
from pathlib import Path

from PIL import Image

FIG = Path(__file__).resolve().parents[1] / "figures"
OUT = FIG / "tiff"
OUT.mkdir(exist_ok=True)
for pdf in sorted(FIG.glob("fig*.pdf")):
    stem = OUT / pdf.stem
    subprocess.run(["pdftoppm", "-r", "600", "-png", "-singlefile", str(pdf), str(stem)], check=True)
    png = stem.with_suffix(".png")
    Image.open(png).convert("RGB").save(stem.with_suffix(".tif"), compression="tiff_lzw", dpi=(600, 600))
    png.unlink()
    print("tiff", stem.with_suffix(".tif").name)
